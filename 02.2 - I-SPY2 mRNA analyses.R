# COHORT I-SPY2 - GSE196096
setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/01 - TNBC HER2-low/I-SPY2")

# Loading data and metadata
load("Environment_ISPY2_GSE196096.Rdata")

# Gene expression -->  ~19,000 genes assayed on Agilent 44K
# A ComBat batch correction process was applied to the raw data on a group of patients (~800). 

library(tidyverse)
library(openxlsx)
library(conflicted)
conflict_prefer_all("dplyr")

myCombat <- data
metadata.filt <- read.xlsx("Metadata TNBC patients with HER2low info.xlsx")
colnames(myCombat) <- gsub("X", "", colnames(myCombat))
metadata.filt$Patient.Identifier <- as.character(metadata.filt$Patient.Identifier)

metadata.filt <- metadata.filt %>%
  mutate(Color = case_when(HER2_group %in% "HER2-zero" ~ "#884B99",
                           HER2_group %in% "HER2-low" ~ "#7CEBB0" ))

# Gene annotations
RNA.anot <- read_tsv("Gene_metadata_annotations.txt")  %>%
  filter(gene_name %in% rownames(myCombat)) %>%
  select(gene_name, gene_id, gene_type)

 myCombat.filt <- myCombat[, metadata.filt$Patient.Identifier]

### Separate between HER2-low and HER2-zero conditions
myCombat.HER2low <- myCombat.filt[,metadata.filt$Patient.Identifier[metadata.filt$HER2_group == "HER2-low"]]
myCombat.HER2zero <- myCombat.filt[,metadata.filt$Patient.Identifier[metadata.filt$HER2_group == "HER2-zero"]]

# Normality test
shapiro.test(as.numeric(myCombat.filt[1,]))

# We will use a t-test to calculate statistical differences, then adjust for FDR

# We will calculate log2FC as mean(log2(condition1)) - mean(log2(condition2))
# The data is already in log2 format https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5826329

# Calculate the mean of each gene
HER2low.mean <- apply(myCombat.HER2low, 1, mean)
HER2zero.mean <- apply(myCombat.HER2zero, 1, mean)
baseMean <- apply(myCombat.filt, 1, mean)

log2FC <- HER2low.mean - HER2zero.mean 

# Compute statistical significance 
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(myCombat.HER2low)) {
  x = myCombat.HER2low[i,]
  y = myCombat.HER2zero[i,]
  
  t = t.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}

p.adj <- p.adjust(pvalue, method = "fdr", n = length(pvalue)) # Correction by False Discovery Rate

res.naive <- as.data.frame(cbind(baseMean, pvalue, p.adj, log2FC)) %>%
  rownames_to_column(var = "gene_name") %>%
  drop_na() %>%
  left_join(RNA.anot, by = "gene_name") %>%
  filter(!duplicated(gene_name))


# 1 - GSEA
# Enrichment analyzes #
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# 2. WITH ALL THE GENES 
# To make the list: we want the log2 fold change 
GSEA.list.df <- res.naive %>%
  drop_na(gene_id)

GSEA.list <- GSEA.list.df$log2FC
names(GSEA.list) <- GSEA.list.df$gene_id

# omit any NA values 
GSEA.list <- na.omit(GSEA.list)

# sort the list in decreasing order (required for clusterProfiler)
GSEA.list <- sort(GSEA.list, decreasing = TRUE)

#We eliminate ... (variants of genes)
names(GSEA.list) <- sub("\\..*", "", names(GSEA.list))

# Gene Set Enrichment Analysis of Gene Ontology
set.seed(231)
GSEA.GO <- gseGO(geneList=GSEA.list, 
               ont ="BP", 
               keyType = "ENSEMBL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = "org.Hs.eg.db",
               nPermSimple = 10000,
               pAdjustMethod = "BH",
               seed = TRUE,
               eps = 0)

GSEA.GO.df <- as.data.frame(GSEA.GO)

# Plot
gse.paper <- GSEA.GO.df[1:10,]

gse.paper <- gse.paper %>%
  mutate(Description = fct_reorder(Description, desc(NES))) %>%
  arrange(NES)

ggplot(gse.paper, aes(y = Description, x = NES)) +
  geom_segment(aes(yend = Description, xend = 0)) +
  geom_point(size = 6, color ="black", fill = alpha("orange", 0.5), shape = 21) + 
  theme_bw() + 
  labs(x = "NES (Normalized Enrichment Score)")+
  coord_cartesian(xlim = c(-2.4, 0))

gseaplot(GSEA.GO,geneSetID=1,title=GSEA.GO$Description[1]) #GSEA


# Cascade plot HLAs
HLAs2 <- res.naive %>%
  filter(grepl("HLA-", gene_name)) %>%
  filter(gene_type %in% "protein_coding") %>%
  mutate(sig = case_when(pvalue < 0.05  ~ "Y",
                         T ~ "N"))

ggplot(HLAs2, aes(x= fct_reorder(gene_name,desc(log2FC)), y = log2FC, fill = sig)) +
  geom_col(color = "black", width = 1) + 
  ylim(-0.5, 0.5) +
  theme_classic() + 
  scale_fill_manual(values = c("N" = "gray", "Y" = "tomato")) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1))