library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(openxlsx)
library(DESeq2)

# PART 1: QUERY AND DATA PREPARATION #
setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/01 - TNBC HER2-low")

# Uploading metadata
metadata <- read.xlsx("TNBC_TCGA_Final.xlsx")
metadata$Sample.ID <- paste0(metadata$Sample.ID, "A")

metadata <- metadata %>%
  mutate(Group = case_when(
    BC_Subtype_HER2 == "4.2_TNBC_HER2low" ~ "HER2-low",
    BC_Subtype_HER2 =="4.1_TNBC_HER2_0" ~ "HER2-zero"))

metadata$Color <- NA
metadata$Color[which(metadata$Group == "HER2-low")] = "orange"
metadata$Color[which(metadata$Group == "HER2-zero")] = "black"

meta_HER2low <- metadata %>%
  filter(Group %in% "HER2-low")
meta_HER2zero <- metadata %>%
  filter(Group %in% "HER2-zero")

# query gene expression data
query_gene_expression <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open",
  barcode = metadata$Sample.ID)

# download data 
GDCdownload(query_gene_expression)

# prepare data
RNA_data <- GDCprepare(query_gene_expression, summarizedExperiment = TRUE)

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(RNA_data))
coldata <- as.data.frame(colData(RNA_data))

# add condition: HER2-low vs HER2-zero
coldata <- coldata %>%
  mutate(condition = case_when(coldata$sample %in% meta_HER2low$Sample.ID ~ "HER2-low",
                               coldata$sample %in% meta_HER2zero$Sample.ID ~ "HER2-zero")) %>%
  mutate(condition = factor(condition, levels = c("HER2-zero", "HER2-low")))

# remove low purity samples #
# Purity data #
data.purity <- read.xlsx("TCGA_cpe_purity.xlsx", startRow = 3)
data.purity <- data.purity %>%
  select(-X8) %>%
  filter(Cancer.type %in% "BRCA")

data.purity <- data.purity %>% # Filter per barcode patients
  filter(data.purity$Sample.ID %in% metadata$Sample.ID)

data.purity <- data.purity %>%
  filter(CPE > 0.6) %>%
  filter(Sample.ID %in% coldata$sample)# Selection of samples high purity

# Filtering metadata
metadata<- metadata %>%
  filter(metadata$Sample.ID %in% data.purity$Sample.ID)
meta_HER2low <- metadata %>%
  filter(Group %in% "HER2-low")
meta_HER2zero <- metadata %>%
  filter(Group %in% "HER2-zero")


# Filtering coldata
coldata <- coldata %>% # Quitamos las muestras del coldata que no tengan una pureza mínima de tumor
  filter(sample %in% data.purity$Sample.ID)

coldata <- coldata %>%
  filter(!is.na(condition))

# Matrix of counts
BRCAMatrix <- assay(RNA_data,"unstranded")

selected_columns <- coldata$barcode
BRCAMatrix <- BRCAMatrix[, selected_columns]


# PART 2: DESEQ2 ANALYSIS - FROM COUNTS MATRIX #
# setting up countData object
dds <- DESeqDataSetFromMatrix(countData = BRCAMatrix,
                              colData = coldata,
                              design = ~ condition)


# Filter out genes with baseMean < 25 identified in the previous analysis and remove them from the matrix
keep <- read.table("BRCA.df.genes.filtered.txt")
keep <- rownames(keep)
dds <- dds[keep,]

# as reference: TNBC
dds$condition <- relevel(dds$condition, ref = "HER2-zero")
dds <- DESeq(dds)

# PART 3: PREPARING RESULTS #
res <- results(dds, alpha = 0.05)
summary(res)

# results without any filters
res.naive <- as.data.frame(res) %>%
  na.omit()

res.naive <- res.naive %>%
  rownames_to_column(var = "gene_id") %>%
  inner_join(gene_metadata, by = "gene_id") %>%
  select(gene_id,gene_name, everything()) %>%
  filter(baseMean > 25) %>%
  mutate(Expression = case_when(padj < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
                         padj < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
                         T ~ "NS")) %>%
  mutate(Expression = factor(Expression, levels = c("Downregulated", "NS", "Upregulated")))

gene_metadata2 <- gene_metadata %>%
  filter(rownames(gene_metadata) %in% res.naive$gene_id)

# Selecting those significant DEGs (p-value < 0.05, baseMean > 25, abs(fold) > 1) 
res.sigs <- res.naive %>%
  filter(padj < 0.05  & abs(log2FoldChange) > 0.5)

res.sigs <- res.sigs %>%
  select(-9:-16) %>%
  mutate(Expression = case_when(padj < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
                              padj < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
                              T ~ "NS"))  %>%
  select(gene_name, gene_id, everything())

table(res.sigs$Expression)


# PART 4: BIOINFORMATIC ANALYZES AND PLOTS #
# Volcano Plot #
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggrepel)

ggplot(res.naive, aes(x = log2FoldChange, y = -log10(padj), color = Expression)) + 
  geom_point(size = 1) + 
  theme_pubr() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  scale_x_continuous(breaks = c(-10, -5, -2.5, 0 , 2.5, 5 , 7.5)) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "dashed", linewidth = 0.8) + 
  geom_vline(xintercept = -0.5, color = "black", linetype = "dashed", linewidth = 0.8) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", linewidth = 0.8)



# Enrichment analyzes ####
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# GSEA (All the genelist)
gene_list <- res.naive$log2FoldChange

# name the vector
names(gene_list) <- gene_metadata2$gene_id

# omit any NA values 
gene_list.2 <- na.omit(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list.2 = sort(gene_list.2, decreasing = TRUE)

#We eliminate ... (variants of genes)
names(gene_list.2) <- sub("\\..*", "", names(gene_list.2))

# Gene Set Enrichment Analysis of Gene Ontology
set.seed(231)
gse.2 <- gseGO(geneList=gene_list.2, 
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

gse.df.2 <- as.data.frame(gse.2)

## Top 10 
gse.paper <- gse.df.2[1:10,]

gse.paper <- gse.paper %>%
  mutate(Description = fct_reorder(Description, desc(NES))) %>%
  arrange(NES)

ggplot(gse.paper, aes(y = Description, x = NES)) +
  geom_segment(aes(yend = Description, xend = 0)) +
  geom_point(size = 6, color ="black", fill = alpha("orange", 0.5), shape = 21) + 
  theme_bw() + 
  labs(x = "NES (Normalized Enrichment Score)")+
  coord_cartesian(xlim = c(-3.5, 0))
  
gseaplot(gse.2,geneSetID=1,title=gse.2$Description[1]) #GSEA



# Enrichment using DEGs #
res.sigs.down <- res.sigs %>%
  filter(Expression %in% "Downregulated")
res.sigs.down <- res.sigs.down$gene_id
res.sigs.down <- sub("\\..*", "", res.sigs.down)


# DOWNREGULATED
# GO
set.seed(46)
enrichment_down <- enrichGO(gene = res.sigs.down,
                            OrgDb = "org.Hs.eg.db",
                            keyType = "ENSEMBL", 
                            ont = "BP",
                            pAdjustMethod = "BH",
                            readable = TRUE)


enrichment_down.df <- as.data.frame(enrichment_down)


enrichment_down.df <- enrichment_down.df %>%
  separate(GeneRatio, into = c("DE", "N"), sep = "/") %>%
  mutate(DE = as.numeric(DE)) %>%
  mutate(N = as.numeric(N)) %>%
  mutate(Ratio = as.numeric(DE/N)*100) %>% # Ratio en tanto por ciento
  arrange((p.adjust))


top.10.down <- enrichment_down.df %>%
  filter(ID %in% c("GO:1902105","GO:0050900", "GO:0002440", "GO:0002768",
         "GO:0051251", "GO:0032943", "GO:0002237", "GO:0002377","GO:0019724", "GO:0045055"))

top.10.down %>%
  mutate(Description = fct_reorder(Description, Ratio)) %>%
  ggplot(aes(x = Description, y = Ratio, fill = p.adjust)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  coord_flip() +
  xlab("") +
  theme_bw() + 
  scale_fill_gradient(low = "red", high = "blue")



# Gene expression of genes controlled by eMQTLs
library(openxlsx)

Bicluster1 <- read.xlsx("Biclusters/Bicluster 1 DEGs.xlsx") 
Bicluster1 <- Bicluster1 %>%
  filter(Gene %in% res.naive$gene_name)
Bicluster1_changes <- inner_join(Bicluster1, res.naive, by = c("Gene" = "gene_name")) %>%
  filter(pvalue < 0.05)

Bicluster2 <- read.xlsx("Biclusters/Bicluster 2 DEGs.xlsx")
Bicluster2 <- Bicluster2 %>%
  filter(Gene %in% res.naive$gene_name)
Bicluster2_changes <- merge(Bicluster2, res.naive, by = c("Gene" = "gene_name")) %>%
  filter(pvalue < 0.05)

Bicluster3 <- read.xlsx("Biclusters/Bicluster 3 DEGs.xlsx")
Bicluster3 <- Bicluster3 %>%
  filter(Gene %in% res.naive$gene_name)
Bicluster3_changes <- merge(Bicluster3, res.naive, by = c("Gene" = "gene_name")) %>%
  filter(pvalue < 0.05)

Bicluster4 <- read.xlsx("Biclusters/Bicluster 4 DEGs.xlsx")
Bicluster4 <- Bicluster4 %>%
  filter(Gene %in% res.naive$gene_name)
Bicluster4_changes <- merge(Bicluster4, res.naive,by = c("Gene" = "gene_name")) %>%
  filter(pvalue < 0.05)

Bicluster5 <- read.xlsx("Biclusters/Bicluster 5 DEGs.xlsx")
Bicluster5 <- Bicluster5 %>%
  filter(Gene %in% res.naive$gene_name)
Bicluster5 <- inner_join(Bicluster5, res.naive, by = c("Gene"= "gene_name"))


# Adjust the p-value here based on the new sample size (N)
library(fuzzySim)
FDR5 <- FDR(Bicluster5$Gene, pvalues = as.data.frame(cbind(Bicluster5$Gene,Bicluster5$pvalue)), simplif = TRUE)
FDR5 <- FDR5 %>%
  rownames_to_column(var = "gene_name")

Bicluster5.2 <- inner_join(Bicluster5, FDR5, by = c("Gene"= "gene_name")) %>%
  arrange(p.adjusted) %>%
  mutate(Expression = case_when(p.adjusted < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
                                p.adjusted < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
                                T ~ "NS"))

Bicluster5_changes <- Bicluster5.2 %>%
  filter(Expression %in% c("Downregulated", "Upregulated"))


# Volcano plot
Bicluster5_full <- inner_join(Bicluster5, FDR5, by = c("Gene"= "gene_name"))

Bicluster5_full <- Bicluster5_full %>%
    mutate(Expression = case_when(p.adjusted < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
                                p.adjusted < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
                                T ~ "NS"))

ggplot(Bicluster5_full, aes(x = log2FoldChange, y = -log10(p.adjusted), color = Expression)) + 
  geom_point(size = 1.2) + 
  theme_pubr() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept = 0.5, color = "black", linetype = "dashed", linewidth = 0.8) + 
  geom_vline(xintercept = -0.5, color = "black", linetype = "dashed", linewidth = 0.8) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", linewidth = 0.8) 


# FC 5
kmplot.5 <- Bicluster5_changes %>%
  arrange(desc(abs(log2FoldChange))) %>%
  filter(gene_type %in% "protein_coding") %>%
  slice_head(n=5) 

write.table(kmplot.5$Gene, file = "Signature KmPlot top 5 genes (Bicluster5).txt", sep = "\t", 
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)


# Top 10
kmplot.10 <- Bicluster5_changes %>%
  arrange(desc(abs(log2FoldChange))) %>%
  filter(gene_type %in% "protein_coding") %>%
  slice_head(n=10)



#### Tumor microenvironment analysis using xCell #### 
# pipeline: Cell-Type Enrichment Analysis of Bulk Transcriptomes Using xCell from Bioinformatics for Cancer
# Immunotherapy book --> https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_19
devtools::install_github("dviraran/xCell")
BiocManager::install(c("GSVA","GSEABase"), version = "3.8")
install.packages("pracma", "utils", " stats ", " MASS ",
                  "digest", "curl", "quadprog")

library(xCell)
xCell.data

# The input for xCell is a gene expression matrix from human mixed samples. 
# It should be read prior to running the xCell functions. 
# The matrix should contain HUGO gene symbols as row names and the
# columns are samples.
library(DESeq2)
DESeq.counts.norm <- counts(dds, normalized = TRUE) # Ojo, si quieres los 60.660 genes ir al principio de todo y darle al dds

gene_metadata_hgcn <- gene_metadata %>%
  select(gene_id, gene_name)

DESeq.counts.norm2 <- DESeq.counts.norm %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  merge(gene_metadata_hgcn, by = "gene_id")

colnames(DESeq.counts.norm2) <- substr(colnames(DESeq.counts.norm2), 1, 16)

DESeq.counts.norm2 <- DESeq.counts.norm2 %>%
  select(-"gene_id") %>%
  select(gene_name, everything())

DESeq.counts.norm2 <- DESeq.counts.norm2 %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  column_to_rownames("gene_name")

genes <- xCell.data$genes
signatures <- xCell.data$signatures

scores <- rawEnrichmentAnalysis(DESeq.counts.norm2, signatures,
                               genes, parallel.type = "SOCK")

fit.vals <- xCell.data$spill$fv # Sequencing based parameters (calibration for raw score)

tscores <- transformScores(scores, fit.vals)

K <- xCell.data$spill$K # compensation

TME <- spillOver(tscores, K, alpha = 0.5) # Alpha equals 0 means no compensation, and 1 means full compensation. In our experiments, a value of 0.5 reduces the dependencies and increases real signals (this is also the default value)

Micro.scores <- microenvironmentScores(TME)
Micro.scores <- Micro.scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "Scores") %>%
  filter(Scores %in% c("ImmuneScore","MicroenvironmentScore", "StromaScore")) %>%
  pivot_longer(-"Scores", names_to = "Sample.ID", values_to = "Values") %>%
  arrange(Scores)

Micro.scores2 <- merge(Micro.scores, metadata, by = "Sample.ID")

Micro.scores2$Group <- factor(Micro.scores2$Group, levels = c("HER2-zero", "HER2-low"))

ggplot(Micro.scores2, aes(x=Group, fill = Group, y=Values)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("HER2-low" = "#7CEBB0", "HER2-zero" = "#884B99")) +
  geom_jitter(color="black", size=1, width = 0.1) +
  facet_wrap(~Scores, scales = "free_y") + 
  theme_bw()

# Statistical test
Micro.scores.immune <- Micro.scores2 %>%
  filter(Scores %in% "ImmuneScore")

Micro.scores.immune.HER20 <- Micro.scores.immune %>%
  filter(Group %in% "HER2-zero")

Micro.scores.immune.HER2low <- Micro.scores.immune %>%
  filter(Group %in% "HER2-low")

# It's not a normal distribution
shapiro.test(Micro.scores.immune.HER20$Values)
shapiro.test(Micro.scores.immune.HER2low$Values)

# Wilcoxon test
wilcox.test(Micro.scores.immune.HER20$Values, Micro.scores.immune.HER2low$Values)


# Statistical test
Micro.scores.Stroma <- Micro.scores2 %>%
  filter(Scores %in% "StromaScore")

Micro.scores.Stroma.HER20 <- Micro.scores.Stroma %>%
  filter(Group %in% "HER2-zero")
Micro.scores.Stroma.HER2low <- Micro.scores.Stroma %>%
  filter(Group %in% "HER2-low")

# It's not a normal distribution
shapiro.test(Micro.scores.Stroma.HER20$Values)
shapiro.test(Micro.scores.Stroma.HER2low$Values)

# Wilcoxon test
wilcox.test(Micro.scores.Stroma.HER20$Values, Micro.scores.Stroma.HER2low$Values)



## Cell deconvolution
TME2 <- TME %>%
  as.data.frame() %>%
  rownames_to_column(var = "Type") %>%
  pivot_longer(-"Type", names_to = "Sample.ID", values_to = "Values") %>%
  arrange(Type)


#
TME2 <- merge(TME2, metadata, by = "Sample.ID")

# Wilcox test
TME3 <- TME2 %>%
  group_by(Type) %>%
  summarise(Wilcox_p_value = wilcox.test(Values ~ Group)$p.value) %>%
  filter(Wilcox_p_value < 0.05) %>%
  ungroup() %>%
  arrange(Type) %>%
  mutate(Wilcox_p_value = round(Wilcox_p_value, 3))


TME2.sig <- TME2 %>%
  filter(Type %in% TME3$Type)

# We remove some redudant type of cells
TME3.sig <- TME2.sig %>%
  filter(!Type %in% c("Tgd cells",
                      "Melanocytes", 
                      "GMP")) %>%
  mutate(Group = factor(Group, levels = c("HER2-zero", "HER2-low")))


ggplot(TME3.sig, aes(x=Group, fill = Group, y=Values)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("HER2-low" = "#7CEBB0", "HER2-zero" = "#884B99")) +
  geom_jitter(color="black", size=1, width = 0.1) +
  facet_wrap(~Type, scales = "free_y") + 
  theme_bw()



# ERBB2 gene expression plot
BRCA.df <- as.data.frame(DESeq.counts.norm)
BRCA.df.genes2 <- BRCA.df %>%
  rownames_to_column(var = "gene_id") %>%
  merge(gene_metadata_hgcn, by = "gene_id") %>%
  select(-gene_id) %>%
  select(gene_name, everything()) %>%
  filter(gene_name %in% c("ERBB2")) %>%
  filter(gene_name %in% res.naive$gene_name) %>%# quitamos aquellos que tengan menos transcritos
  column_to_rownames(var = "gene_name")

BRCA.df.genes2 <- log2(BRCA.df.genes2 + 1)
colnames(BRCA.df.genes2) <- substr(colnames(BRCA.df.genes2), 1, 16)

BRCA.df.genes3 <- BRCA.df.genes2 %>%
  rownames_to_column(var = "gene_name")  %>%
  pivot_longer(-gene_name, names_to = "Sample.ID", values_to = "Expr") %>%
  inner_join(metadata, by = "Sample.ID")

BRCA.df.genes3$Group <- factor(BRCA.df.genes3$Group, levels = c("HER2-zero", "HER2-low"))

ggplot(BRCA.df.genes3, aes(x = Group, y = Expr, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("HER2-low" = "#7CEBB0", "HER2-zero" = "#884B99")) +
  geom_jitter(color="black", size=1, width = 0.1) +
  facet_wrap(~gene_name, scales = "free_y") +
  theme_pubr() + 
  stat_compare_means(method = "t.test") +
  ylab("log2 (RPKMs + 1)")



# Cascade plot
HLAs2 <- res.naive %>%
  filter(grepl("HLA-", gene_name)) %>%
  filter(gene_type %in% "protein_coding") %>%
  mutate(sig = case_when(pvalue < 0.05  ~ "Y",
                         T ~ "N"))

ggplot(HLAs2, aes(x= fct_reorder(gene_name,desc(log2FoldChange)), y = log2FoldChange, fill = sig)) +
  geom_col(color = "black", width = 1) + 
  ylim(-2,2)+
  theme_classic() + 
  scale_fill_manual(values = c("N" = "gray", "Y" = "tomato")) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  labs(y = "log2FC")

