# RNAseq analyis on T CD8 cells purified from the ear and lymph nodes of L .major-infected C57BL/6 mouse, 5w post-infection.
# This document analyze in R environment the transcriptome of CD8 T cells sorted from the ear and lymph nodes of L. major-infected mice (in replicates):
# These samples were sequenced in 2020 in only one run 389. This one had a good quality and quantity. (batch 2)

# Libraries ----
#library(EnsDb.Hsapiens.v86)
#library(corrplot)
#library(reshape2)

#library(patchwork)
library(ggrepel)
library(ggthemes)
#library(Polychrome)
#library(vegan)
#library(FinCal)
#library(ggExtra)
library(gplots)
library(ggpubr)
#library(Hmisc)
#library(plotly)
#library(JLutils)
#library(broom)
#library(reshape2)
library(ggforce)
#library(msigdbr)
#library(gprofiler2)
#library(gt)
#library(UpSetR)

library(tidyverse)

#'%ni%' <- Negate('%in%')

# Color scheme ----
library(scales)
#Nandas <- c("#73B85C","#6A216F","#F7CD6E","#5E807A","#3275B4","#C7598F","#A9C88B","#172866","#D89559","#B46969","#906997","#404F59")
Nandas <- c("#589C48","#7BB662","#FBB149","#F58024","#733381","#994FB2","#1C427F")
Nandas2 <- c("#589C48","#FBB149","#733381","#1C427F")
show_col(Nandas)

# Export
save(Nandas, file = "Color_palette/Nandas")
load("Color_palette/Nandas")

# Importing and process the CD8s dataset ----
# Metadata:
targets <- read_delim("MouseCD8_mapping_studyDesign/studyDesign_CD8s.txt", 
                      "\t", escape_double = FALSE,
                      trim_ws = TRUE)

targets <- targets %>%
  filter(batch == "2")

# Tximports:
library(biomaRt)
library(tximport)
path1 <- file.path("MouseCD8_mapping_studyDesign/Kallisto_mapping_sept2021/",targets$sample_title, "abundance.h5")
anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
Tx <- getBM(attributes=c('ensembl_transcript_id',
                         'external_gene_name'),
            mart = anno)
Tx <- as_tibble(Tx)

Txi_gene1 <- tximport(path1,
                      type = "kallisto", 
                      tx2gene = Tx, 
                      txOut = FALSE, 
                      countsFromAbundance = "lengthScaledTPM",
                      ignoreTxVersion = TRUE)

#write_tsv(as.data.frame(Txi_gene1$counts), "Txitest.txt")
#save(Txi_gene1, file = "RStudio_outputs/data/Txi_gene1")
myDGEList <- DGEList(Txi_gene1$counts)
colnames(myDGEList) <- targets$sample_title
myDGEList <- myDGEList[-1,]
write_tsv(rownames_to_column(as.data.frame(myDGEList$counts), "geneSymbol"),
          "MouseCD8_mapping_studyDesign/CD8s_counts_raw.txt")

#rm(Txi_gene1)

# Data processing:
## Filtering
library(edgeR)
cpm <- cpm(myDGEList)

keepers <- rowSums(cpm>1)>=3 #n of replicates
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList)
dim(myDGEList.filtered)

## Normalization
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
#cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=FALSE)
colnames(log2.cpm.filtered.norm) <- targets$sample_title
#colnames(myDGEList.filtered.norm) <- targets$sample_title
write_tsv(rownames_to_column(as.data.frame(log2.cpm.filtered.norm), "geneSymbol"),
          "MouseCD8_mapping_studyDesign/CD8s_filtered_log2CPM.txt")
write_tsv(rownames_to_column(as.data.frame(log2.cpm.filtered.norm), "geneSymbol"),
          "Submission/CD8s_filtered_log2CPM.txt")

log2.cpm.filtered.norm_caps <- log2.cpm.filtered.norm
rownames(log2.cpm.filtered.norm_caps) <- str_to_upper(rownames(log2.cpm.filtered.norm))
write_tsv(rownames_to_column(as.data.frame(log2.cpm.filtered.norm_caps), "geneSymbol"),
          "MouseCD8_mapping_studyDesign/log2.cpm.filtered.norm_caps.txt")

## Data not filtered
# not filtered, normalized data:
notfilt_norm <- calcNormFactors(myDGEList, method = "TMM")
notfilt_norm_log2cpm <- cpm(notfilt_norm, log=TRUE)
write_tsv(rownames_to_column(as.data.frame(notfilt_norm_log2cpm), "geneSymbol"),
                             "MouseCD8_mapping_studyDesign/CD8s_NOTfiltered_log2CPM.txt")
write_tsv(rownames_to_column(as.data.frame(notfilt_norm_log2cpm), "geneSymbol"),
          "Submission/CD8s_NOTfiltered_log2CPM.txt")

#rownames(notfilt_norm_log2cpm) <- toupper(rownames(notfilt_norm_log2cpm))
#write.table(notfilt_norm_log2cpm, "RStudio_outputs/notfilt_norm_log2cpm_CD8uppper.txt", sep = "\t", quote = FALSE) 
#notfilt_norm_cpm <- cpm(notfilt_norm, log=FALSE)
#save(notfilt_norm_cpm, file = "RStudio_outputs/data/notfilt_norm_cpm")
#save(notfilt_norm_log2cpm, file = "RStudio_outputs/data/notfilt_norm_log2cpm")

# PCA ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

as.data.frame(pca.res$x) %>%
  rownames_to_column("sample_title") %>%
  left_join(targets) %>%
  ggplot(., aes(x=PC1, y=PC2, color=tissue)) +
  theme_classic() + scale_color_manual(values = c(Nandas[5],Nandas[1])) +
  geom_hline(yintercept = 0, color="light gray") + geom_vline(xintercept = 0, color="light gray") +
  geom_point(size=7) +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  #stat_ellipse(level = 0.95) +
  xlab(paste("PC1 -",pc.per[1],"%")) +
  ylab(paste("PC2 -",pc.per[2],"%")) + ylim(-130,130) + xlim(-130,130) +
  coord_fixed()

# DE analysis Lesion vs dLN ----
# design matrix:
design <- model.matrix(~0 + factor(targets$tissue))
colnames(design) <- levels(factor(targets$tissue))

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(tissue = ear - lymphNodeH,
                                 levels=design)
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

# Volcano:
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits <- as_tibble(myTopHits, rownames = "geneSymbol")
write_tsv(myTopHits %>%
          mutate(logFC_EARvsDLN = logFC) %>%
          dplyr::select(geneSymbol, logFC_EARvsDLN, P.Value,adj.P.Val) %>%
          arrange(desc(logFC_EARvsDLN)),
          "Submission/SupplementalTable1_CD8s_DGEs.txt") 

# Column of significant genes:
sig_genes_up <- myTopHits %>%
  filter(adj.P.Val < 0.05) %>% filter(logFC > 0.59) %>%
  arrange(adj.P.Val)
sig_genes_down <- myTopHits %>%
  filter(adj.P.Val < 0.05) %>% filter(logFC < -0.59) %>%
  arrange(adj.P.Val)
sig_genes_up <- sig_genes_up$geneSymbol
sig_genes_down <- sig_genes_down$geneSymbol
write_tsv(as_tibble(sig_genes_up), "GO_CD8s_results/up_ear.txt") # to be used in GO
#write_tsv(as_tibble(sig_genes_down), "RStudio_outputs/up_ln.txt") # to be used in GO

myTopHits$col0 <- myTopHits$geneSymbol
col0_s <- myTopHits$col0 %in% sig_genes_up
col0_s2 <- myTopHits$col0 %in% sig_genes_down
myTopHits$col0[!col0_s & !col0_s2] <- NA

myTopHits$col1 <- myTopHits$geneSymbol
col1_s <- myTopHits$col1 %in% sig_genes_up
col1_s2 <- myTopHits$col1 %in% sig_genes_down
myTopHits$col1[col1_s] <- "sig_up"
myTopHits$col1[col1_s2] <- "sig_down"
myTopHits$col1[!col1_s & !col1_s2] <- "notsig"

myTopHits$col2 <- myTopHits$col0
col2_s <- myTopHits$col2 %in% NA
col2_s2 <- myTopHits$col2 %in% sig_genes_up
col2_s3 <- myTopHits$col2 %in% sig_genes_down
myTopHits$col2[col2_s] <- "notsig"
myTopHits$col2[col2_s2 | col2_s3] <- "sig"

#myTopHits$col3 <- myTopHits$geneSymbol
#col3_v <- myTopHits$col3 %in% c("GZMB","GNLY","PRF1","IL1B")
#myTopHits$col3[!col3_v] <- NA

myTopHits$col4 <- myTopHits$geneSymbol
col4_v <- myTopHits$col4 %in% sig_genes_up
col44_v <- myTopHits$col4 %in% sig_genes_down
myTopHits$col4[!col4_v & !col44_v] <- NA

ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC,
                      color=col1,size=col1)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", Nandas[5], Nandas[3])) +
  scale_size_manual(values = c(0.5,1,1)) +
  geom_text_repel(aes(label = col4), size = 3.5, fontface="bold",color="black",
                  max.overlaps = 200) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = .59) + geom_vline(xintercept = -.59) +
  annotate("text", x=-5, y=3, label=paste(length(sig_genes_down)), size=6, color=Nandas[5], fontface="bold") +
  annotate("text", x=5, y=3, label=paste(length(sig_genes_up)), size=6, color=Nandas[3], fontface="bold") +
  xlab("logFC ear vs. ln (CD44H)")

# Heatmap:
colormat <- colorRampPalette(c(Nandas[1],"white",Nandas[5]))(100)
colormat5 <- colorRampPalette(c(Nandas[1],"white",Nandas[5]))(5)
colormat6 <- colorRampPalette(c("white",Nandas[5]))(5)
colormat2 <- colorRampPalette(c(Nandas[3],"white",Nandas[7]))(100)
colormat3 <- colorRampPalette(c("gray",Nandas[5]))(100)
colormat4 <- colorRampPalette(c("gray","#FF5000"))(100)

heatmap.2(log2.cpm.filtered.norm[sig_genes_up[1:20],],
          dendrogram = "none", Colv = F, Rowv = F,
          col=colormat,
          #RowSideColors = Goear_colors,
          scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=.7, margins=c(5,20))
heatmap.2(log2.cpm.filtered.norm[sig_genes_down[1:20],],
          dendrogram = "none", Colv = F, Rowv = F,
          col=colormat,
          #RowSideColors = Goear_colors,
          scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=.7, margins=c(5,20)) 


# TFs ----
genelist1 <- c("Zeb2","Prdm1","Tbx21","Irf4","Stat4","Batf","Id2")

# Order genes:
TF_order_FC <- myTopHits %>%
  filter(geneSymbol %in% genelist1) %>%
  arrange(desc(logFC))
TF_order_FC <- TF_order_FC$geneSymbol

heatmap.2(log2.cpm.filtered.norm[TF_order_FC,],
          dendrogram = "none",
          Colv = F, Rowv = F,
          col=colormat,
          scale='row', #labCol=F,
          density.info="none", trace="none",
          cexRow=1, cexCol=.7, margins=c(17,20)) 

# Hypoxia sig ----
# Approaches used:
# 1) GSEA using Fernanda's genesets of interest.
# Nanda looked in the literature/internet genesets for investigation and I made a .gmt file manually.
# Description of pathways:
library(GSEABase)
library(GSVA)
genesets <- getGmt("compiledHypoxiaGSEA.gmt", geneIdType=SymbolIdentifier())

# ssGSEA Hypoxia
GSVA.res_CD8 <- gsva(log2.cpm.filtered.norm_caps,
                     genesets,
                     mx.diff=FALSE,
                     method="ssgsea")

heatmap.2(GSVA.res_CD8,
          dendrogram = "none",
          Colv = F, Rowv = F,
          col=colormat,
          scale='row', #labCol=F,
          density.info="none", trace="none",
          cexRow=1, cexCol=.7, margins=c(19,20)) 


# ssGSEA Exhaustion
# From Wherry's paper: MolecularSignatureofCD8+TCellExhaustion duringChronicViralInfection
# DOI10.1016/j.immuni.2007.09.006
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN.html
exh <- getGmt("GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN.v2023.2.Hs.gmt", geneIdType=SymbolIdentifier())

GSVA.res_CD8_exh <- gsva(log2.cpm.filtered.norm_caps,
                     exh,
                     mx.diff=FALSE,
                     method="ssgsea")


