# In this script I investigate whether CD8 T cells from the lesions exibith exhaustion signatures compared to CD8 T cells from the dLN.
# I used the Supplemental table 2 from the paper https://doi.org/10.1016/j.immuni.2019.11.002 that contains:
# the normalized counts of Stem-like, Transitory and Exhausted T cells

# Libraries ----
library(readxl)
library(limma)
library(edgeR)
library(patchwork)
library(tidyverse)

# Import data ----
ggplot_key <- theme(legend.position="right", legend.title = element_text(size = 15),
                    legend.text = element_text(size = 15), axis.title = element_blank(),
                    axis.text = element_text(size = 17))

norm_counts <- read_excel("Exhaustion_signatures/WA_RH2019JI_1-s2.0-S1074761319304601-mmc3.xlsx", 
                                                        col_types = c("skip", "text", "numeric", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric", "numeric"))
norm_counts <- norm_counts %>%
  rename("Gene_name" = "Gene name")

# Importantly, this gene expression matrix has duplicated values. I will remove those from here for downstream analysis.
norm_counts$Gene_name[duplicated(norm_counts$Gene_name)]
norm_counts <- norm_counts %>%
  distinct(Gene_name, .keep_all = TRUE)

norm_counts <- as.matrix(column_to_rownames(norm_counts, "Gene_name"))
myDGEList_Exh <- DGEList(norm_counts)

# Study design ----
WA_RH2019JI_StudyDesign <- read_delim("Exhaustion_signatures/WA_RH2019JI_StudyDesign.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      col_types = cols(Group = col_factor(levels = c("Naive", 
                                                                                     "Stem_like", "Transitory", "Exhausted"))), 
                                      trim_ws = TRUE)

# Preprocessing ----
# Filter for lowly expressed genes and not frequently expressed:
# Mean average of >1cpm in at least 3 samples in the dataset:
cpm_Exh <- cpm(myDGEList_Exh)
keepers_Exh <- rowSums(cpm_Exh>1)>=3
myDGEList.filtered_Exh <- myDGEList_Exh[keepers_Exh,]

# Normalization:
myDGEList.filtered.norm_Exh <- calcNormFactors(myDGEList.filtered_Exh, method = "TMM")
log2.cpm.filtered.norm_Exh <- cpm(myDGEList.filtered.norm_Exh, log=TRUE)

# Function to plot the gene expression distribution patterns of 10 random samples.
distribution_plot <- function(x, step){
  y <- cpm(x, log=TRUE)
  y <- as_tibble(y, rownames = "geneID")
  y.pivot <- pivot_longer(y,
                          cols = !geneID,
                          names_to = "samples",
                          values_to = "expression") 
  y.pivot %>%
    ggplot(., aes(x=samples, y=expression)) +
    geom_violin(trim = FALSE, show.legend = FALSE, adjust = 5) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 6, 
                 color = "black", 
                 show.legend = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_blank()) +
    ggplot_key +
    geom_hline(yintercept = 1, color="red") +
    labs(y="log2 expression", x = "Random 10 samples",
         title=step)
}

p1 <-distribution_plot(myDGEList_Exh, step = "A) Raw data")
p2 <-distribution_plot(myDGEList.filtered_Exh, step = "B) Filtered data")
p3 <-distribution_plot(myDGEList.filtered.norm_Exh, step = "C) Filtered\nNormalized data")

p1 + p2 + p3

# DEGs ----
design_Exh <- model.matrix(~0 + WA_RH2019JI_StudyDesign$Group)
colnames(design_Exh) <- levels(WA_RH2019JI_StudyDesign$Group)
v.DEGList.filtered.norm_Exh <- voom(myDGEList.filtered.norm_Exh, design_Exh, plot = FALSE)
fit_Exh <- lmFit(v.DEGList.filtered.norm_Exh, design_Exh)
contrast.matrix_Exh <- makeContrasts(Stem_like = Stem_like - Naive,
                                 Transitory = Transitory - Naive,
                                 Exhausted = Exhausted - Naive,
                                 levels=design_Exh)
fits_Exh <- contrasts.fit(fit_Exh, contrast.matrix_Exh)
ebFit_Exh <- eBayes(fits_Exh)
myTopHits_Stem_like <- topTable(ebFit_Exh, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits_Transitory <- topTable(ebFit_Exh, adjust ="BH", coef=2, number=40000, sort.by="logFC")
myTopHits_Exhausted <- topTable(ebFit_Exh, adjust ="BH", coef=3, number=40000, sort.by="logFC")

myTopHits_Stem_like <- myTopHits_Stem_like %>% rownames_to_column("geneSymbol")
myTopHits_Transitory <- myTopHits_Transitory %>% rownames_to_column("geneSymbol")
myTopHits_Exhausted <- myTopHits_Exhausted %>% rownames_to_column("geneSymbol")

# Export signatures and add to /Dropbox/CamilaAnalysis/genesets_scRNAseq:
myTopHits_Stem_like %>%
  arrange(desc(logFC), adj.P.Val) -> Stem_like_sig

myTopHits_Transitory %>%
  arrange(desc(logFC), adj.P.Val) -> Transitory_sig

myTopHits_Exhausted %>%
  arrange(desc(logFC), adj.P.Val) -> Exhausted_sig

write_tsv(as.data.frame(Stem_like_sig$geneSymbol[1:100]), "Exhaustion_signatures/Stem_like_sig.txt")
write_tsv(as.data.frame(Transitory_sig$geneSymbol[1:100]), "Exhaustion_signatures/Transitory_sig.txt")
write_tsv(as.data.frame(Exhausted_sig$geneSymbol[1:100]), "Exhaustion_signatures/Exhausted_sig.txt")

# Analyses on the CD8 T cells ----
library(GSEABase)
library(GSVA)
genesets_Exh <- getGmt("../../../CamilaAnalysis/genesets/genesets_RNAseq.gmt", geneIdType=SymbolIdentifier())

# ssGSEA Hypoxia
GSVA.res_CD8_2 <- gsva(log2.cpm.filtered.norm,
                     genesets_Exh,
                     mx.diff=FALSE,
                     method="ssgsea")

heatmap.2(GSVA.res_CD8_2[c(6,7,8),],
          dendrogram = "none",
          Colv = F, Rowv = F,
          col=colormat,
          scale='row', #labCol=F,
          density.info="none", trace="none",
          cexRow=1, cexCol=.7, margins=c(19,20)) 


as.data.frame(GSVA.res_CD8_2) %>%
  rownames_to_column("geneset") %>%
  pivot_longer(!geneset, names_to = "sample_title",
               values_to = "Expression") %>%
  left_join(targets) %>%
  filter(geneset %in% c("Stem_like","Transitory","Exhausted")) %>%
  ggplot(., aes(x=tissue, y=Expression, color=tissue)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=0.5) +
  theme_classic() +
  #scale_color_manual(values = rev(Dark24[1:2])) + 
  ggplot_key + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "t.test",# label.y = 8.2,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 7, color="black") +
  facet_wrap(. ~ geneset)

as.data.frame(GSVA.res_CD8_2[c(6,7,8),]) %>%
  rownames_to_column("geneset") %>%
  mutate(Lesion = (ear_total_rep1_b2 + ear_total_rep2_b2 + ear_total_rep3_b2)/3,
       dLN = (LN_CD44hi_rep1_b2 + LN_CD44hi_rep2_b2 + LN_CD44hi_rep3_b2)/3) %>%
  dplyr::select(geneset, Lesion,dLN) %>%
  pivot_longer(!geneset, names_to = "sample_group",
               values_to = "AVG") %>%
  ggplot(., aes(x=sample_group, y=fct_inorder(geneset), color=AVG)) +
  geom_point(shape=15, size=15) +
  theme_classic() +
  scale_color_gradientn(colours = colormat6) +
  ggplot_key +
  theme()

