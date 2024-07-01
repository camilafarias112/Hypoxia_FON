# In this script I evaluate hypoxia signatures in bulk RNA-seq of biopsies from patients.
# Amorim's 2023 dataset

geneExpression <- read_delim("Human_files/Amorim_LeishMicrobiome_processed_norm_filt_log2cpm.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
geneExpression <- as.matrix(column_to_rownames(geneExpression,"geneSymbol"))

# Study design:
studyDesign <- read_delim("Human_files/Amorim_LeishMicrobiome_SupplemmentalTable1_LeishOmics_StudyDesign.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
#studyDesign <- studyDesign %>%
#  mutate(PatientID = str_replace_all(PatientID, "^", "CL"))


studyDesign <- studyDesign %>% #remove samples that are not from the RNA-seq dataset
  filter(`Bacterial load (qPCR)` != "NA")

# Hypoxia sig ----
GSVA.res <- gsva(geneExpression, genesets, mx.diff=FALSE, method="ssgsea")

# Factoring groups:
group_factor <- as.factor(colnames(geneExpression))
color.mapA <- function(group_factor) { if (group_factor %in% studyDesign$PatientID) Nandas[5]
                                           else "dark gray"}
color.mapA <- unlist(lapply(group_factor, color.mapA))

# Clustering:
distance <- dist(t(GSVA.res), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "ward.D2") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters)

heatmap.2(GSVA.res,
          Rowv = F,
          Colv=as.dendrogram(clusters),
          col=colormat2,
          ColSideColors = color.mapA,
          scale='row', labCol = NA,
          density.info="none", trace="none",
          cexRow=1.6, cexCol=0.5, margins=c(25,15))

group_hc <- heatmap.2(GSVA.res, Rowv = F, Colv=as.dendrogram(clusters),
                      scale='row', labCol = NA, density.info="none", trace="none")

# Extract groups and categorizing samples according to their hypoxia scores:
colnames(GSVA.res)[group_hc$colInd]
cluster_heat_order <- as.data.frame(colnames(GSVA.res)[group_hc$colInd])
colnames(cluster_heat_order)[1] <- "PatientID"
cluster_heat_order <- cluster_heat_order %>%
  mutate(group_hc1 = case_when(
    PatientID %in% cluster_heat_order$PatientID[1:29] ~ "Hypoxic",
    PatientID %in% colnames(geneExpression)[52:57] ~ "Intact skin",
    TRUE ~ "Non-hypoxic")) %>%
  
  mutate(group_hc2 = case_when(
    PatientID %in% cluster_heat_order$PatientID[1:17] ~ "Hypoxic +",
    PatientID %in% cluster_heat_order$PatientID[18:29] ~ "Hypoxic +++",
    PatientID %in% colnames(geneExpression)[52:57] ~ "Intact skin",
    TRUE ~ "Non-hypoxic"))

# merge with studyDesign
studyDesign <- studyDesign %>% full_join(cluster_heat_order)

# Get samples labels:
hypoxia_non <- studyDesign %>% filter(group_hc2 %in% "Non-hypoxic")
hypoxia_non <- hypoxia_non$PatientID

hypoxia_1 <- studyDesign %>% filter(group_hc2 %in% "Hypoxic +")
hypoxia_1 <- hypoxia_1$PatientID

hypoxia_2 <- studyDesign %>% filter(group_hc2 %in% "Hypoxic +++")
hypoxia_2 <- hypoxia_2$PatientID

# Add socres to Study design:
studyDesign <- as.data.frame(GSVA.res) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(!Pathway, names_to = "PatientID", values_to = "Hypoxia_score") %>%
  filter(Pathway == "Hypoxia_Harris") %>%
  full_join(studyDesign) %>%
  mutate(group = case_when(
    group_hc2 == "Intact skin" ~ "HS",
    TRUE ~ "CL"))


studyDesign %>%
ggplot(., aes(x=group, y=Hypoxia_score, color=group_hc2)) +
  #geom_violin() +
  geom_jitter() +
  theme_classic() + scale_fill_manual(values = c(Nandas[1],Nandas[3],Nandas[5],"#C75836")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 18, vjust = 0.5, angle = 90, hjust = 1), 
        axis.text.y = element_text(size = 18), axis.title = element_text(size = 20),
        legend.position="right", legend.text = element_text(size = 20), legend.title = element_blank()) +
  xlab("") + ylab(paste("% explained")) + 
  NULL

# PCA and Permanova ----
#geneExpression_mod <- geneExpression[,1:51]
#geneExpression_mod <- geneExpression[,c(hypoxia_1,hypoxia_2)]
#geneExpression_mod <- geneExpression[,c(hypoxia_non,hypoxia_1)]
geneExpression_mod <- geneExpression[,c(hypoxia_non,hypoxia_2)]

pca.res_hy <- prcomp(t(geneExpression_mod), scale.=F, retx=T)
pc.var_hy<-pca.res_hy$sdev^2
pc.per_hy<-round(pc.var_hy/sum(pc.var_hy)*100, 1)

hypoxia_pca <- as.data.frame(pca.res_hy$x) %>%
  rownames_to_column("PatientID") %>%
  left_join(studyDesign)

hypoxia_pca %>%
  ggplot(., aes(x=PC1, y=PC2, color=group_hc2)) +
  theme_classic() + scale_color_manual(values = c(Nandas[7],Nandas[3])) +
  geom_hline(yintercept = 0, color="light gray") + geom_vline(xintercept = 0, color="light gray") +
  geom_point(size=5) +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  stat_ellipse(level = 0.95) +
  xlab(paste("PC1 -",pc.per_hy[1],"%")) +
  ylab(paste("PC2 -",pc.per_hy[2],"%")) +
  coord_fixed()

#library(vegan)
dist_vegan_hy <- vegdist(t(2^geneExpression_mod), method = "bray")
vegan_hy <- adonis2(dist_vegan_hy ~ hypoxia_pca$group_hc2,
                 data=hypoxia_pca, permutations = 999, method="bray")
vegan_hy

# DEG Hypoxia ----
# Hypoxia vs. baseline:
hypoxia_non_2 <- studyDesign %>% filter(PatientID %in% c(hypoxia_non,hypoxia_2))
hypoxia_non_2x <- hypoxia_non_2$PatientID

# design matrix:
design_hy <- model.matrix(~0 + factor(hypoxia_non_2$group_hc2))
colnames(design_hy) <- c("hypoxia_2","non_hypoxic")

v.DEGList.filtered.norm_hy <- voom(2^geneExpression[,hypoxia_non_2x], design_hy)
fit_hy <- lmFit(v.DEGList.filtered.norm_hy, design_hy)
contrast.matrix_hy <- makeContrasts(Hypoxia = hypoxia_2 - non_hypoxic,
                                 levels=design_hy)
fits_hy <- contrasts.fit(fit_hy, contrast.matrix_hy)
ebFit_hy <- eBayes(fits_hy)

# MyTopHits and Volcano ----
myTopHits_hy <- topTable(ebFit_hy, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits_hy <- as_tibble(myTopHits_hy, rownames = "geneSymbol")

write_tsv(myTopHits_hy %>%
            mutate(logFC_HypoxiavsNormoxia = logFC) %>%
            dplyr::select(geneSymbol, logFC_HypoxiavsNormoxia, P.Value,adj.P.Val) %>%
            arrange(desc(logFC_HypoxiavsNormoxia)),
          "Submission/SupplementalTable3_Human_Hypoxia_DGEs.txt") 

# Colunm of significant genes:
sig_genes_up_hy <- myTopHits_hy %>%
  #arrange(desc(logFC),adj.P.Val)
  filter(adj.P.Val < 0.01) %>% filter(logFC > 2)

sig_genes_up_hy_forGO <- myTopHits_hy %>%
  #arrange(desc(logFC),adj.P.Val)
  filter(adj.P.Val < 0.01) %>% filter(logFC > 0.59)
write_tsv(as.data.frame(sig_genes_up_hy_forGO$geneSymbol),
          "Submission/HumanHypoxicDGEs_forGO.txt") 

myTopHits_hy %>%
  mutate(col_show = case_when(
    geneSymbol %in% sig_genes_up_hy$geneSymbol ~ geneSymbol,
    TRUE ~ NA_character_)) %>%
  mutate(col_show2 = case_when(
    geneSymbol %in% sig_genes_up_hy$geneSymbol ~ "top",
    TRUE ~ "others")) %>%

ggplot(., aes(y=-log10(adj.P.Val), x=logFC,
                      color=col_show2,size=col_show2)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", Nandas[7])) +
  scale_size_manual(values = c(0.5,2,2)) +
  geom_text_repel(aes(label = col_show), size = 5, fontface="bold",
                  max.overlaps = 200, force = 50) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = .59) + geom_vline(xintercept = -.59) +
  xlim(0,4) +
  xlab("logFC Hypoxic vs. Normoxic")

# Cell estimations ----
# MCP counter:
library(immunedeconv)
res_mcp_counter <- deconvolute(2^geneExpression, "mcp_counter")

res_mcp_counter_g <- res_mcp_counter %>%
  gather(PatientID, cell_score, -cell_type) %>%
  filter(cell_type %in% "Neutrophil") %>%
  rename(Neutrophil_score = cell_score) %>%
  dplyr::select(PatientID, Neutrophil_score)

studyDesign <- studyDesign %>%
  left_join(res_mcp_counter_g)
  
# Violins:
studyDesign %>%
  filter(Neutrophil_score != "NA") %>%
  ggplot(., aes(x=group_hc2, y=log10(Neutrophil_score), color=group_hc2)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("dark gray",Nandas[7],Nandas[9],"#C75836")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 20),
        legend.position="right", legend.text = element_text(size = 20), legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 18)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(method = "wilcox.test", #label.y = 0.9,
  #                   aes(label = ..p.signif..), comparisons = my_comparisons,
  #                   label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("MCPcounter score (log10)")) +
  NULL

selected_genes <- as.data.frame(geneExpression) %>%
  rownames_to_column("geneSymbol") %>%
  filter(geneSymbol %in% c("HIF1A","GZMB","PRDM1")) %>%
  pivot_longer(!geneSymbol, names_to = "PatientID", values_to = "Expression") %>%
  pivot_wider(names_from = "geneSymbol", values_from = "Expression")

ExportHuman <- studyDesign %>%
  left_join(selected_genes) %>%
  mutate(Hypoxia_groups = case_when(
    group_hc2 %in% "Hypoxic +++" ~ "Hypoxia",
    group_hc2 %in% "Hypoxic +" ~ "Intermediate",
    group_hc2 %in% "Non-hypoxic" ~ "Baseline",
    group_hc2 %in% "Intact skin" ~ "HS")) %>%
  dplyr::select(PatientID, group,treatment_outcome, healing_time_days,
                Hypoxia_score, Hypoxia_groups, Neutrophil_score,
                GZMB, HIF1A, PRDM1)

write_tsv(as.data.frame(ExportHuman),
          "Submission/SupplementalTable4_HumanClinicalGeneMetadata.txt")
