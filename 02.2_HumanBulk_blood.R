# In this script I evaluate the hypoxia-related gene expression in the blood  of patients with CL
# The dataset that will be analyzed here is the Amorim et al. 2021.

# Libraries ----
library(GSVA)
library(GSEABase)

library(gplots)
library(limma)
library(edgeR)
#library(reshape2)
library(ggrepel)
library(ggpubr)
library(patchwork)

library(tidyverse)

# Import gene expression matrix and study design ----
# Gene expression matrix:
geneExpression_blood <- read_delim("Human_files/GSE162760_Amorim_GEO_filtered_normalized_batchAjusted.txt", 
                                                                 delim = "\t", escape_double = FALSE, 
                                                                 trim_ws = TRUE)
geneExpression_blood <- as.matrix(column_to_rownames(geneExpression_blood,"geneSymbol"))

# Study design:
studyDesign_blood <- read_csv("Human_files/studydesign_blood.csv")

GSVA.res_blood <- gsva(geneExpression_blood, genesets, mx.diff=FALSE, method="ssgsea")

# Harris Hypoxia gene signature ----
# Nat Rev Cancer 2002, PMID:11902584, https://www.gsea-msigdb.org/gsea/msigdb/cards/HARRIS_HYPOXIA

metadata_blood <- as.data.frame(GSVA.res_blood) %>%
  rownames_to_column("geneset") %>%
  pivot_longer(!geneset, names_to = "patient_ID", values_to = "Hypoxia_score") %>%
  full_join(studyDesign_blood) %>%
  mutate(geneset = factor(geneset)) %>%
  mutate(group_type = case_when(
    startsWith(patient_ID, 'CL') ~ "CL blood",
    startsWith(patient_ID, 'HS') ~ "non-Infected")) %>%
  mutate(group_type = factor(group_type, levels = c("non-Infected","CL blood"))) %>%
  filter(geneset %in% "Hypoxia_Harris")
  
# Violins:
metadata_blood %>%
  ggplot(., aes(x=group_type, y=Hypoxia_score)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("dark gray",Nandas[7],Nandas[9],"#C75836")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20, vjust = 0.5), 
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 20),
        legend.position="right", legend.text = element_text(size = 20), legend.title = element_blank()) +
  #stat_summary(fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(method = "t.test", label.y = 0.9,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.5, size = 5, color="black") +
  xlab("") + ylab(paste("Hypoxia score (Harris)")) + 
  NULL

write_tsv(metadata_blood %>%
            select(patient_ID,group_type,Hypoxia_score),
          "Submission/SupplementalTable5_HumanBloodHypoxia.txt")
