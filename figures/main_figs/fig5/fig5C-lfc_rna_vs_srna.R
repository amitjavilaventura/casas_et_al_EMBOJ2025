# R

# This script does analysis to answer reviewer's comments #3 --> are changes in piRNA expression due to changes in RNA expression?

# This script:
# - Takes the log2fold change of all pairwise contrasts, both for RNA-seq and sRNA-seq.
# - Compares log2FC from RNA-seq and sRNA-seq for each contrasts.
# - Plots the result

WORKDIR="~/projects/mouse_pic_var_review/";
setwd(WORKDIR)

# Load packages
library(dplyr)
library(here)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggmitji)
source(("tools/read_files.R"))

### COMPARISON OF LOG2FOLDCHANGES ---------------

# Import small RNA-seq DEA results for Li et al clusters
srna.deg <- read_excel_sheets(("output/03-srna_process/fisher_DE_TEV/lietal_clusters/lietal_clusters.DEGs.xlsx"),"all",T) %>%
  dplyr::bind_rows() %>%
  dplyr::select(piCID,"lfc_srna"=log2FoldChange,Contrast)

# Import RNA-seq DEA results for Li et al clusters
# Mutate contrast 129vsC3H to reverse it, since in sRNA seq was C3H vs 129
rna.deg <- read_excel_sheets(("output/04-rna_process/private/expression/degs/lietal_clusters.degs.xlsx"),"all",T) %>%
  dplyr::bind_rows() %>%
  dplyr::select(piCID,"lfc_rna"=log2FoldChange,Contrast) %>%
  dplyr::mutate(lfc_rna = ifelse(Contrast == "129vsC3H", -lfc_rna, lfc_rna),
                Contrast = ifelse(Contrast == "129vsC3H", "C3Hvs129",Contrast))

# Join log2FoldChange of RNA-seq and smallRNA seq in one data frame
deg <- dplyr::inner_join(srna.deg,rna.deg)%>% dplyr::mutate(piCID = sub("__.*","",piCID))

# Define genes/clusters that have to be highlighted
genes_to_show = c("pi-Ccrn4l", "pi-Zbtb37", "pi-Phf20", "pi-Mrs2", "pi-Ddx19b", "pi-Pou6f1", "17-qC-59")

# Format DEG dataframe
deg <- deg %>%
  dplyr::mutate(highlight = ifelse(piCID %in% genes_to_show, T, F)) 

# Do plot of lfc comparison
plot.lfc.rna_vs_srna.lietal <-
ggplot(deg, aes(lfc_rna, lfc_srna, color=highlight)) +
  # annotate("rect", xmin=1, xmax=12, ymin=1, ymax=7, alpha=0.3, fill="cornflowerblue") +
  # annotate("rect", xmin=1, xmax=12, ymin=-7, ymax=-1, alpha=0.3, fill="khaki1") +
  # annotate("rect", xmin=-12, xmax=-1, ymin=1, ymax=7, alpha=0.3, fill="pink") +
  # annotate("rect", xmin=-12, xmax=-1, ymin=-7, ymax=-1, alpha=0.3, fill="lightgreen") +
  # geom_vline(xintercept = 0, linewidth=.2, color="black") +
  geom_vline(xintercept = 1, linewidth=.3, color="black", linetype=2) +
  geom_vline(xintercept = -1, linewidth=.3, color="black", linetype=2) +
  # geom_hline(yintercept = 0, linewidth=.2, color="black") +
  geom_hline(yintercept = -1, linewidth=.3, color="black", linetype=2) +
  geom_hline(yintercept = 1, linewidth=.3, color="black", linetype=2) +
  geom_point(data=deg%>%dplyr::filter(!highlight), color="gray30", alpha=1) +
  geom_point(data=deg%>%dplyr::filter(highlight), color="darkred", alpha=1) +
  ggrepel::geom_text_repel(aes(label=piCID), data = deg%>%dplyr::filter(highlight),color="darkred",size=2) +
  facet_wrap(~Contrast) +
  ggmitji::theme_custom(title.face = "italic", subtitle.face = "italic", axis.title.face = "plain") +
  ggmitji::add_grid(.1,0,.1,0) +
  labs(y=expression(Log[2]~FoldChange~"(sRNA-seq)"), x=expression(Log[2]~FoldChange~"(RNA-seq)"),
       title="Fold-change comparison - RNA vs small RNA",
       subtitle = "piRNA clusters annotated by Li et al. (2013)",
       caption = "Reads were aligned to the corresponding assembly for each strain.
       Only clusters present in all strains after ENSEMBL Compara are used.") +
  coord_cartesian(ylim = c(-6,6), xlim = c(-10,10))

### SAVE PLOTS TO FILE --------------------------

pdf(("figures/rna_vs_srna/fig5C-lfc_rna_vs_srna"), height = 5, width = 5)
plot.lfc.rna_vs_srna.lietal
dev.off()
