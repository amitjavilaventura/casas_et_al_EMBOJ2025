# R

# This script:
# - Does a barplot showing the percentage of predicted piCs that overlap piRNA clusters from Li et al. (2013)

WORKDIR="~/projects/mouse_pic_var_review/";
setwd(WORKDIR)

# Load packages
library(dplyr)
library(data.table)
library(bedtoolsr)
library(ggplot2)
library(here)

# Import list of piRNA clusters present in ALL strains
lietal.allstrains <- data.table::fread(here::here("output/02-pirna_clusters_orthologs/lietal_clusters/lietal_clusters_present_in_all_strains.tsv"))
proall.allstrains <- data.table::fread(here::here("output/02-pirna_clusters_orthologs/protrac_merged/protrac_merged_present_in_all_strains.tsv"))

# Import piRNA clusters coordinates in mm10 format
lietal.mm10 <- data.table::fread(here::here("output/01-pirna_clusters/lietal_clusters/lietal_clusters.mm10.bed"), header = F)
proall.mm10 <- data.table::fread(here::here("output/01-pirna_clusters/protrac_clusters/protrac_merged/protrac_merged.3rmsk_filt.mm10.bed"), header = F)

# Overlap lietal and protrac
protrac.vs.lietal <- bedtoolsr::bt.intersect(proall.mm10[which(proall.mm10$V4 %in% proall.allstrains$piCid),], lietal.mm10[which(lietal.mm10$V4 %in% lietal.allstrains$piCid),])

# Format protrac all strains to see if they overlap with lietal clusters or not
proall.overlap.lietal <- proall.allstrains %>% dplyr::mutate(known = ifelse(piCid %in% protrac.vs.lietal$V4, "Overlap", "Don't overlap"))

# Do venn diagram of overlap
plotmics::ggVennBed(proall.mm10[which(proall.mm10$V4 %in% proall.allstrains$piCid),], 
                    lietal.mm10[which(lietal.mm10$V4 %in% lietal.allstrains$piCid),], 
                    setnames = c("proTRAC clusters", "Lietal clusters"), stranded = F)

# Do barplot
barplot <- proall.overlap.lietal %>%
  dplyr::count(known) %>%
  ggplot(aes("", n, fill=known)) +
  geom_col(position="fill", color="black") +
  geom_text(aes(label=n), position=position_fill(vjust = .5), angle = 90) +
  theme_classic() +
  ggmitji::remove_x_axis() +
  labs(y="% of predicted clusters overlapping known clusters")


# Write to pdf
pdf("figures/main_figs/fig2/fig2B-pct_predicted_known.barplot.pdf", width = 3, height = 5)
barplot
dev.off()
