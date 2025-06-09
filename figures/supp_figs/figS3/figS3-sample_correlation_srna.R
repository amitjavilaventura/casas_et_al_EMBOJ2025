# R

# This script does the analysis to answer reviewer's comments #2 --> is intra-strain variation in piRNA sequences lower than inter-strain variation?

# This script:
# - Imports sRNA counts (normalized by DESeq2)
# - Does correlation plots so we see pairwise correlations for ALL the samples

# Set workding directory
WORKDIR="~/projects/github_repo_mouse_strains_paper/";
setwd(WORKDIR)

# Load packages
library(dplyr)
library(data.table)
library(reshape2)
library(corrplot)
library(ggplot2)
library(WVPlots)

library(here)

# Import normalize counts for sRNA-seq in Lietal clusters
norm <- data.table::fread(here::here("output/03-srna_process/fisher_DE_TEV/lietal_clusters/lietal_clusters.normcounts.tsv"))

# Do mean of normalized counts for each condition
norm.mean <- norm %>%
  reshape2::melt(variable.name = "sample", value.name = "expr") %>%
  dplyr::mutate(cond = sub("_rep.*","",sample)) %>%
  dplyr::group_by(Geneid, cond) %>%
  dplyr::summarise(expr=mean(expr)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "cond", values_from = "expr") 

# Do spearman correlation of the data
cor <- cor(norm[,2:ncol(norm)], method="spearman")

# Do pairwise scatterplot of samples from testis and spermatogonia separately
pplot.t <- 
  norm %>% 
  magrittr::set_colnames(colnames(.) %>% sub("testis_", "", .) %>% sub("_rep",".",.)) %>%
  dplyr::mutate(across(is.numeric, ~log10(.+0.1))) %>% 
  WVPlots::PairPlot(d=., meas_vars = colnames(.[,10:ncol(.)]), title = "")
pplot.s <- 
  norm %>% 
  magrittr::set_colnames(colnames(.) %>% sub("spgonia_", "", .) %>% sub("_rep",".",.)) %>%
  dplyr::mutate(across(is.numeric, ~log10(.+0.1))) %>% 
  WVPlots::PairPlot(d=., meas_vars = colnames(.[,2:9]), title = "")

# Extract data from pair plots, and remove the data from the lower part of the pair plot
d <- pplot.t$data %>% 
  dplyr::filter(xv!=yv) %>%
  dplyr::filter(!pair_key %in% 
                  c(paste("129.1", c("129.2","129.3","BL6.1","BL6.2","C3H.1","C3H.2","C3H.3","NOD.1","NOD.2")),
                    paste("129.2", c("129.3","BL6.1","BL6.2","C3H.1","C3H.2","C3H.3","NOD.1","NOD.2")),
                    paste("129.3", c("BL6.1","BL6.2","C3H.1","C3H.2","C3H.3","NOD.1","NOD.2")),
                    paste("BL6.1", c("BL6.2","C3H.1","C3H.2","C3H.3","NOD.1","NOD.2")),
                    paste("BL6.2", c("C3H.1","C3H.2","C3H.3","NOD.1","NOD.2")),
                    paste("C3H.1", c("C3H.2","C3H.3","NOD.1","NOD.2")),
                    paste("C3H.2", c("C3H.3","NOD.1","NOD.2")),
                    paste("C3H.3", c("NOD.1","NOD.2")),
                    paste("NOD.1", c("NOD.2")))) %>%
  dplyr::mutate() 
d.s <- pplot.s$data %>%
  dplyr::filter(xv!=yv,
                !pair_key %in% c(paste("BL6.1", c("BL6.2","BL6.3","BL6.4","CAST.1","CAST.2","CAST.3","CAST.4")),
                                 paste("BL6.2", c("BL6.3","BL6.4","CAST.1","CAST.2","CAST.3","CAST.4")),
                                 paste("BL6.3", c("BL6.4","CAST.1","CAST.2","CAST.3","CAST.4")),
                                 paste("BL6.4", c("CAST.1","CAST.2","CAST.3","CAST.4")),
                                 paste("CAST.1", c("CAST.2","CAST.3","CAST.4")),
                                 paste("CAST.2", c("CAST.3","CAST.4")),
                                 paste("CAST.3", c("CAST.4"))))

# Redo pair plots with ggplot
pplot.t <- 
d %>%
  ggplot(aes(x,y)) +
  geom_abline(intercept = 0, slope = 1, linewidth=.2, linetype=2) +
  geom_point(color="gray30", alpha=.5) +
  facet_grid(cols=vars(xv), rows=vars(yv)) +
  ggpubr::stat_cor(mapping=aes(label = ..r.label..), method = "spearman", cor.coef.name = "rho", size = 3, label.y=6, label.x = 1, color="red") +
  scale_y_continuous(position = "right") + scale_x_continuous(position="top") +
  coord_cartesian(ylim = c(0,7), xlim = c(0,7)) +
  ggmitji::theme_custom(legend = "none", margin = T, border = F, 
                        title.face = "italic", subtitle.face = "italic", axis.title.face = "plain", axis.title.size = 7) +
  labs(x=expression(Log[10]~"norm. counts +0.1"), y=expression(Log[10]~"norm. counts +0.1"),
       title = "Pairwise correlation of sRNA expression across testis samples", 
       subtitle = "piRNA clusters annotated by Li et al. (2013)") 

pplot.s <- 
  d.s %>%
  ggplot(aes(x,y)) +
  geom_abline(intercept = 0, slope = 1, linewidth=.2, linetype=2) +
  geom_point(color="gray30", alpha=.5) +
  facet_grid(cols=vars(xv), rows=vars(yv)) +
  ggpubr::stat_cor(mapping=aes(label = ..r.label..), method = "spearman", cor.coef.name = "rho", size = 3, label.y=6, label.x = 1, color="red") +
  scale_y_continuous(position = "right") + scale_x_continuous(position="top") +
  coord_cartesian(ylim = c(0,7), xlim = c(0,7)) +
  ggmitji::theme_custom(legend = "none", margin = T, border = F, 
                        title.face = "italic", subtitle.face = "italic", axis.title.face = "plain", axis.title.size = 7) +
  labs(x=expression(Log[10]~"norm. counts +0.1"), y=expression(Log[10]~"norm. counts +0.1"),
       title = "Pairwise correlation of sRNA expression across spermatogonia samples", 
       subtitle = "piRNA clusters annotated by Li et al. (2013)") 


# Fill strips from facets
pplot.t <- ggmitji::fill_strips(pplot.t, side=c("top","right"), 
                                colors=c(rep("lightyellow2",2), rep("gray",2), rep("orange",3), rep("darkred",2),
                                         rep("lightyellow2",3), rep("gray",2), rep("orange",3), rep("darkred",1)))
pplot.s <- ggmitji::fill_strips(pplot.s, side=c("top","right"), 
                                colors=c(rep("gray",3), rep("khaki",4),
                                         rep("gray",4), rep("khaki",3)))

# Save plot to file
pdf(here::here("figures/individual_var/figS3A-sample_correlation_srna.pdf"), width = 7, height = 7)
pplot.t
dev.off()

pdf(here::here("figures/individual_var/figS3B-sample_correlation_srna.pdf"), width = 7, height = 7)
pplot.t
dev.off()



