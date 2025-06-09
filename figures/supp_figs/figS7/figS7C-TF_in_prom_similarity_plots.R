# R

# This script:
# - Takes the similarity values from pairwise global alignments of piRNA clusters from Li et al. (2013) and lncRNA genes
# - ALL lncRNA genes with orthologs in the 5 strains were used
# - Only piRNA clusters for which I could retrieve orthologs (with ENSEMBL Compara) for all strains are used (192 from 214).
# - For lncRNA genes and protein-coding genes, only 100 random genes are taken!

### SET UP --------------------------------------

# Set workding directory
WORKDIR="~/projects/mouse_pic_var_review/";
setwd(WORKDIR)
# Load packages
library(dplyr)
library(data.table)
library(here)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggmitji)
library(ggpubr)
library(patchwork)
library(bedtoolsr)

# Define annotation
# Set a second name for the annotation, since I put the wrong name in previous analisis for the lietal clusters
annot = "lietal_clusters" # "lietal_clusters", "protrac_merged"
if(annot=="lietal_clusters") { annoname2 = "lietal_orthologs" } else { annoname2 = annot} 

# Import amyb overlapping piC promoters
amyb_in_de_prom.info <- list.files(here::here("output/06-promoters_analysis/lietal_orthologs/TFpeaks/TF_in_DEpiC_prom/"),"bed", full.names = T) %>%
  purrr::set_names(basename(.) %>% sub("amyb_in_promoters.","",.) %>% sub("tcfl5_in_promoters.","",.) %>% sub("_in_bl6.*","",.)) %>%
  purrr::imap(~data.table::fread(.x, header = F, col.names = c("seqnames","start","end","peakid","length","strand")) %>%
                dplyr::mutate(contrast = sub("\\..*","",.y), depic = sub(".*\\.","",.y))) %>%
  dplyr::bind_rows() %>%
  dplyr::transmute(peakid, contrast, depic, tf = peakid %>% sub("_.*","",.))

 


### IMPORT PROMOTER DATA ---------------------------------

# I have calculated similarity from the pairwise alignments with both Ns and without Ns (this is done in scripts in scripts/5similarity_analysis)
# Basically, I read the alignments in FASTA format (output from EMBOSS Stretcher)
# For alignments with Ns:
# - I count the number of positions with identical bases (including N-N), mismatches and gaps in the alignment.
# - Then I calculate the percentage dividing by the full length.
# For the alignments without Ns:
# - I remove the positions where one of the sequences has an N (regardless the other base)
# - I count the number of pisitions with identical bases, mismatches and gaps.
# - I calculate the percentage dividing by the new length.

# Here I import the data, and then I calculate the percentage of matches, mismatches and gaps. 
# Also I add whether the alignments are with or without Ns, as well as the piRNA cluster classes

### List desried files and give them a name
similarity.data <- list.files(here::here("output/06-promoters_analysis/",annoname2,"/TFpeaks/stretcher/similarity/"), "similarity", recursive = T, full.names = T) %>%
  purrr::set_names(~sub(".tsv","",sub(".*similarity.", "", basename(.x)))) %>%
  
  # Read data and add the group (e.g., protrac, coding, lncrna...) and N (with or without N) as a new column
  # Separate the align column into group and N
  purrr::imap(~data.table::fread(.x) %>% dplyr::mutate(n=.y, contrast = paste0(strain2, "vs", strain1))) %>%
  
  # Calculate % of identity, gaps and mismatches.
  purrr::map(~dplyr::mutate(.x, pct_id = identical/len*100, pct_mismatch = mismatch/len*100, pct_gap = gap/len*100)) %>%
  
  # Bind rows
  dplyr::bind_rows()

# Add wether they overlap de promoter or not
similarity.df <- similarity.data %>% dplyr::inner_join(amyb_in_de_prom.info)

# Get classes of lietal clusters that overlap with TF peaks in promoters
# And add it to similarity
if(annot=="lietal_clusters"){
  # Import classes
  cluster_class <- data.table::fread("data/public/lietal2013/lietal2013-piRNA_clusters.classes.csv")
  
  # Import TF peaks in promoters
  amyb_in_prom  <- data.table::fread("output/06-promoters_analysis/lietal_orthologs/TFpeaks/TF_in_prom/lietal_clusters-amyb_in_promoters.bed") %>% 
    dplyr::transmute(peakid=V4, piCid=V7) %>% dplyr::left_join(cluster_class)
  tcfl5_in_prom <- data.table::fread("output/06-promoters_analysis/lietal_orthologs/TFpeaks/TF_in_prom/lietal_clusters-tcfl5_in_promoters.bed") %>%
    dplyr::transmute(peakid=V4, piCid=V7) %>% dplyr::left_join(cluster_class)
  tf_in_prom <- dplyr::bind_rows(amyb_in_prom, tcfl5_in_prom) %>% dplyr::distinct(peakid,class) 
  
  # Add tf class in similarity data
  similarity.df <- similarity.df %>% dplyr::inner_join(tf_in_prom) %>% dplyr::mutate(class=factor(class,c("Pre-pachytene","Hybrid","Pachytene")))
}


## This is just to check that all peaks are present after join
# similarity.df %>% dplyr::distinct(peakid) %>% nrow()                                 
# similarity.data %>% dplyr::distinct(peakid) %>% nrow()                                 
# amyb_in_de_prom.info %>% dplyr::distinct(peakid) %>% nrow()                                 

### Do mean/median of identity percentage
similarity.summary.de <- similarity.df %>%
  dplyr::group_by(n, contrast, strain1, strain2, depic, tf) %>%
  dplyr::summarise(mean_id = mean(pct_id), median_id = median(pct_id), sd_id = sd(pct_id),
                   mean_gap = mean(pct_gap), median_gap = median(pct_gap), sd_gap = sd(pct_gap),
                   mean_mismatch = mean(pct_mismatch), median_mismatch = median(pct_mismatch), sd_mismatch = sd(pct_mismatch),
                   mean_len = mean(len), median_len = median(len), sd_len = sd(len)) %>%
  dplyr::ungroup()

### Do mean/median of identity percentage with class
similarity.summary.de.class <- similarity.df %>%
  dplyr::group_by(n, contrast, strain1, strain2, depic, tf, class) %>%
  dplyr::summarise(mean_id = mean(pct_id), median_id = median(pct_id), sd_id = sd(pct_id),
                   mean_gap = mean(pct_gap), median_gap = median(pct_gap), sd_gap = sd(pct_gap),
                   mean_mismatch = mean(pct_mismatch), median_mismatch = median(pct_mismatch), sd_mismatch = sd(pct_mismatch),
                   mean_len = mean(len), median_len = median(len), sd_len = sd(len)) %>%
  dplyr::ungroup()


### TAKE RANODM NOTDE -------------
dat2 <- list()
for ( comparison in unique(similarity.df$contrast) ){
  
  dat <- similarity.df %>% dplyr::filter(contrast == comparison) %>% dplyr::filter(n == "withoutN")
  
  de <- dat %>% dplyr::filter(depic == "de")
  nd <- dat %>% dplyr::filter(depic == "notde")
  
  num.de <- nrow(de)
  set.seed(666)
  nd <- dplyr::slice_sample(nd, n=num.de)

  dat2[[comparison]] <- dplyr::bind_rows(de,nd) 
    
}

dat2 <- dplyr::bind_rows(dat2)

dat2.summary.de <- dat2 %>%
  dplyr::select(n, contrast, strain1, strain2, depic,tf, matches("pct|len")) %>% 
  unique() %>%
  dplyr::group_by(n, contrast, strain1, strain2, depic,tf) %>%
  dplyr::summarise(mean_id = mean(pct_id), median_id = median(pct_id), sd_id = sd(pct_id),
                   mean_gap = mean(pct_gap), median_gap = median(pct_gap), sd_gap = sd(pct_gap),
                   mean_mismatch = mean(pct_mismatch), median_mismatch = median(pct_mismatch), sd_mismatch = sd(pct_mismatch),
                   mean_len = mean(len), median_len = median(len), sd_len = sd(len)) %>%
  dplyr::ungroup()

### DO PLOTS ------------------------------------

data3 <- similarity.df %>% dplyr::filter(n=="withoutN") %>%
  dplyr::mutate(depic = sub("notde", "No", depic) %>% sub("de", "Yes",.),
                depic = factor(depic, c("Yes", "No"), c("DE", "NotDE")))


dotplot_tf_in_prom_similarity.class <- 
data3 %>% dplyr::filter(strain2!="CAST", n=="withoutN",class=="Pachytene") %>%
  ggplot(aes(depic, pct_id, color=depic)) +
  ggmitji::stat_summary_boxplot() +
  geom_jitter(width=.1) +
  geom_hline(yintercept = 100, linetype=2, linewidth=.2, color="gray30") +
  ggmitji::stat_info_boxplot("n","N=",-1,3,"black") +
  ggmitji::stat_info_boxplot("mean", "µ=", -8, 3, "darkblue") +
  ggmitji::stat_info_boxplot("median", "M=", -15, 3, "darkred") +
  ggpubr::stat_compare_means(method = "wilcox.test", paired = F, comparisons = list(c(1,2)), label = "p.format", tip.length = 0, label.y = 103) +
  facet_grid(rows=vars(tf), cols=vars(contrast)) +
  ggmitji::theme_custom(title.face = "plain", axis.title.face = "plain", subtitle.size = 9) +
  scale_fill_manual(values = c("green","red"), aesthetics = c("fill","colour")) +
  coord_cartesian(ylim = c(-20, 120)) +
  labs(x="", y = "% of identity", 
       title = "Similarity in pairwise global alignments of A-MYB/TCFL5 peaks and their orthologs", 
       subtitle = "Overlapping promoters of DE and not DE pachytene piRNA clusters in BL6",
       caption = "Promoters are TSS +/- 2500bp; piRNA cluters from Li et al. (2013); ChIP-seq peaks from Yu et al. (2023).")

### WRITE PLOTS TO FILE

pdf("figures/supp_figs/figS7/figS7C-tf_in_prom_similarity_plots.dotplot.pdf", height = 11, width = 11)
dotplot_tf_in_prom_similarity.class
dev.off()
