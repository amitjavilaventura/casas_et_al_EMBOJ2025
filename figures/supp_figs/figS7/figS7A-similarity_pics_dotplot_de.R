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
source("helpers/read_files.R")

# Import lietal clusters overlapping protein coding genes
# Import lncRNA genes overlapping protein coding genes
# Import lncRNA genes overlapping piRNA clusters
protcod.lietal <- data.table::fread(here::here("output/05-similarity_analysis/lietal_orthologs/lietal.proteincoding.txt"), header = F, col.names = "pic") %>%
  tidyr::separate("pic", c("piCid", "piCnum"), sep = "__", remove = F)
protcod.lncrna <- data.table::fread(here::here("output/05-similarity_analysis/lncRNA_orthologs/lncRNA.proteincoding.txt"), header = F, col.names = "gene_id")
protcod.protrac <- data.table::fread(here::here("output/05-similarity_analysis/protrac_orthologs/protrac.proteincoding.txt"), header = F, col.names = "pic")
lncrna.pics <- data.table::fread(here::here("output/05-similarity_analysis/lncRNA_orthologs/lncRNA.piRNAclusters.txt"), header = F, col.names = "gene_id")
prot.pics <- data.table::fread(here::here("output/05-similarity_analysis/coding_orthologs/proteincoding.piCs.txt"), header = F, col.names = "gene_id")

# Import lietal degs
lietal.degs <- read_excel_sheets("output/03-srna_process/fisher_DE_TEV/lietal_clusters/lietal_clusters.DEGs.xlsx", "all", T) %>%
  purrr::map(~dplyr::mutate(.x, piCID=sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", piCID)))
lietal.degs.df <- lietal.degs %>% dplyr::bind_rows()

### IMPORT DATA ---------------------------------

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
similarity.data <- list.files(here::here("output/05-similarity_analysis/"), "similarity", recursive = T, full.names = T) %>%
  purrr::discard(~stringr::str_detect(.x, "old")) %>%
  purrr::set_names(~sub(".tsv","",sub("_orthologs.similarity", "", basename(.x)))) %>%
  
  # Read data and add the group (e.g., protrac, coding, lncrna...) and N (with or without N) as a new column
  # Separate the align column into group and N
  purrr::imap(~data.table::fread(.x) %>% dplyr::mutate(align=.y)) %>%
  purrr::map(~tidyr::separate(.x, "align", c("group", "N"), "\\.")) %>%
  
  # Calculate % of identity, gaps and mismatches.
  purrr::map(~dplyr::mutate(.x, pct_id = identical/len*100, pct_mismatch = mismatch/len*100, pct_gap = gap/len*100)) %>%

  # Set class of pic/gene
  purrr::map(~dplyr::mutate(.x, class = ifelse(group == "lietal", class, group))) %>%
  
  # Set whether the genes/piCs overlap with protein-coding genes or not.
  purrr::map(~dplyr::mutate(.x, proteincoding = ifelse(group=="protcoding", "Protein-coding gene", ifelse(geneid %in% c(protcod.lietal$pic, protcod.protrac$pic, protcod.lncrna$gene_id), "Overlap with coding gene", "No overlap with coding gene")))) %>%
  
  # Set whether genes overlap with piRNA cluster or not
  purrr::map(~dplyr::mutate(.x, pirnacluster = ifelse(group %in% c("lietal","protrac"), "piRNA cluster", ifelse(geneid %in% c(prot.pics$gene_id, lncrna.pics$gene_id), "Overlap with piC", "No overlap with piC")))) %>%
  
  # Filter to get only those genes/pics for which I have all the alignments (in case some alignments didn't work correctly)
  purrr::map(~dplyr::group_by(.x, geneid) %>% dplyr::mutate(i=length(geneid)) %>% dplyr::ungroup() %>% dplyr::filter(i==15))


### Join all lists in one dataframe
similarity.df <- dplyr::bind_rows(similarity.data) %>%
  # Factor and relevel class, group, etc
  dplyr::mutate(class = factor(class, c("Pre-pachytene", "Hybrid", "Pachytene", "protrac", "lncRNA", "protcoding"), c("Pre-pachytene", "Hybrid", "Pachytene", "proTRAC", "lncRNA", "Protein-coding")),
                group = factor(group, c("lietal", "protrac", "lncRNA", "protcoding"), c("Li et al., 2013", "proTRAC", "lncRNA", "Protein-coding")),
                proteincoding = factor(proteincoding, c("Overlap with coding gene", "No overlap with coding gene", "Protein-coding gene")),
                pirnacluster = factor(pirnacluster, c( "piRNA cluster" ,"Overlap with piC", "No overlap with piC")),
                strain1 = factor(strain1, c("BL6","NOD","C3H","129","CAST")),
                strain2 = factor(strain2, c("BL6","NOD","C3H","129","CAST")),
                comparison = paste(strain1, "vs", strain2),
                comparison = factor(comparison, c("BL6 vs BL6", "BL6 vs NOD", "BL6 vs C3H", "BL6 vs 129", "BL6 vs CAST",
                                                  "NOD vs NOD", "NOD vs C3H", "NOD vs 129", "NOD vs CAST", 
                                                  "C3H vs C3H", "C3H vs 129", "C3H vs CAST",
                                                  "129 vs 129", "129 vs CAST", 
                                                  "CAST vs CAST")))
                                              

### Do mean/median of identity percentage
similarity.summary.group <- similarity.df %>%
  dplyr::group_by(group, N, proteincoding, pirnacluster, comparison, strain1, strain2) %>%
  dplyr::summarise(mean_id = mean(pct_id), median_id = median(pct_id), sd_id = sd(pct_id),
                   mean_gap = mean(pct_gap), median_gap = median(pct_gap), sd_gap = sd(pct_gap),
                   mean_mismatch = mean(pct_mismatch), median_mismatch = median(pct_mismatch), sd_mismatch = sd(pct_mismatch),
                   mean_len = mean(len), median_len = median(len), sd_len = sd(len)) %>%
  dplyr::ungroup()


similarity.summary.group.all <- similarity.df %>%
  dplyr::group_by(group, N, comparison,pirnacluster, strain1, strain2) %>%
  dplyr::summarise(mean_id = mean(pct_id), median_id = median(pct_id), sd_id = sd(pct_id),
                   mean_gap = mean(pct_gap), median_gap = median(pct_gap), sd_gap = sd(pct_gap),
                   mean_mismatch = mean(pct_mismatch), median_mismatch = median(pct_mismatch), sd_mismatch = sd(pct_mismatch),
                   mean_len = mean(len), median_len = median(len), sd_len = sd(len)) %>%
  dplyr::ungroup()

similarity.summary.class <- similarity.df %>%
  dplyr::group_by(class, N, proteincoding, comparison, strain1, strain2) %>%
  dplyr::summarise(mean_id = mean(pct_id), median_id = median(pct_id), sd_id = sd(pct_id),
                   mean_gap = mean(pct_gap), median_gap = median(pct_gap), sd_gap = sd(pct_gap),
                   mean_mismatch = mean(pct_mismatch), median_mismatch = median(pct_mismatch), sd_mismatch = sd(pct_mismatch),
                   mean_len = mean(len), median_len = median(len), sd_len = sd(len)) %>%
  dplyr::ungroup()


### DO PLOTS ------------------------------------

#### DOTPLOTS and DEA -----
data.similarity.degs <- similarity.df %>%
  dplyr::filter(strain1!=strain2, N=="withoutN",strain1=="BL6",pirnacluster=="piRNA cluster") %>%
  dplyr::mutate(Contrast=paste0(strain2,"vs",strain1)) %>%
  dplyr::inner_join(lietal.degs.df, by = c("piCid2"="piCID","Contrast")) %>% 
  dplyr::mutate(DEG=ifelse(DEG!="NS", "DE", "NotDE")) %>%
  dplyr::select(N, comparison, group, proteincoding, piCid, class, pct_id, DEG, padj, log2FoldChange) %>%
  unique()

dotplot.similarity.degs <- 
data.similarity.degs %>%
  dplyr::filter(class=="Pachytene",
                comparison != "BL6 vs CAST") %>%
  ggplot(aes(DEG, pct_id)) +
  geom_hline(yintercept = 100, linetype=2, linewidth=.2, color="black") +
  # geom_violin(aes(fill=DEG), color="black") +
  geom_boxplot(aes(color=DEG), width = .5) +
  geom_jitter(aes(color=DEG), width = .1) +
  facet_grid(~comparison) +
  ggmitji::stat_info_boxplot("n","N=",-30,3,"black",.5) +
  ggmitji::stat_info_boxplot("mean","µ=",-45,3,"darkblue",.5) +
  ggmitji::stat_info_boxplot("median","M=",-60,3,"darkred",.5) +
  ggpubr::stat_compare_means(method = "wilcox.test", paired = F, comparisons=list(c(1,2)), label = "p.format", label.y = 105, tip.length = 0, size=3) +
  scale_fill_manual(values = c("red","gray"), aesthetics = c("fill", "colour")) +
  labs(title = "Similarity of pachytene piRNA clusters between strains",
       subtitle = "piRNA clusters from Li et al. (2013) separated by differential expression",
       caption = "Names of clusters with less than 75% of identity are shown.
       Only clusters present in all the strains after ENSEMBL Compara are used.",
       y="% of identity", x=NULL) +
  ggmitji::theme_custom(legend = "right",title.face = "italic", axis.title.face = "plain")+
  ggmitji::remove_x_axis() +
  scale_y_continuous(breaks = c(0,25,50,75,100), limits = c(-60, 130))


### WRITE PLOTS TO FILE -------------------------

pdf("figures/similarity/fig7A-similarity_pics_dotplot_de.pdf", height = 11, width = 7)
dotplot.similarity.degs
dev.off()


