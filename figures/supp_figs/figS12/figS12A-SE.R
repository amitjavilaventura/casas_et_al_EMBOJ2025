# R

# This script:
# - Measures splicing efficiency (SE) using split vs total reads mapping to splice sites
# - Does it with several datasets and then merges all
# - Plots the rsults

### ====== ###
### SET UP ###
### ====== ###

WORKDIR="~/projects/mouse_pic_var_review/";
setwd(WORKDIR)


# Load packages: 
library(dplyr)
library(data.table)
library(purrr)
library(magrittr)
library(openxlsx)

library(ggplot2)
library(patchwork)

source("helpers/read_files.R")

# Import info about canonical noct transcript
noct.info <- fread("data/public/ensembl_noct_canonical.csv")

# Import ICR info
icr.sample.info <- data.table::fread("data/private/icr_rnaseq/sample_info.csv") %>%
  magrittr::set_colnames(c("Samplename","altName","sample","diet","iap","species","strain","age","tissue", "x")) %>%
  dplyr::select(-x) %>%
  dplyr::mutate(iap = ifelse(is.na(iap), "IAP ?", iap),
                iap = ifelse(iap=="T", "IAP +", ifelse(iap=="F", "IAP -", NA)))


### =============== ###
### IMPORT EJ READS ###
### =============== ###

# Import reads
# Filter to get noct canonical transcript
# Divide split by total counts (split+notsplit)
# SE is set to 1 if split + not split reads reads are 0

### PRIVATE DATASET ###
noct_ss.priv <- 
list.files("output/04-rna_process/private/expression/bed_counts/splicesites/", "bed", full.names = T, recursive = T) %>% 
  purrr::set_names(~sub(".noct.*", "", basename(.x))) %>% 
  purrr::imap(~data.table::fread(.x, col.names = c("chr", "start", "end", "id", "intron_len", "strand", "counts")) %>%
                dplyr::mutate(cond = .y, sample = sub("\\..*", "", cond), split = sub(".*\\.", "", cond))) %>% 
  dplyr::bind_rows() %>%
  dplyr::transmute(id, counts, sample, split) %>%
  dplyr::filter(stringr::str_detect(id, paste(noct.info$TRANSCRIPT_ID, collapse = "|"))) %>%
  tidyr::separate("id", c("id", "side"), sep = "___") %>%
  dplyr::group_by(id, split,sample) %>%
  dplyr::summarise(counts = mean(counts)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(split,sample) %>%
  dplyr::arrange(id) %>%
  dplyr::mutate(ej = paste0("EJ", row_number())) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "split", values_from = "counts") %>%
  dplyr::transmute(sample=sub("RNA", "_rep", sample), cond = sub("_rep.*","",sample), dataset = "private",
                   tissue = ifelse(stringr::str_detect(sample, "spq"), "Spermatogonia", "Testis"),
                   age = "Adult", strain = sub("spq","",sub("testis","",cond)), ej, iap = ifelse(strain %in% c("BL6", "NOD"), "IAP +", "IAP -"),
                   se = (split)/(split+notsplit), se = ifelse(is.nan(se), 1, se), 
                   zeroreads = ifelse(is.nan(se), "Zero reads", ""))
 

### ICR OUTBRED MICE ###
noct_ss.icr <- list.files("output/04-rna_process/icr/expression/bed_counts/splicesites/", "bed", full.names = T, recursive = T) %>% 
  purrr::set_names(~sub(".noct.*", "", basename(.x))) %>% 
  purrr::imap(~data.table::fread(.x, col.names = c("chr", "start", "end", "id", "intron_len", "strand", "counts")) %>%
                dplyr::mutate(cond = .y, sample = sub("\\..*", "", cond), split = sub(".*\\.", "", cond))) %>% 
  dplyr::bind_rows() %>%
  dplyr::transmute(id, counts, sample, split) %>%
  dplyr::filter(stringr::str_detect(id, paste(noct.info$TRANSCRIPT_ID, collapse = "|"))) %>%
  tidyr::separate("id", c("id", "side"), sep = "___") %>%
  dplyr::group_by(id, split,sample) %>%
  dplyr::summarise(counts = mean(counts)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(split,sample) %>%
  dplyr::arrange(id) %>%
  dplyr::mutate(ej = paste0("EJ", row_number())) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "split", values_from = "counts") %>%
  dplyr::transmute(sample=sub("RNA", "_rep", sample), cond = "ICR", dataset = "ICR",
                   tissue = ifelse(stringr::str_detect(sample, "spq"), "Spermatogonia", "Testis"),
                   age = "Adult", strain = "Outbred", ej, se = split/(split+notsplit), 
                   se = ifelse(is.nan(se), 1, se),
                   zeroreads = ifelse(is.nan(se), "Zero reads", "")) %>%
  dplyr::inner_join(icr.sample.info %>% dplyr::select(sample, iap))

### Yu et al., 2019 ###
noct_ss.yu2019 <- list.files("output/04-rna_process/yu2019/expression/bed_counts/splicesites/", "bed", full.names = T, recursive = T) %>% 
  purrr::set_names(~sub(".noct.*", "", basename(.x))) %>% 
  purrr::imap(~data.table::fread(.x, col.names = c("chr", "start", "end", "id", "intron_len", "strand", "counts")) %>%
                dplyr::mutate(cond = .y, sample = sub("\\..*", "", cond), split = sub(".*\\.", "", cond))) %>% 
  dplyr::bind_rows() %>%
  dplyr::transmute(id, counts, sample, split) %>%
  dplyr::filter(stringr::str_detect(id, paste(noct.info$TRANSCRIPT_ID, collapse = "|"))) %>%
  tidyr::separate("id", c("id", "side"), sep = "___") %>%
  dplyr::group_by(id, split,sample) %>%
  dplyr::summarise(counts = mean(counts)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(split,sample) %>%
  dplyr::arrange(id) %>%
  dplyr::mutate(ej = paste0("EJ", row_number())) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "split", values_from = "counts") %>%
  dplyr::transmute(sample, cond=sample, dataset="yu2019", tissue="Testis",age = "Adult", strain=sample, ej, 
                   iap = ifelse(strain %in% c("C57BL6J", "C57BL6NJ", "AKRJ"), "IAP +", "IAP -"),
                   se = split/(split+notsplit), zeroreads = ifelse(is.nan(se), "Zero reads", ""), se = ifelse(is.nan(se), 1, se))

### ======== ###
### DO PLOTS ###
### ======== ###

### PRIVATE DATASET ###
noct_ss.priv.plot <- noct_ss.priv %>% dplyr::filter(tissue == "Testis") %>%
     ggplot(aes(strain, se, color=strain))+geom_point(size=3)+ggmitji::stat_line_boxplot(color="black",height=.002)+facet_wrap(~ej) +
     ggmitji::theme_custom(legend="right", axis.title.face = "plain", title.face = "italic") +
     ggmitji::remove_x_axis() +
     labs(y = "SE (split reads / total reads)", x = NULL, subtitle = "Private dataset > Testis samples",
          caption = "SE is set to 1 if the sum of split and not split reads is 0") +
      coord_cartesian(ylim = c(.5,1))

### ICR OUTBREED MICE ###
noct_ss.icr.plot <- noct_ss.icr %>%  na.omit() %>%
  ggplot(aes(iap, se, color=iap))+
  geom_point(size=3)+ggmitji::stat_line_boxplot(color="black",height=.001)+facet_wrap(~ej) +
  ggmitji::theme_custom(legend="right", axis.title.face = "plain", title.face = "italic") +
  ggpubr::stat_compare_means(method = "wilcox.test", paired = F, comparisons = list(c("IAP -", "IAP +")), tip.length = 0) +
  ggmitji::remove_x_axis() +
  labs(title="Splicing efficiency (SE) in exon junctions (EJ) of canonical Noct transcripts", 
       y = "SE (split reads / total reads)", x = NULL, subtitle = "ICR outbreed mice",
       caption = "SE is set to 1 if the sum of split and not split reads is 0") +
  coord_cartesian(ylim = c(.5,1.03))


### Yu et al., 2019 ###
noct_ss.yu2019.plot <- 
  noct_ss.yu2019 %>%  
  ggplot(aes(strain, se, color=strain))+
  geom_point(size=3)+ggmitji::stat_line_boxplot(color="black",height=.001)+facet_wrap(~ej) +
  ggmitji::theme_custom(legend="right", axis.title.face = "plain", title.face = "italic") +
  ggmitji::remove_x_axis() +
  labs(title="Splicing efficiency (SE) in exon junctions (EJ) of canonical Noct transcripts", 
       y = "SE (split reads / total reads)", x = NULL, subtitle = "Yu et al., 2019",
       caption = "SE is set to 1 if the sum of split and not split reads is 0") +
  coord_cartesian(ylim = c(.5,1))

### ========================== ###
### DO PLOTS WITH ALL DATASETS ###
### ========================== ###

# Join all datasets
# Add IAP information
noct_ss.all2 <- list(noct_ss.priv, 
                 noct_ss.icr,
                 noct_ss.yu2019) %>% 
  dplyr::bind_rows()

### ICR OUTBREED MICE ###
noct_ss.all.plot<-
  noct_ss.all2 %>% na.omit() %>%
  ggplot(aes(iap, se, color=iap))+
  geom_point(aes(shape=dataset),size=3)+ggmitji::stat_line_boxplot(color="black",height=.001)+facet_wrap(~ej) +
  ggpubr::stat_compare_means(method = "wilcox.test", paired = F, comparisons = list(c("IAP -", "IAP +")), tip.length = 0) +
  ggmitji::theme_custom(legend="right", axis.title.face = "plain", title.face = "italic") +
  ggmitji::remove_x_axis() +
  labs(title="Splicing efficiency (SE) in exon junctions (EJ) of canonical Noct transcripts", 
       y = "SE (split reads / total reads)", x = NULL, subtitle = "All datasets > Adult samples",
       caption = "SE is set to 1 if the sum of split and not split reads is 0") +
  coord_cartesian(ylim = c(0.5, 1.03))


#####################
### WRITE TO FILE ###
#####################

# Join all ier datasets
# Add IAP information
noct_ss.all <- list(noct_ss.priv, 
                noct_ss.icr,
                noct_ss.yu2019, 
                noct_ss.darbellay2020,
                noct_ss.isbel2015) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(iap = ifelse(!is.na(iap), iap, ifelse(stringr::str_detect(cond, "AKR|NOD|BL6"), "IAP +", "IAP -"))) 

pdf("figures/supp_figs/figS12/figS12A-SE.pdf", width = 9, height = 7)
noct_ss.priv.plot
noct_ss.icr.plot
noct_ss.yu2019.plot
noct_ss.all.plot
dev.off()

