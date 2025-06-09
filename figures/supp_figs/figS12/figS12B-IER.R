# R

# This script:
# - Does plots about the inverse intron expression ratio (IER)
#   - IER is calculated as the per-base coverage (counts divided by length) of the intron diveided by the mean per-base coverage of the surrounding exons.
#   - Counts in introns and exons are counted after removing reads mapping to repeats and gaps +/- 150bp (repeats and gaps +/- 150 bp are also removed from annotation)
# - Does it with several datasets, and then merges them and plots the data separating the strains with and without IAP.


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


### ======================== ###
### IMPORT INTRON/EXON READS ###
### ======================== ###

# List fcounts files of introns and exons and set names to them
# Import files and select desired columns, adding a column for sample name and intron/exon
# Filter to get noct canonical transcript
# Bind rows and add intron/exon number (exon1, intron1, exon2, intron 2) 
# Add condition name (e.g., spqBL6RNA2 --> spqBL6), tissue and strain names


### PRIVATE DATASET ###
ei.private <- list.files("output/04-rna_process/private/expression/fcounts_filt_noct/", "counts$", full.names = T, recursive = T) %>%
  purrr::keep(~stringr::str_detect(.x, "exon_id|intron_id")) %>%
  purrr::set_names(~sub(".s0.*", "", basename(.)) %>% sub("-.*","",.)) %>%
  purrr::map(~data.table::fread(.x, col.names = c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Counts"))) %>% 
  purrr::imap(~dplyr::transmute(.x, cond = .y, Geneid, Length, Counts)) %>%
  purrr::map(~dplyr::filter(.x, stringr::str_detect(Geneid, paste(noct.info$TRANSCRIPT_ID, collapse = "|")))) %>%
  dplyr::bind_rows() %>%
  tidyr::separate("cond", c("sample", "region"), "\\.") %>%
  dplyr::arrange(Geneid) %>%
  dplyr::group_by(sample, region) %>% dplyr::mutate(regionNum = paste0(region, row_number())) %>% dplyr::ungroup() %>%
  dplyr::transmute(tissue = ifelse(stringr::str_detect(sample, "spq"), "Spermatogonia", "Testis"),
                   cond = sub("RNA.","",sample), 
                   strain = sub("spq","",sub("testis","",cond)),
                   age = "adult",
                   sample, regionNum, Length, Counts)

### ICR OUTBRED MICE ###
ei.icr <- list.files("output/04-rna_process/icr/expression/fcounts_filt_noct/", "counts$", full.names = T, recursive = T) %>%
  purrr::keep(~stringr::str_detect(.x, "exon_id|intron_id")) %>%
  purrr::set_names(~sub(".s0.*", "", basename(.)) %>% sub("-.*","",.)) %>%
  purrr::map(~data.table::fread(.x, col.names = c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Counts"))) %>% 
  purrr::imap(~dplyr::transmute(.x, cond = .y, Geneid, Length, Counts)) %>%
  purrr::map(~dplyr::filter(.x, stringr::str_detect(Geneid, paste(noct.info$TRANSCRIPT_ID, collapse = "|")))) %>%
  dplyr::bind_rows() %>%
  tidyr::separate("cond", c("sample", "region"), "\\.") %>%
  dplyr::arrange(Geneid) %>%
  dplyr::group_by(sample, region) %>% dplyr::mutate(regionNum = paste0(region, row_number())) %>% dplyr::ungroup() %>%
  dplyr::transmute(tissue = ifelse(stringr::str_detect(sample, "spq"), "Spermatogonia", "Testis"),
                   cond = sub("RNA.","",sample), 
                   strain = sub("spq","",sub("testis","",cond)),
                   age = "adult",
                   sample, regionNum, Length, Counts)


### Yu et al., 2019 ###
ei.yu2019 <- list.files("output/04-rna_process/yu2019/expression/fcounts_filt_noct/", "counts$", full.names = T, recursive = T) %>%
  purrr::keep(~stringr::str_detect(.x, "exon_id|intron_id")) %>%
  purrr::set_names(~sub(".s0.*", "", basename(.)) %>% sub("-.*","",.)) %>%
  purrr::map(~data.table::fread(.x, col.names = c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Counts"))) %>% 
  purrr::imap(~dplyr::transmute(.x, cond = .y, Geneid, Length, Counts)) %>%
  purrr::map(~dplyr::filter(.x, stringr::str_detect(Geneid, paste(noct.info$TRANSCRIPT_ID, collapse = "|")))) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(Geneid) %>%
  tidyr::separate("cond", c("sample", "region"), "\\.") %>%
  dplyr::group_by(sample, region) %>% dplyr::mutate(regionNum = paste0(region, row_number())) %>% dplyr::ungroup() %>%
  dplyr::transmute(tissue = ifelse(stringr::str_detect(sample, "spq"), "Spermatogonia", "Testis"),
                   cond = sub("RNA.","",sample), 
                   strain = sub("spq","",sub("testis","",cond)),
                   sample, regionNum, Length, Counts,
                   age = "adult")


### ====== ###
### DO IER ###
### ====== ###

# Divide counts by length
# Pivot wider to have introns and exons as columns with counts
# divide intron1 by mean of exon1 and exon2, and divide intron2 by mean of exon2 and exon3

### PRIVATE DATASET ###
ier.private <- ei.private %>% 
  dplyr::transmute(sample, cond, strain, tissue, age, regionNum, norm=Counts/Length) %>% 
  tidyr::pivot_wider(names_from = "regionNum", values_from = "norm") %>%
  dplyr::transmute(sample,cond,strain,tissue,age,intron1=intron1/((exon1+exon2)/2),intron2=intron2/((exon2+exon3)/2)) %>%
  reshape2::melt(variable.name="intron", value.name = "ier") %>%
  dplyr::mutate(dataset="private")

### PRIVATE DATASET ###
ier.icr <- ei.icr %>% 
  dplyr::transmute(sample, cond, strain="outbred", tissue, age, regionNum, norm=Counts/Length) %>% 
  tidyr::pivot_wider(names_from = "regionNum", values_from = "norm") %>%
  dplyr::transmute(sample,cond,strain,tissue,age,intron1=intron1/((exon1+exon2)/2),intron2=intron2/((exon2+exon3)/2)) %>%
  reshape2::melt(variable.name="intron", value.name = "ier") %>%
  dplyr::mutate(dataset="icr") %>%
  dplyr::inner_join(icr.sample.info %>% dplyr::select(sample, iap))

### Yu et al., 2019 ###
ier.yu2019 <- ei.yu2019%>% 
  dplyr::transmute(sample, cond, strain, tissue, age, regionNum, norm=Counts/Length) %>% 
  tidyr::pivot_wider(names_from = "regionNum", values_from = "norm") %>%
  dplyr::transmute(sample,cond,strain,tissue,age,intron1=intron1/((exon1+exon2)/2),intron2=intron2/((exon2+exon3)/2)) %>%
  reshape2::melt(variable.name="intron", value.name = "ier") %>%
  dplyr::mutate(dataset="yu2019")

### ======== ###
### DO PLOTS ###
### ======== ###

### PRIVATE DATASET ###
ier.private.plot <- ier.private %>% dplyr::filter(tissue != "Spermatogonia") %>%
     ggplot(aes(strain,ier,color = strain)) + 
     geom_point(size = 3) + ggmitji::stat_line_boxplot(color="black", height=.0005) + facet_wrap(~intron) +
     ggmitji::theme_custom(legend="right", subtitle.face="plain", subtitle.size=6, axis.title.face="plain", axis.title.size=6, axis.text.size = 6) +
     ggmitji::remove_x_axis() +
     labs(y = "IER", x = NULL, subtitle = "Private dataset --> Testis samples",
          caption = "Intron and exon counts were obtained with featureCounts; repeats and gaps +/- 150bp were removed from the annotation.")

### ICR OUTBRED MICE ###
ier.icr.plot <- ier.icr %>%
  ggplot(aes(iap,ier,color = iap)) + 
  geom_point(size = 3) + ggmitji::stat_line_boxplot(color="black", height=.0005) + facet_wrap(~intron) +
  ggpubr::stat_compare_means(method = "wilcox.test", paired = F, comparisons = list(c("IAP +", "IAP -")), tip.length = 0) +
  ggmitji::theme_custom(legend="right", title.face="italic", title.size=9, subtitle.face="plain", subtitle.size=6, axis.title.face="plain", axis.title.size=6, axis.text.size = 6) +
  ggmitji::remove_x_axis() +
  labs(title = "Inverse intron expression ratio (IER) in canonical Noct transcript", y = "IER", x = NULL, subtitle = "ICR outbred mice",
       caption = "Intron and exon counts were obtained with featureCounts; repeats and gaps +/- 150bp were removed from the annotation.")


### Yu et al., 2019 ###
ier.yu2019.plot <- ier.yu2019 %>%
  ggplot(aes(strain,ier,color = strain)) + 
  geom_point(size = 3) + ggmitji::stat_line_boxplot(color="black", height=.0005) + facet_wrap(~intron) +
  ggmitji::theme_custom(legend="right", title.face="italic", title.size=9, subtitle.face="plain", subtitle.size=6, axis.title.face="plain", axis.title.size=6, axis.text.size = 6) +
  ggmitji::remove_x_axis() +
  labs(title = "Inverse intron expression ratio (IER) in canonical Noct transcript", y = "IER", x = NULL, subtitle = "Yu et al., 2019",
       caption = "Intron and exon counts were obtained with featureCounts; repeats and gaps +/- 150bp were removed from the annotation.")


### ========================== ###
### DO PLOTS WITH ALL DATASETS ###
### ========================== ###

# Join all ier datasets
# Add IAP information
ier.all2 <- list(ier.private, 
                ier.icr,
                ier.yu2019) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(iap = ifelse(!is.na(iap), iap, ifelse(stringr::str_detect(cond, "AKR|NOD|BL6"), "IAP +", "IAP -"))) 

# Do plot
ier.all.plot <- ier.all2 %>%
  ggplot(aes(iap,ier,color = iap)) + 
  geom_point(aes(shape=dataset), size = 3) + ggmitji::stat_line_boxplot(color="black", height=.001) + facet_wrap(~intron) +
  ggmitji::theme_custom(legend="right", title.face="italic", title.size=9, subtitle.face="plain", subtitle.size=6, axis.title.face="plain", axis.title.size=6, axis.text.size = 6) +
  ggmitji::remove_x_axis() +
  ggpubr::stat_compare_means(method = "wilcox.test", paired = F, comparisons = list(c(1,2)), tip.length = 0, label.y = .7) +
  labs(title = "Inverse intron expression ratio (IER) in canonical Noct transcript", y = "IER", x = NULL, subtitle = "All datasets - adult samples",
       caption = "Intron and exon counts were obtained with featureCounts; repeats and gaps +/- 150bp were removed from the annotation.")



#####################
### WRITE TO FILE ###
#####################

# Join all ier datasets
# Add IAP information
ier.all <- list(ier.private, 
                ier.icr,
                ier.yu2019) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(iap = ifelse(!is.na(iap), iap, ifelse(stringr::str_detect(cond, "AKR|NOD|BL6"), "IAP +", "IAP -"))) 

### Write data to file
pdf("figures/supp_figs/figS12/figS12B-IER.pdf", width = 9, height = 7)
ier.private.plot
ier.icr.plot
ier.yu2019.plot
ier.all.plot
dev.off()

