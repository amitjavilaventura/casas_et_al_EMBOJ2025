# R

# This script:
# - Imports the counts estimated by featureCounts and creates a counts table
# - Does normalization with DESeq2

workdir <- "~/projects/mouse_pic_var_review/"
setwd(workdir) 

# Load packages
library(dplyr)
library(data.table)
library(purrr)
library(magrittr)
library(tibble)
library(DESeq2)
library(here)

# Define options for later
lowCounts=10

# # Define command line arguments
# args = commandArgs(T)
# anno = args[1] # cluster annotations

annots = c("protrac_merged","lietal_clusters")
for (anno in annots) { ## only if running all together and not one by one

### Import featureCounts files

# List files
fcounts.f <- list.files(here::here("output/mm10_01-srna_process/fcounts/", anno), "counts$", full.names = T, recursive = T)
names(fcounts.f) <- gsub("\\..*", "", basename(fcounts.f)) 

# Keep only s0 counts
fcounts.f <- fcounts.f %>% purrr::keep(~stringr::str_detect(.x, "s0"))

# Read files
fcounts.l <- fcounts.f %>% purrr::map(~data.table::fread(.x, col.names = c("Geneid", "seqnames", "start", "end", "strand", "length", "Counts")))


### Create raw counts table

# Select desired fields
# Join counts of replicates in one table
counts.table <- fcounts.l %>% purrr::imap(~dplyr::select(.x, Geneid, Counts) %>% magrittr::set_colnames(c("Geneid", .y))) %>% plyr::join_all(by = "Geneid", type = "inner")

# Write to file
outdir = here::here("output/mm10_01-srna_process/rawcounts/",anno)
dir.create(outdir, F, T)
write.table(counts.table, paste0(outdir, "/", anno, ".rawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)

}