# R

# This script:
# - Is an adaptation of Pío's script merge_all_clusters_no_intersection.R
# - Takes the mm10 coordinates of all the ENSEMBL Compara orthologs from proTRAC predicted clusters in each sample (they have been converted from mm39 to mm10 first). 
# - Then takes the union of all these clusters.


# Set workding directory
WORKDIR="~/projects/mouse_pic_var_review/"
setwd(WORKDIR)

# Load required packages
library(intervals)
library(dplyr)
library(purrr)

# Define directory where clusters are
CLUSTERSDIR="output/01-pirna_clusters/protrac/protrac_mm10/"

# Define output directory
OUTDIR="output/01-pirna_clusters/protrac/protrac_merged"

# Define files to import
# Retain only those with all (bidir and monodir, and mm10)
l <- list.files(CLUSTERSDIR, "*.bed", full.names = T, recursive = T)
l <- l %>% purrr::discard(~stringr::str_detect(.x, "monodir|bidir")) %>% purrr::keep(~stringr::str_detect(.x, "mm10"))

# Create a list for chromosomes and clusters
list_c <- vector("list", 0)
list_chr <- c()
for (i in l) {
  list_c[[i]] <- read.delim(i, header = FALSE)[, c(1,2,3,6)]  # This differs from Pio because I use BED-like input, Pio uses GTF-like input
  list_chr <- unique(c(list_chr, list_c[[i]]$V1))
}

# Define the number of conditions and the number of samples per condition
type <- 6
sets <- c(4,4,3,2,3,2) # number of replicates per condition, the names are not sample40..., the names are spgonia_bl6..., so the order is not the same as Pio's

# Separate clusters into plus strand, minus strand or bidirectional
lp <- lapply(list_c, function(x){print(x[x[,4]=="+",])})
lm <- lapply(list_c, function(x){print(x[x[,4]=="-",])})
ld <- lapply(list_c, function(x){print(x[x[,4]==".",])})

# Create lists to store the results of each strand
list_res_p <- vector("list", 0)
list_res_m <- vector("list", 0)
list_res_d <- vector("list", 0)

# For each cromosome, take the union of the clusters
#  In the Pio's script, the function Intervals() appends V4 and V5, instead I append V2 and V3 because I have BED file as input, not GTF
for (c in list_chr) {
  
  ### PLUS STRAND
  list_res_p[[c]] <- NULL
  for (sample in lp) {
    sample <- sample[sample$V1 == c, ]
    if (is.null(list_res_p[[c]])) {
      list_res_p[[c]] <- Intervals(matrix(data = append(sample$V2, sample$V3), ncol = 2), closed = c(TRUE, TRUE),  type = "Z")
    } else {
      list_res_p[[c]] <- interval_union(Intervals(matrix(data = append(sample$V2, sample$V3), ncol = 2), closed = c(TRUE, TRUE),  type = "Z"), list_res_p[[c]])
    }
  }
  
  ### MINUS STRAND
  list_res_m[[c]] <- NULL
  for (sample in lm) {
    sample <- sample[sample$V1 == c, ]
    if (is.null(list_res_m[[c]])) {
      list_res_m[[c]] <- Intervals(matrix(data = append(sample$V2, sample$V3), ncol = 2), closed = c(TRUE, TRUE),  type = "Z")
    } else {
      list_res_m[[c]] <- interval_union(Intervals(matrix(data = append(sample$V2, sample$V3), ncol = 2), closed = c(TRUE, TRUE),  type = "Z"), list_res_m[[c]])
    }
  }
  
  ### BIDIRECTIONAL
  list_res_d[[c]] <- NULL
  for (sample in ld) {
    sample <- sample[sample$V1 == c, ]
    if (is.null(list_res_d[[c]])) { list_res_d[[c]] <- Intervals(matrix(data = append(sample$V2, sample$V3), ncol = 2), closed = c(TRUE, TRUE),  type = "Z")
    } else {
      list_res_d[[c]] <- interval_union(Intervals(matrix(data = append(sample$V2, sample$V3), ncol = 2), closed = c(TRUE, TRUE),  type = "Z"), list_res_d[[c]])
    }
  }
  
}

# Create data frame to store the clusters
clusters <- data.frame(chr=character(), start=integer(), end=integer(), strand=character())

# Bind all the clusters in the plus strand to the dataframe
for (n in seq(1,length(list_res_p))) {
  if (length(list_res_p[[n]][,1]) >0) {
    clusters <- rbind(clusters, data.frame(chr=rep(names(list_res_p[n]), length(list_res_p[[n]][,1])), start=list_res_p[[n]][,1], end= list_res_p[[n]][,2], strand = "+"))
  }
}

# Bind all the clusters in the minus strand to the dataframe
for (n in seq(1,length(list_res_m))) {
  if (length(list_res_m[[n]][,1]) >0) {
    clusters <- rbind(clusters, data.frame(chr=rep(names(list_res_m[n]), length(list_res_m[[n]][,1])), start=list_res_m[[n]][,1], end= list_res_m[[n]][,2], strand = "-"))
  }
}

# Bind all the bidirectional clusters to the dataframe
for (n in seq(1,length(list_res_d))) {
  if (length(list_res_d[[n]][,1]) >0) {
    clusters <- rbind(clusters, data.frame(chr=rep(names(list_res_d[n]), length(list_res_d[[n]][,1])), start=list_res_d[[n]][,1], end= list_res_d[[n]][,2], strand = "."))
  }
}

# Order clusteers
clusters$chr <- gsub("chr", "", clusters$chr)
clusters <- clusters[order(clusters$start),]
clusters <- clusters[order(clusters$chr),]
clusters$chr <- paste0("chr", clusters$chr)

# Create fields for GTF and BED
clusters$source     <- "protrac_orthologs_mouse_strains"
clusters$feature    <- "pirna_cluster"
clusters$score      <- "."
clusters$frame      <- "."
row.names(clusters) <- seq(1,length(rownames(clusters)))
clusters$attribute  <- unlist(lapply(row.names(clusters), function(x){paste('gene_id "PTmerged',x,'"', sep="" )}))
clusters$name       <- unlist(lapply(row.names(clusters), function(x){paste('PTmerged',x, sep="" )}))
clusters$length     <- clusters$end-clusters$start

# Create GTF
dir.create(OUTDIR, F, T)
clusters.gtf <- clusters[,c(1,5,6,2,3,7,4,8,9)]
clusters.gtf$start <- clusters.gtf$start+1 # imported files are bed, so I have to add +1 to the start coord
write.table(clusters.gtf, "output/01-pirna_clusters/protrac/protrac_merged/protrac_merged.1nofilt.mm10.gtf", row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)

# Create BED 
clusters.bed        <- clusters[,c(1,2,3,10,11,4)]
clusters.bed$start  <- clusters.bed$start # in Pío's script, here there was a -1, because imported files were GTF, but now are BED.
clusters.bed$length <- clusters.bed$length
write.table(clusters.bed, "output/01-pirna_clusters/protrac/protrac_merged/protrac_merged.1nofilt.mm10.bed", row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
