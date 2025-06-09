# R

# This script:
# - Imports the counts estimated by featureCounts and creates a counts table
# - Imports GTF of the clusters in mm10 and filters by those in counts table
# - Does normalization with DESeq2
# - It does the differential expression analysis

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

# # Define command line arguments
# args = commandArgs(T)
# anno = args[1] # anno = "protrac_merged", anno = "lietal_clusters"

datasets = c("private") 
for (dataname in datasets) {

annots = c("protrac_merged","lietal_clusters")
for (anno in annots) {


# Define options for later
lowCounts=10

### Import featureCounts files

# List files
fcounts.f <- list.files(here::here("strains_analysis/output/04-rna_process/",dataname, "expression/fcounts/", anno), "counts$", full.names = T, recursive = T)
names(fcounts.f) <- gsub("\\..*", "", basename(fcounts.f)) 

# Keep only s0 counts
fcounts.f <- fcounts.f %>% purrr::keep(~stringr::str_detect(.x, "s0"))

# Read files
fcounts.l <- fcounts.f %>% purrr::map(~data.table::fread(.x, col.names = c("Geneid", "seqnames", "start", "end", "strand", "length", "Counts")))


### Create raw counts table

# Select desired fields
# Join counts of replicates in one table (full join)
# Remove NAs to get only clusters present in all strains
counts.table.full <- fcounts.l %>% purrr::imap(~dplyr::select(.x, Geneid, Counts) %>% magrittr::set_colnames(c("Geneid", .y))) %>% plyr::join_all(by = "Geneid", type = "full")
counts.table <- na.omit(counts.table.full)

# Write to file
outdir = here::here("strains_analysis/output/04-rna_process/",dataname,"expression/rawcounts/")
dir.create(outdir, F, T)
write.table(counts.table.full, paste0(outdir, "/", anno, ".fullRawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(counts.table, paste0(outdir, "/", anno, ".rawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)


### Take mm10 coordinates of clusters in counts.table
# Once I have the raw counts.table, I read the mm10 coordinates of the clusters and read it to the clusters present in the rawCounts table

# Read the clusters
clusters <- data.table::fread(here::here(paste0("strains_analysis/output/02-pirna_clusters_orthologs/",anno,"/",anno,".mm10.gtf.gz"))) %>% dplyr::mutate(V9 = sub("original_pic", "original_pic ", V9))

# Filter the clusters with those present in all strains (counts.table)
clusters.filt <- clusters %>% dplyr::mutate(Geneid = V9 %>% gsub("gene_id \"", "", .) %>% gsub("\".*", "", .)) %>% dplyr::filter(Geneid %in% counts.table$Geneid) %>% dplyr::select(-Geneid)

# Write filtered clusters
outdir= here::here("strains_analysis/output/04-rna_process/",dataname,"expression/clusters/")
dir.create(outdir, F, T)
write.table(clusters.filt, paste0(outdir, "/", anno, ".mm10.gtf"), sep = "\t", quote = F, row.names = F, col.names = F)

### Normalize with DESeq2

# Create countData and colData
rownames(counts.table) <- NULL
countData=counts.table %>% tibble::column_to_rownames("Geneid")
colData=data.frame(sample=colnames(countData)) %>%
    dplyr::mutate(condition=sub("RNA.$","",sample), 
                   strain=sub("spq","",condition) %>% sub("testis","",.), 
                   tissue=sub("spq.*", "spgonia", condition) %>% sub("testis.*","testis",.))

# Create DESeq data set
dds <- DESeq2::DESeqDataSetFromMatrix(countData,colData,~strain)
if ( dataname == "icr" ){ dds <- DESeq2::DESeqDataSetFromMatrix(countData,colData,~1) }

# Remove clusters with low counts
dds <- dds[rowSums(counts(dds)) >= lowCounts,]

# Separate testis and spgonia
dds_spgonia <- dds[, dds$tissue != "testis"]
dds_testis  <-  dds[, dds$tissue != "spgonia"]

dds_spgonia$strain <- droplevels(dds_spgonia$strain)
dds_testis$strain  <- droplevels(dds_testis$strain)

# Run DESeq
dds_spgonia <- DESeq2::DESeq(dds_spgonia)
dds_testis  <- DESeq2::DESeq(dds_testis)

# # Save DESeq object to RDS
# outdir = here::here("strains_analysis/output/04-rna_process/",dataname,"expression/dds/")
# dir.create(outdir, F, T)
# saveRDS(dds_spgonia, paste0(outdir,"/", anno, ".dds.spgonia.rds"))
# saveRDS(dds_testis, paste0(outdir,"/", anno, ".dds.testis.rds"))


### Get normalized counts
# Get normalized counts
norm_spgonia <- counts(dds_spgonia, normalized = T) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid")
norm_testis  <- counts(dds_testis, normalized = T) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid")

# Write to file
outdir = here::here("strains_analysis/output/04-rna_process/",dataname,"expression/normcounts/")
dir.create(outdir, F, T)
write.table(norm_spgonia, paste0(outdir, "/", anno, ".normcounts.spgonia.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(norm_testis, paste0(outdir, "/", anno, ".normcounts.testis.tsv"), sep="\t", col.names=T, row.names=F, quote=F)

### Get DEGs

# Get DE results
NODvsBL6  <- DESeq2::results(dds_testis, contrast = c("strain", "NOD", "BL6"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")
C3HvsBL6  <- DESeq2::results(dds_testis, contrast = c("strain", "C3H", "BL6"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")
x129vsBL6 <- DESeq2::results(dds_testis, contrast = c("strain", "129", "BL6"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")
C3HvsNOD  <- DESeq2::results(dds_testis, contrast = c("strain", "C3H", "NOD"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")
x129vsNOD  <- DESeq2::results(dds_testis, contrast = c("strain", "129", "NOD"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")
x129vsC3H  <- DESeq2::results(dds_testis, contrast = c("strain", "129", "C3H"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")

CASTvsBL6  <- DESeq2::results(dds_spgonia, contrast = c("strain", "CAST", "BL6"), lfcThreshold = 1, alpha = .05, altHypothesis = "greaterAbs") %>% as.data.frame() %>% dplyr::arrange(padj) #%>% tibble::rownames_to_column("Geneid")

# Create list
deg_list <- list(NODvsBL6,C3HvsBL6,x129vsBL6,C3HvsNOD,x129vsNOD,x129vsC3H, CASTvsBL6) %>%
  purrr::set_names(c("NODvsBL6","C3HvsBL6","129vsBL6","C3HvsNOD","129vsNOD","129vsC3H","CASTvsBL6")) %>%
  purrr::imap(~tibble::rownames_to_column(.x, "piCID") %>%
                dplyr::mutate(DEG="NS", DEG=ifelse(log2FoldChange>1&padj<0.05, "Upregulated", ifelse(log2FoldChange<-1&padj<0.05, "Downregulated", DEG))) %>%
                dplyr::mutate(Contrast=.y) %>%
                dplyr::mutate(piCID=ifelse(piCID=="pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", piCID)))

# Write results to XLSX
outdir = here::here("strains_analysis/output/04-rna_process/",dataname,"expression/degs/")
dir.create(outdir, F, T)
openxlsx::write.xlsx(deg_list, paste0(outdir,"/", anno, ".degs.xlsx"))

} # end of for (anno in annots)

} # end of for (dataname in datasets)
