# R

# This script:
# - imports counts in exons - transcript_id
# - normalizes them with sizeFactors from piRNA cluster counts normalization.
# - plots counts for canonical Noct transcript

WORKDIR="~/projects/mouse_pic_var_review/";
setwd(WORKDIR)

# Load packages
library(dplyr)
library(data.table)
library(magrittr)
library(stringr)
library(DESeq2)

library(ggplot2)

##################################
## Define parameters
##################################
## Filtering of clusters
removeLowCountClusters <- T
lowCountThreshold <- 10
removeOverlapingClusters <- T

## DESeq parameters
# Adj p value threshold
padj_sig <- 0.05

# Threshold for filtering of DESeq results and test of association
f_sig <- 1 
cflank <- 5000 # Max allowed distance to transposon.
ignoreStrand <- T

###################################



# Import noct canonical info
noct.info <- data.table::fread("data/public/ensembl_noct_canonical.csv")

# Import one to one orthologs in all strains
one2one <- data.table::fread("data/public/genes/mouse_strains_orthologs.csv") %>% na.omit()


########################
### GET SIZE FACTORS ###
########################

# Load exon-gene_id counts
genecounts <- list.files("output/04-rna_process/private/expression/fcounts/exon-gene_id/", "counts$", full.names = T, recursive = T) %>%
  purrr::set_names(sub("RNA","_rep", sub("testis","testis_", sub("spq","spgonia_",sub(".exon.*","",basename(.)))))) %>%
  purrr::imap(~data.table::fread(.x) %>%
                dplyr::select(1,7) %>%
                magrittr::set_colnames(c("Geneid","counts")) %>%
                dplyr::mutate(sample=.y))

# For each sample, take one2one orthologs and maintain the B6 Geneid
b6_one2one <- one2one %>% dplyr::transmute(Gene = mouse_gene, Geneid = mouse_gene)
genecounts.b6 <- genecounts[stringr::str_detect(names(genecounts),"BL6")] %>% purrr::map(~dplyr::inner_join(b6_one2one,.x) %>% dplyr::select(-Geneid)) %>% dplyr::bind_rows()

nod_one2one <- one2one %>% dplyr::transmute(Gene = mouse_gene, Geneid = nod_gene)
genecounts.nod <- genecounts[stringr::str_detect(names(genecounts),"NOD")] %>% purrr::map(~dplyr::inner_join(nod_one2one,.x) %>% dplyr::select(-Geneid)) %>% dplyr::bind_rows()

c3h_one2one <- one2one %>% dplyr::transmute(Gene = mouse_gene, Geneid = c3h_gene)
genecounts.c3h <- genecounts[stringr::str_detect(names(genecounts),"C3H")] %>% purrr::map(~dplyr::inner_join(c3h_one2one,.x) %>% dplyr::select(-Geneid)) %>% dplyr::bind_rows()

x129_one2one <- one2one %>% dplyr::transmute(Gene = mouse_gene, Geneid = `129_gene`)
genecounts.129 <- genecounts[stringr::str_detect(names(genecounts),"129")] %>% purrr::map(~dplyr::inner_join(x129_one2one,.x) %>% dplyr::select(-Geneid)) %>% dplyr::bind_rows()

# Bind all counts from all samples and convert to table
countsDF <- dplyr::bind_rows(list(genecounts.b6, genecounts.nod, genecounts.c3h, genecounts.129)) %>%
  tidyr::pivot_wider(names_from = "sample",values_from = "counts") %>%
  tibble::column_to_rownames("Gene") %>% 
  na.omit()

# Create col data
colDF <- colnames(countsDF) %>% data.frame(sample=., cond=sub("_rep.","",.), Tissue=sub("_.*","",.)) %>% dplyr::mutate(Strain=sub(".*_","",cond)) 

# Do normalization on testis samples
dds <- DESeqDataSetFromMatrix(countData = countsDF, colData = colDF, design = ~ Strain)
dds <- dds[,dds$Strain!="ICR"]
dds <- dds[,dds$Tissue=="testis"]
dds$Strain <- droplevels(dds$Strain)
dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,]
dds <- DESeq(dds)

# Get sizefactors
size_factors <- sizeFactors(dds) 
sizeFactors.df <- sizeFactors(dds) %>% data.frame(sizeFactor=.) %>% tibble::rownames_to_column("sample")

#########################################
### EXPR ON CANONICAL NOCT TRANSCRIPT ###
#########################################

# Read and format counts
# Filter for canonical noct transcript
counts <- list.files("output/04-rna_process/private/expression/fcounts/exon-transcript_id/", "counts$", full.names = T, recursive = T) %>%
  purrr::set_names(sub("RNA","_rep", sub("testis","testis_", sub("spq","spgonia_",sub(".exon.*","",basename(.)))))) %>%
  purrr::imap(~data.table::fread(.x) %>%
               dplyr::filter(stringr::str_detect(Geneid, paste0(noct.info$TRANSCRIPT_ID, collapse = "|"))) %>%
               dplyr::select(1,8) %>%
               magrittr::set_colnames(c("Geneid","counts")) %>%
               dplyr::mutate(sample=.y))

# Create counts table
countsDF <- counts %>% dplyr::bind_rows() %>% dplyr::mutate(Geneid="Noct") %>% tidyr::pivot_wider(names_from = "sample", values_from = "counts") %>% tibble::column_to_rownames("Geneid")

# Create col data
colDF <- colnames(countsDF) %>% data.frame(sample=., cond=sub("_rep.","",.), Tissue=sub("_.*","",.)) %>% dplyr::mutate(Strain=sub(".*_","",cond)) 

# Do normalization on testis samples
dds <- DESeqDataSetFromMatrix(countData = countsDF, colData = colDF, design = ~ Strain)
dds <- dds[,dds$Strain!="ICR"]
dds <- dds[,dds$Tissue=="testis"]
dds$Strain <- droplevels(dds$Strain)
dds$Strain <- relevel(dds$Strain,"BL6")

dds <- dds[rowSums(counts(dds)) >= 10,]
sizeFactors(dds) <- size_factors

###############
### DO PLOT ###
###############

mypch <- 20
pdf("figures/main_figs/fig5/fig5C-Noct_expr_in_exons.pdf", width = 4, height = 4)
plotCounts(dds, gene="Noct", intgroup="Strain", pch=mypch, transform=F, cex=2, ylim=c(0,500))
dev.off()


