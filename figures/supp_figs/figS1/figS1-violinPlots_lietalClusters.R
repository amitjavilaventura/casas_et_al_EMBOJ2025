# R

############################################
### This script makes the figure showing
### the mean normalized expression of the
### 214 piRNA clusters in the testis samples
### from the inbred mouse strains
############################################

# In the script, I do the normalization, but it shouldn't be necessary, since I already have the normalized counts data

## Cleanup before starting
myvars <- ls()
rm(list=myvars)
rm(myvars)

##################################
## load required libraries
library(DESeq2)
library(GenomicRanges)
library(vioplot)

##################################
## set up variables
##################################
# Project home dir
dir <- "~/projects/mouse_pic_var_review/" # Adrià's home

# Big table with all sample identifiers and other info
t_f <- "data/private/sample_info.csv" 

# Just cluster ID and counts for all samples
count_f <- "output/03-srna_process/rawcounts/lietal_clusters.rawcounts.tsv" 

# File containing piRNA cluster classification into prepachytene, pachytene and hybrid 
clustertype_f <- "data/public/lietal2013/lietal2013-piRNA_clusters.classes.csv"

##################################
## Drawing settings
##################################
mycex <- .5
mypch <- 20
mylwd <- .5

######################################################################
## Get piRNA counts
######################################################################
countsDF <-  read.table(paste0(dir, count_f),  header = TRUE, row.names = 1)
rownames(countsDF) <-  sub("pi-Ip6k1,pi-Ip6k1.","pi-Ip6k1",sub("__piC.*","", rownames(countsDF)))

################################################################
##  Get cluster classes
################################################################
clustertypeDF <- read.table(paste0(dir, clustertype_f), header = T, sep = ",")
clustertypeDF$piCid <- sub("__piC.*","",clustertypeDF$piCid)

################################################################
##  Get strain names for samples
################################################################
t_all <- read.table(paste(dir, t_f, sep = ""), header=F, sep = ",") 
colnames(t_all) <- c("sample","condition","genome","sampleNum", "strain","tissue","replicate","batch","other")
rownames(t_all) <- t_all$sample
t_names <- t_all[,c("strain","tissue"), drop=F]
t_names[,"strain"]<-paste0("X",t_names[,"strain"]) #  Adding X to all strain names, because of strain name 129
# Convert to factors
t_names$strain <- factor(t_names$strain)
# Make sure first level is the control level
t_names$strain <- relevel(t_names$strain, "XBL6")

################################################################
## Normalize expression data using DESeq2
################################################################
# Get the subset of counts that are just for testis inbred strains
dds <- DESeqDataSetFromMatrix(countData = countsDF,
                              colData =  t_names[match(names(countsDF),row.names(t_names)),, drop=F],
                              design = ~ strain)
dds <- dds[,dds$strain!="XICR"]
dds <- dds[,dds$tissue=="testis"]
dds$strain <- droplevels(dds$strain)

## Run default DESeq2 - for count normalization
dds <- DESeq(dds)

################################################################
##  Make violin plots
################################################################

# Get counts matrix normalized
countmatrix <- counts(dds, normalized=T)

# Set logical vectors for each class of clusters
prepach <- rownames(countmatrix) %in% clustertypeDF$piCid[clustertypeDF$class=="Pre-pachytene"]
hybrid <- rownames(countmatrix) %in% clustertypeDF$piCid[clustertypeDF$class=="Hybrid"]
pach <- rownames(countmatrix) %in% clustertypeDF$piCid[clustertypeDF$class=="Pachytene"]

## Make violin plots of normalized, log transformed counts for C57BL6 samples
pdf(file=paste0(dir,"figures/supp_figs/figS1/figS1-violinPlots_lietalClusters.pdf"), paper="a4", pointsize=8, onefile=T)
par(mfrow=c(2,2))
vioplot(
  log(rowMeans(countmatrix[prepach, c("testis_BL6_rep1", "testis_BL6_rep2")])+1),
  log(rowMeans(countmatrix[hybrid, c("testis_BL6_rep1", "testis_BL6_rep2")])+1),
  log(rowMeans(countmatrix[pach, c("testis_BL6_rep1", "testis_BL6_rep2")])+1),
  col="grey30", names=c("prepachytene", "hybrid","pachytene"), ylab="Mean normalized counts (log+1)", main="piRNA cluster expression in whole testis (BL6)", pch=mypch, lwd=mylwd, cex=mycex)

vioplot(
  log(rowMeans(countmatrix[prepach, c("testis_NOD_rep1", "testis_NOD_rep2")])+1), 
  log(rowMeans(countmatrix[hybrid, c("testis_NOD_rep1","testis_NOD_rep2")])+1),
  log(rowMeans(countmatrix[pach, c("testis_NOD_rep1","testis_NOD_rep2")])+1),
  col="darkred", names=c("prepachytene", "hybrid","pachytene"), ylab="Mean normalized counts (log+1)", main="piRNA cluster expression in whole testis (NOD)", pch=mypch, lwd=mylwd, cex=mycex)

vioplot(
  log(rowMeans(countmatrix[prepach, c("testis_C3H_rep1", "testis_C3H_rep2","testis_C3H_rep3")])+1), 
  log(rowMeans(countmatrix[hybrid, c("testis_C3H_rep1","testis_C3H_rep2","testis_C3H_rep3")])+1),
  log(rowMeans(countmatrix[pach, c("testis_C3H_rep1","testis_C3H_rep2","testis_C3H_rep3")])+1),
  col="orange", names=c("prepachytene", "hybrid","pachytene"), ylab="Mean normalized counts (log+1)", main="piRNA cluster expression in whole testis (C3H)", pch=mypch, lwd=mylwd, cex=mycex)

vioplot(
  log(rowMeans(countmatrix[prepach, c("testis_129_rep1", "testis_129_rep2","testis_129_rep3")])+1), 
  log(rowMeans(countmatrix[hybrid, c("testis_129_rep1", "testis_129_rep2","testis_129_rep3")])+1),
  log(rowMeans(countmatrix[pach, c("testis_129_rep1", "testis_129_rep2","testis_129_rep3")])+1),
  col="beige", names=c("prepachytene", "hybrid","pachytene"), ylab="Mean normalized counts (log+1)", main="piRNA cluster expression in whole testis (129)", pch=mypch, lwd=mylwd, cex=mycex)

dev.off()

#write.table(round(counts(dds, normalize=T),digits=2),file=paste0(dir, "docs/SupplementaryTables/InbredStrainTestis_214ClusterNormalizedCounts.txt"),quote=F)
#write.table(assay(dds),file=paste0(dir, "docs/SupplementaryTables/InbredStrainTestis_214ClusterCounts.txt"),quote=F)