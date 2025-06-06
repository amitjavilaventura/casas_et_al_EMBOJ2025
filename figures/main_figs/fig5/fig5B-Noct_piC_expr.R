# R

# This script takes the RNA expression of the whole pi-Noct precursor transcript

# Load packages
library(DESeq2)
library(dplyr)
library("RColorBrewer")
library("pheatmap")



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


##################################
## set up variables
##################################
# Project home dir
dir <- "~/projects/mouse_pic_var_review/"

# Big table with all sample identifiers and other info
t_f <- "data/private/sample_info.csv" 
count_f <- "output/04-rna_process/private/expression/rawcounts/lietal_clusters.rawcounts.tsv" 
clustercoords_f <- "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed"


##################################

##################################
## Drawing settings
##################################
mycex <- .5
mypch <- 20
mylwd <- .5
mySigCol <- "purple"
myNonSigCol <- "black"
myStrains <- c("XBL6", "XNOD", "XC3H", "X129")

mycols <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)

###Colours for strains
### BL6 is grey30
### NOD is darkred
### C3H is orange
### 129 is beige
### CAST is lightgoldenrod1 To get rgb values use col2rgb
# Used n the heatmap
myStrainColors = list(Strain = c(XBL6 = "grey30", X129 = "beige", XC3H = "orange", XNOD = "darkred")) 

######################################################################
## Get piRNA counts
######################################################################

countsDF <-  read.table(count_f,  header = TRUE, row.names = 1, fill=TRUE) %>% magrittr::set_colnames(sub("RNA","_rep",sub("testis","testis_",sub("spq","spgonia_",colnames(.))))) 
rownames(countsDF) <- sub("__piC.*","",rownames(countsDF))

## Create a GRanges object with cluster coordinates and names from Li et al Mol Cell 2013
clustersDF <- read.table(paste0(dir,clustercoords_f))
rownames(clustersDF) <- clustersDF$V4
rownames(clustersDF) <- sub("__piC.*","",rownames(clustersDF))

# clustersDF <- countsDF[,1:5]
clustersDF[clustersDF[,6] == ".",6] <- "*" 
clustersGR <- GRanges(
  seqnames=clustersDF[,1],
  #  IRanges(start=clustersDF[,4], end=clustersDF[,5]),
  #  strand=clustersDF[,7],
  #  clusternames=row.names(countsDF))
  #  clusternames=clustersDF[,10]
  IRanges(start=clustersDF[,2], end=clustersDF[,3]),
  strand=clustersDF[,6],
  clusternames=row.names(clustersDF)
)

## Remove bidirectional clusters which overlap other clusters to avoid redundancy
## 
if (removeOverlapingClusters) {
  clustersGR <-
    clustersGR[!(countOverlaps(clustersGR) > 1 &
                   strand(clustersGR) == "*")]
  countsDF <-
    countsDF[rownames(countsDF) %in% clustersGR$clusternames, ]
}

# countsDF <- countsDF[,-c(1:5)]
countsDF <- na.omit(countsDF)



t_all <- read.table(paste(dir, t_f, sep = ""), header=F, sep = ",", col.names = c("sampleID","condition","assembly","sampleNum","Strain","Tissue","rep","batch","NA"))
rownames(t_all) <- t_all$sampleID
t_names <- t_all[,c("Strain","Tissue"), drop=F]

t_names[,"Strain"]<-paste0("X",t_names[,"Strain"]) #  Adding X to all strain names, due to 129

# Convert to factors
t_names$Strain <- factor(t_names$Strain)

# Make sure first level is the control level
t_names$Strain <- relevel(t_names$Strain, "XBL6")



dds <- DESeqDataSetFromMatrix(countData = countsDF,
                              colData =  t_names[match(names(countsDF),row.names(t_names)),, drop=F],
                              design = ~ Strain)
dds <- dds[,dds$Strain!="XICR"]
dds <- dds[,dds$Tissue=="testis"]
dds$Strain <- droplevels(dds$Strain)

if(removeLowCountClusters){
  dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,]
}

dds <- DESeq(dds)

mypch <- 20
pdf("figures/main_figs/fig5/fig5B-Noct_piC_expr.pdf", width = 4, height = 4)
plotCounts(dds, gene="pi-Ccrn4l", intgroup="Strain", pch=mypch, transform=F, cex=2, ylim=c(0,1000))
dev.off()
