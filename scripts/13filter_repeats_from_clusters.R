# R

# This script:
# - Is an addaptation of Pio's script: https://github.com/vavouri-lab/piC-variation/blob/withHISAT/Rscripts/remove_clusters_inside_repeats.R
# - Removes clusters from the protrac_merged clusters which that are inside repeats, from TEVs (Nellaker2012) and repeatMakser
# - It also adds an extra filtering step adapted by Tanya's command in order to filter clusters overlapping repeats by (80%)

# Go to working directory 
# Absolute path in the local machine
WORKDIR="~/projects/mouse_pic_var_review/";
setwd(WORKDIR)

# Load required packages
# Only data.table, because it imports the files superfast and repeatmasker is big
library(data.table)

# Set input/output directory
# Set input file: BED instead of GTF because I use liftover and I prefer BED
# Set intermediate mm10 file after liftover
# Set output file after filtering
# Set output file after second filtering
# Set RMSK file for filtering
INDIR="output/01-pirna_clusters/protrac_clusters/protrac_merged/";
INFILE=paste0(INDIR,"/protrac_merged.1nofilt.mm10.bed");
OUTFILE1=paste0(INDIR,"/protrac_merged.2noTEVs.mm10"); # add extension manually, .bed and .gtf
OUTFILE2=paste0(INDIR,"/protrac_merged.3rmsk_filt.mm10"); # add extension manually, .bed and .gtf
TEVFILE="data/public/nellaker2012/allTEVs.shared.TEVID.mm9.mm10.tab"
RMSKFILE="genomes/GRCm38/rmsk/mm10.ucsc.rmsk.bed"

### ========================== ###
### ADAPTATION OF PIO'S SCRIPT ###
### ========================== ###

# Read merged cluster BED file
# Also define field for chr, start and end
f1 <- read.table(INFILE, skip = 0, header = FALSE, sep = "\t");
chr1 <- 1; start1 <- 2; end1 <- 3

# Read TEV file from Nellaker2012 (see folder data/public/nellaker2012)
# Reorder to have mm10 coords bed file
# Also define field for chr, start and end
### comment nellaker2012 and use repeatmasker
f2 <- read.table(TEVFILE, header = T);
f2$TEVID <- paste0(f2$TEVID, "::", f2$TEtype)
f2 <- f2[,c("chr_mm10", "start_mm10", "end_mm10", "TEVID", "score", "strand_mm10")]
colnames(f2) <- paste0("V", 1:6)
f2$V1 <- sub("chr","", f2$V1)
chr2 <- 1; start2 <- 2; end2 <- 3;

# Format and filter dataframes
f3 <- data.frame(name = f1$V4, inc = FALSE);
f1[,chr1] <- gsub("chr","",f1[,chr1]);
f2[,chr2] <- gsub("chr","",f2[,chr2]);
f1 <- f1[complete.cases(f1),];
f2 <- f2[complete.cases(f2),];

# Remove clusters inside repeats
for (i in seq(1,nrow(f1),1)) {
    small_f2 <- f2[(f2[,chr2] == f1[i,chr1]) & abs(f2[,end2]-f1[i,end1])<5000 & abs(f2[,start2]-f1[i,start1])<5000 ,];
    print(nrow(small_f2));
    for (h in rownames(small_f2)) {
        if (((min(f1[i,end1],f2[h,end2]) - max(f1[i,start1],f2[h,start2])) / (f1[i,end1]-f1[i,start1])) > 0.8) { f3[i,]$inc <- TRUE }
    }
};

# Filter the clusters data.frame (f1)
clusters <- f1[!f3$inc,];
clusters$V1 <- sub("^chrchr", "chr", paste0("chr",clusters$V1));


# Convert to GTF
clusters.gtf           <- clusters[,1:2];
colnames(clusters.gtf) <- c("chr", "source");
clusters.gtf$chr       <- sub("chr","",clusters$V1);
clusters.gtf$source    <- "proTRAC_merged_noTEVs";
clusters.gtf$feauture  <- "pirna_cluster";
clusters.gtf$start     <- clusters$V2+1;
clusters.gtf$end       <- clusters$V3;
clusters.gtf$score     <- ".";
clusters.gtf$strand    <- clusters$V6;
clusters.gtf$frame     <- ".";
clusters.gtf$attribute <- paste0('gene_id "', clusters$V4, '";');
# head(clusters.gtf)

# Write to file
write.table(clusters, paste0(OUTFILE1,".bed"), col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t");
write.table(clusters.gtf, paste0(OUTFILE1,".gtf"), col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t");


### ============================= ###
### ADAPTATION OF TANYA'S COMMAND ###
### ============================= ###

# Define bedtools command
# Run bedtools command
BEDTOOLS_CMD = paste("cat", RMSKFILE, "| bedtools sort -i - | bedtools merge -i - | bedtools intersect -a", paste0(OUTFILE1,".bed"), "-b - -wa -f 0.8 -v >", paste0(OUTFILE2,".bed"));
print(BEDTOOLS_CMD);
system(BEDTOOLS_CMD);

# Import filtered bedfile
# Create GTF-like dataframe
# Write GTF to file
clusters.filt.bed <- read.delim(paste0(OUTFILE2,".bed"), header = F, sep = "\t");
#nrow(clusters.filt.bed)

clusters.filt.gtf           <- clusters.filt.bed[,1:2];
colnames(clusters.filt.gtf) <- c("chr", "source");
clusters.filt.gtf$chr       <- sub("chr","",clusters.filt.bed$V1);
clusters.filt.gtf$source    <- "proTRAC_merged_RMSKfilt";
clusters.filt.gtf$feauture  <- "pirna_cluster";
clusters.filt.gtf$start     <- clusters.filt.bed$V2+1;
clusters.filt.gtf$end       <- clusters.filt.bed$V3;
clusters.filt.gtf$score     <- ".";
clusters.filt.gtf$strand    <- clusters.filt.bed$V6;
clusters.filt.gtf$frame     <- ".";
clusters.filt.gtf$attribute <- paste0('gene_id "', clusters.filt.bed$V4, '";');
# head(clusters.filt.gtf);

write.table(clusters.filt.gtf, paste0(OUTFILE2,".gtf"), col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t");

###########
### END ###
###########
