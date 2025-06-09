# R

# This script:
# - Imports the counts estimated by featureCounts and creates a counts table
# - Imports GTF of the clusters in mm10 and filters by those in counts table
# - Does normalization with DESeq2
# - It does the differential expression analysis
# - It does TEV analysis and Fisher's test
# - It does counts table including multimapping reads

workdir <- "~/projects/mouse_pic_var_review/"
setwd(workdir) 

# Load packages
library(dplyr)
library(data.table)
library(purrr)
library(magrittr)
library(tibble)
library(DESeq2)

# # Define command line arguments
# args = commandArgs(T)
# anno = args[1] # anno = "protrac_merged", anno = "lietal_clusters"

annots = c("protrac_merged","lietal_clusters")
for (anno in annots) {

### ===== DO RAW COUNTS TABLE ===== ###

# Define options for later
lowCounts=10

### Import featureCounts files

# List files
fcounts.f <- list.files(here::here("output/03-srna_process/fcounts/", anno), "counts$", full.names = T, recursive = T)
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
OUT_RAWCOUNTS = here::here("output/03-srna_process/rawcounts/")
dir.create(OUT_RAWCOUNTS, F, T)
write.table(counts.table.full, paste0(OUT_RAWCOUNTS, "/", anno, ".fullRawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(counts.table, paste0(OUT_RAWCOUNTS, "/", anno, ".rawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)


### Take mm10 coordinates of clusters in counts.table
# Once I have the raw counts.table, I read the mm10 coordinates of the clusters and read it to the clusters present in the rawCounts table

# Read the clusters
clusters <- data.table::fread(here::here(paste0("output/02-pirna_clusters_orthologs/",anno,"/",anno,".mm10.gtf.gz"))) %>% dplyr::mutate(V9 = sub("original_pic", "original_pic ", V9))

# Filter the clusters with those present in all strains (counts.table)
clusters.filt <- clusters %>% dplyr::mutate(Geneid = V9 %>% gsub("gene_id \"", "", .) %>% gsub("\".*", "", .)) %>% dplyr::filter(Geneid %in% counts.table$Geneid) %>% dplyr::select(-Geneid)

# Write filtered clusters
OUT_CLUSTERS= here::here("output/03-srna_process/clusters/")
dir.create(OUT_CLUSTERS, F, T)
write.table(clusters.filt, paste0(OUT_CLUSTERS, "/", anno, ".mm10.gtf"), sep = "\t", quote = F, row.names = F, col.names = F)

}

### ===== DO NORMALIZATION, DEA AND FISHER'S TEST ===== ###

for (anno in annots) {

# Define files
tev_f           <- "data/public/nellaker2012/allTEVs.shared.TEVID.mm9.mm10.tab" # TEVs from Nellaker et al., 2012, with mm10 coordinates
count_f         <- here::here(paste0("output/03-srna_process/rawcounts/", anno, ".rawcounts.tsv"))
clustercoords_f <- here::here(paste0("output/03-srna_process/clusters//", anno, ".mm10.gtf"))
t_f             <- "data/private/sample_info.csv"

# Define SOME parameters
tissueIsTestis <- T # to select testis or spermatogonia 
ignoreStrand   <- T # to ignore the strand or not in the overlapping of piCs and TEVs
inverseStrand  <- F # to look whether the overlapping TEV is in the inverse strand or not
cflank         <- 5000 # max allowed distance to transposon.

removeLowCountClusters   <- T
lowCountThreshold        <- 10
removeOverlapingClusters <- T
padj_sig                 <- 0.05 # Adj p value threshold
f_sig                    <- 1 # log2 fold change threshold

ttypes <- c("LINE","SINE","ERV","IAP") # Types of transposons for the analysis


# If the ignoreStrand is T, do not remove transposon with unknown strand.
if(ignoreStrand==T) { removeTEVsWithUnknownStrand <- F } else { removeTEVsWithUnknownStrand <- T }

# Define the strains depending on the selected tissue
myStrains <- c("XBL6", "XNOD", "XC3H", "X129", "XCAST")


### Import data ---------------------------------

# Get piRNA counts
countsDF <- read.table(count_f,  header = TRUE, row.names = 1)

# Get cluster coords 
# Change the bidirectional clusters strand "." to "*"
# Filter undesired chromosomes
# Create a GRanges object with cluster coordinates 
# Remove bidirectional clusters which overlap other clusters.
clustersDF <- read.table(clustercoords_f)
clustersDF[clustersDF[,7] == ".",7] <- "*" 
clustersDF <- clustersDF[clustersDF[,1] != "Un_JH584304" &  clustersDF[,1] != "4_GL456350_random" & clustersDF[,1] != "4_JH584293_random" & clustersDF[,1] != "4_JH584294_random" ,]
clustersGR <- GRanges(seqnames=clustersDF[,1], IRanges(start=clustersDF[,4], end=clustersDF[,5]), strand=clustersDF[,7], clusternames=clustersDF[,10])

if (removeOverlapingClusters) {
  clustersGR <- clustersGR[!(countOverlaps(clustersGR) > 1 & strand(clustersGR) == "*")]
  countsDF   <- countsDF[rownames(countsDF) %in% clustersGR$clusternames, ]
}


# Import TEVs
# Remove TEVs that were not converted to mm10
# Remove clusters with unknown strand (if desired)
# Change clusters with strand N or 0 to "*" (if they are not removed in the previous step)
# Change the strand (if desired)
# Create a GRanges object with TEV coordinates from Nellaker 
# Separate unstranded TEVs
tevDF <- read.table(tev_f, header=T)
tevDF <- tevDF[!is.na(tevDF[,"start_mm10"]),]
if( removeTEVsWithUnknownStrand ){ tevDF <- tevDF[tevDF[,"strand_mm10"] == "+" | tevDF[,"strand_mm10"] == "-",] } 
tevDF[tevDF[,"strand_mm10"] %in% c("N","0"),"strand_mm10"] <- "*"

if(inverseStrand){
  tevDF[tevDF[,"strand_mm10"] == "+","strand_mm10"] <- "1"
  tevDF[tevDF[,"strand_mm10"] == "-","strand_mm10"] <- "+"
  tevDF[tevDF[,"strand_mm10"] == "1","strand_mm10"] <- "-"
}

tevGR <- GRanges(seqnames=gsub('chr','',tevDF[,"chr_mm10"]), IRanges(start=tevDF[,"start_mm10"],end=tevDF[,"end_mm10"]), strand=tevDF[,"strand_mm10"],
                 TEVtype=tevDF[,"TEtype"], XBL6=tevDF[,"C57B6_ref"],XNOD=tevDF[,"NOD"],X129=tevDF[,"X129S1"],XC3H=tevDF[,"C3H"],XCAST=tevDF[,"CAST"])

unstrandedTEVs <- which(strand(tevGR) == "*")


# Get strain names for samples
# Import sample information
# Add column names to the sample information
# Add X to the strain names (e.g., 129 --> X129)
# Convert strain field to factor
# Relevel strain field to have BL6 as the first level
# Add samplenames as rownames
t_names            <- read.csv(t_f, header=F)[,1:8]
colnames(t_names)  <- c("Sample","Condition","Assembly","Number","Strain","Tissue","Replicate","Batch")
t_names[,"Strain"] <- paste0("X",t_names[,"Strain"]) 
t_names$Strain     <- factor(t_names$Strain)
t_names$Strain     <- relevel(t_names$Strain, "XBL6")
rownames(t_names)  <- t_names$Sample 
 

### DIFFERENTIAL EXPRESSION ---------------------

# Get DESeq dataset from matrix
# Remove clusters with low counts
# Separate dds object in testis and spgonia
# Drop the levels of the dds objects
# Run DESeq
dds <- DESeqDataSetFromMatrix(countData = countsDF, colData = t_names[match(names(countsDF),row.names(t_names)),, drop=F], design = ~Strain)
if(removeLowCountClusters){ dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,] }

dds_testis  <- dds[,dds$Tissue=="testis"]
dds_spgonia <- dds[,dds$Tissue!="testis"]

dds_testis$Strain  <- droplevels(dds_testis$Strain)
dds_spgonia$Strain <- droplevels(dds_spgonia$Strain)

dds_testis  <- DESeq(dds_testis)
dds_spgonia <- DESeq(dds_spgonia)

### GET NORMALIZED COUNTS
norm_spgonia <- counts(dds_spgonia, normalized = T) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid")
norm_testis  <- counts(dds_testis, normalized = T) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid")
norm.counts  <- dplyr::inner_join(norm_spgonia, norm_testis)

# Write all to file
OUT_DEA_FISHER=here::here("output/03-srna_process/fisher_DE_TEV/",anno)
dir.create(OUT_DEA_FISHER, F, T)
write.table(norm.counts, paste0(OUT_DEA_FISHER,"/",anno,".normcounts.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)


### Fishers test --------------------------------

# Define function to perform DEG contrast, overlap with TEVs and fisher's tests in a pairwise manner.
deg_tev_fishers <- function(dds = dds_testis, strain1="XC3H", strain2="XBL6", ttypes = c("LINE","SINE","ERV","IAP")) {
  
  # Define empty vectors for Odds ratios, pvalues and strain-pairs.
  ors         <- c() # odds ratios
  pvals       <- c() # p-values
  strainPairs <- c() # strain pairs, not necessary
  
  # If strain is 129, add X
  if(strain1 =="129" ){ strain1 <- "X129" } else if ( strain2 == "129" ){ strain2 <- "X129" }
  
  # Get DESeq2 results based on contrast of strains 1 and 2
  myres <- na.omit(results(dds, contrast = c("Strain", strain1, strain2), alpha=padj_sig, lfcThreshold=f_sig,  altHypothesis = "greaterAbs"))

  # Define the strain pairs
  strainPairs <- c(strainPairs,paste(strain1,"vs",strain2)) 
  # Remove the X from the strain names
  strainPairs <- gsub("X", "", strainPairs)
  
  # Get clusters overlapping TEVs
  ## - LINEs 
  ## - SINEs 
  ## - ERVs 
  ## - IAPs 
  NE_LINE_dif <- subsetByOverlaps( ignore.strand = ignoreStrand, clustersGR, tevGR[(tevGR$TEVtype == "LINE" | tevGR$TEVtype == "LINE_frag") & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]),], maxgap = cflank)$clusternames
  NE_SINE_dif <- subsetByOverlaps( ignore.strand = ignoreStrand, clustersGR, tevGR[(tevGR$TEVtype == "SINE") & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]),], maxgap = cflank)$clusternames
  NE_ERV_dif  <- subsetByOverlaps( ignore.strand = ignoreStrand, clustersGR, tevGR[(tevGR$TEVtype != "SINE" & tevGR$TEVtype != "LINE" & tevGR$TEVtype != "LINE_frag"  ) & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]), ], maxgap = cflank)$clusternames
  NE_IAP_dif  <- subsetByOverlaps( ignore.strand = ignoreStrand, clustersGR, tevGR[(tevGR$TEVtype == "IAP-I") & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]), ], maxgap = cflank)$clusternames
  
  # Create list of clusters overlapping whith each type of TEs  
  NE_TEV_dif_list <- list(LINE=NE_LINE_dif,SINE=NE_SINE_dif,ERV=NE_ERV_dif,IAP=NE_IAP_dif)
  
  # Separate clusters into significant and not significant (diff. expr)
  sig    <- row.names(myres[which(myres$padj < padj_sig & abs(myres$log2FoldChange)>f_sig),])
  notSig <- row.names(myres[(!row.names(myres) %in% row.names(myres[which(myres$padj < padj_sig & abs(myres$log2FoldChange)>f_sig),])),])
  
  # Set all significant and non-significant clusters, USED LATER
  allSig    <- sig
  allNotSig <- notSig
  
  # For all TEV types
  for(k in 1:length(ttypes)){
    
    # Define TEV type
    ttype <- ttypes[k]
    # print(paste0("Perform Fisher's test with differentially expressed clusters overlapping with TEVs: ", ttype) )
    
    # Filter list of DE clusters overlapping TEVs to include only the desired TEV type
    NE_TEV_dif <- unlist(NE_TEV_dif_list[ttype])  
    
    # Remove clusters overlapping TEVs with unknown strand (if desired)
    if(!ignoreStrand){
      # Set significant and not significant clusters
      # Set clusters to skip, overlapping them to TEVs with unknown strand (unknown strand)
      # From the skip clusters, take the ones with strand + or -
      # Filter them from the significant and not significant clusters.
      sig      <- allSig
      notSig   <- allNotSig
      skipClus <- subsetByOverlaps(clustersGR[clustersGR$clusternames %in% NE_TEV_dif,],tevGR[which(strand(tevGR) == "*"),] )
      skipClus <- skipClus[strand(skipClus) == "+" | strand(skipClus) == "-" ,]$clusternames
      sig      <- sig[!sig %in% skipClus]
      notSig   <- notSig[!notSig %in% skipClus]
    }
    
    # Define numbers for the Fisher's test cross-table
    # - significant clusters with TEV
    # - significant clusters without TEV
    # - not significant clusters with TEV
    # - not significant clusters without TEV
    sigWithTEV       <- sig[sig %in% NE_TEV_dif]
    sigWithoutTEV    <- sig[!(sig %in% NE_TEV_dif)]
    notSigWithTEV    <- notSig[notSig %in% NE_TEV_dif]
    notSigWithoutTEV <- notSig[!(notSig %in% NE_TEV_dif)]
    
    # Create the fisher's test matrix with the numbers above
    m <- matrix(c(length(sigWithTEV), length(notSigWithTEV), length(sigWithoutTEV),length(notSigWithoutTEV)),nrow=2, dimnames=list(c("sig","notSig"),c("withTEV","withoutTEV")))
    
    # # Print cross-table
    # print(paste0("Cross-table of DEGs vs ", ttype, " in: ", strain1," vs ",strain2))
    # print(m)
    
    # Run fisher's test
    f     <-  fisher.test(m)
    
    # Add odds ratio to vector
    # Add pvalues to vector
    ors   <- c(ors,round(f$estimate[[1]],2))
    pvals <- c(pvals,round(f$p.value,4))
  
  }
   
  # Add contrast column to ors and pvals
  fishers.df <- data.frame(constrast=strainPairs, TEV=ttypes,odds_ratio=ors, pvalue=pvals )

  # Create dataframe with all the clusters and their TEV status (1 column per TEV)
  predclusterTEVstatus            <- cbind(clustersGR$clusternames, clustersGR$clusternames %in% NE_LINE_dif, clustersGR$clusternames %in% NE_SINE_dif, clustersGR$clusternames %in% NE_ERV_dif, clustersGR$clusternames %in% NE_IAP_dif)
  predclusterTEVstatus            <- as.data.frame(predclusterTEVstatus)
  colnames(predclusterTEVstatus)  <- c("piCID", "LINE", "SINE", "ERV", "IAP")
    
  # Return the results
  return(list("Fishers_res"=fishers.df,"TEVstatus"=predclusterTEVstatus,"DEGs"=myres))
}

# For each contrast, run the function previously defined using all TEVs
fishers_C3HvsBL6 <- deg_tev_fishers(dds = dds_testis, strain1="XC3H", strain2="XBL6", ttypes = c("LINE","SINE","ERV","IAP"))
fishers_NODvsBL6 <- deg_tev_fishers(dds = dds_testis, strain1="XNOD", strain2="XBL6", ttypes = c("LINE","SINE","ERV","IAP"))
fishers_129vsBL6 <- deg_tev_fishers(dds = dds_testis, strain1="X129", strain2="XBL6", ttypes = c("LINE","SINE","ERV","IAP"))
fishers_C3HvsNOD <- deg_tev_fishers(dds = dds_testis, strain1="XC3H", strain2="XNOD", ttypes = c("LINE","SINE","ERV","IAP"))
fishers_129vsNOD <- deg_tev_fishers(dds = dds_testis, strain1="X129", strain2="XNOD", ttypes = c("LINE","SINE","ERV","IAP"))
fishers_C3Hvs129 <- deg_tev_fishers(dds = dds_testis, strain1="XC3H", strain2="X129", ttypes = c("LINE","SINE","ERV","IAP"))

fishers_CASTvsBL6 <- deg_tev_fishers(dds = dds_spgonia, strain1="XCAST", strain2="XBL6", ttypes = c("LINE","SINE","ERV","IAP"))

# Get the fisher's test results for each contrast and bind the different dfs
fishers_res_list <- 
  list("C3HvsBL6" = fishers_C3HvsBL6$Fishers_res, "NODvsBL6" = fishers_NODvsBL6$Fishers_res, "129vsBL6" = fishers_129vsBL6$Fishers_res, 
       "C3HvsNOD" = fishers_C3HvsNOD$Fishers_res, "129vsNOD" = fishers_129vsNOD$Fishers_res, "C3Hvs129" = fishers_C3Hvs129$Fishers_res, 
       "CASTvsBL6" = fishers_CASTvsBL6$Fishers_res)

fishers_pvals <- fishers_res_list %>% purrr::imap(~dplyr::select(.x, 2, 4) %>% magrittr::set_colnames(c("TEV", .y)) ) %>% plyr::join_all(by="TEV")

# Get the TEV status for each contrast and bind the different dfs
TEVstatus_list <- 
  list("C3HvsBL6" = fishers_C3HvsBL6$TEVstatus, "NODvsBL6" = fishers_NODvsBL6$TEVstatus, "129vsBL6" = fishers_129vsBL6$TEVstatus, 
       "C3HvsNOD" = fishers_C3HvsNOD$TEVstatus, "129vsNOD" = fishers_129vsNOD$TEVstatus, "C3Hvs129" = fishers_C3Hvs129$TEVstatus, 
       "CASTvsBL6" = fishers_CASTvsBL6$TEVstatus) %>%
  purrr::imap(~dplyr::mutate(.x, Contrast = .y))

# Get the TEV status for each contrast and bind the different dfs
DEGs_list <- 
  list("C3HvsBL6" = fishers_C3HvsBL6$DEGs, "NODvsBL6" = fishers_NODvsBL6$DEGs, "129vsBL6" = fishers_129vsBL6$DEGs, 
       "C3HvsNOD" = fishers_C3HvsNOD$DEGs, "129vsNOD" = fishers_129vsNOD$DEGs, "C3Hvs129" = fishers_C3Hvs129$DEGs, 
       "CASTvsBL6" = fishers_CASTvsBL6$DEGs) %>%
  purrr::imap(~as.data.frame(.x) %>% 
                dplyr::arrange(padj, desc(abs(log2FoldChange))) %>% 
                tibble::rownames_to_column("piCID") %>%
                dplyr::mutate(DEG = ifelse(padj<=0.05 & log2FoldChange>=1, "Upregulated", ifelse(padj<=.05 & log2FoldChange<=-1, "Downregulated", "NS")), Contrast = .y))

# Merge TEV status and DEG status
TEV_DEG <- list(TEVstatus_list %>% bind_rows(), DEGs_list %>% dplyr::bind_rows()) %>% plyr::join_all(by = c("piCID", "Contrast"), type = "full")
TEV_DEG <- TEV_DEG %>% split(f = as.factor(.$Contrast))

# Write all to file
OUT_DEA_FISHER=here::here("output/03-srna_process/fisher_DE_TEV/",anno)
dir.create(OUT_DEA_FISHER, F, T)

## - Write TEV status to xlsx
## - Write DEG status to xlsx
## - Write Fisher's results to xlsx
## - Write Fisher's pvalues to tsv
openxlsx::write.xlsx(TEVstatus_list, paste0(OUT_DEA_FISHER,"/",anno,".TEVstatus.xlsx"), overwrite = T)
openxlsx::write.xlsx(DEGs_list, paste0(OUT_DEA_FISHER,"/",anno,".DEGs.xlsx"), overwrite = T)
openxlsx::write.xlsx(TEV_DEG, paste0(OUT_DEA_FISHER,"/",anno,".DEGs_TEVs.xlsx"), overwrite = T)
openxlsx::write.xlsx(fishers_res_list, paste0(OUT_DEA_FISHER,"/",anno,".FishersRes.xlsx"), overwrite = T)
write.table(fishers_pvals, paste0(OUT_DEA_FISHER,"/",anno,".FishersPvalues.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

}


### DO RAW COUNTS TABLE WITH MULTI-MAPPING READS
for (anno in annots) {

### Import featureCounts files

# List files
fcounts_multi.f <- list.files(here::here("output/03-srna_process/fcounts_multi/", anno), "counts$", full.names = T, recursive = T)
names(fcounts_multi.f) <- gsub("\\..*", "", basename(fcounts_multi.f)) 

# Keep only s0 counts
fcounts_multi.f <- fcounts_multi.f %>% purrr::keep(~stringr::str_detect(.x, "s0"))

# Read files
fcounts_multi.l <- fcounts_multi.f %>% purrr::map(~data.table::fread(.x, col.names = c("Geneid", "seqnames", "start", "end", "strand", "length", "Counts")))


### Create raw counts table

# Select desired fields
# Join counts of replicates in one table (full join)
# Remove NAs to get only clusters present in all strains
counts_multi.table.full <- fcounts_multi.l %>% purrr::imap(~dplyr::select(.x, Geneid, Counts) %>% magrittr::set_colnames(c("Geneid", .y))) %>% plyr::join_all(by = "Geneid", type = "full")
counts_multi.table <- na.omit(counts_multi.table.full)

# Write to file
outdir = here::here("output/03-srna_process/rawcounts_multi/")
dir.create(outdir, F, T)
write.table(counts_multi.table.full, paste0(outdir, "/", anno, ".fullRawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(counts_multi.table, paste0(outdir, "/", anno, ".rawcounts.tsv"), sep="\t", col.names=T, row.names=F, quote=F)

}
