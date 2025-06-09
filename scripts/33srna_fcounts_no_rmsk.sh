#!/bin/bash
#SBATCH --job-name=fcounts_no_rmsk                          # Job name
#SBATCH --nodes=1                                      #
#SBATCH --ntasks=1                                     # Run on a single CPU
#SBATCH --cpus-per-task=4
#SBATCH --mem=35gb                                     # Job memory request
#SBATCH --array=1-18
#SBATCH --time=72:00:00                                # Time limit hrs:min:sec

set -u; set -e

# This script:
# - Removes reads overlapping to repeatMasker annotations from bam files
# - Runs featureCounts on the specified BAM file and GTF annotation.


# SET UP ------------------------------------------------------------------------------------------

# Load software
TOOLSDIR=$HOME/.conda/envs/mouse_strains_review/bin/

# Define working directory
WORKDIR=$HOME/projects/mouse_pic_var_review
cd $WORKDIR

# Define IN_BAM_DIR
IN_BAM_DIR=$WORKDIR/strain_analysis/output/03-srna_process/align/bowtie

# Define output directories
## BAM no RMSK
OUT_BAM_NO_RMSK=$WORKDIR/strain_analysis/output/03-srna_process/align/bowtie_no_rmsk
mkdir -p $OUT_BAM_NO_RMSK/

# FCOUNTS NO RMSK
OUT_FCOUNTS_NO_RMSK=$WORKDIR/output/03-srna_process/fcounts_no_rmsk/
mkdir -p $OUT_FCOUNTS_NO_RMSK/


# Define temporary directory
TMPDIR=$WORKDIR/tmp
mkdir -p $TMPDIR

# Define FASTQ to align
SAMPLE_INFO=$WORKDIR/environment/sample_info.csv
SAMPLES=$(cat $SAMPLE_INFO | cut -f 1 -d,)
SAMPLE=$(echo ${SAMPLES[@]} | cut -f $SLURM_ARRAY_TASK_ID -d " ") # COMMENT IF NOT RUNNING IN PARALLEL

# Define option for each sample
if [[ "$SAMPLE" == *"BL6"* ]] ; then GENOME="GRCm38"; ASSEMBLY="mm10" ;
elif [[ "$SAMPLE" == *"C3H"* ]] ; then GENOME="C3H_HeJ_v1"; ASSEMBLY="c3h" ;
elif [[ "$SAMPLE" == *"129"* ]] ; then GENOME="129S1_SvImJ_v1"; ASSEMBLY="129" ;
elif [[ "$SAMPLE" == *"NOD"* ]] ; then GENOME="NOD_ShiLtJ_v1"; ASSEMBLY="nod" ;
elif [[ "$SAMPLE" == *"CAST"* ]] ; then GENOME="CAST_EiJ_v1"; ASSEMBLY="cast" ;
fi

# Define RMSK annotation
cat $WORKDIR/genomes/$GENOME/rmsk/*.$GENOME*rmsk*.custom.bed | sed 's/^chr//g' > $TMPDIR/$SAMPLE.rmsk.bed

### =============== REMOVE READS OVERLAPPING REPEATS =============== ### 

# Remove reads overlapping repeats using bedtools intersect
$TOOLSDIR/bedtools intersect -v -a $IN_BAM_DIR/$SAMPLE.bam -b $TMPDIR/$SAMPLE.rmsk.bed > $OUT_BAM_NO_RMSK/$SAMPLE.no_rmsk.bam


### =============== RUN FEATURE COUNTS WITH CLUSTERS =============== ### 

# Define piRNA cluster annotations to analyise
ANNOTATIONS=(lietal_clusters protrac_merged); 

For each of the defined annotations
for ANNO in ${ANNOTATIONS[@]}; do

    # Define annotations file
    if [ $ANNO == "lietal_clusters" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/lietal_clusters/lietal_clusters.$ASSEMBLY.gtf.gz" ; fi
    if [ $ANNO == "protrac_merged" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/protrac_merged/protrac_merged.$ASSEMBLY.gtf.gz" ; fi

    # Create output directory
    mkdir -p $OUT_FCOUNTS_NO_RMSK/$ANNO/

    # Run featureCounts with s0 s1 and s2
    $TOOLSDIR/featureCounts -Q 1 -T 7 -s 0 -a $GTF -F GTF -t pirna_cluster -g gene_id -o $OUT_FCOUNTS_NO_RMSK/$ANNO/$SAMPLE.s0.counts -O --minOverlap 18 $OUT_BAM_NO_RMSK/$SAMPLE.no_rmsk.bam

done

# Run featureCounts again, but allowing for multimapping:
# - Stablish and -M -Q 0, so there won't be filtering by alignment quality and multimapping reads will remain in the count.
for ANNO in ${ANNOTATIONS[@]}; do

    # Define annotations file
    if [ $ANNO == "lietal_clusters" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/lietal_clusters/lietal_clusters.$ASSEMBLY.gtf.gz" ; fi
    if [ $ANNO == "protrac_merged" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/protrac_merged/protrac_merged.$ASSEMBLY.gtf.gz" ; fi

    # Create output directory
    mkdir -p $OUT_FCOUNTS_NO_RMSK/fcounts_no_rmsk_multi/$ANNO/

    # Run featureCounts with s0 s1 and s2
    $TOOLSDIR/featureCounts -M -Q 0 -T 7 -s 0 -a $GTF -F GTF -t pirna_cluster -g gene_id -o $OUT_FCOUNTS_NO_RMSK/fcounts_no_rmsk_multi/$ANNO/$SAMPLE.s0.counts -O --minOverlap 18 $OUT_BAM_NO_RMSK/$SAMPLE.no_rmsk.bam

done


### END