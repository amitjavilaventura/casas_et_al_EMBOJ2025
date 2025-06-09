#!/bin/bash
#SBATCH --job-name=preprocess_sRNA                     # Job name
#SBATCH --nodes=1                                      #
#SBATCH --ntasks=1                                     # Run on a single CPU
#SBATCH --cpus-per-task=7
#SBATCH --mem=35gb                                     # Job memory request
#SBATCH --array=1-18
#SBATCH --time=72:00:00                                # Time limit hrs:min:sec

set -u; set -e

# This script:
# - Aligns the trimmed-filtered sRNA-seq reads to the genomes of the mouse strains
# - Runs featureCounts on different annotations

### ==================== SET UP ==================== ###

# Set location of software
BOWTIE=$HOME/.conda/envs/mouse_strains_review/bin/bowtie # version 1.2.3
SAMTOOLS=$HOME/.conda/envs/mouse_strains_review/bin/samtools # version 1.10
FEATURECOUNTS=$HOME/.conda/envs/mouse_strains_review/bin/featureCounts # version 2.0.1

# Define working directory
WORKDIR=$HOME/projects/mouse_pic_var_review
cd $WORKDIR

# Define temporary directory
TMPDIR=$WORKDIR/tmp
mkdir -p $TMPDIR

# Define FASTQ to align
SAMPLE_INFO=$WORKDIR/data/private/sample_info.csv
SAMPLES=$(cat $SAMPLE_INFO | cut -f 1 -d,)
SAMPLE=$(echo ${SAMPLES[@]} | cut -f $SLURM_ARRAY_TASK_ID -d " ") # COMMENT IF NOT RUNNING IN PARALLEL

# Define option for each sample
if [[ "$SAMPLE" == *"BL6"* ]] ; then GENOME="GRCm38"; ASSEMBLY="mm10" ;
elif [[ "$SAMPLE" == *"C3H"* ]] ; then GENOME="C3H_HeJ_v1"; ASSEMBLY="c3h" ;
elif [[ "$SAMPLE" == *"129"* ]] ; then GENOME="129S1_SvImJ_v1"; ASSEMBLY="129" ;
elif [[ "$SAMPLE" == *"NOD"* ]] ; then GENOME="NOD_ShiLtJ_v1"; ASSEMBLY="nod" ;
elif [[ "$SAMPLE" == *"CAST"* ]] ; then GENOME="CAST_EiJ_v1"; ASSEMBLY="cast" ;
fi

### ==================== ALIGN READS TO GENOME WITH BOWTIE =============== ###

# Define output directory
OUT_ALIGN=$WORKDIR/output/03-srna_process/align/
mkdir -p $OUT_ALIGN/bowtie $OUT_ALIGN/num

# Define index for the reference genome
INDEX=$WORKDIR/genomes/$GENOME/assembly/bowtie_index/$GENOME

# Define FASTQ file
FASTQ="$WORKDIR/data/private/sRNA-seq/filt/$SAMPLE.filt.fastq"

# Uncompress FASTQ
if [ ! -f $TMPDIR/$SAMPLE.fastq  ] ; then 
    if [ -f $FASTQ.gz ] ; then zcat $FASTQ.gz > $TMPDIR/$SAMPLE.fastq
    elif [ -f $FASTQ.gz ] ; then cat $FASTQ > $TMPDIR/$SAMPLE.fastq
    fi
fi

# Align reads to genome with bowtie and convert SAM to BAM
echo "Alignig reads from $FASTQ to genome $GENOME..."
$BOWTIE -t -M 1 --best --strata -v 1 -q -p 7 --seed 666 -S $INDEX $TMPDIR/$SAMPLE.fastq \
 | $SAMTOOLS view -F 4 -bS - -o $OUT_ALIGN/bowtie/$SAMPLE.bam

# Sort BAM file
$SAMTOOLS sort -m 4G -@ 4 -T $OUT_ALIGN/bowtie/$SAMPLE.bam.tmp -o $OUT_ALIGN/bowtie/$SAMPLE.sorted.bam $OUT_ALIGN/bowtie/$SAMPLE.bam
mv $OUT_ALIGN/bowtie/$SAMPLE.sorted.bam $OUT_ALIGN/bowtie/$SAMPLE.bam

# Index BAM file
echo "Idexing resulting BAM file..."
$SAMTOOLS index $OUT_ALIGN/bowtie/$SAMPLE.bam

# Getting the number of aligned reads
$SAMTOOLS view -c -F 260 $OUT_ALIGN/bowtie/$SAMPLE.bam > $OUT_ALIGN/num/$SAMPLE.bowtie.numReads

### =============== RUN FEATURE COUNTS =============== ###

# Set output directory
OUT_FCOUNTS=$WORKDIR/output/03-srna_process/
mkdir -p $OUT_FCOUNTS/fcounts

# Define annotation
ANNOTATIONS=(lietal_clusters protrac_merged); 

# For each of the defined annotations
# Run featureCounts without allowing for multimapping or reads below -Q 1
for ANNO in ${ANNOTATIONS[@]}; do

    # Define annotations file
    if [ $ANNO == "lietal_clusters" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/lietal_clusters/lietal_clusters.$ASSEMBLY.gtf.gz" ; fi
    if [ $ANNO == "protrac_merged" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/protrac_merged/protrac_merged.$ASSEMBLY.gtf.gz" ; fi

    # Create output directory
    mkdir -p $OUT_FCOUNTS/fcounts/$ANNO/

    # Run featureCounts with s0 s1 and s2
    $FEATURECOUNTS -Q 1 -T 7 -s 0 -a $GTF -F GTF -t pirna_cluster -g gene_id -o $OUT_FCOUNTS/fcounts/$ANNO/$SAMPLE.s0.counts -O --minOverlap 18 $OUT_ALIGN/bowtie/$SAMPLE.bam

done

# Run featureCounts again, but allowing for multimapping reads:
# - Stablish and -M -Q 0, so there won't be filtering by alignment quality and multimapping reads will remain in the count.
for ANNO in ${ANNOTATIONS[@]}; do

    # Define annotations file
    if [ $ANNO == "lietal_clusters" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/lietal_clusters/lietal_clusters.$ASSEMBLY.gtf.gz" ; fi
    if [ $ANNO == "protrac_merged" ]; then GTF="$WORKDIR/output/02-pirna_clusters_orthologs/protrac_merged/protrac_merged.$ASSEMBLY.gtf.gz" ; fi

    # Create output directory
    mkdir -p $OUT_FCOUNTS/fcounts_multi/$ANNO/

    # Run featureCounts with s0 s1 and s2
    $FEATURECOUNTS -M -Q 0 -T 7 -s 0 -a $GTF -F GTF -t pirna_cluster -g gene_id -o $OUT_FCOUNTS/fcounts_multi/$ANNO/$SAMPLE.s0.counts -O --minOverlap 18 $OUT_ALIGN/bowtie/$SAMPLE.bam
    
done

