#!/bin/bash
#SBATCH --job-name=bowtie_mm10                          # Job name
#SBATCH --nodes=1                                      #
#SBATCH --ntasks=1                                     # Run on a single CPU
#SBATCH --cpus-per-task=7
#SBATCH --mem=35gb                                     # Job memory request
#SBATCH --array=1-18
#SBATCH --time=72:00:00                                # Time limit hrs:min:sec

set -u; set -e

# This script:
# - Aligns the trimmed-filtered sRNA-seq reads to the mouse genome mm10

# SET UP ------------------------------------------------------------------------------------------

# Set location of software
BOWTIE=$HOME/.conda/envs/mouse_strains_review/bin/bowtie # version 1.2.3
SAMTOOLS=$HOME/.conda/envs/mouse_strains_review/bin/samtools # version 1.10
FEATURECOUNTS=$HOME/.conda/envs/mouse_strains_review/bin/featureCounts # version 2.0.1

# Define working directory
WORKDIR=$HOME/projects/mouse_pic_var_review/
cd $WORKDIR

# Define temp directory
TMPDIR=$WORKDIR/tmp
mkdir -p $TMPDIR

# Output directory for aligning reads to genome
OUT_ALIGN=$WORKDIR/output/mm10_01-srna_process/align/
mkdir -p $OUT_ALIGN/bowtie $OUT_ALIGN/num

# Output directory for feature counts
OUT_FCOUNTS=$WORKDIR/output/mm10_01-srna_process/
mkdir -p $OUT_FCOUNTS/fcounts


# Define index for the reference genome
INDEX=$WORKDIR/genomes/GRCm38/assembly/bowtie_index/$GENOME


# Define SAMPLE to align
SAMPLE_INFO=$WORKDIR/data/private/sample_info_mm10.csv
SAMPLES=$(cat $SAMPLE_INFO | cut -f 1 -d,)
SAMPLE=$(echo ${SAMPLES[@]} | cut -f $SLURM_ARRAY_TASK_ID -d " ") # COMMENT IF NOT RUNNING IN PARALLEL
    
echo "Sample --> $SAMPLE"

# Define FASTQ
FASTQ="$WORKDIR/data/private/sRNA-seq/filt/$SAMPLE.filt.fastq.gz"
    
# Uncompress FASTQ
if [ ! -f $TMPDIR/$SAMPLE.fastq  ] ; then 
    if [ -f $FASTQ.gz ] ; then zcat $FASTQ.gz > $TMPDIR/$SAMPLE.fastq
    elif [ -f $FASTQ.gz ] ; then cat $FASTQ > $TMPDIR/$SAMPLE.fastq
    fi
fi


# =============== ALIGN READS TO REFERENCE GENOME =============== ###

# Align reads to mm10 with bowtie and convert sam to bam with samtools
$BOWTIE -t -M 1 --best --strata -v 1 -q -p 7 --seed 666 -S $INDEX $TMPDIR/$SAMPLE.fastq \
 | $SAMTOOLS view -F 4 -bS - -o $OUT_ALIGN/bowtie/$SAMPLE.bam		

# Sorting bam file
$SAMTOOLS sort -@ 4 -m 4G $OUT_ALIGN/bowtie/$SAMPLE.bam > $OUT_ALIGN/bowtie/$SAMPLE.sorted.bam
mv $OUT_ALIGN/bowtie/$SAMPLE.sorted.bam $OUT_ALIGN/bowtie/$SAMPLE.bam

# Indexing bam file
$SAMTOOLS index $OUT_ALIGN/bowtie/$SAMPLE.bam

# COunting number of reads
$SAMTOOLS view -c -F 260 $OUT_ALIGN/bowtie/$SAMPLE.bam > $OUT_ALIGN/num/$SAMPLE.bowtie.numReads

### =============== RUN FCOUNTS ON PIRNA CLUSTERS =============== ###

# Define annotation
ANNOTATIONS=(lietal_clusters protrac_merged); 

# For each of the defined annotations
for ANNO in ${ANNOTATIONS[@]}; do

    # Define annotations file
    if [ $ANNO == "lietal_clusters" ]; then GTF="$WORKDIR/output/01-pirna_clusters/lietal_clusters/lietal2013-piRNA_clusters.merged.mm10.gtf"; FEAT=gene; ATTR=gene_id ; fi
    if [ $ANNO == "protrac_merged" ]; then GTF="$WORKDIR/output/01-pirna_clusters/protrac_clusters/protrac_merged/protrac_merged.3rmsk_filt.mm10.gtf"; FEAT=pirna_cluster; ATTR=gene_id ; fi

    # Create output directory
    mkdir -p $OUT_FCOUNTS/fcounts/$ANNO/

    # Run featureCounts with s0
    $FEATURECOUNTS -Q 1 -T 7 -s 0 -a $GTF -F GTF -t $FEAT -g $ATTR -o $OUT_FCOUNTS/fcounts/$ANNO/$SAMPLE.s0.counts -O --minOverlap 18 $OUT_ALIGN/bowtie/$SAMPLE.bam

done

### END