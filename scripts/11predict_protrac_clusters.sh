#!/bin/bash
#SBATCH --job-name=run_protrac                         # Job name
#SBATCH --nodes=1                                           #
#SBATCH --ntasks=1                                          # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb                                           # Job memory request
#SBATCH --array=1-18
#SBATCH --time=10:00:00                                     # Time limit hrs:min:sec

set -e
set -u

# This script:
# - Generates "de novo" predicted clusters with proTRAC for each sample
# - It uses the corresponding assembly for each strain

### SET UP ----------------------------------------------------------------------------------------

# Define singularity paramters
module load singularity

## Sing image to predict clusters
sif_pred="$HOME/singularity/220226_picapp.sif"
sing_exec_pred="singularity exec --bind / $sif_pred"

## Sing image to run ensembl compara
sif_compara="$HOME/singularity/220226_picapp_tools.sif"
sing_exec_compara="singularity exec --bind / $sif_compara"


# Define working directory
WORKDIR="$HOME/projects/mouse_pic_var_review/"
cd $WORKDIR

# Define tools directory
PROTRACDIR="$WORKDIR/tools/protrac/"
LIFTDIR="$WORKDIR/tools/liftover/"
COMPARADIR="$WORKDIR/tools/compara/"

# Define temporary directory
TMPDIR=$WORKDIR/tmp
mkdir -p $TMPDIR

# Go to output directory
OUT_PROTRAC=$WORKDIR/output/01-pirna_clusters/protrac_clusters/protrac_samples/
mkdir -p $OUT_PROTRAC

# Define samples to process, fastq, genome assembly...
SAMPLE_INFO=data/private/sample_info.csv
SAMPLES=$(cat $SAMPLE_INFO | cut -f 1 -d,)
SAMPLE=$(echo ${SAMPLES[@]} | cut -f $SLURM_ARRAY_TASK_ID -d " ")
FASTQ="$WORKDIR/data/private/sRNA-seq/filt/$SAMPLE.filt.fastq"
GENOME=$(cat $SAMPLE_INFO | grep $SAMPLE | cut -f 3 -d, )
ASSEMBLY="$WORKDIR/genomes/$GENOME/assembly/Mus_musculus*.ensembl100.*fa"
CONDITION=$(cat $SAMPLE_INFO | grep $SAMPLE | cut -f 2 -d, | uniq)

echo "Predicting clusters with proTRAC in the sample $SAMPLE using the genome $GENOME..."

# Uncompress FASTQ files
if [ ! -f $TMPDIR/strains.$SAMPLE.fastq ]; then
    if [ -f $FASTQ.gz ]; then zcat $FASTQ.gz > $TMPDIR/strains.$SAMPLE.fastq ;
    elif [ -f $FASTQ ]; then cat $FASTQ > $TMPDIR/strains.$SAMPLE.fastq ;
    fi
fi

# Uncompress GENOME file
if [ ! -f $TMPDIR/strains.$SAMPLE.$GENOME.genome.fa ]; then
    if [ -f $ASSEMBLY.gz ]; then zcat $ASSEMBLY.gz | sed 's/ dna.*REF//g' > $TMPDIR/strains.$SAMPLE.$GENOME.genome.fa ;
    elif [ -f $ASSEMBLY ]; then cat $ASSEMBLY  | sed 's/ dna.*REF//g' > $TMPDIR/strains.$SAMPLE.$GENOME.genome.fa ;
    else echo "Genome assembly ($ASSEMBLY) not found."
    fi
fi

# Go to output directory FOR EACH SAMPLE
OUT_PROTRAC_PRED=$WORKDIR/output/01-pirna_clusters/protrac_clusters/protrac_samples/$CONDITION/$SAMPLE
mkdir -p $OUT_PROTRAC_PRED
cd $OUT_PROTRAC_PRED

### =============== RUN PROTRAC PIPELINE =============== ###

if [ ! -f $SAMPLE.seq.collapsed.no-dust.map.weighted-10000-1000-b-0 ]; then

    # Collapse FASTQ files
    if [ ! -f $TMPDIR/strains.$SAMPLE.seq.collapsed ]; then 
        $sing_exec_pred perl $PROTRACDIR/TBr2_collapse.pl -i $TMPDIR/strains.$SAMPLE.fastq -o $TMPDIR/strains.$SAMPLE.seq.collapsed;
    fi

    # Remove low-quality reads
    if [ ! -f $TMPDIR/strains.$SAMPLE.seq.collapsed.no-dust ]; then 
        $sing_exec_pred perl $PROTRACDIR/TBr2_duster.pl -i $TMPDIR/strains.$SAMPLE.seq.collapsed; 
    fi

    # Map reads to genome
    if [ ! -f $TMPDIR/strains.$SAMPLE.seq.collapsed.no-dust.map ]; then 
        $sing_exec_pred perl $PROTRACDIR/sRNAmapper.pl -input $TMPDIR/strains.$SAMPLE.seq.collapsed.no-dust -genome $TMPDIR/strains.$SAMPLE.$GENOME.genome.fa -alignments best
    fi

    # Reallocate multimapping reads based on transcription level of surroundings
    if [ ! -f $TMPDIR/strains.$SAMPLE.seq.collapsed.no-dust.mapweighted-10000-1000-b-0 ]; then 
        $sing_exec_pred perl $PROTRACDIR/reallocate.pl $TMPDIR/strains.$SAMPLE.seq.collapsed.no-dust.map 10000 1000 b 0
    fi

    # Move weighted map file to the output directory
    mv $TMPDIR/strains.$SAMPLE.seq.collapsed.no-dust.map.weighted-10000-1000-b-0 $TMPDIR/strains.$SAMPLE.map.weighted

fi


# Run proTRAC
if [ ! -f $OUT_PROTRAC/pred_pics/results.table ]; then
    $sing_exec_pred  perl $PROTRACDIR/proTRAC_v244.pl -map $TMPDIR/strains.$SAMPLE.map.weighted -genome $TMPDIR/strains.$SAMPLE.$GENOME.genome.fa 
    ls; mv proTRAC_* $OUT_PROTRAC_PRED/pred_pics
fi

# Clean seqname from predicted clusters GTF
if [ -f $OUT_PROTRAC/pred_pics/clusters.gtf ]; then sed -i 's/ dna:chromosome.*REF//' $OUT_PROTRAC/pred_pics/clusters.gtf; fi


### =============== CONVERT ALL TO MM10 WITH ENSEMBL COMPARA ================ ###

# Define output directory
OUT_PROTRAC_MM10="$WORKDIR/output/01-pirna_clusters/protrac/protrac_mm10/"
mkdir -p $OUT_PROTRAC_MM10

# Define assembly/species name for each condition
if [[ "$SAMPLE" == *"BL6"* ]]; then SPECIES="Mus musculus"; STRAIN="Reference"; COMPARA_NAME="mus_musculus";
elif [[ "$SAMPLE" == *"CAST"* ]]; then SPECIES="Mus musculus"; STRAIN="CAST_EiJ"; COMPARA_NAME="mus_musculus_casteij";
elif [[ "$SAMPLE" == *"C3H"* ]]; then SPECIES="Mus musculus"; STRAIN="C3H/HeJ"; COMPARA_NAME="mus_musculus_c3hhej";
elif [[ "$SAMPLE" == *"NOD"* ]]; then SPECIES="Mus musculus"; STRAIN="NOD_ShiLtJ"; COMPARA_NAME="mus_musculus_nodshiltj";
elif [[ "$SAMPLE" == *"129"* ]]; then SPECIES="Mus musculus"; STRAIN="129S1/SvImJ"; COMPARA_NAME="mus_musculus_129s1svimj";
fi

# Go to directory for predicted piCs
cd $OUT_PROTRAC_PRED

# Create a folder to store results from compara
mkdir -p $OUT_PROTRAC_PRED/mm39_coords

# Convert GTF to BED
#   Also remove those that contain strange chromosomes
cat $OUT_PROTRAC_PRED/pred_pics/clusters.gtf | grep -v "dna:scaffold" | sort -k1,1n -k4,4n | awk -v sample=$SAMPLE 'BEGIN{FS=OFS="\t";n=1}{print "chr"$1,$4-1,$5,sample"-piC_"n,$5-$4+1,$7; n++}' > $OUT_PROTRAC_PRED/pred_pics/clusters.bed

# Copy remove chr to use it for compara
cat $OUT_PROTRAC_PRED/pred_pics/clusters.bed | sed 's/chr//g' > $OUT_PROTRAC_PRED/pred_pics/clusters.4compara.bed

# If sample is BL6, convert to mm39 with liftover
if [[ "$SAMPLE" == *"BL6"* ]]; then
    $LIFTDIR/liftOver -minMatch=0.9 -bedPlus=6 $OUT_PROTRAC_PRED/pred_pics/clusters.bed $LIFTDIR/mm10ToMm39.over.chain.gz $OUT_PROTRAC_PRED/pred_pics/clusters.mm39.bed $OUT_PROTRAC_PRED/pred_pics/clusters.mm39.unmapped
    cat $OUT_PROTRAC_PRED/pred_pics/clusters.mm39.bed | sed 's/chr//g' > $OUT_PROTRAC_PREDpred_pics/clusters.4compara.bed
fi

# Convert BED to compara format --> separate mono from bidir
cat $OUT_PROTRAC_PRED/pred_pics/clusters.4compara.bed | awk 'BEGIN{OFS=FS="\t"}{ if ( $6 ~ /\+/ ) print $1,$2,$3,"1",$4; else if ( $6 ~ /-/  ) print $1,$2,$3,"-1",$4 }' | sort -k1,1n -k2,2n > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.inputCompara.mono.tsv 
cat $OUT_PROTRAC_PRED/pred_pics/clusters.4compara.bed | awk 'BEGIN{OFS=FS="\t"}{ if ( $6 ~ /\./ ) print $1,$2,$3,1,$4 }' | sort -k1,1n -k2,2n > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.inputCompara.bi.tsv 


### Call ENSEMBL Compara Perl API ---------------

# Convert Monodirectional clusters
if [ ! -f $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.mono.tsv ]; then
    $sing_exec_compara perl $COMPARADIR/ensembl_compara_api_murinae_v105.pl $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.inputCompara.mono.tsv "$COMPARA_NAME" > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.mono.tsv
fi

# Convert Bidirectional clusters
if [ ! -f $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.bi.tsv ]; then
    $sing_exec_compara perl $COMPARADIR/ensembl_compara_api_murinae_v105.pl $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.inputCompara.bi.tsv "$COMPARA_NAME" > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.bi.tsv
fi


### Format output from ENSEMBL Compara ----------

# Get clusters that are not present in MSA from ENSEMBL Compara
grep "## Query" $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.mono.tsv | cat >  $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.notInMSA.mono.tsv
grep "## Query" $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.bi.tsv | cat >  $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.notInMSA.bi.tsv

# Format to pass coordinates to Rscript
#  Retrieve only these with GRCm39 coordinates
cat $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.mono.tsv | grep -v "##" | grep "musculus" | grep GRCm39 > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.formattedOutputCompara.mono.tsv
cat $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.outputCompara.bi.tsv | grep -v "##" | grep "musculus" | grep GRCm39 > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.formattedOutputCompara.bi.tsv

# Call RScript to format output and get final regions
# Usage: Rscript script.R input.tsv output.bed
$sing_exec_compara Rscript $COMPARADIR/ensembl_compara_get_regions.R $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.formattedOutputCompara.bi.tsv  $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.bi.mm39.bed.tmp
$sing_exec_compara Rscript $COMPARADIR/ensembl_compara_get_regions.R $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.formattedOutputCompara.mono.tsv  $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mono.mm39.bed

# Change strand info for bidirectional clusters
cat $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.bi.mm39.bed.tmp | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,".",$7,$8}' > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.bi.mm39.bed

# Concatenate mono and bidirectional clusters
cat $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mono.mm39.bed $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.bi.mm39.bed | sort -k1,1 -k2,2n > $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mm39.bed
rm $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.bi.mm39.bed* $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mono.mm39.bed*

# Convert final clusters to mm10
$LIFTDIR/liftOver -minMatch=0.9 -bedPlus=6 $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mm39.bed $LIFTDIR/mm39ToMm10.over.chain.gz $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mm10.bed $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mm10.unmapped

# Copy all files to the same directory 
#   This will be useful to merge all the clusters in one cluster list --> script in strains_analysis/scripts/1_pirna_clusters/protrac/3_merge_protrac_mm10.R
cp $OUT_PROTRAC_PRED/mm39_coords/$SAMPLE.protrac.mm10.bed $OUTDIR