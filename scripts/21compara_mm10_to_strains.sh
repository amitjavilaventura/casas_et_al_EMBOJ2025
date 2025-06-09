#!/bin/bash
#SBATCH --job-name=COMPARA_mm10_to_strains                     # Job name
#SBATCH --nodes=1                                           #
#SBATCH --ntasks=1                                          # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb                                           # Job memory request
#SBATCH --time=10:00:00                                     # Time limit hrs:min:sec
#SBATCH --array=0

set -eu
pwd; hostname; date

# This script does:
# - Converts the proTRAC merged clusters or lietal clusters from mm10 to mm39
# - Prepares them as input for compara

### Set up ----------------------------------------

# Define singularity stuff
module load singularity
sing_exec="singularity exec --bind / $HOME/singularity/220226_picapp_tools.sif*"


# Define working directory
WORKDIR="$HOME/projects/mouse_pic_var_review/"
cd $WORKDIR

# Define working directory
TMPDIR="$WORKDIR/tmp"
mkdir -p $TMPDIR

# Define tools directory
LIFTDIR="$WORKDIR/tools/liftover/"
COMPARADIR="$WORKDIR/tools/compara/"

# Define input directory
INDIR=$WORKDIR/output/01-pirna_clusters/

# Define output directory
OUTDIR=$WORKDIR/output/02-pirna_clusters_orthologs/
mkdir -p $OUTDIR

# Define annotations
ANNOTS=(lietal_clusters protrac_merged protrac_pio)
ANNOTS=(lietal_clusters)

### Convert from mm10 to strains -------

# For each annotation
for ANNOT in ${ANNOTS[@]}; do 

    # Define the INPUT files
    if [ "$ANNOT" == "lietal_clusters" ]; then INPUT=$INDIR/lietal_clusters/lietal_clusters.mm10.bed ; 
    elif [ "$ANNOT" == "protrac_merged" ]; then INPUT=$INDIR/protrac_clusters/protrac_merged/protrac_merged.3rmsk_filt.mm10.bed ;
    fi

    # Create output directory for each annotation
    mkdir -p $OUTDIR/$ANNOT/

    # Convert from mm10 to mm39
    echo "Lifting clusters in $ANNOT from mm10 to mm39 and formatting for ENSEMBL Compara Perl API..."
    $LIFTDIR/liftOver $INPUT $LIFTDIR/mm10ToMm39.over.chain.gz $TMPDIR/$ANNOT.mm39.bed $TMPDIR/$ANNOT.mm39.unmapped

    # Prepare input for compara
    #   Separate monodirectional (strand = + or -) from bidirectional (strand = .)
    cat $TMPDIR/$ANNOT.mm39.bed | sed 's/chr//g' | awk 'BEGIN{OFS=FS="\t"}{ if ( $6 ~ /\+/ ) print $1,$2,$3,"1",$4; else if ( $6 ~ /-/  ) print $1,$2,$3,"-1",$4 }' | sort -k1,1n -k2,2n > $OUTDIR/$ANNOT/$ANNOT.inputCompara.mono.tsv 
    cat $TMPDIR/$ANNOT.mm39.bed | sed 's/chr//g' | awk 'BEGIN{OFS=FS="\t"}{ if ( $6 ~ /\./ ) print $1,$2,$3,1,$4 }' | sort -k1,1n -k2,2n > $OUTDIR/$ANNOT/$ANNOT.inputCompara.bi.tsv 

    ### ========================== ###
    ### GET ORTHOLOGS WITH COMPARA ###
    ### ========================== ###


    # Call ENSEMBL Compara Perl API

    echo "Calling ENSEMBL Compara Perl API to convert clusters in $ANNOT from mm39 to each assembly..."

    if [ ! -f $OUTDIR/$ANNOT/$ANNOT.outputCompara.mono.tsv  ]; then 
        $sing_exec perl $COMPARADIR/ensembl_compara_api_murinae_v105.pl $OUTDIR/$ANNOT/$ANNOT.inputCompara.mono.tsv "mus_musculus" > $OUTDIR/$ANNOT/$ANNOT.outputCompara.mono.tsv 
    fi
    if [ ! -f $OUTDIR/$ANNOT/$ANNOT.outputCompara.bi.tsv  ]; then 
        $sing_exec perl $COMPARADIR/ensembl_compara_api_murinae_v105.pl $OUTDIR/$ANNOT/$ANNOT.inputCompara.bi.tsv "mus_musculus" > $OUTDIR/$ANNOT/$ANNOT.outputCompara.bi.tsv 
    fi

    ### FORMAT OUTPUT FROM ENSEMBL COMPARA ------
    echo "Formatting output from ENSEMBL Compara Perl API..."

    # Get clusters that are not present in the MSA
    grep "## Query" $OUTDIR/$ANNOT/$ANNOT.outputCompara.mono.tsv | cat > $OUTDIR/$ANNOT/$ANNOT.notInMSA.mono.tsv
    grep "## Query" $OUTDIR/$ANNOT/$ANNOT.outputCompara.bi.tsv | cat > $OUTDIR/$ANNOT/$ANNOT.notInMSA.bi.tsv


    # Format to pass coordinates to Rscript
    cat $OUTDIR/$ANNOT/$ANNOT.outputCompara.mono.tsv | grep -v "##" | grep "musculus" > $OUTDIR/$ANNOT/$ANNOT.formattedOutputCompara.mono.tsv
    cat $OUTDIR/$ANNOT/$ANNOT.outputCompara.bi.tsv | grep -v "##" | grep "musculus" > $OUTDIR/$ANNOT/$ANNOT.formattedOutputCompara.bi.tsv

    ### GET THE COORDINATES OF THE ORTHOLOGS ------
    echo "Getting the coordinates of the orthologs from $ANNOT..."

    # Call RScript to format output and get final regions
    # - Usage: Rscript script.R input.tsv output.bed
    module load R
    module load R-bundle-Bioconductor
    module load BEDTools
    # $sing_exec \
    Rscript $COMPARADIR/ensembl_compara_get_regions.R $OUTDIR/$ANNOT/$ANNOT.formattedOutputCompara.bi.tsv  $OUTDIR/$ANNOT/$ANNOT.bi.tsv.tmp.tmp
    # $sing_exec \
    Rscript $COMPARADIR/ensembl_compara_get_regions.R $OUTDIR/$ANNOT/$ANNOT.formattedOutputCompara.mono.tsv  $OUTDIR/$ANNOT/$ANNOT.mono.tsv.tmp

    # Change strand info for bidirectional clusters
    cat $OUTDIR/$ANNOT/$ANNOT.bi.tsv.tmp.tmp | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,".",$7,$8}' > $OUTDIR/$ANNOT/$ANNOT.bi.tsv.tmp

    # Concatenate mono and bidirectional clusters
    cat $OUTDIR/$ANNOT/$ANNOT.mono.tsv.tmp $OUTDIR/$ANNOT/$ANNOT.bi.tsv.tmp | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.all.tsv
    rm $OUTDIR/$ANNOT/*tmp

    # Separate the clusters depending on the assembly
    echo "Separate the coordinates for each assembly..."
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "GRCm39" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.mm39.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "129S1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.129.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "CAST_EiJ_v1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.cast.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "NOD_ShiLtJ" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.nod.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "C3H_HeJ_v1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.c3h.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "AKR_J_v1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.akr.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "FVB_NJ_v1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.fvb.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "LP_J_v1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.lp.bed
    cat $OUTDIR/$ANNOT/$ANNOT.all.tsv | grep "C57BL_6NJ_v1" | grep -v NA | sort -k1,1 -k2,2n > $OUTDIR/$ANNOT/$ANNOT.c57bl6nj.bed

    # Convert mm39 to mm10 with liftover
    $LIFTDIR/liftOver -bedPlus=6 $OUTDIR/$ANNOT/$ANNOT.mm39.bed $LIFTDIR/mm39ToMm10.over.chain.gz $OUTDIR/$ANNOT/$ANNOT.mm10.bed $OUTDIR/$ANNOT/$ANNOT.mm39.mm10.unmapped 

    # Create GTF file
    echo "Creating GTF files..."
    for ASSEMBLY in mm39 mm10 129 cast nod c3h; do
        
        # BED 2 GTF
        cat $OUTDIR/$ANNOT/$ANNOT.$ASSEMBLY.bed \
         | sed 's/^chr//' \
         | sort -k1,1n -k2,2n \
         | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1, "compara", "pirna_cluster", $2+1, $3, ".", $6, ".", "gene_id \""$4"\"; original_pic\""$7"\";"}' \
         > $OUTDIR/$ANNOT/$ANNOT.$ASSEMBLY.gtf
    
        # GZip GTF
        gzip --fast --force $OUTDIR/$ANNOT/$ANNOT.$ASSEMBLY.gtf

    done

done

echo "Done!"