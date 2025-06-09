#!/bin/bash

# This script:
# - Copies the clusters from Li et al., 2013 into the output directory

# Define working directory
WORKDIR="$HOME/projects/mouse_pic_var_review/"
cd $WORKDIR

# Define output directory
OUTDIR=$WORKDIR/output/01-pirna_clusters/lietal_clusters/
mkdir -p $OUTDIR

# Copy Lietal clusters from data directory to output
# - mm10 --> original
# - mm39 --> for ENSEMBL Compara
cp data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm39.bed $OUTDIR/lietal_clusters.mm39.bed
cp data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed $OUTDIR/lietal_clusters.mm10.bed
