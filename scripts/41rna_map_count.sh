#!/bin/bash
#SBATCH --job-name=align_RNA_strains                        # Job name
#SBATCH --nodes=1                                           #
#SBATCH --ntasks=1                                          # Run on a single CPU
#SBATCH --cpus-per-task=5
#SBATCH --mem=35gb                                           # Job memory request
#SBATCH --array=0-1
#SBATCH --time=10:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=logs/nod-private-align_RNA_strains_%j.log               # Standard output and error log

set -e
set -u

# This script:
# - Aligns all READS to the corresponding ASSEMBLY for each sample using Hisat2
# - Runs featureCounts in the piC annotations (known and pred) from different assemblies obtained with ENSEMBL Compara Perl API.
# - Runs featureCounts in the genes, introns and exons of each assembly
# - Converts BAM to BED and takes split reads
# - Counts split and not split reads in splice sites

### =============== SET UP =============== ###

# Define name of the study/dataset
DATASETS=( private yu2019 icr)
DATANAME=private # change if want to analize other samples

# Define directory where hisat is
TOOLSDIR=$HOME/.conda/envs/mouse_strains_review/bin/

# Define working directory
WORKDIR=$HOME/projects/mouse_pic_var_review/
cd $WORKDIR

# Define RNAseq out directory
OUTDIR=$WORKDIR/output/04-rna_process/$DATANAME

# Define BAM directory
BAMDIR=$OUTDIR/align
mkdir -p $BAMDIR/hisat2 $BAMDIR/num $BAMDIR/rseqc_bamstat

# Counts DIR
FCOUNTS_DIR=$OUTDIR/expression/fcounts
mkdir -p $FCOUNTS_DIR


# Define directory for featureCounts
BED_COUNTS_DIR=$OUTDIR/expression/bed_counts/
mkdir -p $BED_COUNTS_DIR

# Define temporary directory 
TMPDIR=$WORKDIR/tmp
mkdir -p $TMPDIR

# Define FASTQ directory
if [ "$DATANAME" == "private" ]; then FASTQDIR=$WORKDIR/data/private/RNA-seq/;
elif [ "$DATANAME" == "yu2019" ]; then FASTQDIR=$WORKDIR/data/public/yuetal2019/rnaseq/; 
elif [ "$DATANAME" == "icr" ]; then FASTQDIR=$WORKDIR/data/private/icr_rnaseq;
fi

# Define array of samples
# Define wether it is single end or paired end                                                                     
if [ "$DATANAME" == "private" ]; then
    SAMPLES=( spqBL6RNA2 spqBL6RNA3 spqCASTRNA1 spqCASTRNA3 testis129RNA1 testis129RNA2 testis129RNA3 testisBL6RNA1 testisBL6RNA2 testisC3HRNA1 testisC3HRNA2 testisC3HRNA3 testisNODRNA1 testisNODRNA2 )
    PAIRED="paired"
elif [ "$DATANAME" == "yu2019" ]; then
    SAMPLES=( AKRJ C3H C57BL6J C57BL6NJ LP )
    PAIRED="paired"
elif [ "$DATANAME" == "yu2022" ]; then
    SAMPLES=( BL6_spgonia_rep1 BL6_spgonia_rep2 BL6_spgonia_rep3 BL6_pacSptocytes_rep1 BL6_pacSptocytes_rep2 BL6_pacSptocytes_rep3 BL6_pacSptocytes_rep4 BL6_secSptocytes_rep1 BL6_secSptocytes_rep2 BL6_secSptocytes_rep3 BL6_secSptocytes_rep4 BL6_rSptids_rep1 BL6_rSptids_rep2 BL6_rSptids_rep3 BL6_rSptids_rep4 )
    PAIRED="paired"
elif [ "$DATANAME" == "icr" ]; then
    SAMPLES=( sample01 sample02 sample03 sample04 sample05 sample06 sample07 sample08 sample09 sample10 sample11 sample12 sample13 sample14 sample15 )
    PAIRED="paired"
fi

# Define paired options for featureCounts
# -p to indicate paired end data, -B to count fragments instead of reads
if [ $PAIRED == "paired" ]; then FCOUNTS_PAIRED="-p"; #FCOUNTS_PAIRED="-p -B";
else FCOUNTS_PAIRED=""; 
fi

# Select sample to align
SAMPLE=$(echo ${SAMPLES[@]} | cut -f $SLURM_ARRAY_TASK_ID -d " " ) # COMMENT IF NOT RUN IN PARALLEL

# Define genome for the corresponding sample
echo "Definining the genome depending on the sample $SAMPLE..."
if [[ "$SAMPLE" == *"C57BL6NJ"* ]] ; then GENOME="C57BL_6NJ_v1"; ASSEMBLY=c57bl6nj; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/mm10.ucsc.rmsk.bed;
elif [[ "$SAMPLE" == *"BL6"* ]] ; then GENOME="GRCm38"; ASSEMBLY=mm10; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/mm10.ucsc.rmsk.bed;
elif [[ "$SAMPLE" == *"sample"* ]] ; then GENOME="GRCm38"; ASSEMBLY=mm10; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/mm10.ucsc.rmsk.bed;
elif [[ "$SAMPLE" == *"C3H"* ]] ; then GENOME="C3H_HeJ_v1"; ASSEMBLY=c3h; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
elif [[ "$SAMPLE" == *"129"* ]] ; then GENOME="129S1_SvImJ_v1"; ASSEMBLY=129; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
elif [[ "$SAMPLE" == *"NOD"* ]] ; then GENOME="NOD_ShiLtJ_v1"; ASSEMBLY=nod; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
elif [[ "$SAMPLE" == *"CAST"* ]] ; then GENOME="CAST_EiJ_v1"; ASSEMBLY=cast; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
elif [[ "$SAMPLE" == *"FVB"* ]] ; then GENOME="FVB_NJ_v1"; ASSEMBLY=fvb; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
elif [[ "$SAMPLE" == *"AKR"* ]] ; then GENOME="AKR_J_v1"; ASSEMBLY=akr; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
elif [[ "$SAMPLE" == *"LP"* ]] ; then GENOME="LP_J_v1"; ASSEMBLY=lp; GTF=$WORKDIR/genomes/$GENOME/gene/Mus*$GENOME*100.gtf; RMSK=$WORKDIR/genomes/$GENOME/rmsk/Mus*$GENOME*ucsc.custom.bed;
else echo "The specified sample does not exist. Please check that the name is correct."
fi

echo "Sample $SAMPLE --> Genome $GENOME"

# Define genome index path
GENOMEINDEX=$WORKDIR/genomes/$GENOME/assembly/hisat2_index/$GENOME 
ls $GENOMEINDEX*

### =============== READ MAPPING =============== ###

# Go to output directory
cd $BAMDIR

# Run HISAT2 mapper with the flag --novel-splicesites-outfile
# The index has been created with known exons and splice sites, so giving them now is not necessary
if [ ! -f $BAMDIR/hisat2/$SAMPLE.bam ] ; then 

    # Extract splice-sites from GTF
    if [ ! -f  $TMPDIR/strains.$SAMPLE.gtf  ]; then
        echo "Uncompressing GTF..."
        if [ -f $GTF.gz ]; then zcat $GTF.gz > $TMPDIR/strains.$SAMPLE.gtf
        elif [ -f $GTF ]; then cat $GTF > $TMPDIR/strains.$SAMPLE.gtf
        fi
    fi
    if [ ! -f $TMPDIR/strains.$SAMPLE.hisat2.splicesites.txt ]; then 
            echo "Extracting splice-sites from $GTF..."
            $TOOLSDIR/hisat2_extract_splice_sites.py $TMPDIR/strains.$SAMPLE.gtf > $TMPDIR/strains.$SAMPLE.hisat2.splicesites.txt; 
            fi

    # Run HISAT2
    if [ -f $GENOMEINDEX.1.ht2 ]; then
        if [ ! -f $BAMDIR/hisat2/$SAMPLE.sam ] ; then
            echo "Running Hisat2 to aling reads from sample $SAMPLE to genome $GENOMEINDEX..." 
            # Specify if data is paired or not
            if [ "$PAIRED" == "paired" ]; then $TOOLSDIR/hisat2 --no-softclip -p 5 -x $GENOMEINDEX -1 $FASTQDIR/${SAMPLE}_1.fastq.gz -2 $FASTQDIR/${SAMPLE}_2.fastq.gz -S $BAMDIR/hisat2/$SAMPLE.sam --known-splicesite-infile $TMPDIR/strains.$SAMPLE.hisat2.splicesites.txt;
            else $TOOLSDIR/hisat2 --no-softclip -p 5 -x $GENOMEINDEX -U $FASTQDIR/${SAMPLE}.fastq.gz -S $BAMDIR/hisat2/$SAMPLE.sam --known-splicesite-infile $TMPDIR/strains.$SAMPLE.hisat2.splicesites.txt;
            fi
        fi
        if [ -f $BAMDIR/hisat2/$SAMPLE.sam ] ; then 
            echo "Converting SAM to BAM...."
            $TOOLSDIR/samtools view -F 4 -Sb $BAMDIR/hisat2/$SAMPLE.sam | $TOOLSDIR/samtools sort -m 4G -@ 4 -T $BAMDIR/hisat2/$SAMPLE.bam.tmp -o $BAMDIR/hisat2/$SAMPLE.bam - ; 
        fi
    else
        echo "Genome index $GENOMEINDEX not found. Please create it."
    fi

fi

# Create BAM index
if [ ! -f $BAMDIR/hisat2/$SAMPLE.bam.bai ] ; then if [ -f $BAMDIR/hisat2/$SAMPLE.bam ] ; then echo "Indexing BAM file..."; $TOOLSDIR/samtools index -@ 5 $BAMDIR/hisat2/$SAMPLE.bam $BAMDIR/hisat2/$SAMPLE.bam.bai ; fi; fi

# Remove SAM
if [ -f $BAMDIR/hisat2/$SAMPLE.sam ]; then if [ -f $BAMDIR/hisat2/$SAMPLE.bam ]; then rm $BAMDIR/hisat2/$SAMPLE.sam; fi; fi 

# Getting the number of aligned reads
# Use bam_stat.py with --mapq 1 to remove multimapping reads
$TOOLSDIR/samtools view -c -F 260 $BAMDIR/hisat2/$SAMPLE.bam > $BAMDIR/num/$SAMPLE.hisat2.unique.numReads

### =============== RUN FEATURE COUNTS ON PIRNA CLUSTERS =============== ###

# NOTE THAT FOR ICR OUTBRED MICE, WE USE THE COMPLETE PIRNA CLUSTER ANNOTATIONS, NOT THOSE AFTER ENSEMBL COMPARA

# Define annotations
ANNOTS=(lietal_clusters protrac_merged)
for ANNOT in ${ANNOTS[@]}; do 

    if [[ "$ANNOT" == *"lietal_clusters"* ]] ; then  PIC_GTF="$WORKDIR/output/02-pirna_clusters_orthologs/lietal_clusters/lietal_clusters.$ASSEMBLY.gtf"; FEAT=pirna_cluster; ATTR=gene_id ;
    elif [[ "$ANNOT" == *"protrac_merged"* ]] ; then  PIC_GTF="$WORKDIR/output/02-pirna_clusters_orthologs/protrac_merged/protrac_merged.$ASSEMBLY.gtf"; FEAT=pirna_cluster; ATTR=gene_id ;
    else echo "The specified $ANNOT does not exist";
    fi

    mkdir -p $FCOUNTS_DIR/$ANNOT

    # Uncompress cluster annotation
    if [ ! -f $TMPDIR/$SAMPLE.$ANNOT.$ASSEMBLY.gtf ]; then 
        if [ -f $PIC_GTF.gz ]; then zcat $PIC_GTF.gz > $TMPDIR/$SAMPLE.$ANNOT.$ASSEMBLY.gtf;
        elif [ -f $PIC_GTF ]; then cat $PIC_GTF > $TMPDIR/$SAMPLE.$ANNOT.$ASSEMBLY.gtf;
        fi
    fi
    
    if [ ! -f $FCOUNTS_DIR/$ANNOT/$SAMPLE.$ANNOT.s0.counts ] ; then 
        $TOOLSDIR/featureCounts $FCOUNTS_PAIRED -Q 1 -T 5 -s 0 -O --minOverlap 18 -a $TMPDIR/$SAMPLE.$ANNOT.$ASSEMBLY.gtf -F GTF -t $FEAT -g $ATTR -o $FCOUNTS_DIR/$ANNOT/$SAMPLE.$ANNOT.s0.counts $BAMDIR/$SAMPLE.bam; 
    fi

done


### =============== RUN FEATURE COUNTS for GENE and TRANSCRIPTS =============== ###

# Create output directories
mkdir -p $FCOUNTS_DIR/exon-gene_id/ $FCOUNTS_DIR/exon-transcript_id/ 

# Uncompress annotation
if [ ! -f $TMPDIR/$SAMPLE.$ASSEMBLY.gtf ]; then 
    echo "Uncompressing GTF annotation..."
    if [ -f $GTF.gz ]; then zcat $GTF.gz > $TMPDIR/$SAMPLE.$ASSEMBLY.gtf;
    elif [ -f $GTF ]; then cat $GTF > $TMPDIR/$SAMPLE.$ASSEMBLY.gtf;
    fi
fi

# RUN featureCounts with feature gene and attr gene_id
if [ ! -f $FCOUNTS_DIR/gene-gene_id/fcounts/s0/$SAMPLE.gene-gene_id.s0.counts ] ; then 
   $TOOLSDIR/featureCounts $FCOUNTS_PAIRED -Q 1 -T 5 -s 0 -O --minOverlap 18 -a $TMPDIR/$SAMPLE.$ASSEMBLY.gtf -F GTF -t gene -g gene_id \
    -o $FCOUNTS_DIR/gene-gene_id/$SAMPLE.gene-gene_id.s0.counts $BAMDIR/$SAMPLE.bam; 
fi

# RUN featureCounts with feature exons and attr gene_id
if [ ! -f $FCOUNTS_DIR/exon-gene_id/fcounts/s0/$SAMPLE.exon-gene_id.s0.counts ] ; then 
   $TOOLSDIR/featureCounts $FCOUNTS_PAIRED -Q 1 -T 5 -s 0 -O --minOverlap 18 -a $TMPDIR/$SAMPLE.$ASSEMBLY.gtf -F GTF -t exon -g gene_id \
    -o $FCOUNTS_DIR/exon-gene_id/$SAMPLE.exon-gene_id.s0.counts $BAMDIR/$SAMPLE.bam; 
fi

# RUN featureCounts with feature exons and attr transcript_id   
if [ ! -f $FCOUNTS_DIR/exon-transcript_id/fcounts/s0/$SAMPLE.exon-transcript_id.s0.counts ] ; then 
   $TOOLSDIR/featureCounts $FCOUNTS_PAIRED -Q 1 -T 5 -s 0 -O --minOverlap 18 -a $TMPDIR/$SAMPLE.$ASSEMBLY.gtf -F GTF -t exon -g transcript_id --extraAttributes gene_id \
    -o $FCOUNTS_DIR/exon-transcript_id/$SAMPLE.exon-transcript_id.s0.counts $BAMDIR/$SAMPLE.bam; 
fi



### =============== RUN FEATURE COUNTS for EXONS and INTRONS =============== ###

# Create output directories
mkdir -p $FCOUNTS_DIR/intron-intron_id/ $FCOUNTS_DIR/exon-exon_id/

# Create exon and intron GTF files
if [ ! -f $TMPDIR/$SAMPLE.$ASSEMBLY.intron.gtf ] || [ ! -f $TMPDIR/$SAMPLE.$ASSEMBLY.exon.gtf ]; then

    # Subtract transcripts and exons from GTF and reorder them to have transcripts id as seqname for BED file, add new exon id to exons BED file
    awk 'BEGIN{ OFS=FS="\t" }{ if ($3 == "exon") print $1,$4-1,$5,$7,$9 }' $TMPDIR/$SAMPLE.$ASSEMBLY.gtf | sed 's/gene_id \"//g'  | sed 's/\";.* transcript_id \"/\t/g' | sed 's/\".*//g' | awk 'BEGIN{OFS=FS="\t"}{print $6,$2,$3,$4,$5,$1,$5"::"$6"::"$1":"$2"-"$3"__exon"}' > $TMPDIR/$SAMPLE.$ASSEMBLY.exons.bed.tmp
    awk 'BEGIN{ OFS=FS="\t" }{ if ($3 == "transcript") print $1,$4-1,$5,$7,$9 }' $TMPDIR/$SAMPLE.$ASSEMBLY.gtf | sed 's/gene_id \"//g'  | sed 's/\";.* transcript_id \"/\t/g' | sed 's/\".*//g' | awk 'BEGIN{OFS=FS="\t"}{print $6,$2,$3,$4,$5,$1}' > $TMPDIR/$SAMPLE.$ASSEMBLY.transcript.bed.tmp

    # Subtract exons from transcripts to create introns BED file (chr,start,end,intron_id_coords,width,strand,gene_id,transcriptid)
    $TOOLSDIR/bedtools subtract -a $TMPDIR/$SAMPLE.$ASSEMBLY.transcript.bed.tmp -b $TMPDIR/$SAMPLE.$ASSEMBLY.exons.bed.tmp \
     | awk 'BEGIN{ OFS=FS="\t" }{ print $6,$2,$3,$5"::"$1"::"$6":"$2"-"$3"__intron",$3-$2,$4,$5,$1 }' > $TMPDIR/$SAMPLE.$ASSEMBLY.intron.bed

    # Reorder exon BED file (chr,start,end,exon_id_coords,width,strand,gene_id,transcriptid)
    awk 'BEGIN{ OFS=FS="\t" }{ print $6,$2,$3,$7,$3-$2,$4,$5,$1 }' $TMPDIR/$SAMPLE.$ASSEMBLY.exons.bed.tmp > $TMPDIR/$SAMPLE.$ASSEMBLY.exon.bed

    # Convert BED to GTF
    awk 'BEGIN{ OFS=FS="\t" }{ print $1,"gtf","exon",$2+1,$3,".",$6,".","gene_id \""$7"\"; transcript_id \""$8"\"; exon_id \""$4"\";" }' $TMPDIR/$SAMPLE.$ASSEMBLY.exon.bed > $TMPDIR/$SAMPLE.$ASSEMBLY.exon.gtf
    awk 'BEGIN{ OFS=FS="\t" }{ print $1,"gtf","intron",$2+1,$3,".",$6,".","gene_id \""$7"\"; transcript_id \""$8"\"; intron_id \""$4"\";" }' $TMPDIR/$SAMPLE.$ASSEMBLY.intron.bed > $TMPDIR/$SAMPLE.$ASSEMBLY.intron.gtf

fi

# RUN featureCounts with feature exons and attr exon_id   
if [ ! -f $FCOUNTS_DIR/exon-exon_id/$SAMPLE.exon-exon_id.s0.counts ] ; then 
    $TOOLSDIR/featureCounts $FCOUNTS_PAIRED -Q 1 -T 5 -s 0 -O --minOverlap 18 -a $TMPDIR/$SAMPLE.$ASSEMBLY.exon.gtf -F GTF -t exon -g exon_id \
     -o $FCOUNTS_DIR/exon-exon_id/$SAMPLE.exon-exon_id.s0.counts $BAMDIR/$SAMPLE.bam; 
fi

if [ ! -f $FCOUNTS_DIR/intron-intron_id/$SAMPLE.intron-intron_id.s0.counts ] ; then 
    $TOOLSDIR/featureCounts $FCOUNTS_PAIRED -Q 1 -T 5 -s 0 -O --minOverlap 18 -a $TMPDIR/$SAMPLE.$ASSEMBLY.intron.gtf -F GTF -t intron -g intron_id \
     -o $FCOUNTS_DIR/intron-intron_id/$SAMPLE.intron-intron_id.s0.counts $BAMDIR/$SAMPLE.bam; 
fi


### =============== BAM 2 BED + SEPARATE SPLIT READS =============== ###

# Convert BAM2Bed (not useful because it is paired end and either separates the two mates or does not take the splice reads into account)
if [ ! -f $BAMDIR/hisat2/$SAMPLE.bed ] ; then 
    if [ -f $BAMDIR/hisat2/$SAMPLE.bam ]; then
        mkdir -p bed
        echo "Converting BAM 2 BED for sample $SAMPLE..."
        $TOOLSDIR/bedtools bamtobed -cigar -i $BAMDIR/hisat2/$SAMPLE.bam > $BAMDIR/hisat2/$SAMPLE.bed 
    fi
fi

# Separate split from notsplit reads
if [ -f $BAMDIR/hisat2/$SAMPLE.bed ] ; then 
    if [ ! -f $BAMDIR/hisat2/$SAMPLE.split.bed ]; then  echo  "Getting split reads from $SAMPLE..."; awk 'BEGIN{ FS=OFS="\t" }{ if ($7 ~ /N/) print $1,$2,$3,$4,$5,$6 }' $BAMDIR/hisat2/$SAMPLE.bed > $BAMDIR/hisat2/$SAMPLE.split.bed; fi
    if [ ! -f $BAMDIR/hisat2/$SAMPLE.notsplit.bed ]; then  echo  "Getting split reads from $SAMPLE..."; awk 'BEGIN{ FS=OFS="\t" }{ if ($7 !~ /N/) print $1,$2,$3,$4,$5,$6 }' $BAMDIR/hisat2/$SAMPLE.bed > $BAMDIR/hisat2/$SAMPLE.notsplit.bed; fi
    rm $BAMDIR/hisat2/$SAMPLE.bed
fi


### =============== CREATE SPLICE SITE GTF =============== ###


### Create splicesites BED files
if [ ! -f $TMPDIR/$SAMPLE.$ASSEMBLY.noct.splicesites.bed ]; then
   
    # Subtract transcripts and exons from GTF and reorder them to have transcripts id as seqname for BED file, add new exon id to exons BED file
    awk 'BEGIN{ OFS=FS="\t"}{ if ($3 == "exon") print $1,$4-1,$5,$7,$9 }' $TMPDIR/$SAMPLE.$ASSEMBLY.gtf | grep Noct | sed 's/gene_id \"//g'  | sed 's/\";.* transcript_id \"/\t/g' | sed 's/\".*//g' | awk 'BEGIN{OFS=FS="\t"}{print $6,$2,$3,$4,$5,$1,$5"::"$6"::"$1":"$2"-"$3"__exon"}' > $TMPDIR/$SAMPLE.$ASSEMBLY.exons.bed.tmp
    awk 'BEGIN{ OFS=FS="\t"}{ if ($3 == "transcript") print $1,$4-1,$5,$7,$9 }' $TMPDIR/$SAMPLE.$ASSEMBLY.gtf | grep Noct | sed 's/gene_id \"//g'  | sed 's/\";.* transcript_id \"/\t/g' | sed 's/\".*//g' | awk 'BEGIN{OFS=FS="\t"}{print $6,$2,$3,$4,$5,$1}' > $TMPDIR/$SAMPLE.$ASSEMBLY.transcript.bed.tmp

    # Subtract exons from transcripts to create introns BED file (chr,start,end,intron_id_coords,width,strand,gene_id,transcriptid)
    $TOOLSDIR/bedtools subtract -a $TMPDIR/$SAMPLE.$ASSEMBLY.transcript.bed.tmp -b $TMPDIR/$SAMPLE.$ASSEMBLY.exons.bed.tmp \
    | awk 'BEGIN{OFS=FS="\t"}{print $6,$2,$3,$5"::"$1"::"$6":"$2"-"$3"__intron",$3-$2,$4,$5,$1}' > $TMPDIR/$SAMPLE.$ASSEMBLY.intron.bed

    # Take intron boundaries and add +7 bp inside exon
    awk 'BEGIN{OFS=FS="\t"}{print $1,$2-7,$2+7,$4"___left",$5,$6}' $TMPDIR/$SAMPLE.$ASSEMBLY.intron.bed > $TMPDIR/$SAMPLE.$ASSEMBLY.ss_left.bed.tmp
    awk 'BEGIN{OFS=FS="\t"}{print $1,$3-7,$3+7,$4"___right",$5,$6}' $TMPDIR/$SAMPLE.$ASSEMBLY.intron.bed > $TMPDIR/$SAMPLE.$ASSEMBLY.ss_right.bed.tmp
    cat $TMPDIR/$SAMPLE.$ASSEMBLY.ss_left.bed.tmp $TMPDIR/$SAMPLE.$ASSEMBLY.ss_right.bed.tmp | sort -k1,1n -k2,2n > $TMPDIR/$SAMPLE.$ASSEMBLY.noct.splicesites.bed
    rm $TMPDIR/$SAMPLE.$ASSEMBLY.ss_*.bed.tmp
fi

### =============== BEDTOOLS INTERSECT ON SPLICESITES =============== ###

mkdir -p $BED_COUNTS_DIR/splicesites/
$TOOLSDIR/bedtools intersect -c -f .5 -a $TMPDIR/$SAMPLE.$ASSEMBLY.noct.splicesites.bed -b $BAMDIR/$SAMPLE.noct.split.bed > $BED_COUNTS_DIR/splicesites/$SAMPLE.split.noct.ss.bed
$TOOLSDIR/bedtools intersect -c -f .5 -a $TMPDIR/$SAMPLE.$ASSEMBLY.noct.splicesites.bed -b $BAMDIR/$SAMPLE.noct.notsplit.bed > $BED_COUNTS_DIR/splicesites/$SAMPLE.notsplit.noct.ss.bed


### END 


