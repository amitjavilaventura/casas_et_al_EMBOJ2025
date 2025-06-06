# Mouse piC-variation (review)

Analyses for the review of Casas et al., EMBO J, 2025. 

Casas, Mitjavila et al., **Genetic polymorphisms lead to major, locus-specific, variation in piRNA production in mouse**, *EMBO Journal*, 2025, [DOI: ]().

Pre-print: Casas et al., **Genetic polymorphisms lead to major, locus-specific, variation in piRNA production in mouse**, *bioRxiv*, 2022, [DOI: 10.1101/2022.10.21.513296](https://doi.org/10.1101/2022.10.21.513296).

### TO DO THINGS

0. _**DONE** (can change if we want) **AM** **Change and UPDATE the README file**, specially the links_

1. Add missing scripts for tables:
  + **TV** Table EV1 - sample info + read mapping of ICR samples (inbred samples are already there)
  + **TV** Table EV4 - fraction of variation explained by strain in known piCs
  + **TV** Table EV6 - comparison of genetic and expression distance trees in known piCs
  + **TV** Table EV7 - norm. counts in known piCs on testis samples of ICR
  + **TV** Table EV8 - variance of known piC exression explained by genetics using pedegrees on ICR samples.

2. Add and/or reorganize preprocessing scripts and outputs for the following analyses: 
  + _**DONE** - I think we didn't use SNP data, so we don't need these -**AM** scripts to Download SNP data of mouse strains for review2 analysis_
  + **TV** scripts and output to preprocess data from hybrid mice, both RNA-seq and ChIP-seq H3K4me, including the SNPsplit script if you want. Send them to me and I will adapt them.

3. Add scripts to do figures:
  + **TV** send ALL figure scripts to Adrià, both for preprocessing and generating figures. 
  
4. _**DONE** - **AM** Add scripts to download data in the data folder_

5. **AM** If we consider so, remove all big outputs such as featureCounts, etc (the scripts are there anyways)

6. **TV** Swap A-MYB / TCFL5 labels in Appendix Figure S7C. In the paper, they are swapped in respect to the actual data. The final conclusion does not change and the only reference to figure S7C does not cite A-MYB nor TCFl5 data, so it won't affect the text.

## Directory structure

The organization of the directory is as follows:

```
root ---|--- data/: public data (e.g, lietal clusters) and private data (e.g., sRNA-seq).
        |--- genomes/: genome assemblies and gene/rmsk annotations I use in the analyses, including the scripts to get them.
        |--- tools/: public software (e.g., proTRAC, liftover...), custom scripts (e.g., ENSEMBL Compara Perl API...) used in several analyses or helper functions
        |--- figures/: figures/tables and scripts to create them.
        |--- scripts/: scripts for data processing and analyses.
        |--- outputs/: output files from data processing and analyses.
```

### Where to find...

* [**Union of clusters predicted with proTRAC (mm10 coords)**](https://github.com/amitjavilaventura/mouse-piRNA-variation__review/tree/main/output/01-pirna_clusters/protrac_clusters/protrac_merged). piRNA clusters were predicted in each sample by aligning reads to the genome assembly of the corresponding strain using the proTRAC pipeline. Then predicted piRNA clusters were converted to mm10 with ENSEMBL Compara Perl API (including the clusters in C57BL/6J strain) and then TE variants from *Nellaker et al. (2012)* and repeats from RepeatMasker were filtered. The clusters in this folder are shown in mm10 coordinates. See scripts in directory: [`scripts/1pirna_clusters/protrac/`](https://github.com/amitjavilaventura/mouse-piRNA-variation__review/tree/main/scripts/1pirna_clusters/protrac).

* [**Annotations of piRNA clusters and orthologs in all assemblies (ENSEMBL Compara)**](https://github.com/amitjavilaventura/mouse-piRNA-variation__review/tree/main/output/02-pirna_clusters_orthologs). Clusters in mm10 coordinates were converted to mm39, and these were converted to each genome assembly coordinates (including mm39) using ENSEMBL Compara Perl API. These are the clusters used for further analysis. Have in mind that some clusters are lost after the conversion with ENSEMBL Compara, and not all clusters are present in all strains. See the script: [`scripts/2pirna_cluster_orthologs/1compara_mm10_to_strains.sh`](https://github.com/amitjavilaventura/mouse-piRNA-variation__review/blob/main/scripts/2pirna_cluster_orthologs/1compara_mm10_to_strains.sh)# casas_et_al_EMBOJ2025
