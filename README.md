# Mouse piC-variation 

Analyses for the paper Casas et al., EMBO J, 2025. 

Casas, Mitjavila et al., **Genetic polymorphisms lead to major, locus-specific, variation in piRNA production in mouse**, *EMBO Journal*, 2025, [DOI: 10.1038/s44318-025-00475-4](https://doi.org/10.1038/s44318-025-00475-4).

Pre-print: Casas et al., **Genetic polymorphisms lead to major, locus-specific, variation in piRNA production in mouse**, *bioRxiv*, 2022, [DOI: 10.1101/2022.10.21.513296](https://doi.org/10.1101/2022.10.21.513296).

### Where to find...

* [**Annotations of piRNA clusters and orthologs in all assemblies (ENSEMBL Compara)**](https://github.com/amitjavilaventura/mouse-piRNA-variation__review/tree/main/output/02-pirna_clusters_orthologs). Clusters in mm10 coordinates were converted to mm39, and these were converted to each genome assembly coordinates (including mm39) using ENSEMBL Compara Perl API. These are the clusters used for further analysis. Have in mind that some clusters are lost after the conversion with ENSEMBL Compara, and not all clusters are present in all strains. See the script: [`scripts/2pirna_cluster_orthologs/1compara_mm10_to_strains.sh`](https://github.com/amitjavilaventura/mouse-piRNA-variation__review/blob/main/scripts/2pirna_cluster_orthologs/1compara_mm10_to_strains.sh)# casas_et_al_EMBOJ2025
