# Scripts

This folder contains the scripts to preprocess the data.

The output of these scripts can be found in the `output` folder.

* `1pirna_clusters`: scripts to obtain known and predicted clusters.
* `2pirna_clusters_orthologs`: scripts to convert mm10 piRNA clusters to other assemblies using lifotver (mm10 --> mm39) and ENSEMBL Compara Perl API (mm39 --> others).
* `3srna_process`: scripts to preprocess (i.e., map, fcounts) small RNA-seq data (mapping data to each assembly).
* `4rna_process`: scripts to preprocess (i.e., map, fcounts) RNA-seq data (mapping data to each assembly).
* `5similarity_analysis`: scripts to obtain the sequences, do pairwise sequence alignments and get "identity" scores of piRNA clusters.
* `6promoter_analysis`: scripts to obtain the sequences of promoters of piRNA clusters, perform pairwise sequence alignments and get the 'identity' scores.
* `mm10_1srna_preprocess`: scripts to preprocess small RNA-seq data mapping reads to GRCm38/mm10 assembly. This is used for ICR data and for synteny analysis.
* `mm10_2rna_preprocess`: scripts to preprocess RNA-seq data mapping rads to GRCm38/mm10 assembly. This is used for the ICR data.