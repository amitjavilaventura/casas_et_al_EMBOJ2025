# R

# This script downloads the orthologous lncRNA genes between mouse strains

# Load packages
library(dplyr)
library(data.table)

### lncRNA ORTHOLOGS mm10, 129, C3H, NOD, CAST ###

### Download gene orthologs from ENSEMBL Biomart: see https://www.ensembl.org/biomart/martview/
### Use mouse strains, and dataset for one of the mouse strains
### Then select gene ids and orthology types as attributes and download
### Add gene type and transcript type equal lncRNA in the filters section

# Import downloaded data
orthologs <- data.table::fread("data/public/genes/mart_export_orthologs_lncRNA_mouse_strains.txt", na.strings = "") %>%
  magrittr::set_colnames(c("BL6_geneid", "129_geneid", "BL6_homo", "C3H_geneid", "C3H_homo", "NOD_geneid", "NOD_homo", "CAST_geneid", "CAST_homo"))

# Retain only those with one2one orthologs in all strains
one2one_orth <- orthologs %>% 
  na.omit() %>%
  dplyr::filter(if_all(matches("homo"), ~ . == "ortholog_one2one")) %>%
  dplyr::select(-matches("homo")) 

# Add the mouse geneid as an appendix to the other geneids
one2one_orth <- one2one_orth %>% 
  dplyr::mutate(gene=BL6_geneid) %>%
  dplyr::mutate(across(matches("geneid"), ~paste0(., "::", gene))) %>%
  dplyr::select(-gene)

# Write to file
write.table(one2one_orth, here::here("data/public/genes/mouse_strains_orthologs.lncRNAs.csv"), quote = F, sep = ",", col.names = T, row.names = F)
