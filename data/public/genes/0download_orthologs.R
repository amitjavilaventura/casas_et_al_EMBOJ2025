# R

# This script downloads the orthologous genes between mouse strains

# Load packages
library(dplyr)
library(data.table)

### ORTHOLOGS mm10, 129, C3H, NOD, CAST ###

### Download gene orthologs from ENSEMBL Biomart: see https://www.ensembl.org/biomart/martview/
### Use mouse strains, and dataset for one of the mouse strains
### Then select gene ids and orthology types as attributes and download

# Import downloaded data
orthologs <- data.table::fread("data/public/genes/mart_export.txt", na.strings = "", col.names = c("cast_gene", "cast_homo", "mouse_homo", "129_gene", "nod_gene", "c3h_gene", "c3h_homo", "nod_homo", "mouse_gene"))

# Filter to get only one2one orthologs in all strains
orthologs <- orthologs %>% dplyr::filter(cast_homo == "ortholog_one2one", c3h_homo == "ortholog_one2one", nod_homo == "ortholog_one2one", mouse_homo == "ortholog_one2one") %>% dplyr::select(-matches("homo")) %>% na.omit()

# Write to file
write.table(orthologs, "data/public/genes/mouse_strains_orthologs.csv", quote = F, sep = ",", row.names = F, col.names = T)


### ORTHOLOGS mm10, C57BL6NJ, C3H, AKR, LP ###
# Import downloaded data
orthologs2 <- data.table::fread("data/public/genes/mart_export2.txt", na.strings = "", col.names = c("C57BL6NJ_GENEID", "C57BL6NJ_TRANSCRIPT", "C57BL6J_GENEID", "C57BL6J_HOMO",  "AKR_GENEID", "AKR_HOMO", "C3H_GENEID", "C3H_HOMO",  "LP_GENEID", "LP_HOMO"))

# Filter to get only one2one orthologs in all strains
orthologs2 <- orthologs2 %>% dplyr::filter(C57BL6J_HOMO == "ortholog_one2one", C57BL6J_HOMO == "ortholog_one2one", AKR_HOMO == "ortholog_one2one", LP_HOMO == "ortholog_one2one") %>% dplyr::select(-matches("HOMO")) %>% na.omit()
colnames(orthologs2) <- sub("_GENEID", "", colnames(orthologs2))

# Write to file
write.table(orthologs2, "data/public/genes/mouse_strains_orthologs-yu2019.csv", quote = F, sep = ",", row.names = F, col.names = T)
