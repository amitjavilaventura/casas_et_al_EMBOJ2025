# R

# This script:
# - Takes a cluster annotation and looks which clusters are and are not present in all strains. 

### SET UP --------------------------------------

workdir <- "~/projects/mouse_pic_var_review/"
setwd(workdir)

# Load packages
library(dplyr)
library(data.table)
library(purrr)
library(stringr)
library(here)

# Defina annotation
anno = "lietal_clusters"
anno = "protrac_merged"

### IMPORT DATA ---------------------------------

# Import original clusters
if (anno == "lietal_clusters"){
  pics <- list.files(here::here("output/01-pirna_clusters/", anno), "mm10", full.names = T, recursive = T) %>%
    data.table::fread(.) %>% dplyr::transmute(piCid=sub("pi-Ip6k1,pi-Ip6k1.__piC117","pi-Ip6k1__piC117", V4), piCnum=sub(".*__", "", piCid))
} else if(anno == "protrac_merged"){
  pics <- data.table::fread(here::here("strains_analysis/output/01-pirna_clusters/protrac_clusters/protrac_merged//protrac_merged.3rmsk_filt.mm10.bed")) %>%
    dplyr::transmute(piCid = V4, piCnum = piCid)
}
# List orthologous list
pic_orth <- list.files(here::here("output/02-pirna_clusters_orthologs/", anno), "bed$", full.names = T, recursive = T) %>%
  purrr::discard(stringr::str_detect(., "filt.bed")) %>%
  purrr::set_names(sub(".*\\.","",sub(".bed","",basename(.)))) %>%
  purrr::imap(~data.table::fread(.x) %>% 
                dplyr::transmute(piCid_orth=sub("pi-Ip6k1,pi-Ip6k1.__piC117","pi-Ip6k1__piC117", V4), 
                                 piCid=sub("pi-Ip6k1,pi-Ip6k1.__piC117","pi-Ip6k1__piC117", V7),
                                 piCnum=sub(".*__", "", piCid),
                                 strain=.y)) 

# Create table for orthologs
pic_orth_table <- pic_orth %>%
  dplyr::bind_rows() %>%
  tidyr::pivot_wider(names_from = "strain", values_from = "piCnum")


### JOIN DATA -----------------------------------

# Join pics and pic_orhtologs
# Add T or F is the cluster is present in the corresponding strain
# Add T or F if the clusters is present in ALL the strains
pics_all_strains <- dplyr::left_join(pics, pic_orth_table) %>%
  dplyr::mutate(across(!matches("piC"), ~ifelse(is.na(.), F, T))) %>%
  dplyr::select(piCid, piCnum, piCid_orth, mm10, cast, c3h, "129", nod) %>% 
  dplyr::mutate(AllStrains = if_else(rowSums(across(!matches("piC")))<5, F, T))


### WRITE TO FILE ------------------------------

outdir = here::here("output/02-pirna_clusters_orthologs", anno, "")
write.table(pics_all_strains %>% dplyr::filter(AllStrains), paste0(outdir, anno, "_present_in_all_strains.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(pics_all_strains %>% dplyr::filter(!AllStrains), paste0(outdir, anno, "_notpresent_in_all_strains.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
