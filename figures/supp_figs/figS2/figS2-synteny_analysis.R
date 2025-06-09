# R

# This script:
# - Imports clusters from Li et al., 2013
# - LiftOvers them to rat, macaque and human.
# - Looks wether synteny is broken or not. 

# The idea is to answer part of the reviewer's comment:
## (3) For the differentially expressed loci highlighted in Figure 1, please (1) note whether the locus
##  produces transposon-silencing piRNAs, pre-pachytene piRNAs from the 3′ UTR of an mRNA, or pachytene piRNAs 
##  from one of the annotated pachytene piRNA genes; and (2) indicate which are evolutionarily conserved
##  (i.e., present at the syntenic location) across placental mammals, 
##  only among rodents, and which are mouse-specific (e.g., 14-qA3-284 and 14-qC1-1261).

WORKDIR="~/projects/github_repo_mouse_strains_paper/";
setwd(WORKDIR)

### Load packages
library(dplyr)
library(data.table)
library(ggplot2)

### Import data ---------------------------------

# Import mouse clusters from Li et al., 2013
# Import classes and merge the clusters with the classes.
pics <- data.table::fread("data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed") %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))
pic_class <- data.table::fread("data/public/lietal2013/lietal2013-piRNA_clusters.classes.csv")
pics <- dplyr::inner_join(pics, pic_class, by=c("V4"="piCid"))

# Import all clusters from Ozata et al., 2020
ozata2020 <- list.files("data/public/ozata2022/", "bed", full.names = T) %>%
  purrr::set_names(basename(.) %>% sub(".*clusters.", "", .) %>% sub(".bed$", "", .)) %>%
  purrr::map(~data.table::fread(.x, col.names = c("seqnames","start","end","id","width","strand")) %>%
               dplyr::arrange(seqnames, start) %>%
               dplyr::mutate(width=end-start))


### Synteny ------------------------------------

# Define liftover stuff
liftover="tools/liftover/liftOver"
mouse.rat="tools/liftover/mm10ToRn6.over.chain.gz"
mouse.mac="tools/liftover/mm10ToRheMac8.over.chain.gz"
mouse.hum="tools/liftover/mm10ToHg19.over.chain.gz"
mouse.mar="tools/liftover/mm10ToMarmoset.over.chain.gz"
mouse.opo="tools/liftover/mm10ToOpossum.over.chain.gz"
mouse.pla="tools/liftover/mm10ToPlatypus.over.chain.gz"

# Define outdir
outdir="figures/supp_figs/figS2/"
dir.create(outdir, F, T)

# Define liftover commands
lift_pics_rat = paste(liftover, "-minMatch=0.1 -bedPlus=6", "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed", mouse.rat, paste0(outdir, "/lietal_clusters.rat.bed"), paste0(outdir,"/lietal_clusters.rat.unmapped"))
lift_pics_mac = paste(liftover, "-minMatch=0.1 -bedPlus=6", "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed", mouse.mac, paste0(outdir, "/lietal_clusters.mac.bed"), paste0(outdir,"/lietal_clusters.mac.unmapped"))
lift_pics_hum = paste(liftover, "-minMatch=0.1 -bedPlus=6", "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed", mouse.hum, paste0(outdir, "/lietal_clusters.hum.bed"), paste0(outdir,"/lietal_clusters.hum.unmapped"))
lift_pics_mar = paste(liftover, "-minMatch=0.1 -bedPlus=6", "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed", mouse.mar, paste0(outdir, "/lietal_clusters.mar.bed"), paste0(outdir,"/lietal_clusters.mar.unmapped"))
lift_pics_opo = paste(liftover, "-minMatch=0.1 -bedPlus=6", "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed", mouse.opo, paste0(outdir, "/lietal_clusters.opo.bed"), paste0(outdir,"/lietal_clusters.opo.unmapped"))
lift_pics_pla = paste(liftover, "-minMatch=0.1 -bedPlus=6", "data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed", mouse.pla, paste0(outdir, "/lietal_clusters.pla.bed"), paste0(outdir,"/lietal_clusters.pla.unmapped"))


# Run liftover commands
system(lift_pics_rat)
system(lift_pics_mac)
system(lift_pics_hum)
system(lift_pics_mar)
system(lift_pics_opo)
system(lift_pics_pla)

# Import mouse clusters converted with liftover
mouse_pics.rat <- data.table::fread("figures/supp_figs/figS2/lietal_clusters.rat.bed") %>% dplyr::arrange(V1,V2) %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))
mouse_pics.mac <- data.table::fread("figures/supp_figs/figS2/lietal_clusters.mac.bed") %>% dplyr::arrange(V1,V2) %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))
mouse_pics.hum <- data.table::fread("figures/supp_figs/figS2/lietal_clusters.hum.bed") %>% dplyr::arrange(V1,V2) %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))
mouse_pics.mar <- data.table::fread("figures/supp_figs/figS2/lietal_clusters.mar.bed") %>% dplyr::arrange(V1,V2) %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))
mouse_pics.opo <- data.table::fread("figures/supp_figs/figS2/lietal_clusters.opo.bed") %>% dplyr::arrange(V1,V2) %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))
mouse_pics.pla <- data.table::fread("figures/supp_figs/figS2/lietal_clusters.pla.bed") %>% dplyr::arrange(V1,V2) %>% dplyr::mutate(V4 = sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", V4))

# Intersect clusters from Li et al. 2013 (lifted to other species) with clusters from Ozata et al., 2020
# Merge them into a single df, including the synteny, wether there is overlap or not, and the species
mouse_pics.syntenic <- 
  list(bedtoolsr::bt.intersect(mouse_pics.rat, ozata2020$rat, s=T, loj=T) %>% dplyr::transmute(id=V4,synteny=T,pi=ifelse(V7!=".",T,F),species="Rat"),
       bedtoolsr::bt.intersect(mouse_pics.mac, ozata2020$macaque, s=T, loj=T) %>% dplyr::transmute(id=V4,synteny=T,pi=ifelse(V7!=".",T,F),species="Macaque"), 
       bedtoolsr::bt.intersect(mouse_pics.hum, ozata2020$human, s=T, loj=T) %>% dplyr::transmute(id=V4,synteny=T,pi=ifelse(V7!=".",T,F),species="Human"), 
       bedtoolsr::bt.intersect(mouse_pics.mar, ozata2020$marmoset, s=T, loj=T) %>% dplyr::transmute(id=V4,synteny=T,pi=ifelse(V7!=".",T,F),species="Marmoset"),
       bedtoolsr::bt.intersect(mouse_pics.opo, ozata2020$opossum, s=T, loj=T) %>% dplyr::transmute(id=V4,synteny=T,pi=ifelse(V7!=".",T,F),species="Opossum"),
       bedtoolsr::bt.intersect(mouse_pics.pla, ozata2020$platypus, s=T, loj=T) %>% dplyr::transmute(id=V4,synteny=T,pi=ifelse(V7!=".",T,F),species="Platypus")) %>%
  dplyr::bind_rows()


# Define synteny T or F for rat, mac, hum
pics.synteny <- pics %>% 
  dplyr::transmute(id=V4, class) %>%
  dplyr::mutate(Rat = ifelse(id %in% mouse_pics.rat$V4, T, F),
                Macaque = ifelse(id %in% mouse_pics.mac$V4, T, F),
                Human = ifelse(id %in% mouse_pics.hum$V4, T, F),
                Marmoset = ifelse(id %in% mouse_pics.mar$V4, T, F),
                Opossum = ifelse(id %in% mouse_pics.opo$V4, T, F),
                Platypus = ifelse(id %in% mouse_pics.pla$V4, T, F)) %>%
  tidyr::pivot_longer(cols = !matches("id|class"), names_to = "species", values_to = "synteny") %>%
  dplyr::left_join(mouse_pics.syntenic) %>%
  dplyr::mutate(across(everything(), ~ifelse(is.na(.), F, .))) %>%
  unique() %>%
  dplyr::mutate(pirna_synteny = ifelse(synteny&pi, "Syntenic with piRNA", ifelse(synteny&!pi, "Syntenic without piRNA", "Not syntenic"))) %>%
  dplyr::mutate(class=factor(class, c("Pre-pachytene", "Hybrid", "Pachytene")),
                species=factor(species, rev(c("Rat", "Macaque", "Human", "Marmoset", "Opossum", "Platypus")), rev(c("Rat", "Macaque", "Human", "Marmosete", "Opossum", "Platypus"))),
                pirna_synteny=factor(pirna_synteny, c("Syntenic with piRNA", "Syntenic without piRNA", "Not syntenic")),
                id=sub("__.*","",id))

### Add expression in BL6 ----------------------

# Import expression data 
#  Use data aligned to mm10, so all clusters are present in the analysis
pics.expression <- data.table::fread(here::here("output/mm10_01-srna_process/fisher_DE_TEV/lietal_clusters/lietal_clusters.normcounts.tsv"))

# Do mean expression
pics.expr.mean <- pics.expression %>% 
  reshape2::melt(variable.name = "sample", value.name = "expr") %>%
  dplyr::mutate(Geneid=sub("pi-Ip6k1,pi-Ip6k1.__piC117", "pi-Ip6k1__piC117", Geneid),
                cond = sub("_rep.*", "", sample), 
                #tissue = sub("_.*","",sample)
                ) %>%
  dplyr::group_by(Geneid, cond) %>%
  dplyr::summarise(expr=mean(expr)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "cond", values_from = "expr") %>%
  dplyr::arrange(testis_BL6) %>%
  reshape2::melt(variable.name = "cond", value.name = "expr") %>%
  dplyr::mutate(tissue=sub("_.*", "", cond), 
                strain=sub(".*_","",cond)) %>%
  dplyr::inner_join(pic_class, by=c("Geneid"="piCid")) %>%
  dplyr::mutate(Geneid=sub("__piC.*", "", Geneid),
                Geneid=factor(Geneid, unique(Geneid)))

pics.expr.mean.filt <- pics.expr.mean #%>% dplyr::filter(stringr::str_detect(cond,"BL6"))

### Join expression and synteny

pics.synteny.expr <- dplyr::inner_join(pics.expr.mean.filt, 
                                       pics.synteny %>% dplyr::select(-class), by=c("Geneid"="id")) %>%
  dplyr::mutate(Geneid=factor(Geneid, unique(pics.expr.mean$Geneid)),
                class = factor(class, c("Pre-pachytene","Hybrid","Pachytene"))) %>%
  dplyr::mutate(tissue = ifelse(tissue=="testis", "Testis", "Spermatogonia"),
                condition = paste0(tissue, " (", strain, ")"),
                condition = factor(condition, rev(c("Testis (BL6)", "Testis (NOD)", "Testis (C3H)", "Testis (129)", "Spermatogonia (BL6)", "Spermatogonia (CAST)"))))

### Do plots ------------------------------------

# Do the synteny plot heatmap
pic.synteny.plot <- pics.synteny.expr %>% 
  ggplot(aes(Geneid, species, fill=pirna_synteny)) +
  geom_tile(color="black", linewidth=.1) +
  facet_grid(~class, scales = "free_x", space = "free_x") +
  ggmitji::theme_custom(legend="bottom", x.text.angle=90, x.text.hjust=1, x.text.vjust =.5,
                        title.face = "italic", axis.title.face = "plain", axis.title.size = 9) +
  theme(axis.text.x = element_text(size=3, margin = margin(t=-2)), axis.ticks = element_blank(),
        legend.title = element_blank(), axis.text.y = element_text(size=6))+
  scale_fill_manual(values = c("darkred","pink","gray")) +
  scale_y_discrete(expand = c(0,0)) +
  labs(y=NULL, x=NULL, title="Synteny of piRNA clusters from Li et al. (2013)",
       caption="For synteny information, we used liftover -minMatch=0.1; for piRNA information, 
       clusters were overlapped to clusters from corresponding species annotated by Ozata et al. (2020)")

# Do the expresion heatmap
pic.expr.plot.all <- pics.synteny.expr %>%
  ggplot(aes(Geneid, condition, fill=log10(expr+1))) +
  geom_tile(color="black", linewidth=.1) +
  facet_grid(~class, scales = "free_x", space = "free_x") +
  ggmitji::theme_custom(legend="bottom", x.text.angle=90, x.text.hjust=1, x.text.vjust =.5,
                        title.face = "italic", axis.title.face = "plain", axis.title.size = 9) +
  theme(axis.text.x = element_text(size=3, margin = margin(t=-2)), axis.ticks = element_blank(),
        legend.title = element_text(size=9), axis.text.y = element_text(size=6))+
  scale_fill_gradient2(low = "white", mid = "white", high = "darkblue", midpoint = 2, limits=c(0,7), oob=scales::squish) +
  guides(fill=guide_colourbar(title = "Log10. norm. counts +1", title.position = "bottom",
                              barwidth = unit(3, "in"), barheight = unit(.1, "in"),
                              frame.colour = "black", frame.linewidth = .5, frame.linetype = 1, 
                              ticks = T, ticks.colour = "black", ticks.linewidth = .5)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(y=NULL, x=NULL, title=NULL, subtitle=NULL)


  
# Join both heatmaps
layout="
B
A"

pics_synteny_expr_all.plot <- list(pic.expr.plot.all+ggmitji::rm_strips_x()+labs(y = "Expression")+theme(plot.margin = unit(c(-1,.25,.5,.25), "cm"), legend.title = element_text(size=7)),
                                   pic.synteny.plot+ggmitji::remove_x_axis()+labs(caption="", y = "Synteny")+theme(plot.margin = unit(c(.25,.25,-1,.25), "cm"))) %>%
  patchwork::wrap_plots(design = layout) +
  patchwork::plot_layout(design = layout, heights = c(70,30), guides = "collect") +
  patchwork::plot_annotation(theme=theme(legend.position = "bottom",legend.box = "vertical"),
                             caption="For synteny information, we used liftover -minMatch=0.1; 
                             for piRNA information, clusters were overlapped to clusters from corresponding species annotated by Ozata et al. (2020);
                             for expression, we mapped sRNA-seq data to mm10.") &
  theme(axis.title.y=element_text(size=7)) 


# Write plots to file
pdf("figures/supp_figs/figS2/figS2-lietal_synteny_expr_bl6.pdf", width = 8, height = 5)
pics_synteny_expr_all.plot
dev.off()

