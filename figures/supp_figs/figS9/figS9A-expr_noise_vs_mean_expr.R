
### GET COEFFICIENT OF VARIATION

icr_info <- data.table::fread("data/private/allSampleIdentifiers_mouseSmallRNA_private_enriched.txt") %>% 
  dplyr::filter(Strain=="ICR") %>%
  dplyr::select(SampleID,MouseID) 


# Get cluster names and cluster number
pics_names <- data.table::fread("data/public/lietal2013/lietal2013-piRNA_clusters.merged.mm10.bed") %>%
  dplyr::select("piCID"=V4) %>%
  tidyr::separate("piCID", c("ClusterName","piCID"),"__")

# Import raw counts from ICR mice in lietal clusters
icr_pic_raw <- data.table::fread("data/private/ICR-lietal_clusters.rawcounts.tsv")
colnames(icr_pic_raw) <- sub("V1","piCID",colnames(icr_pic_raw))
icr_pic_raw <- icr_pic_raw %>% dplyr::inner_join(pics_names[,c("piCID","ClusterName")]) %>%
  dplyr::select(piCID,ClusterName, c(icr_info$SampleID)) %>% 
  dplyr::mutate(Geneid=paste0(ClusterName,"__",piCID) %>% sub("\\(pi-Noct\\)","",.)) %>%
  dplyr::select(-piCID,-ClusterName)

# Normalization and variance stabilization transformation
countdata.icr <- icr_pic_raw %>% tibble::column_to_rownames("Geneid")
coldata.icr   <- data.frame(sample=colnames(countdata.icr))
dds.icr <- DESeq2::DESeqDataSetFromMatrix(countdata.icr,coldata.icr,~1)
dds.icr <- DESeq2::DESeq(dds.icr)
vst.icr <- DESeq2::vst(dds.icr, blind = T, nsub = nrow(dds.icr))
icr_pic_expr <- vst.icr@assays@data@listData %>% as.data.frame() %>% tibble::rownames_to_column("Geneid")

# Compute sd and mean of expression
# Also compute CV (sd/mean)
icr_pic_cv <- icr_pic_expr %>% 
  dplyr::transmute(Geneid, 
                   mean_expr = rowMeans(across(matches("sample"))), 
                   sd_expr   = matrixStats::rowSds(as.matrix(across(matches("sample"))))) %>%
  dplyr::mutate(cv = sd_expr/mean_expr, cv_pct=cv*100)

# Arrange and format the dataframe with the expression and the CV
icr_pic_cv <- icr_pic_cv %>% dplyr::arrange(desc(cv))
icr_pic_expr <- icr_pic_expr %>% dplyr::select(Geneid, matches("sample")) 
icr_pic_expr <- icr_pic_expr[match(icr_pic_cv$Geneid, icr_pic_expr$Geneid),]

# Fit a curve with non-linear regression (loess) --> cv (y) in function of mean expression (x)
# and find the predicted cv
loess_fit <- loess(cv ~ mean_expr, data = icr_pic_cv)
predicted_cv <- predict(loess_fit)

# same but with cv squared
loess_fit_sqr <- loess(cv^2 ~ mean_expr, data = icr_pic_cv)
predicted_cv_sqr <- predict(loess_fit_sqr)

# Add predicted cv and residuals to df 
icr_pic_cv.res <- icr_pic_cv %>% dplyr::mutate(predicted_cv = predicted_cv, residuals = cv - predicted_cv)
icr_pic_cv.res <- icr_pic_cv.res %>% dplyr::mutate(predicted_cv_sqr = predicted_cv_sqr, residuals_sqr = cv^2 - predicted_cv_sqr)

# Identify outliers (residuals greater than the 95th percentile)
outlier_threshold <- quantile(icr_pic_cv.res$residuals, 0.95)
outliers <- icr_pic_cv.res %>% dplyr::filter(residuals > outlier_threshold)
outlier_threshold_sqr <- quantile(icr_pic_cv.res$residuals_sqr, 0.95)
outliers_sqr <- icr_pic_cv.res %>% dplyr::filter(residuals_sqr > outlier_threshold_sqr)

#print(outliers$Geneid)

# Join outliers to plot 
icr_pic_cv.res <- icr_pic_cv.res %>% dplyr::mutate(outlier=ifelse(Geneid %in% outliers$Geneid, T, F))
icr_pic_cv.res <- icr_pic_cv.res %>% dplyr::mutate(label=ifelse(outlier, sub("__.*","",Geneid), NA)) 

icr_pic_cv.res <- icr_pic_cv.res %>% dplyr::mutate(outlier_sqr=ifelse(Geneid %in% outliers_sqr$Geneid, T, F))
icr_pic_cv.res <- icr_pic_cv.res %>% dplyr::mutate(label_sqr=ifelse(outlier_sqr, sub("__.*","",Geneid), NA))

# Do plot with highlighted outliers
snp_mean_ouliers_sqr_plot_line <-
icr_pic_cv.res %>% 
  dplyr::arrange(desc(outlier_sqr), desc(predicted_cv_sqr)) %>%
  dplyr::mutate(label_sqr = ifelse(!is.na(label_sqr), sub(".*__","",Geneid), label_sqr)) %>%
  ggplot(aes(x=(mean_expr),y=cv^2,color=outlier_sqr)) + 
  geom_point(alpha=.7) + 
  scale_color_manual(values = c("black","darkred")) + 
  geom_smooth(method = "loess", color="green", fill="gray", se=F) +
  ggrepel::geom_text_repel(aes(label=label_sqr), size=2.5) + 
  ggmitji::theme_clean(clean.y.axis = F, legend = "none") +
  labs(title="CV^2 vs mean expression", x="Mean expression", y="Coefficient of variation squared")

dir.create("figures/supp_figs/figS9/")
ggsave(plot=snp_mean_ouliers_sqr_plot_line, filename = "figures/supp_figs/figS9/figS9A-expr_noise_vs_mean_expr.pdf", width = 7, height = 7, units = "in")
