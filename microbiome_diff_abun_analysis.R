# Microbiome Analysis - Differential Abundance Analysis

# load libraries
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(phyloseq)
library(ALDEx2)
library(Maaslin2)
library(corncob)
library(purrr)
library(DESeq2)
library(apeglm)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name
meta$condition <- factor(meta$condition, levels = c("disease", "healthy")) # factor (but not for ALDEx2)


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

# filter low prevalence taxa (present in less than 10% of samples)
min_prevalence <- ceiling(0.10 * ncol(feat)) # minimum number of samples
feat_filtered <- feat[rowSums(feat > 0) >= min_prevalence, ] 

# make sure names are syntactically valid 
rownames(feat_filtered) <- make.names(rownames(feat_filtered))


###########################################################
#####   ALDEX2 - ANOVA-LIKE DIFFERENTIAL EXPRESSION   #####
###########################################################

all(colnames(feat) == rownames(meta)) # ensure sample names are the same
group <- as.character(meta$condition) # extract condition vector (needs to be a character vector)

# run ALDEx2 (with CLR transformation)
aldex_results <- aldex(reads = feat_filtered, # table of raw counts
                       conditions = group, # specifies group labels (needs to be in the same order as found in feat_filtered)
                       test = "t", # Welch's t and Wilcoxon rank-sum tests
                       effect = TRUE, # calculate effect sizes
                       denom = "all", # denominator for CLR transformation
                       verbose = TRUE)

# arrange differentially abundant features by adjusted Wilcoxon rank sum test
aldex_results %>% arrange(wi.eBH)


### MA plot (log-ratio abundance versus effect size (log-ratio difference between groups (CLR effect size)))
# red dots significant
aldex_plot_data <- aldex_results[aldex_results$rab.all > 0, ] # remove low abundance features for plotting
ggplot(aldex_plot_data, aes(x = rab.all, y = effect, color = wi.eBH < 0.5)) +
  geom_point(alpha = 0.7) + scale_x_log10() + theme_minimal() +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), name = "wi.eBH < 0.05") +
  labs(x = "Mean CLR-transformed abundance (across Monte Carlo instances)", y = "Effect size (log-ratio)", 
       title = "MA Plot (ALDEx2)")


### volcano plot (effect size versus -log10(FDR))
aldex_results$log10p <- -log10(aldex_results$wi.eBH + 1e-10)  # to avoid log(0)
ggplot(aldex_results, aes(x = effect, y = log10p, color = wi.eBH < 0.5)) + 
  geom_point(alpah = 0.7) + theme_minimal() +
  scale_color_manual(values = c("grey", "red"), name = "wi.eBH < 0.05") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  labs(x = "Effect size", y = "-log10(FDR)", 
       title = "Volcano plot (ALDEx2)")


### boxplots of top 5 significant features by effect size
top_features <- aldex_results %>%
  filter(wi.eBH < 0.5) %>% 
  slice_max(abs(effect), n = 5) %>%
  rownames()

# reshape data.frame to long format
feat_long <- feat_filtered[top_features, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = -feature, names_to = "sample_name", values_to = "count") %>%
  left_join(meta, by = "sample_name")

# boxplots of top differential features
ggplot(feat_long, aes(x = condition, y = count + 1, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ feature, scales = "free_y") + scale_y_log10() + theme_minimal() +
  labs(title = "Top differential features (ALDEx2)", y = "CLR-normalized abundance", x = NULL) +
  theme(legend.position = "none")


### heatmap of top 25 features
top_features_aldex <- aldex_results %>%
  arrange(wi.eBH) %>%
  slice_head(n = 25) %>%
  rownames() # get top 25 features by wi.EBH

heat_data <- log10(feat_filtered[top_features_aldex, ] + 1) # subset and transform count data
heat_data_scaled <- t(scale(t(heat_data))) # z-score scale each row to visualize patterns within features

# create annotation for samples
annotation_col <- data.frame(condition = meta$condition)
rownames(annotation_col) <- rownames(meta)

# plot the heatmap (log10(count+1) + z scaled)
pheatmap(heat_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = list(condition = c(disease = "firebrick",healthy = "steelblue")),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 25 ALDEx2 features")


#################################################################################
#####   MAASLIN2 - MICROBIOME MULTIVARIATE ASSOCIATION WITH LINEAR MODELS   #####
#################################################################################

# transpose feature table
feat_maaslin <- as.data.frame(t(feat_filtered))

# ensure sample names are the same
all(rownames(feat_maaslin) == rownames(meta))

# run MaAsLin2
maaslin_results <- Maaslin2(input_data = feat_maaslin,
                            input_metadata = meta,
                            output = "maaslin2_output", # path to the output directory
                            fixed_effects = "condition", # variables to test (or to adjust for)
                            normalization = "CLR", # CLR transformation (MaAsLin2 applies TSS normalization and adds a pseudocount internally)
                            transform = "NONE",
                            analysis_method = "LM", # linear model
                            plot_heatmap = TRUE, # generates heatmap of top associations
                            plot_scatter = TRUE) # generates scatterplots for individual significant associations

# load results
maaslin_all <- read_tsv("maaslin2_output/all_results.tsv", show_col_types = FALSE)

# arrange differentially abundant features by q-value
maaslin_all <- maaslin_all %>% arrange(qval)


### MA plot (log-ratio abundance versus effect size)
# calculate mean abundance per feature (same data used as input for MaAsLin2)
mean_abundance <- feat_maaslin %>%
  colMeans() %>%
  enframe(name = "feature", value = "mean_abundance")

# merge with MaAsLin2 results
ma_plot_data <- left_join(maaslin_all, mean_abundance, by = "feature")

# MA plot (log-ratio abundance versus effect size (regression coefficient (slope) from the linear model))
ggplot(ma_plot_data, aes(x = mean_abundance, y = coef, color = qval < 0.5)) +
  geom_point(alpha = 0.7) + scale_x_log10() + theme_minimal() +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), name = "qval < 0.05") +
  labs(x = "Mean CLR-transformed values", y = "Effect size (regression coefficient)", 
       title = "MA Plot (MaAsLin2)")

### volcano plot (effect size (regression coefficient) versus -log10(FDR))
maaslin_all$log10q <- -log10(maaslin_all$qval + 1e-10)
ggplot(maaslin_all, aes(x = coef, y = log10q, color = qval < 0.55)) +
  geom_point(alpha = 0.7) + theme_minimal() +
  scale_color_manual(values = c("grey", "red"), name = "qval < 0.05") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  labs(x = "Effect size (regression coefficient)", y = "-log10(FDR)", 
       title = "Volcano plot (MaAsLin2)")


### boxplots of top 5 significant features by effect size
top_features <- maaslin_all %>% 
  filter(qval < 0.5) %>% 
  slice_max(abs(coef), n = 5) %>% 
  pull(feature)

# reshape data.frame to long format (feat_maaslin is just feat_filtered transposed)
feat_long <- feat_filtered[top_features, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = -feature, names_to = "sample_name", values_to = "count") %>%
  left_join(meta, by = "sample_name")

# boxplots of top differential features
ggplot(feat_long, aes(x = condition, y = count + 1, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ feature, scales = "free_y") + scale_y_log10() + theme_minimal() +
  labs(title = "Top differential features (MaAsLin2)", y = "CLR-normalized abundance", x = NULL) +
  theme(legend.position = "none")


### heatmap of top 25 features
top_features_maaslin <- maaslin_all %>%
  arrange(qval) %>%
  slice_head(n = 25) %>%
  pull(feature) # get top 25 features by qval

heat_data <- log10(feat_filtered[top_features_maaslin, ] + 1) # subset and transform count data
heat_data_scaled <- t(scale(t(heat_data))) # z-score scale each row to visualize patterns within features

# create annotation for samples
annotation_col <- data.frame(condition = meta$condition)
rownames(annotation_col) <- rownames(meta)

# plot the heatmap (log10(count+1) + z scaled)
pheatmap(heat_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = list(condition = c(disease = "firebrick",healthy = "steelblue")),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 25 MaAsLin2 features")


#############################################################################################
#####   CORNCOB - COUNT REGRESSION FOR CORRELATED OBSERVATIONS WITH THE BETA-BINOMIAL   #####
#############################################################################################

# format feature table and metadata
feat_otu <- as.matrix(feat_filtered) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE)  # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create a phyloseq object
ps <- phyloseq(feat_otu, sampledata)


# run corncob (fits beta-binomial regression model - tests differential abundance and/or differential variability)
corncob_results <- differentialTest(formula = ~ condition, # specifies full model for differential abundance
                                    phi.formula = ~ condition, # set to condition to also test for differential variability (if dispersion varies by condition)
                                    formula_null = ~ 1, # add covariates to adjust for (add to both this and formula, only leave variable of interest out of null formula)
                                    phi.formula_null = ~ 1, # add covariates to adjust for with dispersion (add to both this and formula, only leave variable of interest out of null formula)
                                    data = ps, # phyloseq object
                                    test = "LRT", # statistical test used to compare models (LRT better for small sample sizes than Wald, but more computationally intensive)
                                    boot = TRUE,
                                    nboot = 1000,
                                    fdr_cutoff = 0.05) # cut-off for Benjamini-Hochberg

# name the models using the rownames of feat_filtered
names(corncob_results$all_models) <- rownames(feat_filtered)

# get feature names
features <- names(corncob_results$all_models)

# extract mu and phi coefficients, p-values and fdr from corncob results
coef_df <- purrr::map_dfr(features, function(feature) { # extract model results for each name in features
  model_summary <- corncob_results$all_models[[feature]] # model summary for current feature
  
  if (!inherits(model_summary, "summary.bbdml")) return(NULL) # skip the feature if the model didn't run properly
  
  coefs <- model_summary$coefficients # combined coefficient matrix from summary object
  
  # define the coefficient names
  mu_term <- "mu.conditionhealthy"
  phi_term <- "phi.conditionhealthy"
  
  # skip if either coefficient is missing
  if (!(mu_term %in% rownames(coefs)) || !(phi_term %in% rownames(coefs))) return(NULL)
  
  mu_estimate <- coefs[mu_term, "Estimate"] # extract mu
  phi_estimate <- coefs[phi_term, "Estimate"] # extract phi
  
  data.frame(feature = feature,
             mu_coef = mu_estimate,
             phi_coef = phi_estimate,
             pval = corncob_results$p[feature],
             fdr = corncob_results$p_fdr[feature])
})

# calculate mean counts per feature and join to coef_df
mean_abundance <- rowMeans(feat_filtered) %>%
  enframe(name = "feature", value = "mean_abundance")
plot_data <- left_join(coef_df, mean_abundance, by = "feature")

# add column for dominant effect (mu versus phi vs both)
plot_data <- plot_data %>%
  mutate(effect_type = case_when(abs(mu_coef) > abs(phi_coef) ~ "mu",
                                 abs(phi_coef) > abs(mu_coef) ~ "phi",
                                 TRUE ~ "both"))


### MA plot (mean abundance (log scale) versus effect size (mu coefficient))
ggplot(plot_data, aes(x = mean_abundance, y = mu_coef, color = fdr < 0.05)) +
  geom_point(alpha = 0.7) + scale_x_log10() + theme_minimal() +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), name = "FDR < 0.05") +
  labs(x = "Mean abundance (counts, log scale)", y = "Effect size (mu coefficient)",
       title = "MA Plot (corncob)")


### volcano plot (effect size (mu coefficient) versus -log10FDR)
plot_data <- plot_data %>% mutate(log10_fdr = -log10(fdr + 1e-10)) # add log10(FDR) to plot_data

ggplot(plot_data, aes(x = mu_coef, y = log10_fdr, color = fdr < 0.05)) +
  geom_point(alpha = 0.7) + theme_minimal() +
  scale_color_manual(values = c("grey", "red"), name = "FDR < 0.05") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  labs(x = "Abundance effect size (mu coefficient)", y = "-log10(FDR)", 
       title = "Volcano plot (corncob)")


### boxplots of top significant features (by effect size)
# select top 5 significant features by absolute mu coefficient
top_features <- plot_data %>%
  filter(fdr < 0.05) %>%
  slice_max(order_by = abs(mu_coef), n = 5) %>%
  pull(feature)

# subset count data by top_features and pivot to long format
top_counts <- feat_filtered[top_features, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = -feature,
                      names_to = "sample_name",
                      values_to = "count")

# merge table with metadata
plot_df <- left_join(top_counts, meta, by = "sample_name")

# boxplot with log scale counts
ggplot(plot_df, aes(x = condition, y = count + 1, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ feature, scales = "free_y") + scale_y_log10() + theme_minimal() +
  labs(title = "Top 5 differential features (corncob)", y = "Raw counts (log scale)", x = NULL) +
  theme(legend.position = "none")


### heatmap of top 25 features
top_features_corncob <- plot_data %>%
  arrange(fdr) %>%
  slice_head(n = 25) %>%
  pull(feature) # get top 25 features by fdr

heat_data <- log10(feat_filtered[top_features_corncob, ] + 1) # subset and transform count data
heat_data_scaled <- t(scale(t(heat_data))) # z-score scale each row to visualize patterns within features

# create annotation for samples
annotation_col <- data.frame(condition = meta$condition)
rownames(annotation_col) <- rownames(meta)

# plot the heatmap (log10(count+1) + z scaled)
pheatmap(heat_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = list(condition = c(disease = "firebrick",healthy = "steelblue")),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 25 corncob features")


########################################################
#####   DESEQ2 - DIFFERENTIAL ABUNDANCE ANALYSIS   #####
########################################################

# convert feat_filtered to integer matrix
count_matrix <- round(as.matrix(feat_filtered))
all(colnames(count_matrix) == rownames(meta)) # ensure sample names are the same

# create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta,
                              design = ~ condition)

# run DESeq2 analysis
dds <- DESeq(dds)

# extract differential abundance results (positive log2FoldChange = higher in healthy)
res <- results(dds, contrast = c("condition", "healthy", "disease"))

# shrink the log2fold change estimates using the apeglm method (noise reduction)
res <- lfcShrink(dds, coef = "condition_healthy_vs_disease", type = "apeglm")

# format and view results
# DESeq2 filters/excludes features with very low baseMean from FDR correction (they lack power to detect diff abund)
res_df <- as.data.frame(res) %>% rownames_to_column("feature") %>%
  arrange(padj)


### MA plot (mean of normalized counts versus log2 fold change)
plotMA(res, ylim = c(-5, 5), main = "MA plot (DESeq2)") # DESeq2 MA plot

# MA plot with ggplot to match others
ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = padj < 0.05)) +
  geom_point(alpha = 0.7) + scale_x_log10() + theme_minimal() + 
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"), name = "FDR < 0.05") +
  labs(x = "Mean of normalized counts", y = "Log2 fold change", 
       title = "MA plot (DESeq2)")


### volcano plot (log2 fold change versus -log10 adjusted p-value)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)) +
  geom_point(alpha = 0.7) + theme_minimal() +
  scale_color_manual(values = c("grey", "red"), name = "FDR < 0.05") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  labs(x = "Log2 fold change", y = "-log10 adjusted p-value",
       title = "Volcano plot (DESeq2)")


### boxplots of top significant features (by effect size)
# select top 5 significant features by log2fold change
top_features <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05) %>%
  slice_max(order_by = -abs(log2FoldChange), n = 5) %>%
  pull(feature)

# subset count data by top_features and pivot to long format
plot_data <- feat_filtered[top_features, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = -feature, 
                      names_to = "sample_name", 
                      values_to = "count") 

# merge table with metadata
plot_df <- left_join(plot_data, meta, by = "sample_name")


# boxplot with log scale counts
ggplot(plot_df, aes(x = condition, y = count + 1, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ feature, scales = "free_y") + scale_y_log10() + theme_minimal() +
  labs(title = "Top 5 differential features (DESeq2)", y = "Raw counts (log scale)", x = NULL) +
  theme(legend.position = "none")


### heatmap of top 25 features
top_features_deseq <- res_df %>%
  arrange(padj) %>%
  slice_head(n = 25) %>%
  pull(feature) # get top 25 features by padj

heat_data <- log10(feat_filtered[top_features_deseq, ] + 1) # subset and transform count data
heat_data_scaled <- t(scale(t(heat_data))) # z-score scale each row to visualize patterns within features

# create annotation for samples
annotation_col <- data.frame(condition = meta$condition)
rownames(annotation_col) <- rownames(meta)

# plot the heatmap (log10(count+1) + z scaled)
pheatmap(heat_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = list(condition = c(disease = "firebrick",healthy = "steelblue")),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 25 DESeq2 features")


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Edmonton
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.13               apeglm_1.30.0                 DESeq2_1.48.1                
# [4] SummarizedExperiment_1.38.1   Biobase_2.68.0                MatrixGenerics_1.20.0        
# [7] matrixStats_1.5.0             GenomicRanges_1.60.0          GenomeInfoDb_1.44.1          
# [10] IRanges_2.42.0                S4Vectors_0.46.0              BiocGenerics_0.54.0          
# [13] generics_0.1.4                corncob_0.4.2                 Maaslin2_1.22.0              
# [16] ALDEx2_1.40.0                 latticeExtra_0.6-30           zCompositions_1.5.0-5        
# [19] survival_3.8-3                truncnorm_1.0-9               MASS_7.3-65                  
# [22] phyloseq_1.52.0               doParallel_1.0.17             iterators_1.0.14             
# [25] foreach_1.5.2                 ParBayesianOptimization_1.2.6 Matrix_1.7-3                 
# [28] MLmetrics_1.1.3               pROC_1.19.0.1                 caret_7.0-1                  
# [31] lattice_0.22-7                xgboost_1.7.11.1              compositions_2.0-8           
# [34] lubridate_1.9.4               forcats_1.0.0                 stringr_1.5.1                
# [37] dplyr_1.1.4                   purrr_1.1.0                   readr_2.1.5                  
# [40] tidyr_1.3.1                   tibble_3.3.0                  tidyverse_2.0.0              
# [43] ggplot2_3.5.2                
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.5.0            hardhat_1.4.1            rpart_4.1.24             lifecycle_1.0.4         
# [5] rstatix_0.7.2            processx_3.8.6           globals_0.18.0           vroom_1.6.5             
# [9] DiceKriging_1.6.0        backports_1.5.0          magrittr_2.0.3           remotes_2.5.0           
# [13] pkgbuild_1.4.8           pbapply_1.7-4            bayesm_3.1-6             DBI_1.2.3               
# [17] RColorBrewer_1.1-3       ade4_1.7-23              abind_1.4-8              quadprog_1.5-8          
# [21] hash_2.2.6.3             nnet_7.3-20              tensorA_0.36.2.1         ipred_0.9-15            
# [25] lava_1.8.1               GenomeInfoDbData_1.2.14  listenv_0.9.1            vegan_2.7-1             
# [29] parallelly_1.45.1        permute_0.9-8            codetools_0.2-20         getopt_1.20.4           
# [33] DelayedArray_0.34.1      tidyselect_1.2.1         UCSC.utils_1.4.0         farver_2.1.2            
# [37] jsonlite_2.0.0           multtest_2.64.0          e1071_1.7-16             Formula_1.2-5           
# [41] bbmle_1.0.25.1           dbscan_1.2.2             tools_4.5.0              detectseparation_0.3    
# [45] Rcpp_1.1.0               glue_1.8.0               prodlim_2025.04.28       SparseArray_1.8.1       
# [49] mgcv_1.9-3               numDeriv_2016.8-1.1      withr_3.0.2              BiocManager_1.30.26     
# [53] rhdf5filters_1.20.0      callr_3.7.6              digest_0.6.37            timechange_0.3.0        
# [57] R6_2.6.1                 ROI.plugin.lpsolve_1.0-2 jpeg_0.1-11              dichromat_2.0-0.1       
# [61] utf8_1.2.6               pls_2.8-5                data.table_1.17.8        recipes_1.3.1           
# [65] robustbase_0.99-4-1      class_7.3-23             httr_1.4.7               S4Arrays_1.8.1          
# [69] som_0.3-5.2              ModelMetrics_1.2.2.2     pkgconfig_2.0.3          gtable_0.3.6            
# [73] timeDate_4041.110        registry_0.5-1           XVector_0.48.0           pcaPP_2.0-5             
# [77] carData_3.0-5            biomformat_1.36.0        zigg_0.0.2               scales_1.4.0            
# [81] logging_0.10-108         png_0.1-8                optparse_1.7.5           gower_1.0.2             
# [85] rstudioapi_0.17.1        ROI_1.0-1                tzdb_0.5.0               reshape2_1.4.4          
# [89] coda_0.19-4.1            checkmate_2.3.2          curl_6.4.0               nlme_3.1-168            
# [93] bdsmatrix_1.3-7          proxy_0.4-27             rhdf5_2.52.1             desc_1.4.3              
# [97] pillar_1.11.0            grid_4.5.0               vctrs_0.6.5              VGAM_1.1-13             
# [101] slam_0.1-55              ggpubr_0.6.1             car_3.1-3                chemometrics_1.4.4      
# [105] cluster_2.1.8.1          lpSolveAPI_5.5.2.0-17.14 locfit_1.5-9.12          mvtnorm_1.3-3           
# [109] cli_3.6.5                compiler_4.5.0           rlang_1.1.6              crayon_1.5.3            
# [113] future.apply_1.20.0      ggsignif_0.6.4           labeling_0.4.3           mclust_6.1.1            
# [117] interp_1.1-6             emdbook_1.3.14           ps_1.9.1                 plyr_1.8.9              
# [121] stringi_1.8.7            deldir_2.0-4             BiocParallel_1.42.1      Biostrings_2.76.0       
# [125] lars_1.3                 hms_1.1.3                bit64_4.6.0-1            future_1.67.0           
# [129] Rhdf5lib_1.30.0          trust_0.1-8              Rfast_2.1.5.1            igraph_2.1.4            
# [133] broom_1.0.9              RcppParallel_5.1.10      biglm_0.9-3              DEoptimR_1.1-4          
# [137] directlabels_2025.6.24   bit_4.6.0                ape_5.8-1   

