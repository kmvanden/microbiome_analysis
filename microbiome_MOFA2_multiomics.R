### MOFA2 Workflow - Metagenomics and Metabolomics

# load libraries
library(ggplot2)
library(tidyverse)
library(readxl)
library(compositions)
library(MOFA2)
library(basilisk)
library(pheatmap)
library(umap)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# subset meta to samples present in metabolomics data
sample_key <- read_excel("metabo_batch.xlsx", sheet = "Sample Meta Data") # sample name key
sample_key <- sample_key %>% mutate(CLIENT_SAMPLE_ID = sub("KMV_", "SMS_", CLIENT_SAMPLE_ID)) # rename sample ids
meta <- read.table("metadata.txt", header = TRUE) # metadata (to create key for metabolomics data)
all(sample_key$CLIENT_SAMPLE_ID %in% meta$sample_id)

meta <- meta %>%
  left_join(sample_key %>% dplyr::select(CLIENT_SAMPLE_ID, PARENT_SAMPLE_NAME),
            by = c("sample_id" = "CLIENT_SAMPLE_ID")) %>% # add PARENT_SAMPLE_NAME to meta
  filter(!is.na(PARENT_SAMPLE_NAME)) # subset
rownames(meta) <- meta$sample_name


### feature table for metagenomics data
feat <- read.table("feature_table.txt", header = TRUE)
feat_t <- as.data.frame(t(feat)) # transpose feature table

# subset feature table for metagenomics data to the samples present in the metabolomics analysis
all(rownames(meta) %in% rownames(feat_t))
feat_sub <- feat_t[rownames(feat_t) %in% rownames(meta),]
all(rownames(meta) == rownames(feat_sub))

# filter species present in less than 10% of samples
min_prevalence <- 0.10
feat_filtered <- feat_sub[, colMeans(feat_sub > 0) >= min_prevalence] # subset the feature table to only include features present in at least 10% of samples
dim(feat_filtered)

# convert feature table to relative abundances
feat_rel_abund <- feat_filtered/rowSums(feat_filtered)

# add pseudocount and perform CLR transformation
feat_clr <- clr(feat_rel_abund + 1e-6)
feat <- as.data.frame(feat_clr)


### feature table for metabolomics data
metabo <- read_excel("metabo_batch.xlsx", sheet = "Batch-normalized Data") # batch-normalized feature table
all(metabo$PARENT_SAMPLE_NAME %in% meta$PARENT_SAMPLE_NAME)

# replace names in PARENT_SAMPLE_NAME column in metabo with names from sample_name (from meta)
metabo <- metabo %>%
  left_join(meta %>% dplyr::select(PARENT_SAMPLE_NAME, sample_name), by = "PARENT_SAMPLE_NAME") %>%
  mutate(PARENT_SAMPLE_NAME = sample_name) %>%
  dplyr::select(-sample_name) %>%
  slice(match(meta$sample_name, PARENT_SAMPLE_NAME)) %>% # match order of samples in meta
  column_to_rownames(var = "PARENT_SAMPLE_NAME") # move sample id column to rownames
all(rownames(metabo) == rownames(meta))

# replace CHEM_ID (colnames of metabo) with names from CHEMICAL_NAME (from metabo_key)
metabo_key <- read_excel("metabo_batch.xlsx", sheet = "Chemical Annotation")
all(colnames(metabo) == metabo_key$CHEM_ID)
name_map <- metabo_key %>% dplyr::select(CHEM_ID, CHEMICAL_NAME) %>% deframe() 
colnames(metabo) <- name_map[colnames(metabo)] %>% ifelse(is.na(.), colnames(metabo), .)
all(colnames(metabo) == metabo_key$CHEMICAL_NAME)
colnames(metabo) <- make.names(colnames(metabo)) # make sure names are syntactically valid 


### missingness (zeros listed as NAs in data.frame)
missing_per_feature <- colMeans(is.na(metabo)) * 100 # missingness per feature
summary(missing_per_feature)
hist(missing_per_feature, breaks = 50, main = "Percent missing per metabolite", xlab = "Percent missing")

# remove features with greater than 30% missing values
metabo_filtered <- metabo[, missing_per_feature <= 30]

# log2 transformation (with pseudocount) for non-missing values
metabo_log <- metabo_filtered %>%
  mutate(across(everything(), ~ ifelse(is.na(.), NA, log2(. + 1e-6))))


### samples have to be in the same order for MOFA2
all(rownames(feat) == rownames(metabo_log))

# construct view list (features as rows and samples as columns)
data_list <- list(metagenomics = as.matrix(t(feat)), 
                  metabolomics = as.matrix(t(metabo_log)))
all(colnames(data_list$metagenomics) == colnames(data_list$metabolomics))


##############################################
###   MOFA2 - MULTI-OMIC FACTOR ANALYSIS   ###
##############################################

### create the MOFA object 
mofa_obj <- create_mofa(data_list)


### define options
data_opts <- get_default_data_options(mofa_obj) # data options
data_opts$scale_views <- TRUE # scale metagenomics and metabolomics data to unit variance

model_opts <- get_default_model_options(mofa_obj) # model options
model_opts$num_factors <- 10 # number of latent factors to learn

train_opts <- get_default_training_options(mofa_obj) # training options
train_opts$convergence_mode <- "medium"
train_opts$seed <- 1234 # set seed

# prepare the MOFA object with the defined options
mofa_obj <- prepare_mofa(object = mofa_obj,
                         data_options = data_opts,
                         model_options = model_opts,
                         training_options = train_opts)

# run the model
mofa_trained <- run_mofa(mofa_obj, use_basilisk = TRUE)

# save model
# save_model(mofa_trained, file = "mofa_trained_model.hdf5")
# mofa_trained <- load_model("mofa_trained_model.hdf5")


### elbow plot
# variance explained per factor summed across views
var_exp <- calculate_variance_explained(mofa_trained)
var_exp_per_factor <- rowSums(var_exp$r2_per_factor$group1)

elbow_df <- as.data.frame(var_exp$r2_per_factor$group1)
elbow_df <- elbow_df %>% mutate(both = metagenomics + metabolomics)
elbow_df$factor <- seq_len(nrow(elbow_df))

elbow_long <- elbow_df %>%
  pivot_longer(cols = c("metagenomics", "metabolomics", "both"),
               names_to = "view",
               values_to = "var_exp")

ggplot(elbow_long, aes(x = factor, y = var_exp, color = view)) +
  scale_color_manual(values = c(metagenomics = "steelblue", metabolomics = "indianred3", both = "plum4")) +
  geom_point() + geom_line() + theme_minimal() +
  scale_x_continuous(breaks = elbow_df$factor) +
  labs(title = "Elbow plot (variance explained per latent factor)",
       x = "Latent factor number", y = "Variance explained (%)")


### plot variance explained by each latent factor for each view
plot_variance_explained(mofa_trained)


### add metadata for visualization
all(samples_names(mofa_trained)$group1 == rownames(meta)) # ensure sample names match
meta_mofa <- meta %>% rename(sample = sample_name)
samples_metadata(mofa_trained) <- meta_mofa


### plot latent factor values 
plot_factors(mofa_trained, factors = c(1:2), color_by = "condition")
plot_factors(mofa_trained, factors = c(1:5), color_by = "condition")


### t-test and barplot for latent factors 
all_factor_scores <- as.data.frame(get_factors(mofa_trained)$group1) # latent factor scores

p_value <- sapply(1:ncol(all_factor_scores), function(i) {
  df <- data.frame(factor = all_factor_scores[, i], condition = meta$condition)
  test <- t.test(factor ~ condition, data = df)
  p_value = test$p.value
}) # get p-values

p_value_df = data.frame(factor = colnames(all_factor_scores), p_value = p_value)
p_value_df <- p_value_df %>% arrange(p_value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
p_value_df

all(rownames(all_factor_scores) == rownames(meta)) # check sample order between latent factor scores and metadata
all_factor_scores$sample <- rownames(all_factor_scores)
all_factor_scores$condition = meta$condition

# barplot by latent factor
ggplot(all_factor_scores, aes(x = condition, y = Factor1, color = condition)) +
  geom_boxplot() + geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title = "Factor scores by condition", x = "Condition", y = "Factor score")

# pivot to long format
factor_long <- all_factor_scores %>%
  pivot_longer(cols = starts_with("Factor"), names_to = "Factor", values_to = "Score")

# barplot for all latent factors
ggplot(factor_long, aes(x = condition, y = Score, color = condition)) +
  geom_boxplot() + geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Factor, nrow = 2, ncol = 5) + theme_minimal() +
  labs(title = "Latent factor scores by condition",
       x = "Condition", y = "Factor score")


### loadings per latent factor
# metabolomics
plot_weights(mofa_trained, view = "metabolomics", factor = 1, nfeatures = 20, text_size = 2)

weights <- get_weights(mofa_trained, view = "all", factor = "all")
weights_metabo <- as.data.frame(weights$metabolomics)
weights_metabo$feature <- rownames(weights_metabo)
top_metabo_feat_factor1 <- weights_metabo %>% 
  arrange(desc(abs(Factor1))) %>%
  slice_head(n = 20) %>%
  pull(feature)

heatmap_data_metabo <- metabo_log[, top_metabo_feat_factor1] # subset to top features
heatmap_data_metabo_scaled <- t(scale(heatmap_data_metabo))
pheatmap(heatmap_data_metabo_scaled, annotation_col = meta[, "condition", drop = FALSE],
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)

heatmap_res <- pheatmap(heatmap_data_metabo_scaled, annotation_col = meta[, "condition", drop = FALSE],
                        cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
feature_clusters <- as.data.frame(cutree(heatmap_res$tree_row, k = 4))
feature_clusters <- feature_clusters %>% arrange(cutree(heatmap_res$tree_row, k = 4))


# metagenomics
plot_weights(mofa_trained, view = "metagenomics", factor = 1, nfeatures = 20, text_size = 2)

weights_metagen <- as.data.frame(weights$metagenomics)
weights_metagen$feature <- rownames(weights_metagen)
top_metagen_feat_factor1 <- weights_metagen %>% 
  arrange(desc(abs(Factor1))) %>%
  slice_head(n = 20) %>%
  pull(feature)

heatmap_data_metagen <- feat[, top_metagen_feat_factor1] # subset to top features
heatmap_data_metagen_scaled <- t(scale(heatmap_data_metagen))
pheatmap(heatmap_data_metagen_scaled, annotation_col = meta[, "condition", drop = FALSE],
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)


### latent factor correlation
factor_corr <- cor(all_factor_scores %>% select(starts_with("Factor")), use = "pairwise.complete.obs")
pheatmap(factor_corr, main = "Correlation between latent factors",
         display_numbers = TRUE, 
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(100))


### cluster samples based on latent factors

### hierarchical clustering
dist_matrix <- dist(all_factor_scores %>% select(starts_with("Factor")))
hc <- hclust(dist_matrix)
plot(hc)

hc_clusters <- cutree(hc, k = 5) # cut dendrogram into 4 clusters
meta_mofa$hc_cluster <- as.factor(hc_clusters) # add hierarchical cluster labels
rect.hclust(hc, k = 5, border = "red")


### k-means clustering
factor_scores <- all_factor_scores %>% select(starts_with("Factor")) # get latent factor scores
set.seed(1234) # set seed
wss <- sapply(1:10, function(k) {
  kmeans(factor_scores, centers = k, nstart = 25)$tot.withinss
})
diff(wss) # drops between k clusters

# elbow plot
plot(1:10, wss, type = "b", pch = 19,
     xlab = "Number of clusters", ylab = "Total within-cluster sum of squares",
     main = "Elbow method")

# cluster with selected k value
set.seed(123)
kmeans_res <- kmeans(factor_scores, centers = 2, nstart = 25)

meta_mofa$kmeans_cluster <- as.factor(kmeans_res$cluster) # add k-means cluster labels


### UMAP plot
umap_res <- umap(factor_scores)

# dataframe for plotting
umap_df <- data.frame(UMAP1 = umap_res$layout[,1],
                      UMAP2 = umap_res$layout[,2],
                      condition = meta_mofa$condition,
                      k_cluster = meta_mofa$kmeans_cluster,
                      h_cluster = meta_mofa$hc_cluster)

# plot UMAP colored by condition
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) + theme_minimal() +
  labs(title = "UMAP of latent factor scores (colored by condition)")


# plot UMAP colored by hierarchical clusters
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = h_cluster)) +
  geom_point(size = 2, alpha = 0.8) + theme_minimal() +
  labs(title = "UMAP of latent factor scores (colored by hierarchical cluster)")

# correlation of hierarchical clusters with condition
table(meta_mofa$hc_cluster, meta_mofa$condition)
fisher.test(table(meta_mofa$hc_cluster, meta_mofa$condition)) # fisher's exact test

# plot how hierarchical clusters map to condition
ggplot(meta_mofa, aes(x = hc_cluster, fill = condition)) +
  geom_bar(position = "fill") + theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of condition by hierarchical cluster",
       x = "Hierarchical cluster", y = "Proportion", fill = "Condition")


# plot UMAP colored by k-means clusters
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = k_cluster)) +
  geom_point(size = 2, alpha = 0.8) + theme_minimal() +
  labs(title = "UMAP of latent factor scores (colored by k-means cluster)")

# correlation of k-means clusters with condition
table(meta_mofa$kmeans_cluster, meta_mofa$condition)
fisher.test(table(meta_mofa$kmeans_cluster, meta_mofa$condition)) # fisher's exact test

# plot how k-means clusters map to condition
ggplot(meta_mofa, aes(x = kmeans_cluster, fill = condition)) +
  geom_bar(position = "fill") + theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of condition by k-means cluster",
       x = "K-means cluster", y = "Proportion", fill = "Condition")


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6.1
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] umap_0.2.10.0      pheatmap_1.0.13    basilisk_1.20.0    reticulate_1.43.0  MOFA2_1.18.0      
# [6] compositions_2.0-8 readxl_1.4.5       lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [11] dplyr_1.1.4        purrr_1.1.0        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0      
# [16] tidyverse_2.0.0    ggplot2_3.5.2     
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6          dir.expiry_1.16.0     tensorA_0.36.2.1      ggrepel_0.9.6        
# [5] corrplot_0.95         rhdf5_2.52.1          lattice_0.22-7        tzdb_0.5.0           
# [9] rhdf5filters_1.20.0   vctrs_0.6.5           tools_4.5.0           generics_0.1.4       
# [13] stats4_4.5.0          parallel_4.5.0        DEoptimR_1.1-4        pkgconfig_2.0.3      
# [17] Matrix_1.7-3          RColorBrewer_1.1-3    S4Vectors_0.46.0      lifecycle_1.0.4      
# [21] compiler_4.5.0        farver_2.1.2          pillar_1.11.0         crayon_1.5.3         
# [25] MASS_7.3-65           openssl_2.3.3         uwot_0.2.3            DelayedArray_0.34.1  
# [29] abind_1.4-8           RSpectra_0.16-2       robustbase_0.99-4-1   tidyselect_1.2.1     
# [33] Rtsne_0.17            stringi_1.8.7         reshape2_1.4.4        cowplot_1.2.0        
# [37] grid_4.5.0            colorspace_2.1-1      cli_3.6.5             SparseArray_1.8.1    
# [41] magrittr_2.0.3        S4Arrays_1.8.1        h5mread_1.0.1         dichromat_2.0-0.1    
# [45] withr_3.0.2           filelock_1.0.3        scales_1.4.0          timechange_0.3.0     
# [49] XVector_0.48.0        matrixStats_1.5.0     cellranger_1.1.0      askpass_1.2.1        
# [53] png_0.1-8             hms_1.1.3             HDF5Array_1.36.0      IRanges_2.42.0       
# [57] basilisk.utils_1.20.0 rlang_1.1.6           Rcpp_1.1.0            glue_1.8.0           
# [61] bayesm_3.1-6          BiocGenerics_0.54.0   rstudioapi_0.17.1     jsonlite_2.0.0       
# [65] Rhdf5lib_1.30.0       R6_2.6.1              plyr_1.8.9            MatrixGenerics_1.20.0

