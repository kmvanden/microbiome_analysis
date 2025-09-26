### sPLSDA Workflow - Metabolomics

# load libraries
library(ggplot2)
library(tidyverse)
library(readxl)
library(impute)
library(mixOmics)
library(caret)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
### metadata
meta <- read.table("metadata.txt", header = TRUE)

# subset meta to samples present in metabolomics data
sample_key <- read_excel("metabo_batch.xlsx", sheet = "Sample Meta Data") # sample name key
sample_key <- sample_key %>% mutate(CLIENT_SAMPLE_ID = sub("KMV_", "SMS_", CLIENT_SAMPLE_ID)) # rename sample ids
all(sample_key$CLIENT_SAMPLE_ID %in% meta$sample_id)
meta <- meta %>%
  left_join(sample_key %>% dplyr::select(CLIENT_SAMPLE_ID, PARENT_SAMPLE_NAME),
            by = c("sample_id" = "CLIENT_SAMPLE_ID")) %>% # add PARENT_SAMPLE_NAME to meta
  filter(!is.na(PARENT_SAMPLE_NAME)) # subset

# convert condition column into a factor
rownames(meta) <- meta$sample_name
meta$condition <- as.factor(meta$condition) # needs to be a factor for sPLS-DA
table(meta$condition)


### feature table
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

### feature filtering, zero imputation and log transformation 
# assess missingness (zeros listed as NAs in data.frame)
missing_per_feature <- colMeans(is.na(metabo)) * 100 # missingness per feature
summary(missing_per_feature)
hist(missing_per_feature, breaks = 50, main = "Percent missing per metabolite", xlab = "Percent missing")

# remove features with greater than 30% missing values
metabo_filtered <- metabo[, missing_per_feature <= 30]

# impute remaining zeros with kNN
metabo_imputed <- impute.knn(as.matrix(t(metabo_filtered)))$data # impute.knn expects features as rows
metabo_imputed <- t(metabo_imputed)

# log2 transformation (with pseudocount)
metabo <- log2(metabo_imputed + 1e-6)


##################################################################################
########   MIXOMICS - SPARSE PARTIAL LEAST SQUARES DISCRIMINANT ANALYSIS   #######
##################################################################################

### evaluate performance using all features
# distance not specified = all distances included
metabo_all_feat <- mixOmics::plsda(X = metabo, Y = meta$condition, ncomp = 5)
perf_metabo_all_features <- perf(metabo_all_feat,
                                 validation = "Mfold", 
                                 folds = 5, 
                                 nrepeat = 50, 
                                 progressBar = TRUE)

plot(perf_metabo_all_features) # error rate versus number of components
perf_metabo_all_features$choice.ncomp # optimal number of components
perf_metabo_all_features$error.rate # error rate for each component


### feature selection using sparse method
list_keepX <- c(5, 10, 15, 20, 25, 50, 75, 100) # list of number of features to keep per component

### tune (with repeated cross-validation) number of features (keepX) per component 
tune_metabo_res <- tune.splsda(X = metabo,
                               Y = meta$condition,
                               test.keepX = list_keepX,
                               ncomp = 5, 
                               validation = "Mfold", # k-fold cross-validation (stratified)
                               folds = 5, # 5 folds for cross-validation
                               nrepeat = 50, # number of repeats of cross-validation
                               dist = "max.dist", # distance metric
                               measure = "overall", # classes are balanced
                               progressBar = TRUE)

plot(tune_metabo_res, sd = TRUE) # visualize performance across keepX values
tune_metabo_res$choice.ncomp # optimal number of components
tune_metabo_res$error.rate # error rates for keepX used
tune_metabo_res$choice.keepX # optimal number of features


### function to run tune.splsda for all three distance metrics (from microbiome_sPLSDA_metagen)
# list of number of features to keep per component
list_keepX <- c(5, 10, 15, 20, 25, 50, 75, 100)

### tune (with repeated cross-validation) number of features (keepX) per component for all three distance metrics
tuned_metabo_results <- tune_all_distances(X = metabo,
                                           Y = meta$condition,
                                           list.keepX = list_keepX,
                                           ncomp = 5,
                                           folds = 5,
                                           nrepeat = 50,
                                           measure = "overall")

lapply(tuned_metabo_results, function(dist) dist$error.rate) # overall error rate for each distance metric (for each keepX value)
lapply(tuned_metabo_results, function(dist) dist$choice.keepX) # optimal number of features
lapply(tuned_metabo_results, function(dist) dist$choice.ncomp) # optimal number of components
sapply(tuned_metabo_results, function(res) min(res$error.rate)) # minimum overall error rate (for each distance metric)

# visualize performance across keepX values
plot(tuned_metabo_results$max.dist, sd = TRUE) # max.dist
plot(tuned_metabo_results$centroids.dist, sd = TRUE) # centroids.dist 
plot(tuned_metabo_results$mahalanobis.dist, sd = TRUE) # mahalanobis.dist


### final model using optimal number of features and 1 component for model performance metrics
final_model_metabo <- mixOmics::splsda(X = metabo, Y = meta$condition, ncomp = 1, keepX = tuned_metabo_results$max.dist$choice.keepX)
final_metabo_res <- perf(final_model_metabo,
                         validation = "Mfold", 
                         folds = 5, 
                         nrepeat = 50, 
                         dist = "max.dist", 
                         measure = "overall",
                         progressBar = TRUE)

final_metabo_res$error.rate # error rate for final model
final_metabo_res$error.rate.class # error rate per class for final model


### ROC curve (evaluation of model performance)
auroc(final_model_metabo, roc.comp = 1)


### generate confusion matrix
predictions_df <- as.data.frame(final_metabo_res$class)  # class predictions
# function to get majority vote from the 50 repeats
get_majority_vote <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1] 
}
avg_predictions <- apply(predictions_df, 1, get_majority_vote) # get majority prediction
confusionMatrix(as.factor(avg_predictions), as.factor(meta$condition))


### plot frequency of selection for top selected features for component 1
final_model_features <- final_metabo_res$features$stable$comp1 # features by selection frequency (component 1)
final_model_features_top <- sort(head(final_model_features, 10), decreasing = FALSE) # top 5 features by selection frequency
features_df <- data.frame(feature = names(final_model_features_top),
                          frequency = as.numeric(final_model_features_top)) # create data.frame for ggplot2
features_df$feature <- factor(features_df$feature, levels = features_df$feature) # order features for plot

# ggplot of frequency of feature selection
ggplot(features_df, aes(x = feature, y = frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Feature stability (component 1)", x = "Selected features", y = "Frequency of selection") +
  theme_minimal(base_size = 12) + theme(axis.text.y = element_text(size = 10))


### loading weights
selectVar(final_model_metabo)

# plot loading weights (contribution of features to component)
plotLoadings(final_model_metabo, 
             comp = 1, method = "median", contrib = "max",
             legend.color = c("steelblue", "indianred3"), legend = FALSE,
             title = "Loadings plot")


### final model using optimal number of features and 2 components for plotting
plot_model_metabo <- mixOmics::splsda(X = metabo, Y = meta$condition, ncomp = 2, keepX = tuned_metabo_results$max.dist$choice.keepX)
plot_metabo_res <- perf(plot_model_metabo,
                        validation = "Mfold", 
                        folds = 5, 
                        nrepeat = 50, 
                        dist = "max.dist", 
                        measure = "overall",
                        progressBar = TRUE)

plot_metabo_res$error.rate # error rate for plot model
plot_metabo_res$error.rate.class # error rate per class for plot model


### sample plot (visualization of sample distribution in component space)
plotIndiv(plot_model_metabo, comp = c(1,2),
          ind.names = FALSE, legend = TRUE,
          X.label = "Component 1", Y.label = "Component 2",
          star = TRUE, ellipse = TRUE, ellipse.level = 0.95,
          col = c("steelblue", "indianred3"),
          cex = 3, pch = c(1,1), point.lwd = 0.5,
          title = "Sample plot")


### feature plot/correlation circle plot (visualization of features in the component space)
plotVar(plot_model_metabo, comp = c(1,2),
        var.names = TRUE, legend = TRUE, style = "graphics",
        col = "steelblue", pch = 1, cex = 0.5,
        cutoff = 0.30, # correlation coefficient cut-off
        title = "Variable plot")


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
#   [1] caret_7.0-1     mixOmics_6.32.0 lattice_0.22-7  MASS_7.3-65     impute_1.82.0  
# [6] readxl_1.4.5    lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
# [11] purrr_1.1.0     readr_2.1.5     tidyr_1.3.1     tibble_3.3.0    tidyverse_2.0.0
# [16] ggplot2_3.5.2  
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1     timeDate_4041.110    farver_2.1.2         pROC_1.19.0.1       
# [5] digest_0.6.37        rpart_4.1.24         timechange_0.3.0     lifecycle_1.0.4     
# [9] survival_3.8-3       magrittr_2.0.3       compiler_4.5.0       rlang_1.1.6         
# [13] tools_4.5.0          igraph_2.1.4         data.table_1.17.8    labeling_0.4.3      
# [17] rARPACK_0.11-0       plyr_1.8.9           RColorBrewer_1.1-3   BiocParallel_1.42.1 
# [21] withr_3.0.2          nnet_7.3-20          grid_4.5.0           stats4_4.5.0        
# [25] e1071_1.7-16         future_1.67.0        globals_0.18.0       scales_1.4.0        
# [29] iterators_1.0.14     dichromat_2.0-0.1    cli_3.6.5            ellipse_0.5.0       
# [33] generics_0.1.4       rstudioapi_0.17.1    future.apply_1.20.0  RSpectra_0.16-2     
# [37] reshape2_1.4.4       tzdb_0.5.0           proxy_0.4-27         splines_4.5.0       
# [41] parallel_4.5.0       cellranger_1.1.0     matrixStats_1.5.0    vctrs_0.6.5         
# [45] hardhat_1.4.1        Matrix_1.7-3         hms_1.1.3            ggrepel_0.9.6       
# [49] listenv_0.9.1        foreach_1.5.2        gower_1.0.2          pls_2.8-5           
# [53] recipes_1.3.1        glue_1.8.0           parallelly_1.45.1    codetools_0.2-20    
# [57] stringi_1.8.7        gtable_0.3.6         pillar_1.11.0        ipred_0.9-15        
# [61] lava_1.8.1           R6_2.6.1             corpcor_1.6.10       class_7.3-23        
# [65] Rcpp_1.1.0           gridExtra_2.3        nlme_3.1-168         prodlim_2025.04.28  
# [69] pkgconfig_2.0.3      ModelMetrics_1.2.2.2

