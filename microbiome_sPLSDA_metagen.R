### sPLSDA Workflow - Metagenomics

# load libraries
library(ggplot2)
library(tidyverse)
library(compositions)
library(mixOmics)
library(caret)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name
str(meta)

# convert condition column into a factor
meta$condition <- as.factor(meta$condition) # needs to be a factor for sPLS-DA
table(meta$condition)


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

### filter species present in less than 10% of samples
feat_t <- as.data.frame(t(feat)) # transpose feature table
min_prevalence <- 0.10
feat_filtered <- feat_t[, colMeans(feat_t > 0) >= min_prevalence] # subset the feature table to only include features present in at least 10% of samples
dim(feat_filtered) # 70 935

# convert feature table to relative abundances
feat_rel_abund <- feat_filtered/rowSums(feat_filtered)

# add pseudocount and perform CLR transformation
feat_clr <- clr(feat_rel_abund + 1e-6)
feat <- as.data.frame(feat_clr)

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat))


##################################################################################
########   MIXOMICS - SPARSE PARTIAL LEAST SQUARES DISCRIMINANT ANALYSIS   #######
##################################################################################

### evaluate performance using all features
# distance not specified = all distances included
plsda_all_feat <- mixOmics::plsda(X = feat, Y = meta$condition, ncomp = 5)
perf_res_all <- perf(plsda_all_feat,
                     validation = "Mfold", 
                     folds = 5, 
                     nrepeat = 50, 
                     progressBar = TRUE)

plot(perf_res_all) # error rate versus number of components
perf_res_all$choice.ncomp # optimal number of components
perf_res_all$error.rate # error rate for each component


### feature selection using sparse method
list_keepX <- c(5, 10, 15, 20, 25, 50, 75, 100) # list of number of features to keep per component

### tune (with repeated cross-validation) number of features (keepX) per component 
tune_splsda_res <- tune.splsda(X = feat,
                               Y = meta$condition,
                               test.keepX = list_keepX,
                               ncomp = 5, 
                               validation = "Mfold", # k-fold cross-validation (stratified)
                               folds = 5, # 5 folds for cross-validation
                               nrepeat = 50, # number of repeats of cross-validation
                               dist = "max.dist", # distance metric
                               measure = "BER", # balanced error rate
                               progressBar = TRUE)

plot(tune_splsda_res, sd = TRUE) # visualize performance across keepX values
tune_splsda_res$choice.ncomp # optimal number of components
tune_splsda_res$error.rate # error rates for keepX used
tune_splsda_res$choice.keepX # optimal number of features


### function to run tune.splsda for all three distance metrics
tune_all_distances <- function(X, Y, list.keepX, ncomp = 5, folds = 5, nrepeat = 50, measure = "BER") {

  dist_methods <- c("max.dist", "centroids.dist", "mahalanobis.dist") # define distance methods
  results <- list() # initialize a list to hold results
  
  # loop over each distance method
  for (dm in dist_methods) {
    set.seed(1234)
    
    # tune with keepX using each distance method
    res <- tune.splsda(X = X,
                       Y = Y,
                       test.keepX = list_keepX,
                       ncomp = ncomp,
                       validation = "Mfold",
                       folds = folds,
                       nrepeat = nrepeat,
                       dist = dm,
                       measure = measure,
                       progressBar = TRUE)
    
    results[[dm]] <- res # store results for each distance method
  }
  return(results)
}

# list of number of features to keep per component
list_keepX <- c(5, 10, 15, 20, 25, 50, 75, 100)

### tune (with repeated cross-validation) number of features (keepX) per component for all three distance metrics
tuned_results <- tune_all_distances(X = feat,
                                    Y = meta$condition,
                                    list.keepX = list_keepX,
                                    ncomp = 5,
                                    folds = 5,
                                    nrepeat = 50)

lapply(tuned_results, function(dist) dist$error.rate) # BER for each distance metric (for each keepX value)
lapply(tuned_results, function(dist) dist$choice.keepX) # optimal number of features
lapply(tuned_results, function(dist) dist$choice.ncomp) # optimal number of components
sapply(tuned_results, function(res) min(res$error.rate)) # minimum BER (for each distance metric)

# visualize performance across keepX values
plot(tuned_results$max.dist, sd = TRUE) # max.dist
plot(tuned_results$centroids.dist, sd = TRUE) # centroids.dist 
plot(tuned_results$mahalanobis.dist, sd = TRUE) # mahalanobis.dist


### final model using optimal number of features and 1 component for model performance metrics
final_model_perf <- mixOmics::splsda(X = feat, Y = meta$condition, ncomp = 1, keepX = tuned_results$max.dist$choice.keepX)
final_splsda_perf <- perf(final_model_perf,
                         validation = "Mfold", 
                         folds = 5, 
                         nrepeat = 50, 
                         dist = "max.dist", 
                         measure = "BER",
                         progressBar = TRUE)

final_splsda_perf$error.rate # error rate for final model
final_splsda_perf$error.rate.class # error rate per class for final model


### ROC curve (evaluation of model performance)
auroc(final_model_perf, roc.comp = 1)


### generate confusion matrix
predictions_df <- as.data.frame(final_splsda_perf$class)  # class predictions
# function to get majority vote from the 50 repeats
get_majority_vote <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1] 
}
avg_predictions <- apply(predictions_df, 1, get_majority_vote) # get majority prediction
confusionMatrix(as.factor(avg_predictions), as.factor(meta$condition))


### plot frequency of selection for top selected features for component 1
final_model_features <- final_splsda_perf$features$stable$comp1 # features by selection frequency (component 1)
final_model_features_top <- sort(head(final_model_features, 25), decreasing = FALSE) # top 25 features by selection frequency
features_df <- data.frame(feature = names(final_model_features_top),
                          frequency = as.numeric(final_model_features_top)) # create data.frame for ggplot2
features_df$feature <- factor(features_df$feature, levels = features_df$feature) # order features for plot

# ggplot of frequency of feature selection
ggplot(features_df, aes(x = feature, y = frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Feature stability (component 1)", x = "Selected features", y = "Frequency of selection") +
  theme_minimal(base_size = 12) + theme(axis.text.y = element_text(size = 10))


### loading weights
selectVar(final_model_perf)

### plot loading weights (contribution of features to component)
plotLoadings(final_model_perf, 
             comp = 1, method = "median", contrib = "max",
             legend.color = c("steelblue", "indianred3"), legend = FALSE,
             title = "Loadings plot")


### final model using optimal number of features and 2 components for plotting
plot_model_perf <- mixOmics::splsda(X = feat, Y = meta$condition, ncomp = 2, keepX = tuned_results$max.dist$choice.keepX)
plot_splsda_perf <- perf(plot_model_perf,
                         validation = "Mfold", 
                         folds = 5, 
                         nrepeat = 50, 
                         dist = "max.dist", 
                         measure = "BER",
                         progressBar = TRUE)

plot_splsda_perf$error.rate # error rate for plot model
plot_splsda_perf$error.rate.class # error rate per class for plot model


### sample plot (visualization of sample distribution in component space)
plotIndiv(plot_model_perf, comp = c(1,2),
          ind.names = FALSE, legend = TRUE, 
          X.label = "Component 1", Y.label = "Component 2",
          star = TRUE, ellipse = TRUE, ellipse.level = 0.95,
          col = c("steelblue", "indianred3"),
          cex = 3, pch = c(1,1), point.lwd = 0.5,
          title = "Sample plot")


### feature plot/correlation circle plot (visualization of features in the component space)
plotVar(plot_model_perf, comp = c(1,2),
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
#   [1] caret_7.0-1        mixOmics_6.32.0    lattice_0.22-7     MASS_7.3-65       
# [5] compositions_2.0-8 lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [9] dplyr_1.1.4        purrr_1.1.0        readr_2.1.5        tidyr_1.3.1       
# [13] tibble_3.3.0       tidyverse_2.0.0    ggplot2_3.5.2     
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6         ellipse_0.5.0        tensorA_0.36.2.1    
# [4] recipes_1.3.1        ggrepel_0.9.6        tzdb_0.5.0          
# [7] vctrs_0.6.5          tools_4.5.0          generics_0.1.4      
# [10] stats4_4.5.0         parallel_4.5.0       DEoptimR_1.1-4      
# [13] ModelMetrics_1.2.2.2 rARPACK_0.11-0       pkgconfig_2.0.3     
# [16] Matrix_1.7-3         data.table_1.17.8    RColorBrewer_1.1-3  
# [19] lifecycle_1.0.4      compiler_4.5.0       farver_2.1.2        
# [22] codetools_0.2-20     class_7.3-23         prodlim_2025.04.28  
# [25] pillar_1.11.0        BiocParallel_1.42.1  gower_1.0.2         
# [28] iterators_1.0.14     rpart_4.1.24         foreach_1.5.2       
# [31] parallelly_1.45.1    lava_1.8.1           nlme_3.1-168        
# [34] RSpectra_0.16-2      robustbase_0.99-4-1  digest_0.6.37       
# [37] tidyselect_1.2.1     future_1.67.0        stringi_1.8.7       
# [40] listenv_0.9.1        reshape2_1.4.4       splines_4.5.0       
# [43] grid_4.5.0           cli_3.6.5            magrittr_2.0.3      
# [46] dichromat_2.0-0.1    survival_3.8-3       future.apply_1.20.0 
# [49] corpcor_1.6.10       withr_3.0.2          scales_1.4.0        
# [52] timechange_0.3.0     globals_0.18.0       matrixStats_1.5.0   
# [55] igraph_2.1.4         nnet_7.3-20          timeDate_4041.110   
# [58] gridExtra_2.3        hms_1.1.3            hardhat_1.4.1       
# [61] rlang_1.1.6          Rcpp_1.1.0           glue_1.8.0          
# [64] bayesm_3.1-6         pROC_1.19.0.1        ipred_0.9-15        
# [67] rstudioapi_0.17.1    R6_2.6.1             plyr_1.8.9          

