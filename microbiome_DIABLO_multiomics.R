### DIABLO Workflow - Metagenomics and Metabolomics

# load libraries
library(ggplot2)
library(tidyverse)
library(readxl)
library(compositions)
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

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat))


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

### feature filtering, zero imputation and log transformation 
# assess missingness (zeros listed as NAs in data.frame)
total_missing <- sum(is.na(metabo)) / prod(dim(metabo)) * 100 # overall missingness
missing_per_feature <- colMeans(is.na(metabo)) * 100 # missingness per feature
summary(missing_per_feature)
missing_per_sample <- rowMeans(is.na(metabo)) * 100 # missingness per sample
hist(missing_per_feature, breaks = 50, main = "Percent missing per metabolite", xlab = "Percent missing")

# remove features with greater than 30% missing values
metabo_filtered <- metabo[, missing_per_feature <= 30]

# impute remaining zeros with kNN
metabo_imputed <- impute.knn(as.matrix(t(metabo_filtered)))$data # impute.knn expects features as rows
metabo_imputed <- t(metabo_imputed)

# log2 transformation (with pseudocount)
metabo <- log2(metabo_imputed + 1e-6)


### make sure that samples are in the same order in both the feat and metabo
all(rownames(feat) == rownames(metabo))


##############################################################################################
###   DIABLO - DATA INTEGRATION ANALYSIS FOR BIOMARKER DISCOVERY USING LATENT COMPONENTS   ###
##############################################################################################

# list of data blocks
data_blocks <- list(metagenomics = feat, metabolomics = metabo)

# full design (maximize correlation between blocks)
full_design <- matrix(1, ncol = length(data_blocks), nrow = length(data_blocks),
                 dimnames = list(names(data_blocks), names(data_blocks)))
diag(full_design) <- 0  # no self-correlation

# calculate cross-correlations between blocks using PLS with one component
cor_blocks = pls(data_blocks$metagenomics, data_blocks$metabolomics, ncomp = 1)
cor(cor_blocks$variates$X, cor_blocks$variates$Y)

design <- matrix(c(0, 0.88, 0.88, 0), nrow = 2, byrow = TRUE, dimnames = list(names(data_blocks), names(data_blocks)))


### evaluate performance using all features
# distance not specified = all distances included
full_model_diablo <- block.plsda(X = data_blocks, Y = meta$condition, design = design, ncomp = 5)
full_model_perf <- perf(full_model_diablo,
                        validation = "Mfold", 
                        folds = 5, 
                        nrepeat = 50, 
                        progressBar = TRUE)

plot(full_model_perf) # error rate versus number of components
full_model_perf$choice.ncomp # optimal number of components
full_model_perf$error.rate # error rate for each component


### feature selection 
list_keepX <- list(metagenomics = c(5, 10, 15, 20, 25, 50),
                   metabolomics = c(5, 10, 15, 20, 25, 50))

# tune (with repeated cross-validation) number of features (list_keepX) per component 
tune_diablo <- tune.block.splsda(X = data_blocks,
                                 Y = meta$condition,
                                 test.keepX = list_keepX,
                                 design = design,
                                 ncomp = 2, 
                                 validation = "Mfold", # k-fold cross-validation (stratified)
                                 folds = 5, # 5 folds for cross-validation
                                 nrepeat = 50, # number of repeats of cross-validation
                                 dist = "max.dist", # distance metric
                                 measure = "BER", # balanced error rate
                                 progressBar = TRUE)

plot(tune_diablo, sd = TRUE) # visualize performance across keepX values
tune_diablo$choice.keepX # optimal keepX values
tune_diablo$error.rate # error rates


### final model using optimal number of features and 1 component for model performance metrics
perf_keepX <- list(metagenomics = tune_diablo$choice.keepX$metagenomics[[1]], metabolomics = tune_diablo$choice.keepX$metabolomics[[1]])
final_model_diablo <- block.splsda(X = data_blocks, Y = meta$condition, keepX = perf_keepX, design = design, ncomp = 1)
final_diablo_perf <- perf(final_model_diablo,
                          validation = "Mfold", 
                          folds = 5, 
                          nrepeat = 50, 
                          dist = "max.dist", 
                          measure = "overall",
                          progressBar = TRUE)

final_diablo_perf$error.rate # error rate for final model
final_diablo_perf$error.rate.per.class # error rate per class for final model
plotDiablo(final_model_diablo) # correlation between the two blocks


### ROC curve (evaluation of model performance)
auroc(final_model_diablo, roc.block = 1, roc.comp = 1) # metagenomics
auroc(final_model_diablo, roc.block = 2, roc.comp = 1) # metabolomics


### generate confusion matrix
predictions_df <- as.data.frame(final_diablo_perf$AveragedPredict.class)  # class predictions
# function to get majority vote from the 50 repeats
get_majority_vote <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1] 
}
avg_predictions <- apply(predictions_df, 1, get_majority_vote) # get majority prediction
confusionMatrix(as.factor(avg_predictions), as.factor(meta$condition))


### feature stability across repeats (metagenomics)
# extract feature names (with non-zero stability) from each repeat
feature_list_metagen <- lapply(final_diablo_perf$features$stable, function(rep) {
  comp_vec <- rep$metagenomics$comp1
  names(comp_vec[comp_vec > 0])
})
feature_stability_metagen <- as.data.frame(table(unlist(feature_list_metagen)))
colnames(feature_stability_metagen) <- c("feature", "count")
feature_stability_metagen <- feature_stability_metagen[order(-feature_stability_metagen$count), ]

# plot feature stability across repeats
features_df <- as.data.frame(head(feature_stability_metagen, 50)) # top features by selection frequency (metagenomics)
features_df$feature <- factor(features_df$feature, levels = features_df$feature) # order features for plot
ggplot(features_df, aes(x = feature, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Feature stability (metagenomics)", x = "Selected features", y = "Frequency of selection") +
  theme_minimal(base_size = 12) + theme(axis.text.y = element_text(size = 10))


### feature stability across repeats (metabolomics)
# extract feature names (with non-zero stability) from each repeat
feature_list_metabo <- lapply(final_diablo_perf$features$stable, function(rep) {
  comp_vec <- rep$metabolomics$comp1
  names(comp_vec[comp_vec > 0])
})
feature_stability_metabo <- as.data.frame(table(unlist(feature_list_metabo)))
colnames(feature_stability_metabo) <- c("feature", "count")
feature_stability_metabo <- feature_stability_metabo[order(-feature_stability_metabo$count), ]

# plot feature stability across repeats
features_df <- as.data.frame(head(feature_stability_metabo, 10)) # top features by selection frequency (metabolomics)
features_df$feature <- factor(features_df$feature, levels = features_df$feature) # order features for plot
ggplot(features_df, aes(x = feature, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Feature stability (metabolomics)", x = "Selected features", y = "Frequency of selection") +
  theme_minimal(base_size = 12) + theme(axis.text.y = element_text(size = 10))


### loading weights
selectVar(final_model_diablo, block = "metagenomics")
selectVar(final_model_diablo, block = "metabolomics")

# plot loading weights (metagenomics)
plotLoadings(final_model_diablo, 
             comp = 1, block = "metagenomics",
             method = "median", contrib = "max", ndisplay = 50,
             legend.color = c("steelblue", "indianred3"), legend = FALSE,
             title = "Loadings plot") 

# plot loading weights (metabolomics)
plotLoadings(final_model_diablo, 
             comp = 1, block = "metabolomics",
             method = "median", contrib = "max", ndisplay = 5,
             legend.color = c("steelblue", "indianred3"), legend = FALSE,
             title = "Loadings plot") 


### circos plot (relationship between features in different blocks)
circosPlot(final_model_diablo, cutoff = 0.50,
           line = TRUE, showIntraLinks = TRUE,
           color.blocks = c("steelblue", "indianred3"),
           title = "DIABLO - Circos Plot")


### plot relevance network
network(final_model_diablo, blocks = c(1,2), cutoff = 0.50, 
        color.node = c("steelblue", "indianred3"),
        color.edge = c("red4","red3","red2","green2","green3","green4"),
        cex.node.name = 0.3, lwd.edge = 1.5, block.var.names = TRUE, 
        interactive = TRUE)


### final model using optimal number of features and 2 component for plotting
plot_model_diablo <- block.splsda(X = data_blocks, Y = meta$condition, keepX = tune_diablo$choice.keepX, design = design, ncomp = 2)

### sample plot (visualization of sample distribution in component space)
plotIndiv(plot_model_diablo, comp = c(1,2),
          ind.names = FALSE, legend = TRUE, 
          X.label = "Component 1", Y.label = "Component 2",
          star = TRUE, ellipse = TRUE, ellipse.level = 0.95,
          col = c("steelblue", "indianred3"),
          cex = 3, pch = c(1,1), point.lwd = 0.5,
          title = "Sample plot")

### feature plot/correlation circle plot (visualization of features in the component space)
plotVar(plot_model_diablo, comp = c(1,2),
        var.names = TRUE, legend = TRUE, style = "graphics",
        col = c("steelblue", "indianred3"), pch = c(1,1), cex = c(0.5, 0.5),
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
# [5] impute_1.82.0      compositions_2.0-8 readxl_1.4.5       lubridate_1.9.4   
# [9] forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.1.0       
# [13] readr_2.1.5        tidyr_1.3.1        tibble_3.3.0       tidyverse_2.0.0   
# [17] ggplot2_3.5.2     
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1     timeDate_4041.110    farver_2.1.2         tensorA_0.36.2.1    
# [5] pROC_1.19.0.1        digest_0.6.37        rpart_4.1.24         timechange_0.3.0    
# [9] lifecycle_1.0.4      survival_3.8-3       magrittr_2.0.3       compiler_4.5.0      
# [13] rlang_1.1.6          tools_4.5.0          igraph_2.1.4         data.table_1.17.8   
# [17] labeling_0.4.3       rARPACK_0.11-0       plyr_1.8.9           RColorBrewer_1.1-3  
# [21] BiocParallel_1.42.1  withr_3.0.2          stats4_4.5.0         nnet_7.3-20         
# [25] grid_4.5.0           e1071_1.7-16         future_1.67.0        globals_0.18.0      
# [29] scales_1.4.0         iterators_1.0.14     dichromat_2.0-0.1    cli_3.6.5           
# [33] ellipse_0.5.0        generics_0.1.4       rstudioapi_0.17.1    future.apply_1.20.0 
# [37] robustbase_0.99-4-1  RSpectra_0.16-2      reshape2_1.4.4       tzdb_0.5.0          
# [41] proxy_0.4-27         bayesm_3.1-6         splines_4.5.0        parallel_4.5.0      
# [45] cellranger_1.1.0     matrixStats_1.5.0    vctrs_0.6.5          hardhat_1.4.1       
# [49] Matrix_1.7-3         hms_1.1.3            ggrepel_0.9.6        listenv_0.9.1       
# [53] foreach_1.5.2        gower_1.0.2          recipes_1.3.1        glue_1.8.0          
# [57] parallelly_1.45.1    DEoptimR_1.1-4       codetools_0.2-20     stringi_1.8.7       
# [61] gtable_0.3.6         pillar_1.11.0        ipred_0.9-15         lava_1.8.1          
# [65] R6_2.6.1             corpcor_1.6.10       class_7.3-23         Rcpp_1.1.0          
# [69] gridExtra_2.3        nlme_3.1-168         prodlim_2025.04.28   pkgconfig_2.0.3     
# [73] ModelMetrics_1.2.2.2

