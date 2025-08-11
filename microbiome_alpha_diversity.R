# Microbiome Analysis - Alpha Diversity

# load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(vegan)
library(breakaway)
library(DivNet)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name

# feature table
feat <- read.table("feature_table.txt", header = TRUE)


##################################################################################
#####   ALPHA DIVERSITY ANALYSIS - TRADITONAL ECOLOGICAL DIVERSITY METRICS   #####
##################################################################################

### alpha diversity analysis with phyloseq
# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)

# rarefy the phyloseq object (using the minimum sequencing depth of the samples)
sample_sums(ps) %>% summary() # check sequencing depth (rarefaction depth based on sequencing depth)
ps_rarefied <- rarefy_even_depth(ps, rngseed = 1234, sample.size = min(sample_sums(ps)), verbose = FALSE)

# calculate alpha diversity using phyloseq
alpha_phyloseq <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon", "Simpson", "Chao1"))
alpha_phyloseq$sample_name <- rownames(alpha_phyloseq) # add column with sample names
alpha_phyloseq <- left_join(alpha_phyloseq, meta, by = "sample_name") # merge phyloseq object with metadata (for plotting and statistical tests) 


### alpha diversity analysis with vegan
otu_vegan <- t(feat) # transpose feature table
all(rownames(otu_vegan) == meta$sample_name) # ensure sample names are the same

# rarefy feature table (using the minimum sequencing depth of the samples)
rowSums(otu_vegan) %>% summary() # check sequencing depth (rarefaction depth based on sequencing depth)
rarecurve(otu_vegan, step = 1000, col = "blue", label = TRUE) # visualize rarefaction
set.seed(1234)  # set seed
otu_vegan_rarefied <- rrarefy(otu_vegan, sample = min(rowSums(otu_vegan)))

# calculate alpha diversity using vegan
richness <- specnumber(otu_vegan_rarefied)
shannon <- diversity(otu_vegan_rarefied, index = "shannon")
simpson <- diversity(otu_vegan_rarefied, index = "simpson")

# combine alpha diversity metrics into a data.frame
alpha_vegan <- data.frame(sample_name = rownames(otu_vegan_rarefied),
                          richness_vegan = richness,
                          shannon_vegan = shannon,
                          simpson_vegan = simpson)

# merge with metadata
alpha_vegan <- left_join(alpha_vegan, meta, by = "sample_name")

### combine alpha diversity metrics from phyloseq and vegan analyses
alpha_combined <- left_join(alpha_phyloseq, alpha_vegan, by = "sample_name", suffix = c("_phyloseq", "_vegan"))
alpha_combined$condition_phyloseq <- NULL
alpha_combined <- alpha_combined %>% rename(condition = condition_vegan)
alpha_combined$condition_vegan <- "condition"


### plot alpha diversity metrics and compute statistics
# observed richness
ggplot(alpha_combined, aes(x = condition, y = Observed, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Observed richness - phyloseq", y = "Number of unique taxa")
wilcox.test(Observed ~ condition, data = alpha_combined)

# Shannon diversity
ggplot(alpha_combined, aes(x = condition, y = Shannon, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Shannon diversity - phyloseq", y = "Shannon index")
wilcox.test(Shannon ~ condition, data = alpha_combined)

# Simpson diversity
ggplot(alpha_combined, aes(x = condition, y = Simpson, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Simpson diversity - phyloseq", y = "Simpson index")
wilcox.test(Simpson ~ condition, data = alpha_combined)

# Chao1 diversity
ggplot(alpha_combined, aes(x = condition, y = Chao1, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Chao1 Estimator - phyloseq", y = "Simpson index")
wilcox.test(Chao1 ~ condition, data = alpha_combined)


#############################################################################
#####   ALPHA DIVERSITY ANALYSIS - STATISTICAL MODEL-BASED APPROACHES   #####
#############################################################################

### create a phyloseq object 
# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)


### alpha diversity analysis with breakaway
# feature table must include singletons (features with 1 read) and doubletons (features with 2 reads)
# breakaway uses these low frequency taxa ti estimate unseen richness
# many preprocesing steps (including bracken with the parameter --threshold/-t) remove taxa with very low reads
feat_otu_trans <- t(feat_otu) # transpose count matrix 

# run breakway per sample
breakaway_results <- lapply(1:nrow(feat_otu_trans), function(i) {
  tryCatch(breakaway(feat_otu_trans[i, ]),
           error = function(e) return(NULL))
})

# extract richness estimates and sample names
rich_est <- data.frame(sample_name = rownames(feat_otu_trans),
                       richness = sapply(breakaway_results, function(x) if (!is.null(x)) x$estimate else NA),
                       error = sapply(breakaway_results, function(x) if (!is.null(x)) x$error else NA))

# merge estimates with metadata
rich_est <- left_join(rich_est, meta, by = "sample_name")

# design matrix with condition as grouping variable
design_matrix_rich <- model.matrix(~ condition, data = rich_est)

# plot richness
ggplot(rich_est, aes(x = condition, y = richness, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Richness estimate - breakaway", y = "Richness estimate", x = "")

# model-based hypothesis testing for richness
betta_rich <- betta(rich_est$richness,
                    rich_est$error,
                    X = design_matrix_rich)
betta_rich$table


### alpha diversity analysis with DivNet
# run DivNet
# most abundant feature in each sample
most_abundant_taxa <- apply(feat, 2, function(x) {
  feature_index <- which.max(x) # get index of max value
  feature_name <- rownames(feat)[feature_index] # get feature name
  return(feature_name)
})
divnet_results <- divnet(ps, B = 5, base = NULL)
# saveRDS(divnet_results, file = "divnet_results.rds")
# divnet_results <- readRDS("divnet_results.rds")


# extract Shannon and Simpson diversity estimates
alpha_diver_est <- data.frame(sample_name = names(divnet_results$shannon),
                              shannon_estimate = sapply(divnet_results$shannon, function(x) x$estimate),
                              shannon_error = sapply(divnet_results$shannon, function(x) x$error),
                              simpson_estimate = sapply(divnet_results$simpson, function(x) x$estimate),
                              simpson_error = sapply(divnet_results$simpson, function(x) x$error))

# merge with metadata
alpha_diver_est <- left_join(alpha_diver_est, meta, by = "sample_name")

# design matrix with condition as grouping variable
design_matrix_alpha <- model.matrix(~ condition, data = alpha_diver_est)

# plot Shannon diversity
ggplot(alpha_diver_est, aes(x = condition, y = shannon_estimate, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Shannon diversity - DivNet", y = "Shannon estimate", x = "")

# plot Simpson diversity
ggplot(alpha_diver_est, aes(x = condition, y = simpson_estimate, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Simpson diversity - DivNet", y = "Simpson estimate", x = "")

# model-based hypothesis testing for shannon diversity 
betta_shan <- betta(alpha_diver_est$shannon_estimate, 
                    alpha_diver_est$shannon_error, 
                    X = design_matrix_alpha)
betta_shan$table

# model-based hypothesis testing for simpson diversity 
betta_simp <- betta(alpha_diver_est$simpson_estimate, 
                    alpha_diver_est$simpson_error, 
                    X = design_matrix_alpha)
betta_simp$table


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BiocParallel_1.42.1 doParallel_1.0.17   iterators_1.0.14    foreach_1.5.2       DivNet_0.4.1       
# [6] breakaway_4.8.4     vegan_2.7-1         permute_0.9-7       phyloseq_1.52.0     lubridate_1.9.4    
# [11] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4         purrr_1.0.4         readr_2.1.5        
# [16] tidyr_1.3.1         tibble_3.3.0        tidyverse_2.0.0     ggplot2_3.5.2      
# 
# loaded via a namespace (and not attached):
#   [1] ade4_1.7-23             tidyselect_1.2.1        farver_2.1.2            Biostrings_2.76.0      
# [5] digest_0.6.37           timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.8.1        
# [9] multcompView_0.1-10     processx_3.8.6          survival_3.8-3          magrittr_2.0.3         
# [13] compiler_4.5.0          rlang_1.1.6             tools_4.5.0             igraph_2.1.4           
# [17] data.table_1.17.4       labeling_0.4.3          pkgbuild_1.4.8          curl_6.3.0             
# [21] plyr_1.8.9              RColorBrewer_1.1-3      abind_1.4-8             withr_3.0.2            
# [25] desc_1.4.3              BiocGenerics_0.54.0     grid_4.5.0              stats4_4.5.0           
# [29] multtest_2.64.0         biomformat_1.36.0       Rhdf5lib_1.30.0         scales_1.4.0           
# [33] MASS_7.3-65             dichromat_2.0-0.1       cli_3.6.5               crayon_1.5.3           
# [37] reformulas_0.4.1        generics_0.1.4          remotes_2.5.0           rstudioapi_0.17.1      
# [41] httr_1.4.7              reshape2_1.4.4          tzdb_0.5.0              minqa_1.2.8            
# [45] ape_5.8-1               rhdf5_2.52.1            splines_4.5.0           BiocManager_1.30.26    
# [49] XVector_0.48.0          vctrs_0.6.5             boot_1.3-31             Matrix_1.7-3           
# [53] jsonlite_2.0.0          callr_3.7.6             IRanges_2.42.0          hms_1.1.3              
# [57] S4Vectors_0.46.0        glue_1.8.0              nloptr_2.2.1            ps_1.9.1               
# [61] codetools_0.2-20        mvnfast_0.2.8           stringi_1.8.7           gtable_0.3.6           
# [65] GenomeInfoDb_1.44.0     UCSC.utils_1.4.0        lme4_1.1-37             pillar_1.10.2          
# [69] rhdf5filters_1.20.0     GenomeInfoDbData_1.2.14 R6_2.6.1                Rdpack_2.6.4           
# [73] lattice_0.22-7          Biobase_2.68.0          rbibutils_2.3           Rcpp_1.0.14            
# [77] nlme_3.1-168            mgcv_1.9-3              pkgconfig_2.0.3 

