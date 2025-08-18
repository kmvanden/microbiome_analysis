# Microbiome Analysis - Beta Diversity

# load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(vegan)
library(microbiome)
library(DivNet)
library(breakaway)
library(ape)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name

# feature table
feat <- read.table("feature_table.txt", header = TRUE)


#########################################################################
#####   BETA DIVERSITY ANALYSIS - NO PREVALENCE FILTERING OF TAXA   #####
#########################################################################

### beta diversity analysis with phyloseq (mainly)
# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x)) # relative abundance for bray-curtis and canberra
ps_pa <- transform_sample_counts(ps, function(x) as.numeric(x > 0)) # presence/absence for jaccard
ps_log <- transform_sample_counts(ps, function(x) log1p(x))  # log-transformation for euclidean
ps_clr <- transform(ps, "clr")  # clr transformation (microbiome package adds pseudocount) for aitchison

# compute distance matrices 
bray_dist <- distance(ps_rel, method = "bray") # bray-curtis
jaccard_dist <- distance(ps_pa, method = "jaccard", binary = TRUE) # jaccard
euc_dist <- distance(ps_log, method = "euclidean") # euclidean
canberra_dist <- distance(ps_rel, method = "canberra") # canberra
aitchison_dist <- distance(ps_clr, method = "euclidean")  # aitchison (euclidean in clr space)

# statistical testing with PERMANOVA and betadispar
adonis2(bray_dist ~ condition, data = meta) # bray-curtis
betadisper(bray_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(bray_dist, group = meta$condition))

adonis2(jaccard_dist ~ condition, data = meta) # jaccard
betadisper(jaccard_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(jaccard_dist, group = meta$condition))

adonis2(euc_dist ~ condition, data = meta) # euclidean
betadisper(euc_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(euc_dist, group = meta$condition))

adonis2(canberra_dist ~ condition, data = meta) # canberra
betadisper(canberra_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(canberra_dist, group = meta$condition))

adonis2(aitchison_dist ~ condition, data = meta) # aitchison
betadisper(aitchison_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(canberra_dist, group = meta$condition))


# PCoA ordination
ordination_pcoa_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist) # bray-curtis
ordination_pcoa_jaccard <- ordinate(ps_pa, method = "PCoA", distance = jaccard_dist) # jaccard
ordination_pcoa_euc <- ordinate(ps_log, method = "PCoA", distance = euc_dist) # euclidean
ordination_pcoa_canberra <- ordinate(ps_rel, method = "PCoA", distance = canberra_dist) # canberra
ordination_pcoa_aitchison <- ordinate(ps_clr, method = "PCoA", distance = aitchison_dist) # aitchison

# PCoA plots
plot_ordination(ps_rel, ordination_pcoa_bray, color = "condition") +
  geom_point(size = 2) + ggtitle("PCoA - Bray-Curtis") + theme_minimal() # bray-curtis
plot_ordination(ps_pa, ordination_pcoa_jaccard, color = "condition") +
  geom_point(size = 2) + ggtitle("PCoA - Jaccard") + theme_minimal() # jaccard
plot_ordination(ps_log, ordination_pcoa_euc, color = "condition") +
  geom_point(size = 2) + ggtitle("PCoA - Euclidean") + theme_minimal() # euclidean
plot_ordination(ps_rel, ordination_pcoa_canberra, color = "condition") +
  geom_point(size = 2) + ggtitle("PCoA - Canberra") + theme_minimal() # canberra
plot_ordination(ps_clr, ordination_pcoa_aitchison, color = "condition") +
  geom_point(size = 2) + ggtitle("PCoA - Aitchison") + theme_minimal() # aitchison


# NMDS ordination
ordination_nmds_bray <- ordinate(ps_rel, method = "NMDS", distance = bray_dist) # bray-curtis
ordination_nmds_jaccard <- ordinate(ps_pa, method = "NMDS", distance = jaccard_dist) # jaccard
# ordination_nmds_euc <- ordinate(ps_log, method = "NMDS", distance = euc_dist) # euclidean
ordination_nmds_canberra <- ordinate(ps_rel, method = "NMDS", distance = canberra_dist) # canberra
# ordination_nmds_aitchison <- ordinate(ps_clr, method = "NMDS", distance = aitchison_dist) # aitchison

# NMDS plots
plot_ordination(ps_rel, ordination_nmds_bray, color = "condition") +
  geom_point(size = 2) + ggtitle("NMDS - Bray-Curtis") + theme_minimal() # bray-curtis
plot_ordination(ps_pa, ordination_nmds_jaccard, color = "condition") +
  geom_point(size = 2) + ggtitle("NMDS - Jaccard") + theme_minimal() # jaccard
plot_ordination(ps_log, ordination_nmds_euc, color = "condition") +
  geom_point(size = 2) + ggtitle("NMDS - Euclidean") + theme_minimal() # euclidean
plot_ordination(ps_rel, ordination_nmds_canberra, color = "condition") +
  geom_point(size = 2) + ggtitle("NMDS - Canberra") + theme_minimal() # canberra
plot_ordination(ps_clr, ordination_nmds_aitchison, color = "condition") +
  geom_point(size = 2) + ggtitle("NMDS - Aitchison") + theme_minimal() # aitchison


### PCA on CLR-transformed data
# prepare data
clr_otu <- t(otu_table(ps_clr)) # extract CLR-transformed data and transpose
pca_clr <- prcomp(clr_otu, center = TRUE, scale. = FALSE) # run PCA
pca_df <- as.data.frame(pca_clr$x) # extract PCA scores (coordinates)
pca_df$sample_name <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = "sample_name")  # join with metadata
var_explained <- round(100 * summary(pca_clr)$importance[2, 1:2], 1) # extract percentage of variance explained

# PCA plot of CLR-transformed data
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) + 
  geom_point(size = 2) + theme_minimal() + stat_ellipse(type = "norm", level = 0.95) +
  labs(title = "PCA of CLR-transformed data", 
       x = paste0("PC1 (", var_explained[1], "% variance)"), 
       y = paste0("PC2 (", var_explained[2], "% variance)"))


######################################################################
#####   BETA DIVERSITY ANALYSIS - PREVALENCE FILTERING OF TAXA   #####
######################################################################

### beta diversity analysis with phyloseq (mainly)
# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object and filter low prevalence taxa (present in less than 10% of samples)
ps <- phyloseq(feat_otu, sampledata)
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= 0.1 * nsamples(ps), prune = TRUE)

ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x)) # relative abundance for bray-curtis and canberra
ps_pa <- transform_sample_counts(ps_filt, function(x) as.numeric(x > 0)) # presence/absence for jaccard
ps_log <- transform_sample_counts(ps_filt, function(x) log1p(x))  # log-transformation for euclidean
ps_clr <- transform(ps_filt, "clr")  # clr transformation (microbiome package adds pseudocount) for aitchison

# compute distance matrices 
bray_dist <- distance(ps_rel, method = "bray") # bray-curtis
jaccard_dist <- distance(ps_pa, method = "jaccard", binary = TRUE) # jaccard
euc_dist <- distance(ps_log, method = "euclidean") # euclidean
canberra_dist <- distance(ps_rel, method = "canberra") # canberra
aitchison_dist <- distance(ps_clr, method = "euclidean")  # aitchison (euclidean in clr space)

# statistical testing with PERMANOVA and betadisper
adonis2(bray_dist ~ condition, data = meta) # bray-curtis
betadisper(bray_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(bray_dist, group = meta$condition))

adonis2(jaccard_dist ~ condition, data = meta) # jaccard
betadisper(jaccard_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(jaccard_dist, group = meta$condition))

adonis2(euc_dist ~ condition, data = meta) # euclidean
betadisper(euc_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(euc_dist, group = meta$condition))

adonis2(canberra_dist ~ condition, data = meta) # canberra
betadisper(canberra_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(canberra_dist, group = meta$condition))

adonis2(aitchison_dist ~ condition, data = meta) # aitchison
betadisper(aitchison_dist, group = meta$condition) %>% permutest()
boxplot(betadisper(aitchison_dist, group = meta$condition))


# PCoA ordination
ordination_pcoa_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist) # bray-curtis
ordination_pcoa_jaccard <- ordinate(ps_pa, method = "PCoA", distance = jaccard_dist) # jaccard
ordination_pcoa_euc <- ordinate(ps_log, method = "PCoA", distance = euc_dist) # euclidean
ordination_pcoa_canberra <- ordinate(ps_rel, method = "PCoA", distance = canberra_dist) # canberra
ordination_pcoa_aitchison <- ordinate(ps_clr, method = "PCoA", distance = aitchison_dist) # aitchison

# PCoA plots
plot_ordination(ps_rel, ordination_pcoa_bray, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("PCoA - Bray-Curtis") + theme_minimal() # bray-curtis
plot_ordination(ps_pa, ordination_pcoa_jaccard, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("PCoA - Jaccard") + theme_minimal() # jaccard
plot_ordination(ps_log, ordination_pcoa_euc, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("PCoA - Euclidean") + theme_minimal() # euclidean
plot_ordination(ps_rel, ordination_pcoa_canberra, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("PCoA - Canberra") + theme_minimal() # canberra
plot_ordination(ps_clr, ordination_pcoa_aitchison, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("PCoA - Aitchison") + theme_minimal() # aitchison


# NMDS ordination
ordination_nmds_bray <- ordinate(ps_rel, method = "NMDS", distance = bray_dist) # bray-curtis
ordination_nmds_jaccard <- ordinate(ps_pa, method = "NMDS", distance = jaccard_dist) # jaccard
ordination_nmds_euc <- ordinate(ps_log, method = "NMDS", distance = euc_dist) # euclidean
ordination_nmds_canberra <- ordinate(ps_rel, method = "NMDS", distance = canberra_dist) # canberra
ordination_nmds_aitchison <- ordinate(ps_clr, method = "NMDS", distance = aitchison_dist) # aitchison

# NMDS plots
plot_ordination(ps_rel, ordination_nmds_bray, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("NMDS - Bray-Curtis") + theme_minimal() # bray-curtis
plot_ordination(ps_pa, ordination_nmds_jaccard, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("NMDS - Jaccard") + theme_minimal() # jaccard
plot_ordination(ps_log, ordination_nmds_euc, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("NMDS - Euclidean") + theme_minimal() # euclidean
plot_ordination(ps_rel, ordination_nmds_canberra, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("NMDS - Canberra") + theme_minimal() # canberra
plot_ordination(ps_clr, ordination_nmds_aitchison, color = "condition") + stat_ellipse(type = "norm", level = 0.95) +
  geom_point(size = 2) + ggtitle("NMDS - Aitchison") + theme_minimal() # aitchison


### PCA on CLR-transformed data
# prepare data
clr_otu <- t(otu_table(ps_clr)) # extract CLR-transformed data and transpose
pca_clr <- prcomp(clr_otu, center = TRUE, scale. = FALSE) # run PCA
pca_df <- as.data.frame(pca_clr$x) # extract PCA scores (coordinates)
pca_df$sample_name <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = "sample_name")  # join with metadata
var_explained <- round(100 * summary(pca_clr)$importance[2, 1:2], 1) # extract percentage of variance explained

# PCA plot of CLR-transformed data
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) + 
  geom_point(size = 2) + theme_minimal() + stat_ellipse(type = "norm", level = 0.95) +
  labs(title = "PCA of CLR-transformed data", 
       x = paste0("PC1 (", var_explained[1], "% variance)"), 
       y = paste0("PC2 (", var_explained[2], "% variance)"))


###################################################################################
#####   BETA DIVERSITY ANALYSIS - STATISTICAL MODEL-BASED APPROACH - DIVNET   #####
###################################################################################

### create a phyloseq object 
# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)


### beta diversity analysis with DivNet
# most abundant feature in each sample
most_abundant_taxa <- apply(feat, 2, function(x) {
  feature_index <- which.max(x) # get index of max value
  feature_name <- rownames(feat)[feature_index] # get feature name
  return(feature_name)
})

# run DivNet
divnet_results <- divnet(ps, B = 5, base = NULL)
# saveRDS(divnet_results, file = "divnet_results.rds")
# divnet_results <- readRDS("divnet_results.rds")


# create sample specimen matrix
sample_specimen_matrix <- diag(nrow(sample_data(ps)))
rownames(sample_specimen_matrix) <- rownames(meta)
colnames(sample_specimen_matrix) <- rownames(meta)


test_bray <- testBetaDiversity(dv = divnet_results, 
                               h0 = "bray-curtis",
                               groups = as.numeric(factor(meta$condition)),
                               sample_specimen_matrix = sample_specimen_matrix,
                               n_boot = 1000)
test_bray$p_value

test_euc <- testBetaDiversity(dv = divnet_results, 
                              h0 = "euclidean",
                              groups = as.numeric(factor(meta$condition)),
                              sample_specimen_matrix = sample_specimen_matrix,
                              n_boot = 1000)
test_euc$p_value


# extract the bray-curtis and euclidean matrices
bray_divnet <- divnet_results$`bray-curtis`
euc_divnet <- divnet_results$euclidean
# isSymmetric(bray_divnet)
# all(diag(bray_divnet) == 0) 

# PCoA ordination
# use pcoa() from ape which works on generic, externally computed distance matrices (not on phyloseq objects)
bray_pcoa <- pcoa(bray_divnet)
ordination_coords_bray <- bray_pcoa$vectors[, 1:2]  # extract first two PCoA axes
var_explained_bray <- 100 * bray_pcoa$values$Relative_eig[1:2] # extract percentage of variance explained

euc_pcoa <- pcoa(euc_divnet)
ordination_coords_euc <- euc_pcoa$vectors[, 1:2]  # extract first two PCoA axes
var_explained_euc <- 100 * euc_pcoa$values$Relative_eig[1:2] # extract percentage of variance explained

# PCoA plot of bray-curtis (divnet)
ordination_df_bray <- as.data.frame(ordination_coords_bray)
ordination_df_bray$condition <- meta$condition  # add metadata to data.frame
ggplot(ordination_df_bray, aes(x = Axis.1, y = Axis.2, color = condition)) +
  geom_point(size = 2) + theme_minimal() + stat_ellipse(type = "norm", level = 0.95) +
  labs(title = "PCoA - Bray-Curtis (DivNet)", 
       x = paste0("PCoA 1 (", round(var_explained_bray[1], 1), "%)"), 
       y = paste0("PCoA 2 (", round(var_explained_bray[2], 1), "%)"))
  
# PCoA plot of euclidean (divnet)
ordination_df_euc <- as.data.frame(ordination_coords_euc)
ordination_df_euc$condition <- meta$condition  # add metadata to data.frame
ggplot(ordination_df_euc, aes(x = Axis.1, y = Axis.2, color = condition)) +
  geom_point(size = 2) + theme_minimal() + stat_ellipse(type = "norm", level = 0.95) +
  labs(title = "PCoA - Euclidean (DivNet)", 
       x = paste0("PCoA 1 (", round(var_explained_euc[1], 1), "%)"), 
       y = paste0("PCoA 2 (", round(var_explained_euc[2], 1), "%)"))


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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ape_5.8-1         microbiome_1.30.0 DivNet_0.4.1      breakaway_4.8.4   vegan_2.7-1       permute_0.9-7    
# [7] phyloseq_1.52.0   lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.4      
# [13] readr_2.1.5       tidyr_1.3.1       tibble_3.3.0      tidyverse_2.0.0   ggplot2_3.5.2    
# 
# loaded via a namespace (and not attached):
#   [1] ade4_1.7-23             tidyselect_1.2.1        farver_2.1.2            Biostrings_2.76.0      
# [5] digest_0.6.37           timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.8.1        
# [9] multcompView_0.1-10     survival_3.8-3          magrittr_2.0.3          compiler_4.5.0         
# [13] rlang_1.1.6             tools_4.5.0             igraph_2.1.4            data.table_1.17.4      
# [17] labeling_0.4.3          plyr_1.8.9              RColorBrewer_1.1-3      abind_1.4-8            
# [21] Rtsne_0.17              withr_3.0.2             BiocGenerics_0.54.0     grid_4.5.0             
# [25] stats4_4.5.0            multtest_2.64.0         biomformat_1.36.0       Rhdf5lib_1.30.0        
# [29] scales_1.4.0            iterators_1.0.14        MASS_7.3-65             dichromat_2.0-0.1      
# [33] cli_3.6.5               crayon_1.5.3            reformulas_0.4.1        generics_0.1.4         
# [37] rstudioapi_0.17.1       httr_1.4.7              reshape2_1.4.4          tzdb_0.5.0             
# [41] minqa_1.2.8             rhdf5_2.52.1            splines_4.5.0           parallel_4.5.0         
# [45] XVector_0.48.0          vctrs_0.6.5             boot_1.3-31             Matrix_1.7-3           
# [49] jsonlite_2.0.0          IRanges_2.42.0          hms_1.1.3               S4Vectors_0.46.0       
# [53] foreach_1.5.2           glue_1.8.0              nloptr_2.2.1            codetools_0.2-20       
# [57] mvnfast_0.2.8           stringi_1.8.7           gtable_0.3.6            GenomeInfoDb_1.44.0    
# [61] UCSC.utils_1.4.0        lme4_1.1-37             pillar_1.10.2           rhdf5filters_1.20.0    
# [65] GenomeInfoDbData_1.2.14 R6_2.6.1                Rdpack_2.6.4            doParallel_1.0.17      
# [69] lattice_0.22-7          Biobase_2.68.0          rbibutils_2.3           Rcpp_1.0.14            
# [73] nlme_3.1-168            mgcv_1.9-3              pkgconfig_2.0.3

