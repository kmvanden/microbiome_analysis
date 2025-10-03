# Enterotyping Workflow

# load libraries
library(ggplot2)
library(tidyverse)
library(DirichletMultinomial)
library(parallel)
library(vegan)
library(phyloseq)
library(cluster)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ashr)


# set.seed
set.seed(1234)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

### filter species present in less than 30% of samples
feat_t <- as.data.frame(t(feat)) # transpose feature table
min_prevalence <- 0.30
feat_filtered <- feat_t[, colMeans(feat_t > 0) >= min_prevalence] # subset the feature table to only include features present in at least 10% of samples
# feat_species <- as.matrix(feat_filtered) # feature table using species

### collapse to genus
genus_names <- sub("_.*", "", colnames(feat_filtered)) # collapse feature names to genus
feat_filtered_t <- as.data.frame(t(feat_filtered)) # transpose
feat_filtered_t$genus <- genus_names # add genus column
feat_long_genus <- feat_filtered_t %>%
  group_by(genus) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") # sum rows by genus

feat_genus <- as.matrix(feat_long_genus[,-1])
rownames(feat_genus) <- feat_long_genus$genus
feat_genus <- t(feat_genus) # feature table using genus

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat_genus))


##############################################
#####   DIRICHLET-MULTINOMIAL MODELING   #####
##############################################

# fit DMM models for k = 1 to k = 5
fit_dmm <- mclapply(1:5, function(k) dmn(feat_genus, k = k), mc.cores = parallel::detectCores() - 1)

# compare fit of DMM models using BIC values
BIC_values <- sapply(fit_dmm, BIC)
best_k_val <- which.min(BIC_values) # index of best model

# plot BIC values versus number of clusters
plot(1:5, BIC_values, type = "b", pch = 16,
     xlab = "Number of clusters", ylab = "BIC values",
     main = "Model selection using BIC")

# choose best model and assign clusters
# best_model <- fit_dmm[[best_k_val]]
best_model <- fit_dmm[[3]]

cluster_probs <- best_model@group
cluster_assign <- apply(cluster_probs, 1, which.max)
cluster_ids <- data.frame(sample = rownames(feat_genus), cluster = cluster_assign)


### plot proportions of clusters by condition and condition by cluster 
# convert to long format for plotting
cluster_probs_df <- as.data.frame(cluster_probs)
cluster_probs_df$sample_id <- rownames(cluster_probs_df) # sample_id
cluster_probs_df$condition <- meta$condition # condition
cluster_probs_long <- cluster_probs_df %>%
  pivot_longer(cols = -c(sample_id, condition), names_to = "cluster", values_to = "abundance")

# plot proportion of cluster by condition
ggplot(cluster_probs_long, aes(x = cluster, y = abundance, fill = condition)) +
  geom_bar(stat = "identity", position = "fill") + theme_minimal() +
  labs(title = "Enterotype composition", y = "Proportion", x = "Cluster")

# plot proportion of condition by cluster
ggplot(cluster_probs_long, aes(x = condition, y = abundance, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") + theme_minimal() +
  labs(title = "Condition composition", y = "Proportion", x = "Condition")

# enterotype clusters associated with condition
table(cluster_ids$cluster, meta$condition)
fisher.test(table(cluster_ids$cluster, meta$condition)) # Fisher’s exact test


#####################################
#####   PCOA WITH BRAY-CURTIS   #####
#####################################

# format feature table and metadata
feat_otu <- t(feat_genus) # transpose genus feature table
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
meta$cluster <- as.factor(cluster_assign) # add cluster assignments
sampledata <- sample_data(meta) # convert to sample_data

ps <- phyloseq(feat_otu, sampledata) # create phyloseq object

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x)) # relative abundance for bray-curtis
bray_dist <- phyloseq::distance(ps_rel, method = "bray") # bray-curtis

# PCoA ordination
ordination_pcoa_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)

# PCoA plot by condition
plot_ordination(ps_rel, ordination_pcoa_bray, color = "condition") +
  geom_point(size = 2) + ggtitle("PCoA - Bray-Curtis") + theme_minimal()
adonis2(bray_dist ~ condition, data = meta) # PERMANOVA condition
betadisper(bray_dist, group = meta$condition) %>% permutest() # betadispar condition
boxplot(betadisper(bray_dist, group = meta$condition)) # boxplot condition
sil_con <- silhouette(as.numeric(as.factor(meta$condition)), dist = bray_dist)
plot(sil_con, border = NA, main = "Silhouette Plot (Bray-Curtis)") # silhouette plot condition

# PCoA plot by cluster
plot_ordination(ps_rel, ordination_pcoa_bray, color = "cluster") +
  geom_point(size = 2) + ggtitle("PCoA - Bray-Curtis") + theme_minimal()
adonis2(bray_dist ~ cluster, data = meta) # PERMANOVA cluster
betadisper(bray_dist, group = meta$cluster) %>% permutest() # betadispar cluster
boxplot(betadisper(bray_dist, group = meta$cluster)) # boxplot cluster
sil_clus <- silhouette(cluster_assign, dist = bray_dist)
plot(sil_clus, border = NA, main = "Silhouette Plot (Bray-Curtis)") # silhouette plot cluster


##########################################
#####   HEATMAP OF GENUS ABUNDANCE   #####
##########################################

### heatmap of genus abundances by cluster
feat_genus_df <- as.data.frame(feat_genus)
feat_genus_df$sample <- rownames(feat_genus_df)
feat_genus_df$cluster <- cluster_ids$cluster # add cluster ids to genus feature table

# mean abundance of genus per cluster
mean_abundance <- feat_genus_df %>%
  group_by(cluster) %>%
  summarise(across(-sample, mean), .groups = "drop") %>%
  column_to_rownames("cluster")
mean_abundance <- as.data.frame(t(mean_abundance))

# subset to top 25 most variable genera
taxa_variance <- apply(mean_abundance, 1, var) # calculate variance for each genus across clusters
top25_taxa <- names(sort(taxa_variance, decreasing = TRUE))[1:25] # top 25 most variable genera
mean_abundance_top25 <- mean_abundance[top25_taxa, ] # subset mean_abundance to top 25 genera
mean_abundance_scaled <- t(scale(t(mean_abundance_top25))) # scale data (z-score scaling)

# plot heatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
pheatmap(mean_abundance_scaled, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Top 25 most variable taxa",
         fontsize_row = 8,
         fontsize_col = 10,
         border_color = NA,
         color = my_palette,
         clustering_method = "ward.D2")


#######################################################
#####   DESEQ2 - DIFFERENTIAL ABUNDANCE TESTING   #####
#######################################################

# transpose feat_genus and convert to integer matrix
count_matrix <- round(as.matrix(t(feat_genus)))

# add cluster assignment to metadata
meta$cluster <- as.factor(cluster_assign)
table(meta$cluster)

# ensure sample names are the same
all(colnames(count_matrix) == rownames(meta))


# create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta,
                              design = ~ cluster)
  
# run DESeq2 analysis
dds <- DESeq(dds)

# pairwise comparisons and shrink log2fold change estimates using the ashr method
res_1_2 <- lfcShrink(dds, contrast = c("cluster", "1", "2"), type = "ashr")
res_1_3 <- lfcShrink(dds, contrast = c("cluster", "1", "3"), type = "ashr")
res_2_3 <- lfcShrink(dds, contrast = c("cluster", "2", "3"), type = "ashr")

# filter for significant results
res_1_2_sig <- res_1_2[which(res_1_2$padj < 0.05), ] %>% as.data.frame() %>% arrange(log2FoldChange)
res_1_3_sig <- res_1_3[which(res_1_3$padj < 0.05), ] %>% as.data.frame() %>% arrange(log2FoldChange)
res_2_3_sig <- res_2_3[which(res_2_3$padj < 0.05), ] %>% as.data.frame() %>% arrange(log2FoldChange)

# increased abundance in cluster one
inc_cluster_one <- intersect(rownames(res_1_2_sig[which(res_1_2_sig$log2FoldChange > 0), ]), 
                             rownames(res_1_3_sig[which(res_1_3_sig$log2FoldChange > 0), ]))

# increased abundance in cluster two
inc_cluster_two <- intersect(rownames(res_1_2_sig[which(res_1_2_sig$log2FoldChange < 0), ]), 
                             rownames(res_2_3_sig[which(res_2_3_sig$log2FoldChange > 0), ]))

# increased abundance in cluster three
inc_cluster_three <- intersect(rownames(res_1_3_sig[which(res_1_3_sig$log2FoldChange < 0), ]), 
                               rownames(res_2_3_sig[which(res_2_3_sig$log2FoldChange < 0), ]))


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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ashr_2.2-67                 DESeq2_1.48.2               SummarizedExperiment_1.38.1
# [4] Biobase_2.68.0              MatrixGenerics_1.20.0       matrixStats_1.5.0          
# [7] GenomicRanges_1.60.0        GenomeInfoDb_1.44.3         RColorBrewer_1.1-3         
# [10] pheatmap_1.0.13             cluster_2.1.8.1             phyloseq_1.52.0            
# [13] vegan_2.7-1                 permute_0.9-8               DirichletMultinomial_1.50.0
# [16] IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0        
# [19] generics_0.1.4              lubridate_1.9.4             forcats_1.0.1              
# [22] stringr_1.5.2               dplyr_1.1.4                 purrr_1.1.0                
# [25] readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0               
# [28] tidyverse_2.0.0             ggplot2_4.0.0              
# 
# loaded via a namespace (and not attached):
#   [1] ade4_1.7-23             tidyselect_1.2.1        farver_2.1.2            Biostrings_2.76.0      
# [5] S7_0.2.0                digest_0.6.37           timechange_0.3.0        lifecycle_1.0.4        
# [9] survival_3.8-3          invgamma_1.2            magrittr_2.0.4          compiler_4.5.0         
# [13] rlang_1.1.6             tools_4.5.0             igraph_2.1.4            data.table_1.17.8      
# [17] labeling_0.4.3          S4Arrays_1.8.1          DelayedArray_0.34.1     plyr_1.8.9             
# [21] BiocParallel_1.42.2     abind_1.4-8             withr_3.0.2             grid_4.5.0             
# [25] multtest_2.64.0         biomformat_1.36.0       Rhdf5lib_1.30.0         scales_1.4.0           
# [29] iterators_1.0.14        MASS_7.3-65             dichromat_2.0-0.1       cli_3.6.5              
# [33] crayon_1.5.3            rstudioapi_0.17.1       httr_1.4.7              reshape2_1.4.4         
# [37] tzdb_0.5.0              ape_5.8-1               rhdf5_2.52.1            splines_4.5.0          
# [41] XVector_0.48.0          vctrs_0.6.5             Matrix_1.7-4            jsonlite_2.0.0         
# [45] hms_1.1.3               mixsqp_0.3-54           irlba_2.3.5.1           locfit_1.5-9.12        
# [49] foreach_1.5.2           glue_1.8.0              codetools_0.2-20        stringi_1.8.7          
# [53] gtable_0.3.6            UCSC.utils_1.4.0        pillar_1.11.1           rhdf5filters_1.20.0    
# [57] truncnorm_1.0-9         GenomeInfoDbData_1.2.14 R6_2.6.1                lattice_0.22-7         
# [61] SQUAREM_2021.1          Rcpp_1.1.0              SparseArray_1.8.1       nlme_3.1-168           
# [65] mgcv_1.9-3              pkgconfig_2.0.3 

