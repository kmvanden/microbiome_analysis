### Longitudinal Analysis of Gavage Data

# load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(nlme)
library(microbiome)
library(pairwiseAdonis)
library(scales)
library(reshape2)
library(pheatmap)
library(grid)
library(colorspace)
library(Maaslin2)
library(mgcv)
library(MetaLonDA)
library(LinDA)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load data
# metadata
meta <- read.table("metadata_gavage.txt", header = TRUE)
rownames(meta) <- meta$sample_id
meta$day_factor <- factor(meta$day) # add day as a factor for plotting

# feature table
feat <- read.table("otu_table_gavage.txt", header = TRUE)


########################################
#####   ALPHA DIVERSITY ANALYSIS   #####
########################################

# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)

# rarefy the phyloseq object (for Observed richness)
sample_sums(ps) %>% summary() # check sequencing depth (rarefaction depth based on sequencing depth)
ps_rarefied <- rarefy_even_depth(ps, rngseed = 1234, sample.size = min(sample_sums(ps)), verbose = FALSE)

# relative abundance for Shannon and Simpson
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x)) 

# calculate alpha diversity using phyloseq
alpha_phyloseq <- estimate_richness(ps_rel, measures = c("Shannon", "Simpson"))
chao_df <- estimate_richness(ps, measures = "Chao1")
observed_df <- estimate_richness(ps_rarefied, measures = "Observed")
alpha_phyloseq$Chao1 <- chao_df$Chao1 # add Chao1 to alpha_phyloseq
alpha_phyloseq$Observed <- observed_df$Observed # add Observed to alpha_phyloseq
alpha_phyloseq$sample_name <- rownames(alpha_phyloseq) # add column with sample names
alpha_phyloseq <- left_join(alpha_phyloseq, meta, by = c("sample_name" = "sample_id")) # merge phyloseq object with metadata (for plotting and statistical tests) 


### plot alpha diversity metrics
# Observed richness
ggplot(alpha_phyloseq, aes(x = day, y = Observed, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Observed richness", y = "Number of unique taxa", x = "Day")

# Shannon diversity
ggplot(alpha_phyloseq, aes(x = day, y = Shannon, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Shannon diversity", y = "Shannon index", x = "Day")

# Simpson diversity
ggplot(alpha_phyloseq, aes(x = day, y = Simpson, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Simpson diversity", y = "Simpson index", x = "Day")

# Chao1 diversity
ggplot(alpha_phyloseq, aes(x = day, y = Chao1, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Chao1 estimator", y = "Chao1 estimator", x = "Day")


### compute statistics (linear mixed-effects model + post-hoc pairwise comparisons using t-tests on estimated marginal means)
# Observed richness
lme_Obs_f <- lmer(Observed ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Obs_f)
anova(lme_Obs_f)
emm_Obs_f <- emmeans(lme_Obs_f, ~ gavage | day_factor)
pairs(emm_Obs_f, adjust = "tukey")

lme_Obs_n <- lmer(Observed ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Obs_n)
anova(lme_Obs_n)

# Shannon diversity
lme_Shan_f <- lmer(Shannon ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Shan_f)
anova(lme_Shan_f)
emm_Shan_f <- emmeans(lme_Shan_f, ~ gavage | day_factor)
pairs(emm_Shan_f, adjust = "tukey")

lme_Shan_n <- lmer(Shannon ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Shan_n)
anova(lme_Shan_n)

# Simpson diversity
lme_Simp_f <- lmer(Simpson ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Simp_f)
anova(lme_Simp_f)
emm_Simp_f <- emmeans(lme_Simp_f, ~ gavage | day_factor)
pairs(emm_Simp_f, adjust = "tukey")

lme_Simp_n <- lmer(Simpson ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Simp_n)
anova(lme_Simp_n)

# Chao1 diversity
lme_Chao_f <- lmer(Chao1 ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Chao_f)
anova(lme_Chao_f)
emm_Chao_f <- emmeans(lme_Chao_f, ~ gavage | day_factor)
pairs(emm_Chao_f, adjust = "tukey")

lme_Chao_n <- lmer(Chao1 ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq)
summary(lme_Chao_n)
anova(lme_Chao_n)


######################################################################
#####   BETA DIVERSITY ANALYSIS - PREVALENCE FILTERING OF TAXA   #####
######################################################################

# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)

# filter low prevalence taxa (present in less than 10% of samples)
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


### statistical testing with PERMANOVA
# function for pairwise by timepoint
pairwise_by_timepoint <- function(dist_matrix, meta, group_var, time_var, 
                                  perm = 999, p_adjust_method = "fdr") {
  
  # ensure the distance matrix is dist obj
  if(!inherits(dist_matrix, "dist")) {
    dist_matrix <- as.dist(dist_matrix)
  }
  
  time_points <- unique(meta[[time_var]]) # get unique timepoints
  
  results_list <- list() # list to store results
  
  for(t in time_points){
    
    # subset by time
    idx <- meta[[time_var]] == t
    dist_sub <- as.dist(as.matrix(dist_matrix)[idx, idx])
    meta_sub <- meta[idx, ]
    
    # run pairwise.adonis
    res <- pairwise.adonis(dist_sub, factors = meta_sub[[group_var]], perm = perm)
    res$Time <- t  # add time information
    
    # store results
    results_list[[as.character(t)]] <- res
  }
  
  combined <- do.call(rbind, results_list) # combine into a single data.frame
  combined$p.adjusted_global <- p.adjust(combined$p.value, method = p_adjust_method) # apply global FDR across all comparisons
  
  return(combined)
}

# bray-curtis
adonis2(bray_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
pairwise_by_timepoint(dist_matrix = bray_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") # pairwise by group per time point

# jaccard
adonis2(jaccard_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
pairwise_by_timepoint(dist_matrix = jaccard_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") # pairwise by group per time point

# euclidean
adonis2(euc_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
pairwise_by_timepoint(dist_matrix = euc_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") # pairwise by group per time point

# canberra
adonis2(canberra_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
pairwise_by_timepoint(dist_matrix = canberra_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") # pairwise by group per time point

# aitchison
adonis2(aitchison_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
pairwise_by_timepoint(dist_matrix = aitchison_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") # pairwise by group per time point


### PCoA ordination
ordination_pcoa_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist) # bray-curtis
ordination_pcoa_jaccard <- ordinate(ps_pa, method = "PCoA", distance = jaccard_dist) # jaccard
ordination_pcoa_euc <- ordinate(ps_log, method = "PCoA", distance = euc_dist) # euclidean
ordination_pcoa_canberra <- ordinate(ps_rel, method = "PCoA", distance = canberra_dist) # canberra
ordination_pcoa_aitchison <- ordinate(ps_clr, method = "PCoA", distance = aitchison_dist) # aitchison


### PCoA plots
plot_ordination(ps_rel, ordination_pcoa_bray, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("PCoA - Bray-Curtis") + theme_minimal() # bray-curtis

plot_ordination(ps_rel, ordination_pcoa_canberra, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("PCoA - Canberra") + theme_minimal() # canberra

plot_ordination(ps_pa, ordination_pcoa_jaccard, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("PCoA - Jaccard") + theme_minimal() # jaccard

plot_ordination(ps_log, ordination_pcoa_euc, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("PCoA - Euclidean") + theme_minimal() # euclidean

plot_ordination(ps_clr, ordination_pcoa_aitchison, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("PCoA - Aitchison") + theme_minimal() # aitchison


### PCoA plot using ggplot for more control over plotting
df_bray <- as.data.frame(ordination_pcoa_bray$vectors) # extract coordinates
df_bray$sample_id <- rownames(df_bray)
df_bray <- left_join(df_bray, meta, by = "sample_id") # join with metadata

var_explained <- ordination_pcoa_bray$values$Relative_eig * 100 # get variance explained

# compute centroids per gavage x day
df_centroids <- df_bray %>%
  group_by(gavage, day_factor) %>%
  summarize(mean_Axis1 = mean(Axis.1),
            mean_Axis2 = mean(Axis.2),
            .groups = "drop") %>%
  arrange(gavage, as.numeric(as.character(day_factor)))

ggplot(df_bray, aes(x = Axis.1, y = Axis.2)) +
  geom_path(data = df_centroids, 
            aes(x = mean_Axis1, y = mean_Axis2, group = gavage, color = gavage), 
            arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, color = gavage, size = day_factor), 
             shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = gavage, size = day_factor), 
             shape = 21, stroke = 0, alpha = 0.35) +
  stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_fill_discrete(guide = "none") + 
  scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3,"28" = 4)) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  ggtitle("PCoA – Bray-Curtis") + theme_minimal()


### distance-based redundancy analysis (db-RDA) - constrained PCoA
# fits multivariate linear model in distance space

# set permutation scheme within mice
perm_mouse <- how(blocks = meta$mouse_id, nperm = 999) 

# bray-curtis
cap_bray <- capscale(bray_dist ~ gavage + day, data = meta) # db-RDA model
anova(cap_bray, permutations = perm_mouse) # gavage and day together
anova(cap_bray, by = "terms", permutations = perm_mouse) # gavage or day
anova(cap_bray, by = "axis", permutations = perm_mouse) # axis tests

cap_scores <- scores(cap_bray, display = "sites") # get CAP coordinates
cap_df <- as.data.frame(cap_scores)
cap_df$sample_id <- rownames(cap_df)
cap_df <- left_join(cap_df, meta, by = "sample_id")

var_explained <- cap_bray$CCA$eig/sum(cap_bray$CCA$eig) * 100 # get variance explained

# compute centroids
df_centroids <- cap_df %>%
  group_by(gavage, day_factor) %>%
  summarize(mean_CAP1 = mean(CAP1),
            mean_CAP2 = mean(CAP2),
            .groups = "drop") %>%
  arrange(gavage, as.numeric(as.character(day_factor)))

ggplot(cap_df, aes(x = CAP1, y = CAP2)) +
  geom_path(data = df_centroids, 
            aes(x = mean_CAP1, y = mean_CAP2, group = gavage, color = gavage),
            arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
  geom_point(aes(x = CAP1, y = CAP2, color = gavage, size = day_factor), 
             shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
  geom_point(aes(x = CAP1, y = CAP2, fill = gavage, size = day_factor), 
             shape = 21, stroke = 0, alpha = 0.35) +
  stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_fill_discrete(guide = "none") +
  scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3, "28" = 4)) +
  xlab(paste0("CAP1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("CAP2 (", round(var_explained[2], 1), "%)")) +
  ggtitle("db-RDA - Bray-Curtis") + theme_minimal()

# jaccard
cap_jaccard <- capscale(jaccard_dist ~ gavage + day, data = meta) # db-RDA model
anova(cap_jaccard, permutations = perm_mouse) # gavage and day together
anova(cap_jaccard, by = "terms", permutations = perm_mouse) # gavage or day
anova(cap_jaccard, by = "axis", permutations = perm_mouse) # axis tests

# canberra
cap_canberra <- capscale(canberra_dist ~ gavage + day, data = meta) # db-RDA model
anova(cap_canberra, permutations = perm_mouse) # gavage and day together
anova(cap_canberra, by = "terms", permutations = perm_mouse) # gavage or day
anova(cap_canberra, by = "axis", permutations = perm_mouse) # axis tests

# euclidean/aitchison
clr_mat <- t(otu_table(ps_clr))  # clr-transformed feature table
rda_model <- rda(clr_mat ~ gavage + day, data = meta) # RDA model

anova(rda_model, permutations = perm_mouse) # gavage and day together
anova(rda_model, by = "terms", permutations = perm_mouse) # gavage or day
anova(rda_model, by = "axis", permutations = perm_mouse) # axis tests


### statistical testing using a distance-based linear mixed model (dbLMM)
# extract coordinates for Axis.1
dbLMM_df <- data.frame(pcoa_bray = ordination_pcoa_bray$vectors[,1],
                       pcoa_canberra = ordination_pcoa_canberra$vectors[,1],
                       pcoa_jaccard = ordination_pcoa_jaccard$vectors[,1],
                       pcoa_euc = ordination_pcoa_euc$vectors[,1],
                       pcoa_aitchison = ordination_pcoa_aitchison$vectors[,1])
dbLMM_df$sample_id <- rownames(ordination_pcoa_bray$vectors)
dbLMM_df <- left_join(dbLMM_df, meta, by = "sample_id")  # join metadata

# bray-curtis
dblmm_bray_f <- lmer(pcoa_bray ~ gavage * day_factor + (1|mouse_id), data = dbLMM_df)
summary(dblmm_bray_f)
anova(dblmm_bray_f)
emm_bray_f <- emmeans(dblmm_bray_f, ~ gavage | day_factor)
pairs(emm_bray_f, adjust = "tukey")

dblmm_bray_n <- lmer(pcoa_bray ~ gavage * day + (1|mouse_id), data = dbLMM_df)
summary(dblmm_bray_n)
anova(dblmm_bray_n)

# canberra
dblmm_can_f <- lmer(pcoa_canberra ~ gavage * day_factor + (1|mouse_id), data = dbLMM_df)
summary(dblmm_can_f)
anova(dblmm_can_f)
emm_can_f <- emmeans(dblmm_can_f, ~ gavage | day_factor)
pairs(emm_can_f, adjust = "tukey")

dblmm_can_n <- lmer(pcoa_canberra ~ gavage * day + (1|mouse_id), data = dbLMM_df)
summary(dblmm_can_n)
anova(dblmm_can_n)

# jaccard
dblmm_jacc_f <- lmer(pcoa_jaccard ~ gavage * day_factor + (1|mouse_id), data = dbLMM_df)
summary(dblmm_jacc_f)
anova(dblmm_jacc_f)
emm_jacc_f <- emmeans(dblmm_jacc_f, ~ gavage | day_factor)
pairs(emm_jacc_f, adjust = "tukey")

dblmm_jacc_n <- lmer(pcoa_jaccard ~ gavage * day + (1|mouse_id), data = dbLMM_df)
summary(dblmm_jacc_n)
anova(dblmm_jacc_n)

# euclidean
dblmm_euc_f <- lmer(pcoa_euc ~ gavage * day_factor + (1|mouse_id), data = dbLMM_df)
summary(dblmm_euc_f)
anova(dblmm_euc_f)
emm_euc_f <- emmeans(dblmm_euc_f, ~ gavage | day_factor)
pairs(emm_euc_f, adjust = "tukey")

dblmm_euc_n <- lmer(pcoa_euc ~ gavage * day + (1|mouse_id), data = dbLMM_df)
summary(dblmm_euc_n)
anova(dblmm_euc_n)

# aitchison
dblmm_aitch_f <- lmer(pcoa_aitchison ~ gavage * day_factor + (1|mouse_id), data = dbLMM_df)
summary(dblmm_aitch_f)
anova(dblmm_aitch_f)
emm_aitch_f <- emmeans(dblmm_aitch_f, ~ gavage | day_factor)
pairs(emm_aitch_f, adjust = "tukey")

dblmm_aitch_n <- lmer(pcoa_aitchison ~ gavage * day + (1|mouse_id), data = dbLMM_df)
summary(dblmm_aitch_n)
anova(dblmm_aitch_n)


### NMDS ordination
ordination_nmds_bray <- ordinate(ps_rel, method = "NMDS", distance = bray_dist) # bray-curtis
ordination_nmds_jaccard <- ordinate(ps_pa, method = "NMDS", distance = jaccard_dist) # jaccard
ordination_nmds_canberra <- ordinate(ps_rel, method = "NMDS", distance = canberra_dist) # canberra


### NMDS plots
plot_ordination(ps_rel, ordination_nmds_bray, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("NMDS - Bray-Curtis") + theme_minimal() # bray-curtis

plot_ordination(ps_rel, ordination_nmds_jaccard, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("NMDS - Jaccard") + theme_minimal() # jaccard

plot_ordination(ps_rel, ordination_nmds_canberra, color = "gavage", shape = "day_factor") + 
  stat_ellipse(aes(group = gavage), type = "norm", level = 0.95) + geom_point(size = 1.5) + 
  ggtitle("NMDS - Canberra") + theme_minimal() # canberra


### NMDS plot using ggplot for more control over plotting
df_nmds <- as.data.frame(ordination_nmds_bray$points) # extract coordinates
df_nmds$sample_id <- rownames(df_nmds)
df_nmds <- left_join(df_nmds, meta, by = "sample_id")  # join with metadata

# compute centroids per gavage x day
df_centroids <- df_nmds %>%
  group_by(gavage, day_factor) %>%
  summarize(mean_MDS1 = mean(MDS1),
            mean_MDS2 = mean(MDS2),
            .groups = "drop") %>%
  arrange(gavage, as.numeric(as.character(day_factor)))

ggplot(df_nmds, aes(x = MDS1, y = MDS2)) +
  geom_path(data = df_centroids, 
            aes(x = mean_MDS1, y = mean_MDS2, group = gavage, color = gavage), 
            arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
  geom_point(aes(color = gavage, size = day_factor), 
             shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
  geom_point(aes(fill = gavage, size = day_factor), 
             shape = 21, stroke = 0, alpha = 0.35) +
  stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_fill_discrete(guide = "none") + 
  scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3,"28" = 4)) +
  xlab("NMDS1") +
  ylab("NMSD2") +
  ggtitle("NMDS - Bray-Curtis") + theme_minimal()


### extract PCA data and plot - euclidean
log_otu <- t(otu_table(ps_log)) # extract CLR-transformed data and transpose
pca_euc <- prcomp(log_otu, center = TRUE, scale. = FALSE) # run PCA
pca_euc_df <- as.data.frame(pca_euc$x) # extract PCA scores (coordinates)
pca_euc_df$sample_name <- rownames(pca_euc_df)
pca_euc_df <- left_join(pca_euc_df, meta, by = c("sample_name" = "sample_id"))  # join with metadata
euc_var_explained <- round(100 * summary(pca_euc)$importance[2, 1:2], 1) # extract percentage of variance explained

# compute centroids per gavage x day
df_centroids <- pca_euc_df %>%
  group_by(gavage, day_factor) %>%
  summarize(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            .groups = "drop") %>%
  arrange(gavage, as.numeric(as.character(day_factor)))

ggplot(pca_euc_df, aes(x = PC1, y = PC2)) +
  geom_path(data = df_centroids, 
            aes(x = mean_PC1, y = mean_PC2, group = gavage, color = gavage), 
            arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
  geom_point(aes(color = gavage, size = day_factor), 
             shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
  geom_point(aes(fill = gavage, size = day_factor), 
             shape = 21, stroke = 0, alpha = 0.35) +
  stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_fill_discrete(guide = "none") + 
  scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3,"28" = 4)) +
  xlab(paste0("PC1 (", euc_var_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", euc_var_explained[2], "% variance)")) +
  ggtitle("PCA of log-transformed data (euclidean)") + theme_minimal()


### extract PCA data and plot
# aitchison
clr_otu <- t(otu_table(ps_clr)) # extract CLR-transformed data and transpose
pca_ait <- prcomp(clr_otu, center = TRUE, scale. = FALSE) # run PCA
pca_ait_df <- as.data.frame(pca_ait$x) # extract PCA scores (coordinates)
pca_ait_df$sample_name <- rownames(pca_ait_df)
pca_ait_df <- left_join(pca_ait_df, meta, c("sample_name" = "sample_id"))  # join with metadata
ait_var_explained <- round(100 * summary(pca_ait)$importance[2, 1:2], 1) # extract percentage of variance explained

# compute centroids per gavage x day
df_centroids <- pca_ait_df %>%
  group_by(gavage, day_factor) %>%
  summarize(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            .groups = "drop") %>%
  arrange(gavage, as.numeric(as.character(day_factor)))

ggplot(pca_ait_df, aes(x = PC1, y = PC2)) +
  geom_path(data = df_centroids, 
            aes(x = mean_PC1, y = mean_PC2, group = gavage, color = gavage), 
            arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
  geom_point(aes(color = gavage, size = day_factor), 
             shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
  geom_point(aes(fill = gavage, size = day_factor), 
             shape = 21, stroke = 0, alpha = 0.35) +
  stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_fill_discrete(guide = "none") + 
  scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3,"28" = 4)) +
  xlab(paste0("PC1 (", ait_var_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", ait_var_explained[2], "% variance)")) +
  ggtitle("PCA of CLR-transformed data (aitchison)") + theme_minimal()


##########################################
#####   HEATMAPS OF TAXA ABUNDANCE   #####
##########################################

### plot heatmap with ggplot
# get top 25 taxa by abundance
top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:25]
ps_top <- prune_taxa(top_taxa, ps_rel)

df_heat <- psmelt(ps_top)  # long format data.frame for plotting

df_heat <- df_heat %>%
  group_by(OTU) %>%
  mutate(Abundance_scaled = scales::rescale(Abundance)) %>%
  ungroup() # scale abundance

# plot heatmap
ggplot(df_heat, aes(x = gavage_day, y = OTU, fill = Abundance_scaled)) +
  geom_tile() + scale_fill_viridis_c() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Heatmap of top 25 taxa")


### plot heatmap with pheatmap
mat <- as.data.frame(otu_table(ps_top))   # otu table
meta <- as.data.frame(sample_data(ps_top))  # metadata
mat_t <- as.data.frame(t(mat))  # transpose metadata 
mat_t$gavage_day <- meta$gavage_day # add metadata to otu table

mat_day <- mat_t %>%
  group_by(gavage_day) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame() # aggregate by gavage data

mat_day_t <- t(as.matrix(mat_day[,-1]))   # drop first column and transpose
colnames(mat_day_t) <- mat_day$gavage_day  # set to gavage_day

pheatmap(mat_day_t, scale = "row", # scale taxa abundance across samples
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         main = "Top 25 taxa (aggregated by gavage_day)")


##########################################################
#####   TRAJECTORY PLOT - TAXA ABUNDANCE OVER TIME   #####
##########################################################

# long format data.frame for plotting
df_taxa <- psmelt(ps_rel)  

# summarize relative abundance of taxa over time and gavage group
df_traj <- df_taxa %>%
  group_by(Day = day, Gavage = gavage, OTU) %>%
  summarize(mean_abundance = mean(Abundance), 
            sem = sd(Abundance)/sqrt(n()),
            .groups = "drop")

# get top 6 taxa by mean abundance
top6 <- df_traj %>%
  group_by(OTU) %>%
  summarize(total = sum(mean_abundance)) %>%
  top_n(6, total) %>%
  pull(OTU)

# plot abundance trajectory for top 6 taxa
ggplot(df_traj %>% filter(OTU %in% top6), 
       aes(x = Day, y = mean_abundance, color = Gavage, group = Gavage, fill = Gavage)) +
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  facet_wrap(~OTU, scales = "free_y") + theme_minimal() +
  labs(title = "Trajectory of top 6 species over time", 
       x = "Day", y = "Mean relative abundance")

# plot abundance over time for specific taxa
ggplot(df_traj %>% filter(OTU == "Akkermansia_muciniphila"), 
       aes(x = Day, y = mean_abundance, color = Gavage, group = Gavage, fill = Gavage)) +
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  facet_wrap(~OTU, scales = "free_y") + theme_minimal() +
  labs(title = "Trajectory of top 5 species over time", 
       x = "Day", y = "Mean relative abundance")

ggplot(df_traj %>% filter(OTU == "Dorea_formicigenerans"), 
       aes(x = Day, y = mean_abundance, color = Gavage, group = Gavage, fill = Gavage)) +
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  facet_wrap(~OTU, scales = "free_y") + theme_minimal() +
  labs(title = "Trajectory of top 5 species over time", 
       x = "Day", y = "Mean relative abundance")


################################
#####   STACKED BARPLOTS   #####
################################

# get top 15 taxa by abundance
top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:15]
ps_top <- prune_taxa(top_taxa, ps_rel)
df_stack <- psmelt(ps_rel)  # long format data.frame for plotting

# assign all none top_taxa as Other and order so Other is last
df_stack$OTU_mod <- ifelse(df_stack$OTU %in% top_taxa, as.character(df_stack$OTU), "Other")
df_stack$OTU_mod <- factor(df_stack$OTU_mod,
                           levels = c(setdiff(unique(df_stack$OTU_mod), "Other"), "Other"))

# color palette
colors <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#999999", "#AA4499", "#44AA99")
my_colors <- rep(colors, length.out = length(unique(df_stack$OTU_mod)))

ggplot(df_stack, aes(x = gavage_day, y = Abundance, fill = OTU_mod)) +
  geom_bar(stat = "identity") + theme_minimal() +
  scale_fill_manual(values = my_colors)+
  labs(title = "Stacked barplot of top 20 taxa over time", 
       x = "Gavage group and day", y = "Mean relative abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] LinDA_0.2.0          MetaLonDA_1.1.8      mgcv_1.9-3           Maaslin2_1.22.0      colorspace_2.1-2    
# [6] pheatmap_1.0.13      reshape2_1.4.4       scales_1.4.0         pairwiseAdonis_0.4.1 cluster_2.1.8.1     
# [11] microbiome_1.30.0    nlme_3.1-168         emmeans_1.11.2-8     lmerTest_3.1-3       lme4_1.1-37         
# [16] Matrix_1.7-4         vegan_2.7-1          permute_0.9-8        phyloseq_1.52.0      lubridate_1.9.4     
# [21] forcats_1.0.1        stringr_1.5.2        dplyr_1.1.4          purrr_1.1.0          readr_2.1.5         
# [26] tidyr_1.3.1          tibble_3.3.0         tidyverse_2.0.0      ggplot2_4.0.0       
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4            DBI_1.2.3               sandwich_3.1-1          rlang_1.1.6            
# [5] magrittr_2.0.4          clue_0.3-66             multcomp_1.4-28         ade4_1.7-23            
# [9] compiler_4.5.0          vctrs_0.6.5             rmutil_1.1.10           pkgconfig_2.0.3        
# [13] crayon_1.5.3            backports_1.5.0         XVector_0.48.0          labeling_0.4.3         
# [17] modeest_2.4.0           gss_2.2-9               pracma_2.4.4            tzdb_0.5.0             
# [21] UCSC.utils_1.4.0        nloptr_2.2.1            GenomeInfoDb_1.44.3     jsonlite_2.0.0         
# [25] biomformat_1.36.0       rhdf5filters_1.20.0     Rhdf5lib_1.30.0         broom_1.0.10           
# [29] parallel_4.5.0          R6_2.6.1                biglm_0.9-3             stringi_1.8.7          
# [33] RColorBrewer_1.1-3      rpart_4.1.24            boot_1.3-32             numDeriv_2016.8-1.1    
# [37] estimability_1.5.1      Rcpp_1.1.0              iterators_1.0.14        zoo_1.8-14             
# [41] IRanges_2.42.0          splines_4.5.0           igraph_2.1.4            timechange_0.3.0       
# [45] tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1       timeDate_4041.110      
# [49] doParallel_1.0.17       codetools_0.2-20        lattice_0.22-7          plyr_1.8.9             
# [53] Biobase_2.68.0          withr_3.0.2             S7_0.2.0                stable_1.1.6           
# [57] coda_0.19-4.1           Rtsne_0.17              survival_3.8-3          getopt_1.20.4          
# [61] Biostrings_2.76.0       pillar_1.11.1           foreach_1.5.2           stats4_4.5.0           
# [65] reformulas_0.4.1        pcaPP_2.0-5             generics_0.1.4          S4Vectors_0.48.0       
# [69] hms_1.1.3               timeSeries_4041.111     minqa_1.2.8             xtable_1.8-4           
# [73] glue_1.8.0              statip_0.2.3            tools_4.5.0             robustbase_0.99-6      
# [77] data.table_1.17.8       spatial_7.3-18          fBasics_4041.97         mvtnorm_1.3-3          
# [81] rhdf5_2.52.1            optparse_1.7.5          ape_5.8-1               rbibutils_2.3          
# [85] GenomeInfoDbData_1.2.14 cli_3.6.5               viridisLite_0.4.2       gtable_0.3.6           
# [89] DEoptimR_1.1-4          stabledist_0.7-2        digest_0.6.37           BiocGenerics_0.54.0    
# [93] pbkrtest_0.5.5          ggrepel_0.9.6           TH.data_1.1-4           farver_2.1.2           
# [97] multtest_2.64.0         lifecycle_1.0.4         httr_1.4.7              MASS_7.3-65

