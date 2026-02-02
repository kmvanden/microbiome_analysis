### Longitudinal Analysis of Gavage Data

# load libraries
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(purrr)
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(breakaway)
library(DivNet)
library(doParallel)
library(nlme)
library(microbiome)
library(pairwiseAdonis)
library(scales)
library(reshape2)
library(pheatmap)
library(grid)
library(colorspace)
library(Maaslin2)
library(corncob)
library(compositions)
library(mgcv)
library(LinDA)
library(glmmTMB)
library(DHARMa)
library(gratia)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load data
# metadata
meta <- read.table("metadata_gavage.txt", header = TRUE)
rownames(meta) <- meta$sample_id
meta$day_factor <- factor(meta$day) # add day as a factor for plotting

# feature table
feat <- read.table("otu_table_gavage.txt", header = TRUE)


################################################################
#####   ALPHA DIVERSITY ANALYSIS - COMPUTATION AND PLOTS   #####
################################################################

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


### calculate traditional alpha diversity metrics using phyloseq
alpha_phyloseq <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
chao_df <- estimate_richness(ps, measures = "Chao1")
observed_df <- estimate_richness(ps_rarefied, measures = "Observed")
alpha_phyloseq$Chao1 <- chao_df$Chao1 # add Chao1 to alpha_phyloseq
alpha_phyloseq$Observed <- observed_df$Observed # add Observed to alpha_phyloseq
alpha_phyloseq$sample_name <- rownames(alpha_phyloseq) # add column with sample names
alpha_phyloseq <- left_join(alpha_phyloseq, meta, by = c("sample_name" = "sample_id")) # merge phyloseq object with metadata (for plotting and statistical tests) 


### calculate breakaway richness
# transpose feature table
feat_t <- as.data.frame(t(feat))

# breakaway richness per sample (breakaway expects a count vector per sample)
breakaway_richness_df <- purrr::map_dfr(rownames(feat_t), function(s) {
  ba <- breakaway(as.numeric(feat_t[s, ]))
  data.frame(sample_id = s,
             estimate = ba$estimate,
             se = ba$se)
}) %>% left_join(meta, by = "sample_id") # merge with metadata


### calculate Shannon and Simpson diversity using DivNet
divnet_diversity <- divnet(ps, ncores = parallel::detectCores() - 1)

# convert to data.frame
divnet_alpha_diversity <- data.frame(sample_id = names(divnet_diversity$shannon),
                                     shannon_estimate = sapply(divnet_diversity$shannon, function(x) x$estimate),
                                     shannon_error = sapply(divnet_diversity$shannon, function(x) x$error),
                                     simpson_estimate = sapply(divnet_diversity$simpson, function(x) x$estimate),
                                     simpson_error = sapply(divnet_diversity$simpson, function(x) x$error))

# merge with metadata
divnet_alpha_diversity_df <- divnet_alpha_diversity %>%
  left_join(meta, by = "sample_id")


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
  labs(title = "Chao1 estimator", y = "Chao1 richness estimate", x = "Day")

# breakaway richness
ggplot(breakaway_richness_df, aes(x = day, y = estimate, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Richness (breakaway)", y = "Estimated total taxa", x = "Day")

# DivNet Shannon diversity
ggplot(divnet_alpha_diversity_df, aes(x = day, y = shannon_estimate, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Shannon diversity (DivNet)", y = "Shannon index", x = "Day")

# DivNet Simpson diversity
ggplot(divnet_alpha_diversity_df, aes(x = day, y = simpson_estimate, color = gavage)) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
  stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
  labs(title = "Simpson diversity (DivNet)", y = "Simpson index", x = "Day")


######################################################################
#####   ALPHA DIVERSITY ANALYSIS - LINEAR MIXED-EFFECTS MODELS   #####
######################################################################

### compute statistics (linear mixed-effects models)
### Observed richness
lme_Obs_n <- lmer(Observed ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq) # numeric (linear dynamics)
summary(lme_Obs_n)
anova(lme_Obs_n)

lme_Obs_f <- lmer(Observed ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq) # factor
summary(lme_Obs_f)
anova(lme_Obs_f)

plot(resid(lme_Obs_f) ~ fitted(lme_Obs_f)) # check homoscedasticity 
qqnorm(resid(lme_Obs_f)); qqline(resid(lme_Obs_f)) # check normality

# pairwise contrasts within each day
emm_Obs_f <- emmeans(lme_Obs_f, ~ gavage | day_factor)
pairs(emm_Obs_f, adjust = "tukey")


### Shannon diversity
lme_Shan_n <- lmer(Shannon ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq) # numeric (linear dynamics)
anova(lme_Shan_n)
summary(lme_Shan_n)

lme_Shan_f <- lmer(Shannon ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq) # factor
anova(lme_Shan_f)
summary(lme_Shan_f)

# pairwise contrasts within each day
emm_Shan_f <- emmeans(lme_Shan_f, ~ gavage | day_factor) 
pairs(emm_Shan_f, adjust = "tukey")


### Simpson diversity
lme_Simp_n <- lmer(Simpson ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq) # numeric (linear dynamics)
summary(lme_Simp_n)
anova(lme_Simp_n)

lme_Simp_f <- lmer(Simpson ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq) # factor
summary(lme_Simp_f)
anova(lme_Simp_f)

# pairwise contrasts within each day
emm_Simp_f <- emmeans(lme_Simp_f, ~ gavage | day_factor)
pairs(emm_Simp_f, adjust = "tukey")


### Chao1 diversity
lme_Chao_n <- lmer(Chao1 ~ gavage * day + (1 | mouse_id), data = alpha_phyloseq) # numeric
summary(lme_Chao_n)
anova(lme_Chao_n)

lme_Chao_f <- lmer(Chao1 ~ gavage * day_factor + (1 | mouse_id), data = alpha_phyloseq) # factor
summary(lme_Chao_f)
anova(lme_Chao_f)

plot(resid(lme_Chao_f) ~ fitted(lme_Chao_f)) # check homoscedasticity 
qqnorm(resid(lme_Chao_f)); qqline(resid(lme_Chao_f)) # check normality

# pairwise contrasts within each day
emm_Chao_f <- emmeans(lme_Chao_f, ~ gavage | day_factor)
pairs(emm_Chao_f, adjust = "tukey")


### richness (breakaway)
lme_b_Obs_n <- lmer(estimate ~ gavage * day + (1 | mouse_id), data = breakaway_richness_df, # numeric (linear dynamics)
                    weights = 1 / se^2) # inverse-variance weighting
anova(lme_b_Obs_n)
summary(lme_b_Obs_n)

lme_b_Obs_f <- lmer(estimate ~ gavage * day_factor + (1 | mouse_id), data = breakaway_richness_df, # factor
                    weights = 1 / se^2) # inverse-variance weighting
anova(lme_b_Obs_f)
summary(lme_b_Obs_f)

plot(resid(lme_b_Obs_f) ~ fitted(lme_b_Obs_f)) # check homoscedasticity 
qqnorm(resid(lme_b_Obs_f)); qqline(resid(lme_b_Obs_f)) # check normality

# pairwise contrasts within each day
emm_b_Obs_f <- emmeans(lme_b_Obs_f, ~ gavage | day_factor)
pairs(emm_b_Obs_f, adjust = "tukey")


### Shannon diversity (DivNet)
lme_dn_Shan_n <- lmer(shannon_estimate ~ gavage * day + (1 | mouse_id), data = divnet_alpha_diversity_df, # numeric (linear dynamics)
                      weights = 1 / (shannon_error^2)) # inverse-variance weighting
anova(lme_dn_Shan_n)
summary(lme_dn_Shan_n)

lme_dn_Shan_f <- lmer(shannon_estimate ~ gavage * day_factor + (1 | mouse_id), data = divnet_alpha_diversity_df,# factor
                      weights = 1 / (shannon_error^2)) # inverse-variance weighting
anova(lme_dn_Shan_f)
summary(lme_dn_Shan_f)

# pairwise contrasts within each day
emm_dn_Shan_f <- emmeans(lme_dn_Shan_f, ~ gavage | day_factor)
pairs(emm_dn_Shan_f, adjust = "tukey")


### Simpson diversity (DivNet)
lme_dn_Simp_n <- lmer(simpson_estimate ~ gavage * day + (1 | mouse_id), data = divnet_alpha_diversity_df, # numeric (linear dynamics)
                   weights = 1 / (simpson_error^2)) # inverse-variance weighting
anova(lme_dn_Simp_n)
summary(lme_dn_Simp_n)

lme_dn_Simp_f <- lmer(simpson_estimate ~ gavage * day_factor + (1 | mouse_id), data = divnet_alpha_diversity_df, # factor
                   weights = 1 / (simpson_error^2)) # inverse-variance weighting
anova(lme_dn_Simp_f)
summary(lme_dn_Simp_f)

# pairwise contrasts within each day
emm_dn_Simp_f <- emmeans(lme_dn_Simp_f, ~ gavage | day_factor)
pairs(emm_dn_Simp_f, adjust = "tukey")


############################################################################
#####   ALPHA DIVERSITY ANALYSIS - GENERALIZED ADDITIVE MIXED MODELS   #####
############################################################################

### compute statistics (generalized additive mixed models)
alpha_phyloseq$gavage <- as.factor(alpha_phyloseq$gavage) 
alpha_phyloseq <- alpha_phyloseq %>% mutate(day_c = day - mean(day)) # center day (day zero does not exist)

breakaway_richness_df$gavage <- as.factor(breakaway_richness_df$gavage)
divnet_alpha_diversity_df$gavage <- as.factor(divnet_alpha_diversity_df$gavage)


### Observed richness
gamm_Obs <- gamm(Observed ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                 random = list(mouse_id = ~1),
                 data = alpha_phyloseq)
summary(gamm_Obs$gam)
gam.check(gamm_Obs$gam) # check residuals
mgcv::concurvity(gamm_Obs$gam) # check concurvity
draw(gamm_Obs$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_Obs$gam, select = "s(day_c)")
draw(pair_comp)


### Shannon diversity
gamm_Shan <- gamm(Shannon ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                  random = list(mouse_id = ~1),
                  data = alpha_phyloseq)
summary(gamm_Shan$gam)
gam.check(gamm_Shan$gam) # check residuals
mgcv::concurvity(gamm_Shan$gam) # check concurvity
draw(gamm_Shan$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_Shan$gam, select = "s(day_c)")
draw(pair_comp)


### Simpson diversity
gamm_Simp <- gamm(Simpson ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                  random = list(mouse_id = ~1),
                  data = alpha_phyloseq)
summary(gamm_Simp$gam)
gam.check(gamm_Simp$gam) # check residuals
mgcv::concurvity(gamm_Simp$gam) # check concurvity
draw(gamm_Simp$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_Simp$gam, select = "s(day_c)")
draw(pair_comp)


### Chao1 diversity
gamm_Chao <- gamm(Chao1 ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                  random = list(mouse_id = ~1),
                  data = alpha_phyloseq)
summary(gamm_Chao$gam)
gam.check(gamm_Chao$gam) # check residuals
mgcv::concurvity(gamm_Chao$gam) # check concurvity
draw(gamm_Chao$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_Chao$gam, select = "s(day_c)")
draw(pair_comp)


### richness (breakaway)
gamm_b_Obs <- gamm(estimate ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                   random = list(mouse_id = ~1),
                   data = breakaway_richness_df,
                   weights = 1 / se^2) # inverse-variance weighting
summary(gamm_b_Obs$gam)
gam.check(gamm_b_Obs$gam) # check residuals
mgcv::concurvity(gamm_b_Obs$gam) # check concurvity
draw(gamm_b_Obs$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_b_Obs$gam, select = "s(day_c)")
draw(pair_comp)


### Shannon diversity (DivNet)
gamm_dn_Shan <- gamm(shannon_estimate ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                     random = list(mouse_id = ~1),
                     data = divnet_alpha_diversity_df,
                     weights = 1 / (shannon_error^2)) # inverse-variance weighting
summary(gamm_dn_Shan$gam)
gam.check(gamm_dn_Shan$gam) # check residuals
mgcv::concurvity(gamm_dn_Shan$gam) # check concurvity
draw(gamm_dn_Shan$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_dn_Shan$gam, select = "s(day_c)")
draw(pair_comp)


### Simpson diversity (DivNet)
gamm_dn_Simp <- gamm(simpson_estimate ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                     random = list(mouse_id = ~1),
                     data = divnet_alpha_diversity_df,
                     weights = 1 / (simpson_error^2)) # inverse-variance weighting
summary(gamm_dn_Simp$gam)
gam.check(gamm_dn_Simp$gam) # check residuals
mgcv::concurvity(gamm_dn_Simp$gam) # check concurvity
draw(gamm_dn_Simp$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_dn_Simp$gam, select = "s(day_c)")
draw(pair_comp)


### plot smooths for Shannon
# create data.frame of smooths for all gavage x centered day combinations
smooth_df <- expand.grid(day_c = seq(min(alpha_phyloseq$day_c), max(alpha_phyloseq$day_c),
                                   length = 200), gavage = levels(alpha_phyloseq$gavage))

# predict smooths
pred <- predict(gamm_Shan$gam, newdata = smooth_df, se.fit = TRUE, type = "response")

# add predictions to smooth_df
smooth_df$fit   <- pred$fit
smooth_df$upper <- pred$fit + 2 * pred$se.fit
smooth_df$lower <- pred$fit - 2 * pred$se.fit

# add day to smooth_df
smooth_df$day <- smooth_df$day_c + mean(alpha_phyloseq$day)

ggplot() + theme_minimal() +
  geom_point(data = alpha_phyloseq, aes(day, Shannon, color = gavage), alpha = 0.4) +
  geom_line(data = smooth_df, aes(day, fit, color = gavage), linewidth = 1) +
  geom_ribbon(data = smooth_df, aes(day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
  labs(title = "GAMM Smooths for Shannon Diversity", y = "Shannon diversity", x = "Day")


### function to plot smooths for alpha diversity metrics
plot_alpha_gamm <- function(data, gamm_obj, alpha_div_metric) {
  
  # create data.frame of smooths for all gavage x centered day combinations
  smooth_df <- expand.grid(day_c = seq(min(data$day_c), max(data$day_c), length = 200),
                           gavage = levels(data$gavage))
  
  # predict smooths
  pred <- predict(gamm_obj$gam, newdata = smooth_df, se.fit = TRUE, type = "response")
  
  # add predictions to smooth_df
  smooth_df$fit   <- pred$fit
  smooth_df$upper <- pred$fit + 2 * pred$se.fit
  smooth_df$lower <- pred$fit - 2 * pred$se.fit
  
  # add day to smooth_df
  smooth_df$day <- smooth_df$day_c + mean(alpha_phyloseq$day)
  
  ggplot() + theme_minimal() +
    geom_point(data = data, aes(x = day, y = .data[[alpha_div_metric]], color = gavage), alpha = 0.4) +
    geom_line(data = smooth_df, aes(day, fit, color = gavage), linewidth = 1) +
    geom_ribbon(data = smooth_df, aes(day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
    labs(title = paste0("GAMM Smooths for ", alpha_div_metric), y = alpha_div_metric, x = "Day")
}

plot_alpha_gamm(data = alpha_phyloseq, gamm_obj = gamm_Obs, alpha_div_metric = "Observed")
plot_alpha_gamm(data = alpha_phyloseq, gamm_obj = gamm_Shan, alpha_div_metric = "Shannon")
plot_alpha_gamm(data = alpha_phyloseq, gamm_obj = gamm_Simp, alpha_div_metric = "Simpson")
plot_alpha_gamm(data = alpha_phyloseq, gamm_obj = gamm_Chao, alpha_div_metric = "Chao1")


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


############################################################################
#####   BETA DIVERSITY ANALYSIS - STATISTICAL TESTING WITH PERMANOVA   #####
############################################################################

### tests whether group centroids differ in multivariate space

# function for pairwise by time point
pairwise_by_timepoint <- function(dist_matrix, meta, group_var, time_var, 
                                  perm = 999, p_adjust_method = "BH") {
  
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

### bray-curtis
adonis2(bray_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
betadisper(bray_dist, group = meta$gavage) %>% permutest()

# pairwise by group per time point
pairwise_by_timepoint(dist_matrix = bray_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") 

lapply(unique(meta$day), function(d) {
  idx <- meta$day == d
  bd <- betadisper(as.dist(as.matrix(bray_dist)[idx, idx]), group = meta$gavage[idx])
  return(bd)
}) %>% lapply(permutest)


### jaccard
adonis2(jaccard_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
betadisper(jaccard_dist, group = meta$gavage) %>% permutest()

# pairwise by group per time point
pairwise_by_timepoint(dist_matrix = jaccard_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") 

lapply(unique(meta$day), function(d) {
  idx <- meta$day == d
  bd <- betadisper(as.dist(as.matrix(jaccard_dist)[idx, idx]), group = meta$gavage[idx])
  return(bd)
}) %>% lapply(permutest)


### euclidean
adonis2(euc_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
betadisper(euc_dist, group = meta$gavage) %>% permutest()

# pairwise by group per time point
pairwise_by_timepoint(dist_matrix = euc_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day") 

lapply(unique(meta$day), function(d) {
  idx <- meta$day == d
  bd <- betadisper(as.dist(as.matrix(euc_dist)[idx, idx]), group = meta$gavage[idx])
  return(bd)
}) %>% lapply(permutest)


### canberra
adonis2(canberra_dist ~ gavage * day, data = meta, strata = meta$mouse_id) 
betadisper(canberra_dist, group = meta$gavage) %>% permutest() 

# pairwise by group per time point
pairwise_by_timepoint(dist_matrix = canberra_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day")

lapply(unique(meta$day), function(d) {
  idx <- meta$day == d
  bd <- betadisper(as.dist(as.matrix(canberra_dist)[idx, idx]), group = meta$gavage[idx])
  return(bd)
}) %>% lapply(permutest)


### aitchison
adonis2(aitchison_dist ~ gavage * day, data = meta, strata = meta$mouse_id)
betadisper(aitchison_dist, group = meta$gavage) %>% permutest()

# pairwise by group per time point
pairwise_by_timepoint(dist_matrix = aitchison_dist,
                      meta = meta, 
                      group_var = "gavage", 
                      time_var = "day")

lapply(unique(meta$day), function(d) {
  idx <- meta$day == d
  bd <- betadisper(as.dist(as.matrix(aitchison_dist)[idx, idx]), group = meta$gavage[idx])
  return(bd)
}) %>% lapply(permutest)


###################################################################
#####   BETA DIVERSITY ANALYSIS - PCOA ORDINATION AND PLOTS   #####
###################################################################

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


############################################################################
#####   BETA DIVERSITY ANALYSIS - DISTANCE-BASED REDUNDANCY ANALYSIS   #####
############################################################################

# variation that can be attributed to explanatory variables (constrained PCoA) 
# capscale() for non-euclidean (bray-curtis, canberra, jaccard) and rda() for euclidean (euclidean and aitchison)

# stratified permutation scheme
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


###################################################################
#####   BETA DIVERSITY ANALYSIS - NMDS ORDINATION AND PLOTS   #####
###################################################################

### preserves rank order of dissmilarities

### NMDS ordination
ordination_nmds_bray <- ordinate(ps_rel, method = "NMDS", distance = bray_dist) # bray-curtis
ordination_nmds_bray$stress
ordination_nmds_jaccard <- ordinate(ps_pa, method = "NMDS", distance = jaccard_dist) # jaccard
ordination_nmds_jaccard$stress 
ordination_nmds_canberra <- ordinate(ps_rel, method = "NMDS", distance = canberra_dist) # canberra
ordination_nmds_canberra$stress

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


#####################################################################
#####   BETA DIVERSITY ANALYSIS - EXTRACT PCA SCORES AND PLOT   #####
#####################################################################

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


########################################################################################################################
#####   LONGITUDINAL DIFFERENTIAL ABUNDANCE - MICROBIOME MULTIVARIABLE ASSOCIATION WITH LINEAR MODELS (MAASLIN2)   #####
########################################################################################################################

# filter low prevalence taxa (present in less than 10% of samples)
feat_filtered <- feat[rowSums(feat > 0) >= ceiling(0.10 * ncol(feat)), ] 

feat_t <- t(feat_filtered) # transpose feature table
all(rownames(feat_t) == meta$sample_id)  # check that sample ids match

meta$gavage <- factor(meta$gavage, levels = c("G_1DMD", "G_4DMD", "G_4W7C", "G_4WMD")) # gavage as factor

# create interaction columns
X <- model.matrix(~ gavage * day, data = meta) # design matrix
int_cols <- X[, grep(":", colnames(X)), drop = FALSE] # interaction terms
colnames(int_cols) <- gsub(":", ".", colnames(int_cols)) # rename column names
meta <- cbind(meta, int_cols) # add interaction columns to metadata

# run MaAsLin2
mas_long_inter <- Maaslin2(input_data = feat_t,
                           input_metadata = meta,
                           output = "maaslin2_longitudinal",
                           fixed_effects = c("gavage", "day", colnames(int_cols)),
                           random_effects = "mouse_id",
                           normalization = "CLR",
                           transform = "NONE",
                           analysis_method = "LM")

maaslin2_res <- read.table("maaslin2_longitudinal/all_results.tsv", header = TRUE)
maaslin2_res <- maaslin2_res %>%
  filter(qval < 0.05) %>%
  arrange(pval)


###################################################################################################################################
#####   LONGITUDINAL DIFFERENTIAL ABUNDANCE - COUNT REGRESSION FOR CORRELATED OBSERVATIONS WITH THE BETA-BINOMIAL (CORNCOB)   #####
###################################################################################################################################

# filter low prevalence taxa (present in less than 10% of samples)
feat_filtered <- feat[rowSums(feat > 0) >= ceiling(0.10 * ncol(feat)), ] 

# format feature table and metadata
feat_otu <- as.matrix(feat_filtered) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE)  # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create a phyloseq object
ps <- phyloseq(feat_otu, sampledata)

# run corncob
corncob_res <- differentialTest(formula = ~ gavage * day,
                                phi.formula = ~ 1,
                                formula_null = ~ gavage + day,
                                phi.formula_null = ~ 1,
                                data = ps,
                                test = "LRT",
                                random = ~ mouse_id,
                                fdr_cutoff = 0.05)


taxa_names <- rownames(feat_filtered) # get taxa names
# iterate function over all significant models and bind data.frames together
corncob_coef_df <- map_dfr(seq_along(corncob_res$all_models), function(i) {
  model <- corncob_res$all_models[[i]] # get i-th model
  taxa_name <- taxa_names[i] # get taxa corresponding to i-th model
  
  # extract coeeficients matrix from model
  coefs <- model$coefficients
  
  # keep only mu rows (abundance)
  mu_coefs <- coefs[grep("^mu\\.", rownames(coefs)), , drop = FALSE]
  
  # convert to data.frame and remove mu. prefix
  df <- as.data.frame(mu_coefs) %>%
    rownames_to_column(var = "term") %>%
    mutate(term = gsub("^mu\\.", "", term),
           taxa = taxa_name)
  
  return(df)
})

# calculate adjusted p-values and filter
corncob_coef_df <- corncob_coef_df %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, method = "BH")) %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(Estimate)) %>%
  filter(term != "(Intercept)")


#####################################################################################################
#####   LONGITUDINAL DIFFERENTIAL ABUNDANCE - SPLINE-BASED - GENERALIZED ADDITIVE MIXED MODEL   #####
#####################################################################################################

# allows non-linear changes temporal changes 

# filter low prevalence taxa (present in less than 10% of samples)
feat_filtered <- feat[rowSums(feat > 0) >= ceiling(0.10 * ncol(feat)), ] 

# CLR transformation
feat_filtered <- t(feat_filtered) # transpose
feat_rel_abund <- feat_filtered/rowSums(feat_filtered) # convert feature table to relative abundances
feat_clr <- clr(feat_rel_abund + 1e-6) # add pseudocount and perform CLR transformation
feat_clr <- t(feat_clr)
feat_clr <- as.data.frame(feat_clr)

# ensure sample names are the same
all(colnames(feat_clr) == rownames(meta))

# pivot feature table to long form
feat_long <- feat_clr %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(cols = -taxon,
               names_to = "sample_id",
               values_to = "abundance")

# merge with metadata
df_long <- feat_long %>%
  left_join(meta, by = "sample_id")

# convert gavage to factor
df_long$gavage <- as.factor(df_long$gavage)

# center day (day zero does not exist)
df_long <- df_long %>% mutate(day_c = day - mean(day))


### fit model with one feature
taxon_name <- "Dorea_longicatena"
df_taxon <- df_long %>% filter(taxon == taxon_name)
gamm_model <- gamm(abundance ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                   random = list(mouse_id = ~1),
                   data = df_taxon)
summary(gamm_model$gam)
gam.check(gamm_model$gam) # check residuals
mgcv::concurvity(gamm_model$gam) # check concurvity (collinearity for smooth terms)
draw(gamm_model$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_model$gam, select = "s(day_c)")
draw(pair_comp)


taxon_name <- "Faecalibacterium_prausnitzii"
df_taxon <- df_long %>% filter(taxon == taxon_name)
gamm_model <- gamm(abundance ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                   random = list(mouse_id = ~1),
                   data = df_taxon)
summary(gamm_model$gam)
gam.check(gamm_model$gam) # check residuals
mgcv::concurvity(gamm_model$gam) # check concurvity (collinearity for smooth terms)
draw(gamm_model$gam) # smooth plots

# pairwise comparisons
pair_comp <- difference_smooths(gamm_model$gam, select = "s(day_c)")
draw(pair_comp)


### plot smooths for one feature
# create data.frame of smooths for all gavage x centered day combinations
smooth_df <- expand.grid(day_c = seq(min(df_taxon$day_c), max(df_taxon$day_c), length = 200),
                         gavage = levels(df_taxon$gavage))

# predict smooths
pred <- predict(gamm_model$gam, newdata = smooth_df, se.fit = TRUE, type = "response")

# add predictions to smooth_df
smooth_df$fit   = pred$fit
smooth_df$upper = pred$fit + 2 * pred$se.fit
smooth_df$lower = pred$fit - 2 * pred$se.fit

# add day to smooth_df
smooth_df$day <- smooth_df$day_c + mean(df_taxon$day)

ggplot() + theme_minimal() +
  geom_point(data = df_taxon, aes(day, abundance, color = gavage), alpha = 0.4) +
  geom_line(data = smooth_df, aes(day, fit, color = gavage), linewidth = 1) +
  geom_ribbon(data = smooth_df, aes(x = day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
  labs(title = paste0("GAMM Smooths for ", taxon_name), y = "CLR abundance", x = "Day") # + facet_wrap(~gavage) # include for faceted plots


### fit model to all features
fit_gamm_taxon <- function(df, taxon_name) {
  df_taxon <- df %>% filter(taxon == taxon_name)
  
  # try-catch wrapper
  out <- tryCatch({
    
    fit <- gamm(abundance ~ gavage + s(day_c, by = gavage, bs = "fs", k = 4),
                random = list(mouse_id = ~1),
                data = df_taxon)
    
    sm <- summary(fit$gam)
    
    ### extract parameteric and smooth terms separately and rename columns
    # parameteric terms (gavage)
    param_df <- as.data.frame(sm$p.table) %>%
      mutate(term = rownames(sm$p.table),
             taxon = taxon_name,
             .before = 1) %>%
      rename(estimate = Estimate, std_error = `Std. Error`,
             t_value = `t value`, p_value = `Pr(>|t|)`)
    
    # smooth terms (s(day):gavage)
    smooth_df <- as.data.frame(sm$s.table) %>%
      mutate(smooth = rownames(sm$s.table),
             taxon = taxon_name,
             .before = 1) %>%
      rename(edf = edf, ref_df = Ref.df, 
             F_value = F, p_value = `p-value`)
    
    list(param = param_df,
         smooth = smooth_df,
         model_gam = fit$gam,
         success = TRUE,
         error = NA)
    
    }, error = function(e) {
      
    # return empty rows, but keep track of taxon and error
    list(param = tibble(taxon = taxon_name, term = NA, estimate = NA, 
                        std_error = NA, t_value = NA, p_value = NA),
         smooth = tibble(taxon = taxon_name, smooth = NA, edf = NA, 
                         ref_df = NA, F_value = NA, p_value = NA),
         model_gam = NULL,
         success = FALSE, error = as.character(e))
  })
  
  return(out)
}

taxa_list <- unique(df_long$taxon) # list of all taxa
gamm_results <- map(taxa_list, ~ fit_gamm_taxon(df_long, .x)) 

# combine parametric coefficients (gavage) into a data.frame
param_results <- purrr::map_df(gamm_results, "param") 
param_results$p_adj <- p.adjust(param_results$p_value, method = "BH")
param_results <- param_results %>%
  filter(p_adj < 0.05) %>%
  arrange(p_value)

# combine smooth term results (s(day):gavage) into a data.frame and keep only significant values
smooth_results <- purrr::map_df(gamm_results, "smooth") 
smooth_results$p_adj <- p.adjust(smooth_results$p_value, method = "BH")
smooth_results <- smooth_results %>%
  filter(p_adj < 0.05) %>%
  arrange(p_value)

# plot pairwise differences
taxon_name <- "Faecalibacterium_prausnitzii"
smooth_diff <- gamm_results[[which(taxa_list == taxon_name)]]
pair_comp <- difference_smooths(smooth_diff$model_gam, select = "s(day_c)")
draw(pair_comp)

# extract success and failure flags
fit_status <- tibble(taxon = names(gamm_results),
                     success = purrr::map_lgl(gamm_results, "success"),
                     error = purrr::map_chr(gamm_results, "error"))

top_taxa_gamm <- c("Faecalibacterium_prausnitzii", "Bacteroides_nordii", "Parabacteroides_distasonis", 
                   "Turicibacter_sp_LA61", "Dorea_longicatena", "Parabacteroides_goldsteinii", 
                   "Bifidobacterium_pseudolongum", "Lactobacillus_johnsonii", "Blautia_sp_MarseilleP2398",
                   "Akkermansia_muciniphila", "Bacteroides_fragilis", "Anaerostipes_caccae",
                   "Bacteroides_thetaiotaomicron", "Bacteroides_vulgatus", "Roseburia_intestinalis")


#######################################################################
#####   LONGITUDINAL DIFFERENTIAL ABUNDANCE - LME-BASED - LINDA   #####
#######################################################################

# linear model with correction for compositional effects (estimates biased by heavy tails and dense signal)

# filter low prevalence taxa (present in less than 10% of samples)
feat_filtered <- feat[rowSums(feat > 0) >= ceiling(0.10 * ncol(feat)), ] 

# ensure sample names are the same
all(colnames(feat_filtered) == rownames(meta))

### main effects (group + time)
linda_main <- linda(otu.tab = feat_filtered, meta = meta,
                    formula = "~ gavage + day + (1 | mouse_id)",
                    pseudo.cnt = 0.5) # add pseudocount

# combine all linda_main$output data.frames into a single data.frame
linda_main_df <- imap_dfr(linda_main$output, ~ {
  df <- .x
  df <- df %>% 
    rownames_to_column("taxon")
  df %>%
    mutate(effect = .y)
})

linda_main_df <- linda_main_df %>%
  filter(padj < 0.05) %>%
  arrange(padj)


### interaction (group x time)
linda_inter <- linda(otu.tab = feat_filtered, meta = meta,
                     formula = "~ gavage * day + (1 | mouse_id)",
                     pseudo.cnt = 0.5) # add pseudocount

# combine all linda_inter$output data.frames into a single data.frame
linda_inter_df <- imap_dfr(linda_inter$output, ~ {
  df <- .x
  df <- df %>% 
    rownames_to_column("taxon")
  df %>%
    mutate(effect = .y)
})

linda_inter_df <- linda_inter_df %>%
  filter(padj < 0.05) %>%
  arrange(padj)


##########################################################################
#####   LONGITUDINAL DIFFERENTIAL ABUNDANCE - GLMM-BASED - GLMMTMB   #####
##########################################################################

# negative binomial and zero-inflated, negative binomial model

# filter low prevalence taxa (present in less than 10% of samples)
feat_filtered <- feat[rowSums(feat > 0) >= ceiling(0.10 * ncol(feat)), ] 

# ensure sample names are the same
all(colnames(feat_filtered) == rownames(meta))

# pivot feature table to long form
feat_count_long <- feat_filtered %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(cols = -taxon,
               names_to = "sample_id",
               values_to = "count")

# merge with metadata
df_count_long <- feat_count_long %>%
  left_join(meta, by = "sample_id")

# convert gavage to factor
df_count_long$gavage <- as.factor(df_count_long$gavage)


### fit model with one feature - negative binomial model
taxon_name <- "Faecalibacterium_prausnitzii"
df_count_taxon <- df_count_long %>% filter(taxon == taxon_name)

model_nb <- glmmTMB(count ~ gavage * day + (1 | mouse_id),
                    family = nbinom2, 
                    data = df_count_taxon)
summary(model_nb)$coefficients$cond # Wald z-tests
testUniformity(simulateResiduals(fittedModel = model_nb, n = 1000)) # check residuals


### fit model with one feature - zero-inflated, negative binomial model
taxon_name <- "Bacteroides_intestinalis"
df_count_taxon <- df_count_long %>% filter(taxon == taxon_name)

model_zinb <- glmmTMB(count ~ gavage * day + (1 | mouse_id),
                      ziformula = ~ gavage, # probability of zero inflation dependent on gavage group
                      family = nbinom2,
                      data = df_count_taxon)
summary(model_zinb)$coefficients$cond # Wald z-tests
testUniformity(simulateResiduals(fittedModel = model_zinb, n = 1000)) # check residuals


### fit (zero-inflated) negative binomial model to all features
fit_glmmTMB_taxon <- function(df, taxon_name, zero_inflated = FALSE) {
  df_taxon <- df %>% filter(taxon == taxon_name)
  
  # try-catch wrapper
  out <- tryCatch({
    
    # choose model type
    if (zero_inflated) {
      fit <- glmmTMB(count ~ gavage * day + (1 | mouse_id),
                     ziformula = ~ gavage, 
                     family = nbinom2, 
                     data = df_taxon)
    } else {
      fit <- glmmTMB(count ~ gavage * day + (1 | mouse_id),
                     family = nbinom2,
                     data = df_taxon)
    }

    # extract coefficients
    coefs <- summary(fit)$coefficients$cond
    
    # create data.frame with results
    param_df <- as.data.frame(coefs) %>%
      mutate(term = rownames(coefs),
             taxon = taxon_name,
             .before = 1)
    
    list(param = param_df, 
         success = TRUE, 
         error = NA)
    
  }, error = function(e) {
    
    # return empty rows, but keep track of taxon and error
    list(param = tibble(taxon = taxon_name, term = NA, Estimate = NA, 
                        Std.Error = NA, z_value = NA, p_value = NA),
         success = FALSE, error = as.character(e))
  })
  
  return(out)
}


### fit model to multiple features - negative binomial model
taxa_list <- unique(df_count_long$taxon) # list of all taxa
glmmTMB_nb_results_list <- map(taxa_list, ~ fit_glmmTMB_taxon(df_count_long, .x))

# combine results into a data.frame
glmmTMB_nb_results <- purrr::map_df(glmmTMB_nb_results_list, "param") 

glmmTMB_nb_results <- glmmTMB_nb_results %>%
  mutate(p_adj = p.adjust(`Pr(>|z|)`, method = "BH")) %>%
  filter(p_adj < 0.05) %>%
  arrange(`Pr(>|z|)`)

# extract success and failure flags
fit_status_nb <- tibble(taxon = taxa_list,
                     success = purrr::map_lgl(glmmTMB_nb_results_list, "success"),
                     error = purrr::map_chr(glmmTMB_nb_results_list, "error"))

### fit model to multiple features - zero-inflated, negative binomial model
df_zero_inflated_taxa <- df_count_long %>%
  group_by(taxon) %>%
  summarise(prop_zero = mean(count == 0)) %>%
  filter(prop_zero > 0.60) # subset data.frame where > 60% of counts are zero 

zero_inflated_taxa <- unique(df_zero_inflated_taxa$taxon) # list of all zero-inflated taxa 

glmmTMB_zinb_results_list <- map(zero_inflated_taxa, ~ fit_glmmTMB_taxon(df_count_long, .x, zero_inflated = TRUE))

# combine results into a data.frame
glmmTMB_zinb_results <- purrr::map_df(glmmTMB_zinb_results_list, "param") 

glmmTMB_zinb_results <- glmmTMB_zinb_results %>%
  mutate(p_adj = p.adjust(`Pr(>|z|)`, method = "BH")) %>%
  filter(p_adj < 0.05) %>%
  arrange(`Pr(>|z|)`)

# extract success and failure flags
fit_status_zinb <- tibble(taxon = zero_inflated_taxa,
                          success = purrr::map_lgl(glmmTMB_zinb_results_list, "success"),
                          error = purrr::map_chr(glmmTMB_zinb_results_list, "error"))


############################################################################################
#####   RELATIVE ABUNDANCE PHYLOSEQ OBJECT FOR DIFFERENTIAL ABUNDANCE VISUALIZATIONS   #####
############################################################################################

# format feature table and metadata
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

all(colnames(feat_otu) == rownames(meta)) # ensure sample names are the same
sampledata <- sample_data(meta) # convert to sample_data

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)

# filter low prevalence taxa (present in less than 10% of samples)
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= 0.1 * nsamples(ps), prune = TRUE)

ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x)) 


##########################################
#####   HEATMAPS OF TAXA ABUNDANCE   #####
##########################################

### subset phyloseq object to desired taxa
top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:25] # get top 25 taxa by abundance
top_taxa_gamm # top taxa from gamm analysis
ps_top <- prune_taxa(top_taxa_gamm, ps_rel) # filter by top_taxa, top_taxa_gamm, etc


### plot heatmap with ggplot
df_heat <- psmelt(ps_top)  # long format data.frame for plotting

df_heat <- df_heat %>%
  group_by(OTU, gavage_day) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(OTU) %>%
  mutate(Abundance_scaled = scales::rescale(Abundance)) %>%
  ungroup() # scale each taxon

# plot heatmap
ggplot(df_heat, aes(x = gavage_day, y = OTU, fill = Abundance_scaled)) +
  geom_tile() + scale_fill_viridis_c() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Heatmap of top taxa")


### plot heatmap with pheatmap
mat <- as.data.frame(otu_table(ps_top)) # otu table
mat_t <- as.data.frame(t(mat)) # transpose metadata 
mat_t$gavage_day <- meta$gavage_day # add metadata to otu table

mat_day <- mat_t %>%
  group_by(gavage_day) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame() # aggregate by gavage data

mat_day_t <- t(as.matrix(mat_day[,-1])) # drop first column and transpose
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
  group_by(day, gavage, OTU) %>%
  summarize(mean_abundance = mean(Abundance), 
            sem = sd(Abundance)/sqrt(n()),
            .groups = "drop")


### plot abundance over time for specific taxa
ggplot(df_traj %>% filter(OTU == "Faecalibacterium_prausnitzii"), 
       aes(x = day, y = mean_abundance, color = gavage, group = gavage, fill = gavage)) +
  geom_line() + geom_point() + theme_minimal() +
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  labs(title = "Trajectory of species over time", 
       x = "Day", y = "Mean relative abundance")

ggplot(df_traj %>% filter(OTU == "Turicibacter_sp_LA61"), 
       aes(x = day, y = mean_abundance, color = gavage, group = gavage, fill = gavage)) +
  geom_line() + geom_point() + theme_minimal() +
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  labs(title = "Trajectory of species over time", 
       x = "Day", y = "Mean relative abundance")


### plot abundance trajectory for top 6 taxa
# get top 6 taxa by mean abundance
top6 <- df_traj %>%
  group_by(OTU) %>%
  summarize(total = sum(mean_abundance)) %>%
  slice_max(total, n = 6) %>%
  pull(OTU)

ggplot(df_traj %>% filter(OTU %in% top6), 
       aes(x = day, y = mean_abundance, color = gavage, group = gavage, fill = gavage)) +
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  facet_wrap(~OTU, scales = "free_y") + theme_minimal() +
  labs(title = "Trajectory of top 6 species over time", 
       x = "Day", y = "Mean relative abundance")


### plot abundance trajectory for top taxa from gamm analysis
ggplot(df_traj %>% filter(OTU %in% top_taxa_gamm), 
       aes(x = day, y = mean_abundance, color = gavage, group = gavage, fill = gavage)) +
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
  facet_wrap(~OTU, scales = "free_y") + theme_minimal() +
  labs(title = "Trajectory of top species from gamm analysis over time", 
       x = "Day", y = "Mean relative abundance")


################################
#####   STACKED BARPLOTS   #####
################################

### color palette
colors <- c("#862185", "#009E73", "#88CCEE", "#CC6677", "#D55E00", 
            "#44AA99", "#332288", "#E69F00", "#0072B2", "#AA4499",
            "#F0E442", "#117733", "#DB72FB", "#619CFF", "#882255",
            "#999999")

### long format data.frame for plotting
df_stack <- psmelt(ps_rel)

df_stack <- df_stack %>%
  group_by(gavage_day, OTU) %>%
  summarise(abundance = mean(Abundance),
            gavage = first(gavage),
            day = first(day),
            day_factor = first(day_factor),
            .groups = "drop") %>%
  arrange(desc(abundance))


### add columns for plotting to df_stack
top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:15] # get top 15 taxa by abundance
df_stack$top_taxa <- ifelse(df_stack$OTU %in% top_taxa, as.character(df_stack$OTU), "Other") # assign all non-top_taxa as Other
df_stack$top_taxa <- factor(df_stack$top_taxa, levels = c(setdiff(unique(df_stack$top_taxa), "Other"), "Other"))

top_taxa_gamm # top taxa from gamm analysis
df_stack$top_taxa_gamm <- ifelse(df_stack$OTU %in% top_taxa_gamm, as.character(df_stack$OTU), "Other") # assign all non-top_taxa_gamm as Other
df_stack$top_taxa_gamm <- factor(df_stack$top_taxa_gamm, levels = c(setdiff(unique(df_stack$top_taxa_gamm), "Other"), "Other"))


### plot using top_taxa or top_taxa_gamm
top <- "top_taxa"
my_colors <- rep(colors, length.out = length(unique(df_stack[[top]])))

# plot all gavage groups on one graph
ggplot(df_stack, aes(x = gavage_day, y = abundance, fill = .data[[top]])) +
  geom_bar(stat = "identity") + theme_minimal() +
  scale_fill_manual(values = my_colors) +
  labs(title = "Stacked barplot of important taxa over time", 
       x = "Gavage group and day", y = "Mean relative abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# facet gavage groups
ggplot(df_stack, aes(x = day_factor, y = abundance, fill = .data[[top]])) +
  geom_bar(stat = "identity") + theme_minimal() + 
  scale_fill_manual(values = my_colors) +
  facet_wrap(~gavage, scales = "free_x") +
  labs(title = "Stacked barplot of important taxa over time", 
       x = "Day", y = "Mean relative abundance") + 
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
#   [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] glmmTMB_1.1.12       LinDA_0.2.0          mgcv_1.9-3           compositions_2.0-9  
# [5] Maaslin2_1.22.0      colorspace_2.1-2     pheatmap_1.0.13      reshape2_1.4.4      
# [9] scales_1.4.0         pairwiseAdonis_0.4.1 cluster_2.1.8.1      microbiome_1.30.0   
# [13] nlme_3.1-168         doParallel_1.0.17    iterators_1.0.14     foreach_1.5.2       
# [17] DivNet_0.4.1         breakaway_4.8.4      emmeans_1.11.2-8     lmerTest_3.1-3      
# [21] lme4_1.1-37          Matrix_1.7-4         vegan_2.7-1          permute_0.9-8       
# [25] phyloseq_1.52.0      lubridate_1.9.4      forcats_1.0.1        stringr_1.5.2       
# [29] dplyr_1.1.4          purrr_1.1.0          readr_2.1.5          tidyr_1.3.1         
# [33] tibble_3.3.0         tidyverse_2.0.0      ggplot2_4.0.0       
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3               Rdpack_2.6.4            sandwich_3.1-1          rlang_1.1.6            
# [5] magrittr_2.0.4          clue_0.3-66             multcomp_1.4-28         ade4_1.7-23            
# [9] compiler_4.5.0          vctrs_0.6.5             rmutil_1.1.10           pkgconfig_2.0.3        
# [13] crayon_1.5.3            XVector_0.48.0          modeest_2.4.0           tzdb_0.5.0             
# [17] UCSC.utils_1.4.0        nloptr_2.2.1            GenomeInfoDb_1.44.3     jsonlite_2.0.0         
# [21] biomformat_1.36.0       rhdf5filters_1.20.0     Rhdf5lib_1.30.0         biglm_0.9-3            
# [25] R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3      rpart_4.1.24           
# [29] boot_1.3-32             numDeriv_2016.8-1.1     estimability_1.5.1      Rcpp_1.1.0             
# [33] zoo_1.8-14              mvnfast_0.2.8           IRanges_2.42.0          splines_4.5.0          
# [37] igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.17.1      
# [41] dichromat_2.0-0.1       abind_1.4-8             timeDate_4041.110       TMB_1.9.17             
# [45] codetools_0.2-20        lattice_0.22-7          plyr_1.8.9              Biobase_2.68.0         
# [49] withr_3.0.2             S7_0.2.0                stable_1.1.6            coda_0.19-4.1          
# [53] Rtsne_0.17              survival_3.8-3          bayesm_3.1-6            getopt_1.20.4          
# [57] Biostrings_2.76.0       pillar_1.11.1           tensorA_0.36.2.1        stats4_4.5.0           
# [61] pcaPP_2.0-5             reformulas_0.4.1        generics_0.1.4          S4Vectors_0.48.0       
# [65] hms_1.1.3               timeSeries_4041.111     minqa_1.2.8             xtable_1.8-4           
# [69] glue_1.8.0              statip_0.2.3            tools_4.5.0             robustbase_0.99-6      
# [73] data.table_1.17.8       spatial_7.3-18          fBasics_4041.97         mvtnorm_1.3-3          
# [77] rhdf5_2.52.1            optparse_1.7.5          ape_5.8-1               rbibutils_2.3          
# [81] GenomeInfoDbData_1.2.14 cli_3.6.5               DEoptimR_1.1-4          gtable_0.3.6           
# [85] stabledist_0.7-2        digest_0.6.37           BiocGenerics_0.54.0     ggrepel_0.9.6          
# [89] TH.data_1.1-4           farver_2.1.2            multtest_2.64.0         lifecycle_1.0.4        
# [93] httr_1.4.7              multcompView_0.1-10     MASS_7.3-65 

