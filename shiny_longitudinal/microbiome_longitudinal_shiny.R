### Longitudinal Analysis of Gavage Data - Shiny App

## Data Analysis/Precomputation

# load libraries
library(phyloseq)
library(tidyverse)
library(vegan)
library(microbiome)
library(lme4)
library(lmerTest)
library(emmeans)
library(pairwiseAdonis)
library(mgcv)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio")

### load and format data
# metadata
meta <- read.table("metadata_gavage.txt", header = TRUE)
rownames(meta) <- meta$sample_id
meta$day_factor <- factor(meta$day) # add day as a factor for plotting
sampledata <- sample_data(meta) # convert to sample_data

# save metadata as RDS file for Shiny app
saveRDS(meta, "microbiome_analysis/shiny_longitudinal/metadata.rds")

# feature table
feat <- read.table("otu_table_gavage.txt", header = TRUE)
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

# ensure sample names are the same
stopifnot(all(colnames(feat_otu) == rownames(meta)))

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)


################################################################
#####     ALPHA DIVERSITY PRECOMPUTATION FOR SHINY APP     #####
################################################################

# rarefy the phyloseq object (for Observed richness)
ps_rarefied <- rarefy_even_depth(ps, rngseed = 1234, sample.size = min(sample_sums(ps)), verbose = FALSE)


### calculate traditional alpha diversity metrics using phyloseq
alpha_diversity_metrics <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_diversity_metrics$Chao1 <- estimate_richness(ps, measures = "Chao1")$Chao1 # add Chao1 to alpha_diversity_metrics
alpha_diversity_metrics$Observed <- estimate_richness(ps_rarefied, measures = "Observed")$Observed # add Observed to alpha_diversity_metrics
alpha_diversity_metrics$sample_name <- rownames(alpha_diversity_metrics) # add column with sample names
alpha_div_df <- left_join(alpha_diversity_metrics, meta, by = c("sample_name" = "sample_id")) # merge phyloseq object with metadata (for plotting and statistical tests) 


### linear mixed effects models
fit_lmm <- function(metric, data = data) {
  formula <- as.formula(paste0(metric, " ~ gavage * day_factor + (1 | mouse_id)")) # build formula string
  model <- lmer(formula, data = data) # fit linear-mixed effects model
  
  # list to store
  list(model = model, 
       anova = anova(model),
       emmeans = emmeans(model, ~ gavage | day_factor) %>% # commute estimated marginal means
         pairs(adjust = "tukey") %>% # tukey pairwise contrasts within each day
         as.data.frame())
}


### generalized additive mixed models
fit_gamm <- function(metric, data) {
  data$gavage <- as.factor(data$gavage)
  formula <- as.formula(paste0(metric, " ~ s(day, k = 4, by = gavage) + gavage"))
  gamm_obj <- mgcv::gamm(formula, random = list(mouse_id = ~1), data = data)
  
  # precompute residuals
  res <- resid(gamm_obj$gam)
  
  # precompute concurvity
  conc <- mgcv::concurvity(gamm_obj$gam, full = TRUE)
  
  # precompute smooth predictions for plotting
  smooth_df <- expand.grid(day = seq(min(data$day), max(data$day), length = 200),
                           gavage = levels(data$gavage))
  
  # predict smooths
  pred <- predict(gamm_obj$gam, newdata = smooth_df, se.fit = TRUE, type = "response") 
  
  # add predictions to smooth_df
  smooth_df$fit <- pred$fit
  smooth_df$upper <- pred$fit + 2 * pred$se.fit
  smooth_df$lower <- pred$fit - 2 * pred$se.fit
  
  # list to store
  list(model = gamm_obj,
       summary = summary(gamm_obj$gam),
       residuals = res,
       concurvity = conc,
       smooth_df = smooth_df)
}

# fit models for all metrics
metrics <- c("Observed", "Shannon", "Simpson", "Chao1")
models <- set_names(metrics) %>% # assigns names to vector elements (resulting list will have named elements)
  map(function(metric){
    list(lmm = fit_lmm(metric, alpha_div_df),
         gamm = fit_gamm(metric, alpha_div_df))
  })


# save alpha diversity outputs to avoid recalculation in the Shiny app
alpha_diversity <- list(metrics = alpha_diversity_metrics,
                        models = models)
saveRDS(alpha_diversity, "microbiome_analysis/shiny_longitudinal/alpha_diversity.rds")



###############################################################
#####     BETA DIVERSITY PRECOMPUTATION FOR SHINY APP     #####
###############################################################

# filter low prevalence taxa (present in less than 10% of samples)
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= 0.1 * nsamples(ps), prune = TRUE)

ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x)) # relative abundance for bray-curtis and canberra
ps_pa <- transform_sample_counts(ps_filt, function(x) as.numeric(x > 0)) # presence/absence for jaccard
ps_log <- transform_sample_counts(ps_filt, function(x) log1p(x))  # log-transformation for euclidean
ps_clr <- transform(ps_filt, "clr")  # clr transformation (microbiome package adds pseudocount) for aitchison


### distance matrices 
distances <- list(bray_curtis = distance(ps_rel, method = "bray"), # bray-curtis
                  jaccard = distance(ps_pa, method = "jaccard", binary = TRUE), # jaccard
                  euclidean = distance(ps_log, method = "euclidean"), # euclidean
                  canberra = distance(ps_rel, method = "canberra"), # canberra
                  aitchison = distance(ps_clr, method = "euclidean"))  # aitchison (euclidean in clr space)


### PCoA ordination
ordination <- list(bray_curtis = ordinate(ps_rel, method = "PCoA", distance = distances$bray_curtis), # bray-curtis
                   jaccard = ordinate(ps_pa, method = "PCoA", distance = distances$jaccard), # jaccard
                   euclidean = ordinate(ps_log, method = "PCoA", distance = distances$euclidean), # euclidean
                   canberra = ordinate(ps_rel, method = "PCoA", distance = distances$canberra), # canberra
                   aitchison = ordinate(ps_clr, method = "PCoA", distance = distances$aitchison)) # aitchison


### overall PERMANOVA
overall_permanova <- function(dist_list, meta) {
  lapply(dist_list, function(d) {
    adonis2(d ~ gavage * day, data = meta, strata = meta$mouse_id)
  })
}
permanova_overall <- overall_permanova(distances, meta)


### overall beta dispersion
overall_betadisper <- function(dist_list, meta) {
  lapply(dist_list, function(d) {
    betadisper(d, group = meta$gavage) %>% permutest()
  })
}
dispersion_overall <- overall_betadisper(distances, meta)


### time-resolved PERMANOVA
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

# pairwise by group per time point
permanova_by_time <- function(dist_list, meta) {
  lapply(dist_list, function(d) {
    df <- pairwise_by_timepoint(dist_matrix = d,
                                meta = meta,
                                group_var = "gavage",
                                time_var = "day")
    # remove sig and Df columns
    df <- df[, !names(df) %in% "sig", drop = FALSE]
    df <- df[, !names(df) %in% "Df", drop = FALSE]
    return(df)
    
  })
}

permanova_time <- permanova_by_time(distances, meta)


### time-resolved beta-dispersion
betadisper_by_time <- function(dist_list, meta) {
  lapply(dist_list, function(d) {
    lapply(unique(meta$day), function(t) {
      idx <- meta$day == t
      bd <- betadisper(as.dist(as.matrix(d)[idx, idx]),
                       group = meta$gavage[idx])
      permutest(bd)
    })
  })
}

dispersion_time <- betadisper_by_time(distances, meta)


### save beta-diversity outputs to avoid recalculation in the Shiny app
beta_diversity <- list()

for (metric in names(distances)) {
  beta_diversity[[metric]] <- list(distance = distances[[metric]],
                                   permanova = list(overall = permanova_overall[[metric]],
                                                    by_time = permanova_time[[metric]]),
                                   dispersion = list(overall = dispersion_overall[[metric]],
                                                     by_time = dispersion_time[[metric]]),
                                   ordination = list(pcoa = ordination[[metric]]))
}

saveRDS(beta_diversity, "microbiome_analysis/shiny_longitudinal/beta_diversity.rds")

