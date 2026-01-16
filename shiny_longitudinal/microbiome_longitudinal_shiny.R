### Longitudinal Analysis of Gavage Data - Shiny App

## Data Analysis/Precomputation

# load libraries
library(phyloseq)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(mgcv)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio")

### load and format data
# metadata
meta <- read.table("metadata_gavage.txt", header = TRUE)
rownames(meta) <- meta$sample_id
meta$day_factor <- factor(meta$day) # add day as a factor for plotting
sampledata <- sample_data(meta) # convert to sample_data

# feature table
feat <- read.table("otu_table_gavage.txt", header = TRUE)
feat_otu <- as.matrix(feat) # convert feature table to matrix
feat_otu <- otu_table(feat_otu, taxa_are_rows = TRUE) # convert to otu_table

# ensure sample names are the same
stopifnot(all(colnames(feat_otu) == rownames(meta)))

# create phyloseq object
ps <- phyloseq(feat_otu, sampledata)

# rarefy the phyloseq object (for Observed richness)
ps_rarefied <- rarefy_even_depth(ps, rngseed = 1234, sample.size = min(sample_sums(ps)), verbose = FALSE)


### calculate traditional alpha diversity metrics using phyloseq
alpha_div_df <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_div_df$Chao1 <- estimate_richness(ps, measures = "Chao1")$Chao1 # add Chao1 to alpha_div_df
alpha_div_df$Observed <- estimate_richness(ps_rarefied, measures = "Observed")$Observed # add Observed to alpha_div_df
alpha_div_df$sample_name <- rownames(alpha_div_df) # add column with sample names
alpha_div_df <- left_join(alpha_div_df, meta, by = c("sample_name" = "sample_id")) # merge phyloseq object with metadata (for plotting and statistical tests) 


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


# save outputs to avoid recalculation in the Shiny app
saveRDS(alpha_div_df, "microbiome_analysis/shiny_longitudinal/alpha_diversity_df.rds")
saveRDS(models, "microbiome_analysis/shiny_longitudinal/alpha_diversity_models.rds")

