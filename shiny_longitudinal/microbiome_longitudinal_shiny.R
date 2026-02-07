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
library(compositions)
library(purrr)
library(gratia)


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
  formula <- as.formula(paste0(metric, " ~ gavage + s(day, by = gavage, bs = 'fs', k = 4)"))
  gamm_obj <- mgcv::gamm(formula, random = list(mouse_id = ~1), data = data)
  
  # precompute residuals
  res <- resid(gamm_obj$gam)
  
  # precompute concurvity
  conc <- mgcv::concurvity(gamm_obj$gam, full = TRUE)
  
  # precompute pairwise smooth differences
  smooth_diff <- difference_smooths(gamm_obj$gam, select = "s(day)") %>%
    mutate(metric = metric, .before = 1) 
  
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
       smooth_df = smooth_df,
       smooth_diff = smooth_diff)
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
ps_clr <- transform(ps_filt, "clr")  # clr transformation (microbiome package adds pseudocount) for aitchison


### distance matrices 
distances <- list(bray_curtis = distance(ps_rel, method = "bray"), # bray-curtis
                  jaccard = distance(ps_pa, method = "jaccard", binary = TRUE), # jaccard
                  canberra = distance(ps_rel, method = "canberra"), # canberra
                  aitchison = distance(ps_clr, method = "euclidean"))  # aitchison (euclidean in clr space)


### PCoA ordination
ordination <- list(bray_curtis = ordinate(ps_rel, method = "PCoA", distance = distances$bray_curtis), # bray-curtis
                   jaccard = ordinate(ps_pa, method = "PCoA", distance = distances$jaccard), # jaccard
                   canberra = ordinate(ps_rel, method = "PCoA", distance = distances$canberra), # canberra
                   aitchison = ordinate(ps_clr, method = "PCoA", distance = distances$aitchison)) # aitchison


### NMDS ordination (non-Euclidean distances)
set.seed(1234) # NMDS is stochastic
nmds_ordination <- list(bray_curtis = ordinate(ps_rel, method = "NMDS", distance = distances$bray_curtis), # bray-curtis
                        jaccard = ordinate(ps_pa, method = "NMDS", distance = distances$jaccard), # jaccard
                        canberra = ordinate(ps_rel, method = "NMDS", distance = distances$canberra)) # canberra


### PCA ordination (Euclidean distances)
pca_ordination <- list(
  aitchison = {
    clr_otu <- t(otu_table(ps_clr)) # extract CLR-transformed data and transpose
    prcomp(clr_otu, center = TRUE, scale. = FALSE) # run PCA
    }
)

### constrained ordination (db-RDA | RDA)
# capscale (distance-based RDA) used on non-Euclidean distances
con_ordination <- list(bray_curtis = capscale(distances$bray_curtis ~ gavage + day, data = meta), # bray-curtis
                       jaccard = capscale(distances$jaccard ~ gavage + day, data = meta), # jaccard
                       canberra = capscale(distances$canberra ~ gavage + day, data = meta), # canberra
                       aitchison = rda(t(otu_table(ps_clr)) ~ gavage + day, data = meta)) # aitchison


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


### constrained ordination permutation tests
perm_mouse <- how(blocks = meta$mouse_id, nperm = 999) # define permutation scheme

anova_perm <- function(con_ord_list, perm) {
  lapply(con_ord_list, function(c) {
    list(overall = anova(c, permutations = perm),
         terms = anova(c, by = "terms", permutations = perm),
         axis = anova(c, by = "axis", permutations = perm))
  })
}
anova_permutation <- anova_perm(con_ordination, perm_mouse)


### save beta-diversity outputs to avoid recalculation in the Shiny app
beta_diversity <- list()

for (metric in names(distances)) {
  beta_diversity[[metric]] <- list(distance = distances[[metric]],
                                   permanova = list(overall = permanova_overall[[metric]],
                                                    by_time = permanova_time[[metric]]),
                                   anova_perm = list(overall = anova_permutation[[metric]]$overall,
                                                     terms = anova_permutation[[metric]]$terms,
                                                     axis = anova_permutation[[metric]]$axis),
                                   dispersion = list(overall = dispersion_overall[[metric]],
                                                     by_time = dispersion_time[[metric]]),
                                   ordination = list(pcoa = ordination[[metric]],
                                                     nmds = if (metric %in% names(nmds_ordination)) nmds_ordination[[metric]] else NULL, # only calculated for non-Euclidean distances
                                                     pca = if (metric %in% names(pca_ordination)) pca_ordination[[metric]] else NULL, # only calculated for Euclidean distances
                                                     con_ord = con_ordination[[metric]]))
}

saveRDS(beta_diversity, "microbiome_analysis/shiny_longitudinal/beta_diversity.rds")


################################################################################
#####     DIFFERENTIAL ABUNDANCE ANALYSIS PRECOMPUTATION FOR SHINY APP     #####
################################################################################

### relative abundance table
# use the created phyloseq object

# filter low prevalence taxa (present in less than 10% of samples)
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= 0.1 * nsamples(ps), prune = TRUE)

# convert to relative abundance 
ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x)) 

# long format data.frame for plotting
df_taxa <- psmelt(ps_rel)  

# summarize relative abundance of taxa over time and gavage group
rel_abun_df <- df_taxa %>%
  group_by(day, gavage, OTU) %>%
  summarize(mean_abundance = mean(Abundance), 
            sem = sd(Abundance)/sqrt(n()),
            .groups = "drop") %>%
  rename(taxon = OTU) %>%
  as.data.frame()


### data.frame for heatmap and stacked barplot
# use ps_rel from the relative abundance table 

# melt to long format
comm_over_df <- psmelt(ps_rel)

# aggregate by gavage_day (mean across samples)
comm_over_df <- comm_over_df %>%
  group_by(gavage_day, OTU, gavage, day, day_factor) %>%
  summarise(abundance = mean(Abundance), .groups = "drop")


### generalized additive mixed models

# filter low prevalence taxa (present in less than 10% of samples)
feat_filtered <- feat[rowSums(feat > 0) >= ceiling(0.10 * ncol(feat)), ] 

# CLR transformation
feat_filtered <- t(feat_filtered) # transpose
feat_rel_abund <- feat_filtered/rowSums(feat_filtered) # convert feature table to relative abundances
feat_clr <- clr(feat_rel_abund + 1e-6) # add pseudocount and perform CLR transformation
feat_clr <- t(feat_clr)
feat_clr <- as.data.frame(feat_clr)

# ensure sample names are the same
stopifnot(all(colnames(feat_clr) == rownames(meta)))

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


### function for generalized additive mixed model
fit_gamm_taxon <- function(taxon_name, data) {
  
  # filter long data.frame by taxon
  df_taxon <- data %>% filter(taxon == taxon_name)
  
  out <- tryCatch({
  gamm_obj <- mgcv::gamm(abundance ~ gavage + s(day, by = gavage, bs = "fs", k = 4),
                         random = list(mouse_id = ~1),
                         data = df_taxon)
  
  # precompute residuals
  res <- resid(gamm_obj$gam)
  
  # precompute concurvity
  conc <- mgcv::concurvity(gamm_obj$gam, full = TRUE)
  
  # precompute pairwise smooth differences
  smooth_diff <- difference_smooths(gamm_obj$gam, select = "s(day)") %>%
    mutate(taxon = taxon_name, .before = 1) 
  
  # precompute summary
  sum <- summary(gamm_obj$gam)
  
  
  ### extract parameteric and smooth terms separately and rename columns
  # parameteric terms (gavage)
  param_df <- as.data.frame(sum$p.table) %>%
    mutate(term = rownames(sum$p.table),
           taxon = taxon_name,
           .before = 1) %>%
    rename(estimate = Estimate, std_error = `Std. Error`,
           t_value = `t value`, p_value = `Pr(>|t|)`)
  
  # smooth terms (s(day):gavage)
  smooth_df <- as.data.frame(sum$s.table) %>%
    mutate(smooth = rownames(sum$s.table),
           taxon = taxon_name,
           .before = 1) %>%
    rename(edf = edf, ref_df = Ref.df, 
           F_value = F, p_value = `p-value`)
  
    
  ### precompute smooth predictions for plotting
  smooth_pred <- expand.grid(day = seq(min(df_taxon$day), max(df_taxon$day), length = 200),
                             gavage = levels(df_taxon$gavage))
  
  # predict smooths
  pred <- predict(gamm_obj$gam, newdata = smooth_pred, se.fit = TRUE, type = "response") 
  
  # add predictions to smooth_pred
  smooth_pred$fit <- pred$fit
  smooth_pred$upper <- pred$fit + 2 * pred$se.fit
  smooth_pred$lower <- pred$fit - 2 * pred$se.fit
  smooth_pred$taxon <- taxon_name
 
  # list to store
  list(model = gamm_obj,
       data = df_taxon,
       summary_param = param_df,
       summary_smooth = smooth_df,
       residuals = res,
       concurvity = conc,
       smooth_pred = smooth_pred,
       smooth_diff = smooth_diff,
       success = TRUE,
       error = NA) 
  
  }, error = function(e) {
    
    list(model = NULL,
         data = NULL,
         summary_param = tibble(term = NA, taxon = taxon_name),
         summary_smooth = tibble(smooth = NA, taxon = taxon_name),
         residuals = NULL,
         concurvity = NULL,
         smooth_pred = NULL,
         smooth_diff = NULL,
         success = FALSE,
         error = as.character(e)) 
    
  })
  out
}


# fit models for all taxa
taxa_list <- unique(df_long$taxon) # list of all taxa
gamm_results_taxa <- purrr::map(taxa_list, ~ fit_gamm_taxon(.x, df_long))
names(gamm_results_taxa) <- taxa_list

# calculate adjusted p-value for parametric coefficients and smooths
param_all <- purrr::map_df(gamm_results_taxa, "summary_param")
smooth_all <- purrr::map_df(gamm_results_taxa, "summary_smooth")

param_all$p_adj <- p.adjust(param_all$p_value, method = "BH")
smooth_all$p_adj <- p.adjust(smooth_all$p_value, method = "BH")

# add p_adj values to summany_param and summary_smooth individual taxon models
gamm_results_taxa <- lapply(gamm_results_taxa, function(x) {
  x$summary_smooth <- x$summary_smooth %>%
    left_join(smooth_all %>% select(smooth, taxon, p_adj), by = c("smooth", "taxon"))
  
  x$summary_param <- x$summary_param %>%
    left_join(param_all %>% select(term, taxon, p_adj), by = c("term", "taxon"))
  
  return(x)
})

# pairwise smooth differences
all_smooth_diffs <- purrr::map_df(gamm_results_taxa, "smooth_diff")


# save differential abundance outputs to avoid recalculation in the Shiny app
differential_abundance <- list(comm_over = comm_over_df,
                               rel_abun = rel_abun_df,
                               gamm = list(data = df_long,
                                           models = gamm_results_taxa,
                                           smooth_diff = all_smooth_diffs))

saveRDS(differential_abundance, "microbiome_analysis/shiny_longitudinal/differential_abundance.rds")

