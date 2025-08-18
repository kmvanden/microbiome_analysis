# Microbiome Analysis
## :microbe::chart_with_upwards_trend: Why Microbiome Data Are Statistically Challenging
**Compositionality**: microbiome data represent relative abundances, not absolute counts, which means that the data are constrained. Since all values are non-negative and sum to a constraint, the data exist in a simplex, not in standard Euclidean space. Furthermore, since the total abundance is fixed, spurious correlations and dependencies are created between taxa (i.e., an increase in the proportion of one taxa must be accompanied by decreases in the others, even if their absolute abundances did not change). Many statistical models assume that variables are independent and unconstrained in Euclidean space, and microbiome data violate these assumptions.
  - **Solution**: Log-ratio transformations move compositional data from the simplex into a standard Euclidean (unconstrained) space suitable for standard statistics. Alternatively, some methods model compositionality explicitly, like ALDEx2, which uses CLR transformation and Monte Carlo sampling from a Dirichlet-multinomial distribution.

**High dimensionality**: microbiome datasets typically contain many more features (taxa) than samples (p >> n), whereas most standard statistical tests assume that the number of samples exceeds the number of features. In small datasets, chance fluctuations in relative abundances can appear as significant differences and high numbers of features relative to the number of samples increases the risk of this happening. Additionally, with high dimensional data, the regression problem becomes undetermined: you have more unknowns (coefficients) than equations (observations), so there are an infinite number of solutions for the coefficients. In logistic regression, this creates a flat optimization surface (ill-conditioned), which can result in non-convergence, numerical instability (very large coefficients) and overfitting.
  - **Solution**: Dimensionality reduction techniques (e.g., PCA or PCoA) reduce the number of features to a smaller set of composite (latent) features. Feature selection by filtering low-abundance or low-prevalence taxa, regularization techniques (e.g., Lasso) or recursive feature elimination can help to identify a smaller subset of features to be used in subsequent analyses.

**Overdispersion and zero-inflation**: in microbiome data, most taxa are absent (zero counts) in most samples, due to both rare taxa being below detection thresholds (rounded or sampling zeros) and due to true biological absence of certain taxa (true or structural zeros). Sparsity is a major contributor to the overdispersion observed in microbiome count data: a species may be absent in most samples, but have very large values in other samples when it is present. Overdispersion violates the assumptions made by models used with count data, like the Poisson distribution, which assume that the variance does not exceed the mean. Additionally, the sparsity of microbiome data also results in zero-inflation (i.e., more zeros than expected under typical count distributions), due to the presence of both sampling and structural zeros.
  - **Solution**: negative binomial models (e.g., in DESeq2) have a dispersion parameter that allows them to model variance independently from the mean and to therefore capture the overdispersion present in microbiome data. Zero-inflation can also be controlled for by using zero-inflated negative binomial (ZINB) models, which explicitly model the zero counts as separate processes (a structural zero process and a negative binomial count process). Additionally, feature filtering before analysis can help to reduce the impact of zero-inflation. 

## :moneybag::balance_scale: Alpha Diversity
Alpha diversity measures the diversity within a single sample. Traditional alpha diversity metrics come from classical ecology and include: 
  - **Richness**: total number of unique taxa
  - **Shannon index**: richness and evenness (how many taxa are present and how evenly their abundances are distributed; penalizes dominance)
  - **Simpson index**: the probability that two individuals randomly selected from a sample belong to the same species, thus emphasizing dominant taxa

These traditional diversity metrics assume an equal sampling effort, but micobiome data typically have uneven sequencing depths.

**Rarefaction**: subsampling of each sample to the same sequencing depth (typically that of the shallowest sample).
  - Rarefaction removes differences in library size (sequencing depth), with the goal of ensuring fair comparison between the samples (i.e., deeper samples won’t appear more diverse).
  - However, the use of rarefaction is disputed, because it discards data (statistical power is decreased by shallower samples), introduces random variation (the subsampling is stochastic), can bias diversity estimation (especially for rare taxa) and does not correct for unobserved species.

Alpha diversity values often exhibit heteroscedasticity (unequal variances across groups) and non-normal distributions. As a result, non-parametric statistical tests (**Wilcoxon Rank-Sum test** for comparing two groups and the **Kruskal-Wallis test** for comparing more than two groups) are typically used. These tests are based on ranked data rather than raw values, and do not assume normality. However, they do assume that the distributions being compared have a similar shape (e.g., similar variance and skewness).

To avoid rarefaction and to better model the underlying microbial diversity, newer methods infer latent (unseen) diversity using statistical methods.

**Chao1**: estimates total species richness. Chao1 infers the number of unseen species by using the observed counts of singletons and doubletons and the assumption that individuals are randomly and independently sampled from the community.

**Breakaway**: estimates true species richness by fitting a rational function model (a ratio of polynomials) to the frequency count data (frequency of the frequencies) using weighted nonlinear least squares to model the unobserved portion of the community (unseen taxa). Weights are assigned based on the estimated variance of each frequency count (i.e., high variance points like singletons receive less weight).

**DivNet**: estimates Shannon and Simpson diversity by modeling latent relative abundances (true but unobserved proportions) of the taxa in the community. It assumes that the additive log-ratio (alr)-transformed latent proportions follow a multivariate normal distribution and that the observed counts result from a multinomial sampling process of the inverse alr-transformed latent proportions. DivNet uses maximum likelihood estimation to jointly estimate the mean vector (μ) and the covariance matrix (Σ) that define the latent distribution, and from that, derives the diversity estimates (the latent proportions most likely to have produced the observed counts) and their associated uncertainty.

**Betta**: a regression framework that performs statistical testing on alpha diversity estimates produced by breakaway or DivNet, while explicitly accounting for sampling variability in those estimates. It uses the alpha diversity values and their standard errors (breakaway) or variance estimates (DivNet) to fit a linear model, where the data points are weighed with respect to the precision (inverse variance) of each estimate. This allows for hypothesis testing or covariate modeling that properly incorporate uncertainty due to unobserved taxa, sequencing depth, and variability in taxon detection.

## :dotted_line_face::corn: Beta-Diversity
Beta diversity measures between-sample differences in whole microbiome communities (it asks are samples from different groups more different than samples from within a given group).

**Bray-Curtis distance**: measures the compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. It is computed as the weighted sum of absolute differences, where weights are the abundances, thus the metric is dominated by abundant taxa.
  - Values range from 0 (samples are identical) to 1 (completely dissimilar/ no shared taxa)
  - **Pre-processing**: low prevalence taxa should be removed (sensitive to sparsity) and counts should be transformed to relative abundances to avoid bias from sequencing depth.
  - **Limitations**: does not explicitly account for compositionality and tends to be biased toward abundant taxa.

**Canberra distance**: measures compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. It sums the ratios of the absolute differences to the sums of the abundances for each taxon, thereby giving relatively equal weight to rare and abundant taxa (i.e., rare taxa contribute proportionally more than they would in Bray-Curtis).
  - Values range from 0 (samples are identical) to 1 (completely dissimilar/ no shared taxa)
  - **Pre-processing**: low prevalence taxa should be removed (sensitive to sparsity) and counts should be transformed to relative abundances to avoid bias from sequencing depth.
  **Limitations**: does not explicitly account for compositionality

**Jaccard distance**: measures the dissimilarity between two samples based on presence/absence of taxa, ignoring their abundances. It quantifies how many taxa are shared versus unique to each sample.
  - Values range from 0 (samples have identical sets of taxa) to 1 (no shared taxa)
  - **Pre-processing**: low prevalence taxa should be removed (sensitive to sparsity) and counts should be converted to presence/absence data (any taxon with a count > 0 becomes 1 (present), otherwise it becomes 0 (absent)). No normalization is needed since it does not depend on abundance. 
  - **Limitations**: ignores abundance information and presence/absence data can be heavily influenced by sampling noise and sequencing depth.

**Euclidean distanc**e: the geometric (straight-line) distance between two samples in multidimensional space, where each dimension corresponds to a feature and and the coordinate of a sample along each axis is determined by the abundance of that feature.
  - Values range from zero (samples are identical) to larger positive values (greater overall differences in abundance)
  - **Pre-processing**: low prevalence taxa should be removed (sensitive to sparsity) and counts should be log transformed (reduces the influence of very abundant taxa)
  - **Limitations**: sensitive to magnitude differences (differences are squared) and thus influenced by highly abundant taxa and sequencing depths.

**Aitchison distance**: the Euclidean distance between samples after centered log-ratio (CLR) transformation (which projects the compositional data from the simplex to Euclidean space).
  - Values range from zero (log-ratios are identical) to large values (greater differences in log-ratio structure)
  - **Pre-processing**: counts need to be CLR-transformed (handles the compositional constraint of microbiome data).
  - **Limitations**: interpretation is less intuitive (represents how much log-ratios between taxa, rather than absolute or relative abundances, differ between samples).

**Principal Coordinate Analysis (PCoA)**: an ordination method that represents samples in a low-dimensional, Euclidean space, while preserving the pairwise distances (or dissimilarities) between them as faithfully as possible. It assumes that the input distance matrix (not the raw data) can be faithfully represented in Euclidean space. The resulting axes represent the amount of variation in the distance matrix explained by each coordinate (interpretable when using metric distances).

**Non-metric Multidimensional Scaling (NMDS)**: an ordination method that represents samples in a low-dimensional space based on rank orders of pairwise distances (not the distances themselves). NMDS does not assume that the data can be projected into Euclidean space and calculates a stress value to quantify how well the configuration preserves the rank-order of the original dissimilarities. It works with any distance metric, but the axes have no direct meaning (no percentage of variance explained).

**Principal Component Analysis (PCA)**: an ordination method that represents samples in a low-dimensional, Euclidean space by preserving the variance in the original feature table (e.g., taxa abundances, not pairwise distances). PCA assumes the data lie in Euclidean space and that linear combinations of features meaningfully capture variation. PCA is valid only when data have been transformed appropriately (e.g., CLR-transformed data). The resulting axes are ordered by the amount of variance they explain in the data.

**DivNet**: estimates Bray-Curtis and Euclidean distances by modeling latent relative abundances (true but unobserved proportions) of the taxa in each community. It assumes that the additive log-ratio (alr)-transformed latent proportions follow a multivariate normal distribution and that the observed counts result from a multinomial sampling process of the inverse alr-transformed latent proportions. DivNet uses maximum likelihood estimation to jointly estimate the mean vector (μ) and the covariance matrix (Σ) that define the latent distribution. Bray-Curtis and Euclidean distances are then calculated from the inferred latent relative abundances, capturing both observed data and sampling variability of each sample. 

**testBetaDiversity**: simulates multiple plausible versions of the latent relative compositions for each sample from the estimated multivariate normal distribution. For each simulation, it calculates Bray-Curtis and Euclidean distances between samples, reflecting both observed and unobserved diversity. To test whether within-group versus between-group distances differ, it permutes the sample labels across the simulations to construct a null distribution. This approach explicitly accounts for sampling uncertainty and compositional structure. 

**PERMANOVA**: uses a single, fixed distance matrix computed from observed data (doesn’t account for sampling variability or latent uncertainty). In order to test whether within-group versus between-group distances differ, a null distribution is generated by permuting sample labels. PERMANOVA assumes equal dispersion (variance) across groups. If this assumption is violated (e.g., one group is more dispersed), it can yield false positives even when group centroids are not meaningfully different.

