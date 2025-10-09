#!/usr/bin/env python
# coding: utf-8



### load librarires
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from skbio.diversity import beta_diversity
from skbio.stats.composition import clr
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from sklearn.manifold import MDS
from sklearn.decomposition import PCA

### set working directory
import os
os.chdir("/Users/kristinvandenham/kmvanden/Jupyter")




### load data
# load metadata
meta = pd.read_csv("metadata.txt", sep=r'\s+', quotechar='"', index_col=0)
print(meta.head(10))

# load feature table
feat = pd.read_csv("feature_table.txt", sep=r'\s+', quotechar='"', index_col=0)
feat = feat.T  # transpose
print(feat.iloc[0:3, 0:3])

# ensure sample names are the same 
print((feat.index == meta.index).all())




### data transformations
feat_rel = feat.div(feat.sum(axis=1), axis=0) # relative abundance for bray-curtis and canberra
feat_pa = feat.gt(0).astype(int) # # presence/absence for jaccard
feat_log = np.log1p(feat) # log transforamtion (with pseudocount) for euclidean
feat_clr = pd.DataFrame(clr(feat + 1), index=feat.index, columns=feat.columns) # clr transformation (with pseudocount) for aitchison




### compute distance matrices
bray_dist = beta_diversity("braycurtis", feat_rel.values, ids=feat_rel.index) # bray-curtis
jaccard_dist = beta_diversity("jaccard", feat_pa.values, ids=feat_pa.index) # jaccard
euc_dist = beta_diversity("euclidean", feat_log.values, ids=feat_log.index) # euclidean
canberra_dist = beta_diversity("canberra", feat_rel.values, ids=feat_rel.index) # canberra

# aitchison (euclidean on CLR)
clr_matrix = feat_clr.values  # extract numeric matrix (NumPy array)

# use Euclidean distance on CLR space
dists = pdist(clr_matrix, metric='euclidean') # euclidean on clr
aitchison_dist = squareform(dists) # convert to 2D square matrix
aitchison_dist = DistanceMatrix(aitchison_dist, ids=feat_clr.index) # wrap as a DistanceMatrix for downstream analyses




### PERMANOVA tests
distance_matrices = {
    'bray-curtis': bray_dist,
    'jaccard': jaccard_dist,
    'euclidean': euc_dist,
    'canberra': canberra_dist,
    'aitchison': aitchison_dist
} # define dictionary of distance matrices

# Run PERMANOVA for each
for name, dist in distance_matrices.items():
    print(f"{name} (PERMANOVA)")
    result = permanova(distance_matrix=dist, grouping=meta['condition'])
    print(result, "\n")




### PCOA ordination
# define dictionary of distance matrices
distance_matrices = {
    'bray': bray_dist,
    'jaccard': jaccard_dist,
    'euclidean': euc_dist,
    'canberra': canberra_dist,
    'aitchison': aitchison_dist
}

# output: ordination result and dataframe per distance type
pcoa_results = {} # ordination object (eigenvalues)
pcoa_dfs = {} # ordination dataframe

for name, dist in distance_matrices.items():
    ord_res = pcoa(dist) # run PCoA
    df = ord_res.samples.copy() # get coordinates
    df['condition'] = meta.loc[df.index, 'condition']  # add condition
    pcoa_results[name] = ord_res # store ordination object
    pcoa_dfs[name] = df # store dataframe




### PCoA plots

# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(2, 3, figsize=(10, 8))

# bray curtis
sns.scatterplot(data=pcoa_dfs['bray'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[0, 0])
axes[0, 0].set_title("Bray-Curtis (PCoA)")
axes[0, 0].set_xlabel(f"PC1 ({pcoa_results['bray'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[0, 0].set_ylabel(f"PC2 ({pcoa_results['bray'].proportion_explained.iloc[1]*100:.1f}% variance)")

# jaccard
sns.scatterplot(data=pcoa_dfs['jaccard'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[0, 1])
axes[0, 1].set_title("Jaccard (PCoA)")
axes[0, 1].set_xlabel(f"PC1 ({pcoa_results['jaccard'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[0, 1].set_ylabel(f"PC2 ({pcoa_results['jaccard'].proportion_explained.iloc[1]*100:.1f}% variance)")

# euclidean
sns.scatterplot(data=pcoa_dfs['euclidean'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[0, 2])
axes[0, 2].set_title("Euclidean (PCoA)")
axes[0, 2].set_xlabel(f"PC1 ({pcoa_results['euclidean'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[0, 2].set_ylabel(f"PC2 ({pcoa_results['euclidean'].proportion_explained.iloc[1]*100:.1f}% variance)")

# canberra
sns.scatterplot(data=pcoa_dfs['canberra'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[1, 0])
axes[1, 0].set_title("Canberra (PCoA)")
axes[1, 0].set_xlabel(f"PC1 ({pcoa_results['canberra'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[1, 0].set_ylabel(f"PC2 ({pcoa_results['canberra'].proportion_explained.iloc[1]*100:.1f}% variance)")

# aitchison
sns.scatterplot(data=pcoa_dfs['aitchison'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[1, 1])
axes[1, 1].set_title("Aitchison (PCoA)")
axes[1, 1].set_xlabel(f"PC1 ({pcoa_results['aitchison'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[1, 1].set_ylabel(f"PC2 ({pcoa_results['aitchison'].proportion_explained.iloc[1]*100:.1f}% variance)")

# hide empty subplot
axes[1, 2].axis('off')

# remove default legends in first five subplots
for ax in axes.flat:
    if ax.legend_ is not None:
        ax.legend_.remove()  

# legend in empty subplot
handles, labels = axes[0, 0].get_legend_handles_labels()
axes[1, 2].legend(handles, labels, title="Condition", loc='center')

plt.tight_layout()
plt.show()




### NMDS ordination
nmds_results = {}
nmds_dfs = {} # NMDS cooridinate dataframes

for name, dist in distance_matrices.items():
    if name in ['bray', 'jaccard', 'canberra']: # only use bray, jaccard and canberra distance matrices
        mds = MDS(n_components=2, dissimilarity="precomputed", metric=False, random_state=1234, n_init=20) # non-metric MND model
        coords = mds.fit_transform(dist.data) # fit NMDS model to distance matrix
        df = pd.DataFrame(coords, index=dist.ids, columns=['NMDS1', 'NMDS2'])
        df['condition'] = meta.loc[df.index, 'condition']
        nmds_results[name] = mds
        nmds_dfs[name] = df # store coordinates and condition




### NMDS plots
# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# bray curtis
sns.scatterplot(data=nmds_dfs['bray'], x='NMDS1', y='NMDS2', hue='condition', s=100, ax=axes[0])
axes[0].set_title("Bray-Curtis (NMDS)")
axes[0].set_xlabel("NMDS1")
axes[0].set_ylabel("NMDS2")

# jaccard
sns.scatterplot(data=nmds_dfs['jaccard'], x='NMDS1', y='NMDS2', hue='condition', s=100, ax=axes[1])
axes[1].set_title("Jaccard (NMDS)")
axes[1].set_xlabel("NMDS1")
axes[1].set_ylabel("NMDS2")

# canberra
sns.scatterplot(data=nmds_dfs['canberra'], x='NMDS1', y='NMDS2', hue='condition', s=100, ax=axes[2])
axes[2].set_title("Canberra (NMDS)")
axes[2].set_xlabel("NMDS1")
axes[2].set_ylabel("NMDS2")

plt.tight_layout()
plt.show()


# In[ ]:


# in R, NMDS is done with metaMDS() from vegan
# it uses a non-metric MDS with monotonic regression and multiple random starts 
# it tries to optimize stress and improves the starting configuration (like Shepard diagrams and runs multiple iterations)




### extract PCA data
# euclidean
pca_log = PCA(n_components=2) 
pca_log_scores = pca_log.fit_transform(feat_log.values) # run PCA
pca_log_df = pd.DataFrame(pca_log_scores, index=feat_log.index, columns=['PC1', 'PC2'])
pca_log_df = pca_log_df.join(meta)
var_explained_log = pca_log.explained_variance_ratio_ * 100 # variance explained

# aitchison
pca_clr = PCA(n_components=2)
pca_clr_scores = pca_clr.fit_transform(feat_clr.values) # run PCA
pca_clr_df = pd.DataFrame(pca_clr_scores, index=feat_clr.index, columns=['PC1', 'PC2'])
pca_clr_df = pca_clr_df.join(meta)
var_explained_clr = pca_clr.explained_variance_ratio_ * 100 # variance explained




### PCA plots
# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# euclidean
sns.scatterplot(data=pca_log_df, x='PC1', y='PC2', hue='condition', s=80, ax=axes[0])
axes[0].set_title("PCA of log-transformed data (Euclidean)")
axes[0].set_xlabel(f"PC1 ({var_explained_log[0]:.1f}% variance)")
axes[0].set_xlabel(f"PC2 ({var_explained_log[1]:.1f}% variance)")

# aitchison
sns.scatterplot(data=pca_clr_df, x='PC1', y='PC2', hue='condition', s=80, ax=axes[1])
axes[1].set_title("PCA of CLR-transformed data (Aitchison)")
axes[1].set_xlabel(f"PC1 ({var_explained_clr[0]:.1f}% variance)")
axes[1].set_xlabel(f"PC2 ({var_explained_clr[1]:.1f}% variance)")

plt.tight_layout()
plt.show()




### prevalence filtering of taxa
# check metadata
print(meta.head(10))

# check feature table
print(feat.iloc[0:3, 0:3])
print(feat.shape)

# ensure sample names are the same 
print((feat.index == meta.index).all())

# filter any features present in 10% or fewer samples
feat_filt = feat.loc[:, (feat > 0).sum(axis=0) >= int(0.1 * feat.shape[0])]
print(feat_filt.shape)




### data transformations
feat_rel = feat_filt.div(feat_filt.sum(axis=1), axis=0) # relative abundance for bray-curtis and canberra
feat_pa = feat_filt.gt(0).astype(int) # # presence/absence for jaccard
feat_log = np.log1p(feat_filt) # log transforamtion (with pseudocount) for euclidean
feat_clr = pd.DataFrame(clr(feat_filt + 1), index=feat_filt.index, columns=feat_filt.columns) # clr transformation (with pseudocount) for aitchison




### compute distance matrices
bray_dist = beta_diversity("braycurtis", feat_rel.values, ids=feat_rel.index) # bray-curtis
jaccard_dist = beta_diversity("jaccard", feat_pa.values, ids=feat_pa.index) # jaccard
euc_dist = beta_diversity("euclidean", feat_log.values, ids=feat_log.index) # euclidean
canberra_dist = beta_diversity("canberra", feat_rel.values, ids=feat_rel.index) # canberra

# aitchison (euclidean on CLR)
clr_matrix = feat_clr.values  # extract numeric matrix (NumPy array)

# use Euclidean distance on CLR space
dists = pdist(clr_matrix, metric='euclidean') # euclidean on clr
aitchison_dist = squareform(dists) # convert to 2D square matrix
aitchison_dist = DistanceMatrix(aitchison_dist, ids=feat_clr.index) # wrap as a DistanceMatrix for downstream analyses




### PERMANOVA tests
distance_matrices = {
    'bray-curtis': bray_dist,
    'jaccard': jaccard_dist,
    'euclidean': euc_dist,
    'canberra': canberra_dist,
    'aitchison': aitchison_dist
} # define dictionary of distance matrices

# Run PERMANOVA for each
for name, dist in distance_matrices.items():
    print(f"{name} (PERMANOVA)")
    result = permanova(distance_matrix=dist, grouping=meta['condition'])
    print(result, "\n")




### PCOA ordination
# define dictionary of distance matrices
distance_matrices = {
    'bray': bray_dist,
    'jaccard': jaccard_dist,
    'euclidean': euc_dist,
    'canberra': canberra_dist,
    'aitchison': aitchison_dist
}

# output: ordination result and dataframe per distance type
pcoa_results = {} # ordination object (eigenvalues)
pcoa_dfs = {} # ordination dataframe

for name, dist in distance_matrices.items():
    ord_res = pcoa(dist) # run PCoA
    df = ord_res.samples.copy() # get coordinates
    df['condition'] = meta.loc[df.index, 'condition']  # add condition
    pcoa_results[name] = ord_res # store ordination object
    pcoa_dfs[name] = df # store dataframe




### PCoA plots
# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(2, 3, figsize=(10, 8))

# bray curtis
sns.scatterplot(data=pcoa_dfs['bray'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[0, 0])
axes[0, 0].set_title("Bray-Curtis (PCoA)")
axes[0, 0].set_xlabel(f"PC1 ({pcoa_results['bray'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[0, 0].set_ylabel(f"PC2 ({pcoa_results['bray'].proportion_explained.iloc[1]*100:.1f}% variance)")

# jaccard
sns.scatterplot(data=pcoa_dfs['jaccard'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[0, 1])
axes[0, 1].set_title("Jaccard (PCoA)")
axes[0, 1].set_xlabel(f"PC1 ({pcoa_results['jaccard'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[0, 1].set_ylabel(f"PC2 ({pcoa_results['jaccard'].proportion_explained.iloc[1]*100:.1f}% variance)")

# euclidean
sns.scatterplot(data=pcoa_dfs['euclidean'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[0, 2])
axes[0, 2].set_title("Euclidean (PCoA)")
axes[0, 2].set_xlabel(f"PC1 ({pcoa_results['euclidean'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[0, 2].set_ylabel(f"PC2 ({pcoa_results['euclidean'].proportion_explained.iloc[1]*100:.1f}% variance)")

# canberra
sns.scatterplot(data=pcoa_dfs['canberra'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[1, 0])
axes[1, 0].set_title("Canberra (PCoA)")
axes[1, 0].set_xlabel(f"PC1 ({pcoa_results['canberra'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[1, 0].set_ylabel(f"PC2 ({pcoa_results['canberra'].proportion_explained.iloc[1]*100:.1f}% variance)")

# aitchison
sns.scatterplot(data=pcoa_dfs['aitchison'], x='PC1', y='PC2', hue='condition', s=50, ax=axes[1, 1])
axes[1, 1].set_title("Aitchison (PCoA)")
axes[1, 1].set_xlabel(f"PC1 ({pcoa_results['aitchison'].proportion_explained.iloc[0]*100:.1f}% variance)")
axes[1, 1].set_ylabel(f"PC2 ({pcoa_results['aitchison'].proportion_explained.iloc[1]*100:.1f}% variance)")

# hide empty subplot
axes[1, 2].axis('off')

# remove default legends in first five subplots
for ax in axes.flat:
    if ax.legend_ is not None:
        ax.legend_.remove()  

# legend in empty subplot
handles, labels = axes[0, 0].get_legend_handles_labels()
axes[1, 2].legend(handles, labels, title="Condition", loc='center')

plt.tight_layout()
plt.show()




### NMDS ordination
nmds_results = {}
nmds_dfs = {} # NMDS cooridinate dataframes

for name, dist in distance_matrices.items():
    if name in ['bray', 'jaccard', 'canberra']: # only use bray, jaccard and canberra distance matrices
        mds = MDS(n_components=2, dissimilarity="precomputed", metric=False, random_state=1234, n_init=20) # non-metric MND model
        coords = mds.fit_transform(dist.data) # fit NMDS model to distance matrix
        df = pd.DataFrame(coords, index=dist.ids, columns=['NMDS1', 'NMDS2'])
        df['condition'] = meta.loc[df.index, 'condition']
        nmds_results[name] = mds
        nmds_dfs[name] = df # store coordinates and condition




### NMDS plots
# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# bray curtis
sns.scatterplot(data=nmds_dfs['bray'], x='NMDS1', y='NMDS2', hue='condition', s=100, ax=axes[0])
axes[0].set_title("Bray-Curtis (NMDS)")
axes[0].set_xlabel("NMDS1")
axes[0].set_ylabel("NMDS2")

# jaccard
sns.scatterplot(data=nmds_dfs['jaccard'], x='NMDS1', y='NMDS2', hue='condition', s=100, ax=axes[1])
axes[1].set_title("Jaccard (NMDS)")
axes[1].set_xlabel("NMDS1")
axes[1].set_ylabel("NMDS2")

# canberra
sns.scatterplot(data=nmds_dfs['canberra'], x='NMDS1', y='NMDS2', hue='condition', s=100, ax=axes[2])
axes[2].set_title("Canberra (NMDS)")
axes[2].set_xlabel("NMDS1")
axes[2].set_ylabel("NMDS2")

plt.tight_layout()
plt.show()




### extract PCA data
# euclidean
pca_log = PCA(n_components=2) 
pca_log_scores = pca_log.fit_transform(feat_log.values) # run PCA
pca_log_df = pd.DataFrame(pca_log_scores, index=feat_log.index, columns=['PC1', 'PC2'])
pca_log_df = pca_log_df.join(meta)
var_explained_log = pca_log.explained_variance_ratio_ * 100 # variance explained

# aitchison
pca_clr = PCA(n_components=2)
pca_clr_scores = pca_clr.fit_transform(feat_clr.values) # run PCA
pca_clr_df = pd.DataFrame(pca_clr_scores, index=feat_clr.index, columns=['PC1', 'PC2'])
pca_clr_df = pca_clr_df.join(meta)
var_explained_clr = pca_clr.explained_variance_ratio_ * 100 # variance explained




### PCA plots
# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# euclidean
sns.scatterplot(data=pca_log_df, x='PC1', y='PC2', hue='condition', s=80, ax=axes[0])
axes[0].set_title("PCA of log-transformed data (Euclidean)")
axes[0].set_xlabel(f"PC1 ({var_explained_log[0]:.1f}% variance)")
axes[0].set_xlabel(f"PC2 ({var_explained_log[1]:.1f}% variance)")

# aitchison
sns.scatterplot(data=pca_clr_df, x='PC1', y='PC2', hue='condition', s=80, ax=axes[1])
axes[1].set_title("PCA of CLR-transformed data (Aitchison)")
axes[1].set_xlabel(f"PC1 ({var_explained_clr[0]:.1f}% variance)")
axes[1].set_xlabel(f"PC2 ({var_explained_clr[1]:.1f}% variance)")

plt.tight_layout()
plt.show()

