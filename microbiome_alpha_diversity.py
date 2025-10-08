#!/usr/bin/env python
# coding: utf-8



### load libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.diversity.alpha import shannon, simpson, chao1, sobs
from skbio.stats import subsample_counts
from scipy.stats import mannwhitneyu

### set working directory
import os
os.chdir("/Users/kristinvandenham/kmvanden/Jupyter")




### load data
# load metadata
meta = pd.read_csv("metadata.txt", sep=r'\s+', quotechar='"', index_col=0)
print(meta.head(10))
print(meta.shape)

# load feature table
feat = pd.read_csv("feature_table.txt", sep=r'\s+', quotechar='"', index_col=0)
feat = feat.T  # transpose
print(feat.iloc[0:3, 0:3])
print(feat.shape)

# ensure sample names are the same 
print((feat.index == meta.index).all())




### rarefaction 
min_depth = feat.sum(axis=1).min() # min sequencing depth

def rarefy_row(row, depth):
    return pd.Series(subsample_counts(row.values.astype(int), int(depth)), index=row.index)

feat_rarefied = feat.apply(rarefy_row, axis=1, args=(min_depth,)) # rarefy feature table to minimum depth without replacement
print(feat_rarefied.sum(axis=1))




### alpha diversity analysis
# richness, shannon and simpson using rarefied data
alpha_df = pd.DataFrame({
    'observed': feat_rarefied.apply(sobs, axis=1),
    'shannon': feat_rarefied.apply(shannon, axis=1),
    'simpson': feat_rarefied.apply(simpson, axis=1),
}) # rarefied data

# chao1 using unrarefied data
alpha_df['chao1'] = feat.apply(chao1, axis=1)

print(alpha_df.head(5))

# combine with metadata
alpha_df['sample_name'] = alpha_df.index
alpha_combined = alpha_df.merge(meta[['condition']], left_on='sample_name', right_index=True)
print(alpha_combined.head(5))




### plot alpha diversity measures

# plotting style
sns.set(style="whitegrid")

# plot arrangement
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# observed richness
sns.boxplot(data=alpha_combined, x='condition', y='observed', ax=axes[0, 0])
sns.stripplot(data=alpha_combined, x='condition', y='observed', color='black', alpha=0.5, ax=axes[0, 0])
axes[0, 0].set_title('Observed richness')
axes[0, 0].set_xlabel('Condition')
axes[0, 0].set_ylabel('Number of unique taxa')

# Shannon diversity
sns.boxplot(data=alpha_combined, x='condition', y='shannon', ax=axes[0, 1])
sns.stripplot(data=alpha_combined, x='condition', y='shannon', color='black', alpha=0.5, ax=axes[0, 1])
axes[0, 1].set_title('Shannon Diversity')
axes[0, 1].set_xlabel('Condition')
axes[0, 1].set_ylabel('Shannon index')

# Simpson diversity
sns.boxplot(data=alpha_combined, x='condition', y='simpson', ax=axes[1, 0])
sns.stripplot(data=alpha_combined, x='condition', y='simpson', color='black', alpha=0.5, ax=axes[1, 0])
axes[1, 0].set_title('Simpson Diversity')
axes[1, 0].set_xlabel('Condition')
axes[1, 0].set_ylabel('Simpson index')

# Chao1 estimator
sns.boxplot(data=alpha_combined, x='condition', y='chao1', ax=axes[1, 1])
sns.stripplot(data=alpha_combined, x='condition', y='chao1', color='black', alpha=0.5, ax=axes[1, 1])
axes[1, 1].set_title('Chao1 estimator')
axes[1, 1].set_xlabel('Condition')
axes[1, 1].set_ylabel('Chao1 estimator')

plt.tight_layout()
plt.show()




### Wilcoxon tests

# observed otus
group1 = alpha_combined[alpha_combined['condition'] == 'healthy']['observed']
group2 = alpha_combined[alpha_combined['condition'] == 'disease']['observed']
stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
print(f"Wilcoxon test p-value (observed otus): {p}")

# shannon diversity
group1 = alpha_combined[alpha_combined['condition'] == 'healthy']['shannon']
group2 = alpha_combined[alpha_combined['condition'] == 'disease']['shannon']
stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
print(f"Wilcoxon test p-value (shannon diversity): {p}")

# simpson diversity otus
group1 = alpha_combined[alpha_combined['condition'] == 'healthy']['simpson']
group2 = alpha_combined[alpha_combined['condition'] == 'disease']['simpson']
stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
print(f"Wilcoxon test p-value (simpson diversity): {p}")

# chao1 estimator
group1 = alpha_combined[alpha_combined['condition'] == 'healthy']['chao1']
group2 = alpha_combined[alpha_combined['condition'] == 'disease']['chao1']
stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
print(f"Wilcoxon test p-value (chao1 estimator): {p}")

