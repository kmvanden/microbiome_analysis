# NetCoMi Workflow - same filter as other analyses

# load libraries
library(ggplot2)
library(tidyverse)
library(NetCoMi)

# set.seed
set.seed(1234)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name
group = as.factor(meta$condition) # label as a factor


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

### filter species present in less than 10% of samples
feat_t <- as.data.frame(t(feat)) # transpose feature table
min_prevalence <- 0.10
feat_filtered <- feat_t[, colMeans(feat_t > 0) >= min_prevalence] # subset the feature table to only include features present in at least 10% of samples

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat_filtered))

# make sure names are syntactically valid 
colnames(feat_filtered) <- make.names(colnames(feat_filtered))

# convert feature table to matrix
feat <- as.matrix(feat_filtered)


###############################################################
###   SPARCC - SPARSE CORRELATIONS FOR COMPOSITIONAL DATA   ###
###############################################################

# normal feature table
dim(feat)

### network construction
sparcc_net_tot <- netConstruct(data = feat,
                               group = group,
                               measure = "sparcc", 
                               zeroMethod = "pseudo", # add pseudocount
                               normMethod = "none", # no norm method used for SparCC
                               sparsMethod = "threshold",
                               thresh = 0.3,
                               cores = parallel::detectCores() - 1,
                               verbose = 2,
                               seed = 1234)

# saveRDS(sparcc_net_tot, file = "sparcc_net_tot.rds")
# sparcc_net_tot <- readRDS("sparcc_net_tot.rds")


### degree centrality
#network analysis 
sparcc_analysis_tot_deg <- netAnalyze(sparcc_net_tot, 
                                      clustMethod = "cluster_fast_greedy", 
                                      hubPar = "degree", 
                                      hubQuant = 0.95)
   
# summary of degree centrality analysis
summary(sparcc_analysis_tot_deg, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by degree centrality)
pdf("sparcc_plot_tot_deg_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_tot_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestDegree", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_tot_deg = netCompare(sparcc_analysis_tot_deg,
                                    permTest = TRUE, nPerm = 1000,
                                    adjust = "adaptBH", trueNullMethod = "convest",
                                    cores = parallel::detectCores() - 1, 
                                    storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_tot_deg",
                                    verbose = TRUE, seed = 1234)
summary(sparcc_compare_tot_deg, groupNames = c("resistant", "susceptible"))


### betweenness centrality
# network analysis
sparcc_analysis_tot_bet <- netAnalyze(sparcc_net_tot, 
                                      clustMethod = "cluster_fast_greedy", 
                                      hubPar = "betweenness", 
                                      hubQuant = 0.95)

# summary of betweenness centrality analysis
summary(sparcc_analysis_tot_bet, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by betweenness centrality)
pdf("sparcc_plot_tot_bet_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_tot_bet,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestBetween", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_tot_bet = netCompare(sparcc_analysis_tot_bet,
                                    permTest = TRUE, nPerm = 1000,
                                    adjust = "adaptBH", trueNullMethod = "convest",
                                    cores = parallel::detectCores() - 1,
                                    storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_tot_bet",
                                    verbose = TRUE, seed = 1234)
summary(sparcc_compare_tot_bet, groupNames = c("resistant", "susceptible"))


### closeness centrality
# network analysis
sparcc_analysis_tot_clo <- netAnalyze(sparcc_net_tot, 
                                      clustMethod = "cluster_fast_greedy", 
                                      hubPar = "closeness", 
                                      hubQuant = 0.95)

# summary of closeness centrality analysis
summary(sparcc_analysis_tot_clo, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by closeness centrality)
pdf("sparcc_plot_tot_clo_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_tot_clo,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestClose", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_tot_clo = netCompare(sparcc_analysis_tot_clo,
                                    permTest = TRUE, nPerm = 1000,
                                    adjust = "adaptBH", trueNullMethod = "convest",
                                    cores = parallel::detectCores() - 1,
                                    storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_tot_clo",
                                    verbose = TRUE, seed = 1234)
summary(sparcc_compare_tot_clo, groupNames = c("resistant", "susceptible"))


### eigenvector centrality
# network analysis
sparcc_analysis_tot_eig <- netAnalyze(sparcc_net_tot, 
                                      clustMethod = "cluster_fast_greedy", 
                                      hubPar = "eigenvector", 
                                      hubQuant = 0.95)

# summary of eigenvector centrality analysis
summary(sparcc_analysis_tot_eig, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by eigenvector centrality)
pdf("sparcc_plot_tot_eig_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_tot_eig,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestEigen", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_tot_eig = netCompare(sparcc_analysis_tot_eig,
                                    permTest = TRUE, nPerm = 1000,
                                    adjust = "adaptBH", trueNullMethod = "convest",
                                    cores = parallel::detectCores() - 1,
                                    storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_tot_eig",
                                    verbose = TRUE, seed = 1234)
summary(sparcc_compare_tot_eig, groupNames = c("resistant", "susceptible"))


### centrality measures
sparcc_tot_centralities <- data.frame(feature = names(sparcc_analysis_tot_deg$centralities$degree1),
                                      deg_res = as.numeric(sparcc_analysis_tot_deg$centralities$degree1),
                                      deg_sus = as.numeric(sparcc_analysis_tot_deg$centralities$degree2),
                                      bet_res = as.numeric(sparcc_analysis_tot_bet$centralities$between1),
                                      bet_sus = as.numeric(sparcc_analysis_tot_bet$centralities$between2),
                                      clo_res = as.numeric(sparcc_analysis_tot_clo$centralities$close1),
                                      clo_sus = as.numeric(sparcc_analysis_tot_clo$centralities$close2),
                                      eig_res = as.numeric(sparcc_analysis_tot_eig$centralities$eigenv1),
                                      eig_sus = as.numeric(sparcc_analysis_tot_eig$centralities$eigenv2),
                                      row.names = NULL)

### cluster ids
sparcc_tot_cluster_ids <- data.frame(feature = names(sparcc_analysis_deg$clustering$clust1),
                                     deg_res = as.numeric(sparcc_analysis_tot_deg$clustering$clust1),
                                     deg_sus = as.numeric(sparcc_analysis_tot_deg$clustering$clust2),
                                     bet_res = as.numeric(sparcc_analysis_tot_bet$clustering$clust1),
                                     bet_sus = as.numeric(sparcc_analysis_tot_bet$clustering$clust2),
                                     clo_res = as.numeric(sparcc_analysis_tot_clo$clustering$clust1),
                                     clo_sus = as.numeric(sparcc_analysis_tot_clo$clustering$clust2),
                                     eig_res = as.numeric(sparcc_analysis_tot_eig$clustering$clust1),
                                     eig_sus = as.numeric(sparcc_analysis_tot_eig$clustering$clust2),
                                     row.names = NULL)


###########################################################################
###   NETWORK ANALYSIS WITH MORE STRINGENT FILTERING OF FEATURE TABLE   ###
###########################################################################

# NetCoMi Workflow - more stringent filter than other analyses

# load libraries
library(ggplot2)
library(tidyverse)
library(NetCoMi)

# set.seed
set.seed(1234)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name
group = as.factor(meta$condition) # label as a factor


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

### filter species present in less than 30% of samples
feat_t <- as.data.frame(t(feat)) # transpose feature table
min_prevalence <- 0.30
feat_filtered <- feat_t[, colMeans(feat_t > 0) >= min_prevalence] # subset the feature table to only include features present in at least 10% of samples

### keep only top species by total relative abundance
feat_rel_abun <- feat_filtered/rowSums(feat_filtered) # convert to relative abundance
top_abun <- 200
total_rel_abun <- colSums(feat_rel_abun)
top_taxa <- names(sort(total_rel_abun, decreasing = TRUE))[1:top_abun]
feat_filtered <- feat_filtered[, top_taxa] # filter raw count feature table by top taxa

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat_filtered))

# make sure names are syntactically valid 
colnames(feat_filtered) <- make.names(colnames(feat_filtered))

# convert feature table to matrix
feat <- as.matrix(feat_filtered)


######################################################################
###   FILT | SPARCC - SPARSE CORRELATIONS FOR COMPOSITIONAL DATA   ###
######################################################################

# reduced feature table
dim(feat)

### network construction
sparcc_net <- netConstruct(data = feat,
                           group = group,
                           measure = "sparcc", 
                           zeroMethod = "pseudo", # add pseudocount
                           normMethod = "none", # no norm method used for SparCC
                           sparsMethod = "threshold",
                           thresh = 0.3,
                           cores = parallel::detectCores() - 1,
                           verbose = 2,
                           seed = 1234)

# saveRDS(sparcc_net, file = "sparcc_net.rds")
# sparcc_net <- readRDS("sparcc_net.rds")


### degree centrality
#network analysis 
sparcc_analysis_deg <- netAnalyze(sparcc_net, 
                                  clustMethod = "cluster_fast_greedy", 
                                  hubPar = "degree", 
                                  hubQuant = 0.95)

# summary of degree centrality analysis
summary(sparcc_analysis_deg, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by degree centrality)
pdf("sparcc_plot_deg_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestDegree", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_deg = netCompare(sparcc_analysis_deg,
                                permTest = TRUE, nPerm = 1000,
                                adjust = "adaptBH", trueNullMethod = "convest",
                                cores = parallel::detectCores() - 1, 
                                storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_deg",
                                verbose = TRUE, seed = 1234)
summary(sparcc_compare_deg, groupNames = c("resistant", "susceptible"))


### betweenness centrality
# network analysis
sparcc_analysis_bet <- netAnalyze(sparcc_net, 
                                  clustMethod = "cluster_fast_greedy", 
                                  hubPar = "betweenness", 
                                  hubQuant = 0.95)

# summary of betweenness centrality analysis
summary(sparcc_analysis_bet, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by betweenness centrality)
pdf("sparcc_plot_bet_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_bet,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestBetween", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_bet = netCompare(sparcc_analysis_bet,
                                permTest = TRUE, nPerm = 1000,
                                adjust = "adaptBH", trueNullMethod = "convest",
                                cores = parallel::detectCores() - 1,
                                storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_bet",
                                verbose = TRUE, seed = 1234)
summary(sparcc_compare_bet, groupNames = c("resistant", "susceptible"))


### closeness centrality
# network analysis
sparcc_analysis_clo <- netAnalyze(sparcc_net, 
                                  clustMethod = "cluster_fast_greedy", 
                                  hubPar = "closeness", 
                                  hubQuant = 0.95)

# summary of closeness centrality analysis
summary(sparcc_analysis_clo, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by closeness centrality)
pdf("sparcc_plot_clo_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_clo,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestClose", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_clo = netCompare(sparcc_analysis_clo,
                                permTest = TRUE, nPerm = 1000,
                                adjust = "adaptBH", trueNullMethod = "convest",
                                cores = parallel::detectCores() - 1,
                                storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_clo",
                                verbose = TRUE, seed = 1234)
summary(sparcc_compare_clo, groupNames = c("resistant", "susceptible"))


### eigenvector centrality
# network analysis
sparcc_analysis_eig <- netAnalyze(sparcc_net, 
                                  clustMethod = "cluster_fast_greedy", 
                                  hubPar = "eigenvector", 
                                  hubQuant = 0.95)

# summary of eigenvector centrality analysis
summary(sparcc_analysis_eig, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by eigenvector centrality)
pdf("sparcc_plot_eig_node50.pdf", width = 60, height = 30)
plot(sparcc_analysis_eig,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestEigen", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
sparcc_compare_eig = netCompare(sparcc_analysis_eig,
                                permTest = TRUE, nPerm = 1000,
                                adjust = "adaptBH", trueNullMethod = "convest",
                                cores = parallel::detectCores() - 1,
                                storeAssoPerm = TRUE, fileStoreAssoPerm = "sparcc_compare_eig",
                                verbose = TRUE, seed = 1234)
summary(sparcc_compare_eig, groupNames = c("resistant", "susceptible"))


### centrality measures
sparcc_centralities <- data.frame(feature = names(sparcc_analysis_deg$centralities$degree1),
                                  deg_res = as.numeric(sparcc_analysis_deg$centralities$degree1),
                                  deg_sus = as.numeric(sparcc_analysis_deg$centralities$degree2),
                                  bet_res = as.numeric(sparcc_analysis_bet$centralities$between1),
                                  bet_sus = as.numeric(sparcc_analysis_bet$centralities$between2),
                                  clo_res = as.numeric(sparcc_analysis_clo$centralities$close1),
                                  clo_sus = as.numeric(sparcc_analysis_clo$centralities$close2),
                                  eig_res = as.numeric(sparcc_analysis_eig$centralities$eigenv1),
                                  eig_sus = as.numeric(sparcc_analysis_eig$centralities$eigenv2),
                                  row.names = NULL)

### cluster ids
sparcc_cluster_ids <- data.frame(feature = names(sparcc_analysis_deg$clustering$clust1),
                                 deg_res = as.numeric(sparcc_analysis_deg$clustering$clust1),
                                 deg_sus = as.numeric(sparcc_analysis_deg$clustering$clust2),
                                 bet_res = as.numeric(sparcc_analysis_bet$clustering$clust1),
                                 bet_sus = as.numeric(sparcc_analysis_bet$clustering$clust2),
                                 clo_res = as.numeric(sparcc_analysis_clo$clustering$clust1),
                                 clo_sus = as.numeric(sparcc_analysis_clo$clustering$clust2),
                                 eig_res = as.numeric(sparcc_analysis_eig$clustering$clust1),
                                 eig_sus = as.numeric(sparcc_analysis_eig$clustering$clust2),
                                 row.names = NULL)


#######################################################################################
###   FILT | CCLASSO - CORRELATION INFERENCE FOR COMPOSITIONAL DATA THROUGH LASSO   ###
#######################################################################################

# reduced feature table
dim(feat)

### network construction
cclasso_net <- netConstruct(data = feat,
                            group = group,
                            measure = "cclasso", 
                            zeroMethod = "pseudo", # adds pseudocount of 1
                            normMethod = "none", # no norm method used for CCLasso
                            sparsMethod = "none", # automatically produces sparse networks
                            cores = parallel::detectCores() - 1,
                            verbose = 2,
                            seed = 1234)

# saveRDS(cclasso_net, file = "cclasso_net.rds")
# cclasso_net <- readRDS("cclasso_net.rds")


### degree centrality
#network analysis 
cclasso_analysis_deg <- netAnalyze(cclasso_net, 
                                   clustMethod = "cluster_fast_greedy", 
                                   hubPar = "degree", 
                                   hubQuant = 0.95)

# summary of degree centrality analysis
summary(cclasso_analysis_deg, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by degree centrality)
pdf("cclasso_plot_deg_node50.pdf", width = 60, height = 30)
plot(cclasso_analysis_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestDegree", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
cclasso_compare_deg = netCompare(cclasso_analysis_deg,
                                 permTest = TRUE, nPerm = 1000,
                                 adjust = "adaptBH", trueNullMethod = "convest",
                                 cores = parallel::detectCores() - 1, 
                                 storeAssoPerm = TRUE, fileStoreAssoPerm = "cclasso_compare_deg",
                                 verbose = TRUE, seed = 1234)
summary(cclasso_compare_deg, groupNames = c("resistant", "susceptible"))


### betweenness centrality
# network analysis
cclasso_analysis_bet <- netAnalyze(cclasso_net, 
                                   clustMethod = "cluster_fast_greedy", 
                                   hubPar = "betweenness", 
                                   hubQuant = 0.95)

# summary of betweenness centrality analysis
summary(cclasso_analysis_bet, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by betweenness centrality)
pdf("cclasso_plot_bet_node50.pdf", width = 60, height = 30)
plot(cclasso_analysis_bet,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestBetween", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
cclasso_compare_bet = netCompare(cclasso_analysis_bet,
                                 permTest = TRUE, nPerm = 1000,
                                 adjust = "adaptBH", trueNullMethod = "convest",
                                 cores = parallel::detectCores() - 1,
                                 storeAssoPerm = TRUE, fileStoreAssoPerm = "cclasso_compare_bet",
                                 verbose = TRUE, seed = 1234)
summary(cclasso_compare_bet, groupNames = c("resistant", "susceptible"))


### closeness centrality
# network analysis
cclasso_analysis_clo <- netAnalyze(cclasso_net, 
                                   clustMethod = "cluster_fast_greedy", 
                                   hubPar = "closeness", 
                                   hubQuant = 0.95)

# summary of closeness centrality analysis
summary(cclasso_analysis_clo, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by closeness centrality)
pdf("cclasso_plot_clo_node50.pdf", width = 60, height = 30)
plot(cclasso_analysis_clo,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestClose", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
cclasso_compare_clo = netCompare(cclasso_analysis_clo,
                                 permTest = TRUE, nPerm = 1000,
                                 adjust = "adaptBH", trueNullMethod = "convest",
                                 cores = parallel::detectCores() - 1,
                                 storeAssoPerm = TRUE, fileStoreAssoPerm = "cclasso_compare_clo",
                                 verbose = TRUE, seed = 1234)
summary(cclasso_compare_clo, groupNames = c("resistant", "susceptible"))


### eigenvector centrality
# network analysis
cclasso_analysis_eig <- netAnalyze(cclasso_net, 
                                   clustMethod = "cluster_fast_greedy", 
                                   hubPar = "eigenvector", 
                                   hubQuant = 0.95)

# summary of eigenvector centrality analysis
summary(cclasso_analysis_eig, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by eigenvector centrality)
pdf("cclasso_plot_eig_node50.pdf", width = 60, height = 30)
plot(cclasso_analysis_eig,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestEigen", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
cclasso_compare_eig = netCompare(cclasso_analysis_eig,
                                 permTest = TRUE, nPerm = 1000,
                                 adjust = "adaptBH", trueNullMethod = "convest",
                                 cores = parallel::detectCores() - 1,
                                 storeAssoPerm = TRUE, fileStoreAssoPerm = "cclasso_compare_eig",
                                 verbose = TRUE, seed = 1234)
summary(cclasso_compare_eig, groupNames = c("resistant", "susceptible"))


### centrality measures
cclasso_centralities <- data.frame(feature = names(cclasso_analysis_deg$centralities$degree1),
                                   deg_res = as.numeric(cclasso_analysis_deg$centralities$degree1),
                                   deg_sus = as.numeric(cclasso_analysis_deg$centralities$degree2),
                                   bet_res = as.numeric(cclasso_analysis_bet$centralities$between1),
                                   bet_sus = as.numeric(cclasso_analysis_bet$centralities$between2),
                                   clo_res = as.numeric(cclasso_analysis_clo$centralities$close1),
                                   clo_sus = as.numeric(cclasso_analysis_clo$centralities$close2),
                                   eig_res = as.numeric(cclasso_analysis_eig$centralities$eigenv1),
                                   eig_sus = as.numeric(cclasso_analysis_eig$centralities$eigenv2),
                                   row.names = NULL)

### cluster ids
cclasso_cluster_ids <- data.frame(feature = names(cclasso_analysis_deg$clustering$clust1),
                                  deg_res = as.numeric(cclasso_analysis_deg$clustering$clust1),
                                  deg_sus = as.numeric(cclasso_analysis_deg$clustering$clust2),
                                  bet_res = as.numeric(cclasso_analysis_bet$clustering$clust1),
                                  bet_sus = as.numeric(cclasso_analysis_bet$clustering$clust2),
                                  clo_res = as.numeric(cclasso_analysis_clo$clustering$clust1),
                                  clo_sus = as.numeric(cclasso_analysis_clo$clustering$clust2),
                                  eig_res = as.numeric(cclasso_analysis_eig$clustering$clust1),
                                  eig_sus = as.numeric(cclasso_analysis_eig$clustering$clust2),
                                 row.names = NULL)


########################################################################################################################
###   FILT | SPIECEASI - SPARSE INVERSE COVARIANCE ESTIMATION FOR ECOLOGICAL ASSOCIATION AND STATISTICAL INFERENCE   ###
########################################################################################################################

# reduced feature table
dim(feat)

### network construction
spieceasi_net <- netConstruct(data = feat,
                              group = group,
                              measure = "spieceasi", 
                              zeroMethod = "none", # no pseudocount added for SpiecEasi
                              normMethod = "none", # SpiecEasi applies CLR internally
                              sparsMethod = "none", # automatically produces sparse networks
                              cores = parallel::detectCores() - 1,
                              verbose = 2,
                              seed = 1234)

# saveRDS(spieceasi_net, file = "spieceasi_net.rds")
# spieceasi_net <- readRDS("spieceasi_net.rds")


### degree centrality
# network analysis 
spieceasi_analysis_deg <- netAnalyze(spieceasi_net, 
                                     clustMethod = "cluster_fast_greedy", 
                                     hubPar = "degree", 
                                     hubQuant = 0.95)

# summary of degree centrality analysis
summary(spieceasi_analysis_deg, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by degree centrality)
pdf("spieceasi_plot_deg_node50.pdf", width = 60, height = 30)
plot(spieceasi_analysis_deg,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestDegree", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
spieceasi_compare_deg = netCompare(spieceasi_analysis_deg,
                                   permTest = TRUE, nPerm = 1000,
                                   adjust = "adaptBH", trueNullMethod = "convest",
                                   cores = parallel::detectCores() - 1, 
                                   storeAssoPerm = TRUE, fileStoreAssoPerm = "spieceasi_compare_deg",
                                   verbose = TRUE, seed = 1234)
summary(spieceasi_compare_deg, groupNames = c("resistant", "susceptible"))


### betweenness centrality
# network analysis
spieceasi_analysis_bet <- netAnalyze(spieceasi_net, 
                                     clustMethod = "cluster_fast_greedy", 
                                     hubPar = "betweenness", 
                                     hubQuant = 0.95)

# summary of betweenness centrality analysis
summary(spieceasi_analysis_bet, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by betweenness centrality)
pdf("spieceasi_plot_bet_node50.pdf", width = 60, height = 30)
plot(spieceasi_analysis_bet,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestBetween", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
spieceasi_compare_bet = netCompare(spieceasi_analysis_bet,
                                   permTest = TRUE, nPerm = 1000,
                                   adjust = "adaptBH", trueNullMethod = "convest",
                                   cores = parallel::detectCores() - 1,
                                   storeAssoPerm = TRUE, fileStoreAssoPerm = "spieceasi_compare_bet",
                                   verbose = TRUE, seed = 1234)
summary(spieceasi_compare_bet, groupNames = c("resistant", "susceptible"))


### closeness centrality
# network analysis
spieceasi_analysis_clo <- netAnalyze(spieceasi_net, 
                                     clustMethod = "cluster_fast_greedy", 
                                     hubPar = "closeness", 
                                     hubQuant = 0.95)

# summary of closeness centrality analysis
summary(spieceasi_analysis_clo, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by closeness centrality)
pdf("spieceasi_plot_clo_node50.pdf", width = 60, height = 30)
plot(spieceasi_analysis_clo,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestClose", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
spieceasi_compare_clo = netCompare(spieceasi_analysis_clo,
                                   permTest = TRUE, nPerm = 1000,
                                   adjust = "adaptBH", trueNullMethod = "convest",
                                   cores = parallel::detectCores() - 1,
                                   storeAssoPerm = TRUE, fileStoreAssoPerm = "spieceasi_compare_clo",
                                   verbose = TRUE, seed = 1234)
summary(spieceasi_compare_clo, groupNames = c("resistant", "susceptible"))


### eigenvector centrality
# network analysis
spieceasi_analysis_eig <- netAnalyze(spieceasi_net, 
                                     clustMethod = "cluster_fast_greedy", 
                                     hubPar = "eigenvector", 
                                     hubQuant = 0.95)

# summary of eigenvector centrality analysis
summary(spieceasi_analysis_eig, groupNames = c("resistant", "susceptible"))

# network plot of top 50 nodes (by eigenvector centrality)
pdf("spieceasi_plot_eig_node50.pdf", width = 60, height = 30)
plot(spieceasi_analysis_eig,
     layout = "spring", repulsion = 0.8, 
     sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth",
     nodeFilter = "highestEigen", nodeFilterPar = 50,
     # edgeFilter = "highestWeight", edgeFilterPar = 3000,
     nodeSize = "degree", nodeSizeSpread = 2, nodeColor = "cluster", nodeTransp = 50,
     borderCol = "gray40", borderWidth = 0.5, 
     hubBorderCol = "gray20", hubBorderWidth = 0.5,
     edgeWidth = 0.5, edgeTranspLow = 50, edgeTranspHigh = 30,
     groupNames = c("resistant", "susceptible"),
     labelScale = FALSE,cexLabels = 1, cexTitle = 4,
     labelFont = 1, hubLabelFont = 1)
dev.off()

# network comparison with permutation
spieceasi_compare_eig = netCompare(spieceasi_analysis_eig,
                                   permTest = TRUE, nPerm = 1000,
                                   adjust = "adaptBH", trueNullMethod = "convest",
                                   cores = parallel::detectCores() - 1,
                                   storeAssoPerm = TRUE, fileStoreAssoPerm = "spieceasi_compare_eig",
                                   verbose = TRUE, seed = 1234)
summary(spieceasi_compare_eig, groupNames = c("resistant", "susceptible"))


### centrality measures
spieceasi_centralities <- data.frame(feature = names(spieceasi_analysis_deg$centralities$degree1),
                                     deg_res = as.numeric(spieceasi_analysis_deg$centralities$degree1),
                                     deg_sus = as.numeric(spieceasi_analysis_deg$centralities$degree2),
                                     bet_res = as.numeric(spieceasi_analysis_bet$centralities$between1),
                                     bet_sus = as.numeric(spieceasi_analysis_bet$centralities$between2),
                                     clo_res = as.numeric(spieceasi_analysis_clo$centralities$close1),
                                     clo_sus = as.numeric(spieceasi_analysis_clo$centralities$close2),
                                     eig_res = as.numeric(spieceasi_analysis_eig$centralities$eigenv1),
                                     eig_sus = as.numeric(spieceasi_analysis_eig$centralities$eigenv2),
                                     row.names = NULL)

### cluster ids
spieceasi_cluster_ids <- data.frame(feature = names(spieceasi_analysis_deg$clustering$clust1),
                                    deg_res = as.numeric(spieceasi_analysis_deg$clustering$clust1),
                                    deg_sus = as.numeric(spieceasi_analysis_deg$clustering$clust2),
                                    bet_res = as.numeric(spieceasi_analysis_bet$clustering$clust1),
                                    bet_sus = as.numeric(spieceasi_analysis_bet$clustering$clust2),
                                    clo_res = as.numeric(spieceasi_analysis_clo$clustering$clust1),
                                    clo_sus = as.numeric(spieceasi_analysis_clo$clustering$clust2),
                                    eig_res = as.numeric(spieceasi_analysis_eig$clustering$clust1),
                                    eig_sus = as.numeric(spieceasi_analysis_eig$clustering$clust2),
                                    row.names = NULL)

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] NetCoMi_1.2.0    SpiecEasi_1.99.0 lubridate_1.9.4  forcats_1.0.0    stringr_1.5.2   
# [6] dplyr_1.1.4      purrr_1.1.0      readr_2.1.5      tidyr_1.3.1      tibble_3.3.0    
# [11] tidyverse_2.0.0  ggplot2_4.0.0   
# 
# loaded via a namespace (and not attached):
#   [1] ggtext_0.1.2                    fs_1.6.6                        matrixStats_1.5.0              
# [4] DirichletMultinomial_1.50.0     httr_1.4.7                      RColorBrewer_1.1-3             
# [7] doParallel_1.0.17               dynamicTreeCut_1.63-1           tools_4.5.0                    
# [10] backports_1.5.0                 R6_2.6.1                        vegan_2.7-1                    
# [13] lazyeval_0.2.2                  mgcv_1.9-3                      rhdf5filters_1.20.0            
# [16] permute_0.9-8                   withr_3.0.2                     gridExtra_2.3                  
# [19] preprocessCore_1.70.0           fdrtool_1.2.18                  qgraph_1.9.8                   
# [22] WGCNA_1.73                      cli_3.6.5                       Biobase_2.68.0                 
# [25] textshaping_1.0.1               sandwich_3.1-1                  slam_0.1-55                    
# [28] mvtnorm_1.3-3                   S7_0.2.0                        pbapply_1.7-4                  
# [31] pbivnorm_0.6.0                  systemfonts_1.2.3               yulab.utils_0.2.0              
# [34] SPRING_1.0.4                    foreign_0.8-90                  dichromat_2.0-0.1              
# [37] scater_1.35.0                   decontam_1.28.0                 parallelly_1.45.1              
# [40] limma_3.64.3                    readxl_1.4.5                    filematrix_1.3                 
# [43] fillpattern_1.0.2               huge_1.3.5                      VGAM_1.1-13                    
# [46] rstudioapi_0.17.1               impute_1.82.0                   RSQLite_2.4.2                  
# [49] generics_0.1.4                  shape_1.4.6.1                   gtools_3.9.5                   
# [52] rbiom_2.2.1                     GO.db_3.21.0                    Matrix_1.7-3                   
# [55] biomformat_1.36.0               ggbeeswarm_0.7.2                S4Vectors_0.46.0               
# [58] DECIPHER_3.4.0                  abind_1.4-8                     lifecycle_1.0.4                
# [61] multcomp_1.4-28                 SummarizedExperiment_1.38.1     rhdf5_2.52.1                   
# [64] SparseArray_1.8.1               lavaan_0.6-20                   grid_4.5.0                     
# [67] blob_1.2.4                      crayon_1.5.3                    lattice_0.22-7                 
# [70] beachmat_2.24.0                 KEGGREST_1.48.1                 pillar_1.11.0                  
# [73] knitr_1.50                      GenomicRanges_1.60.0            pulsar_0.3.11                  
# [76] estimability_1.5.1              corpcor_1.6.10                  codetools_0.2-20               
# [79] glue_1.8.0                      data.table_1.17.8               MultiAssayExperiment_1.34.0    
# [82] vctrs_0.6.5                     png_0.1-8                       treeio_1.32.0                  
# [85] Rdpack_2.6.4                    cellranger_1.1.0                gtable_0.3.6                   
# [88] cachem_1.1.0                    xfun_0.52                       rbibutils_2.3                  
# [91] S4Arrays_1.8.1                  coda_0.19-4.1                   pcaPP_2.0-5                    
# [94] survival_3.8-3                  SingleCellExperiment_1.30.1     iterators_1.0.14               
# [97] statmod_1.5.0                   bluster_1.18.0                  TH.data_1.1-3                  
# [100] nlme_3.1-168                    phyloseq_1.52.0                 bit64_4.6.0-1                  
# [103] GenomeInfoDb_1.44.3             irlba_2.3.5.1                   vipor_0.4.7                    
# [106] rpart_4.1.24                    mixedCCA_1.6.2                  colorspace_2.1-1               
# [109] BiocGenerics_0.54.0             DBI_1.2.3                       Hmisc_5.2-3                    
# [112] nnet_7.3-20                     ade4_1.7-23                     mnormt_2.1.1                   
# [115] tidyselect_1.2.1                emmeans_1.11.2                  bit_4.6.0                      
# [118] compiler_4.5.0                  glmnet_4.1-10                   htmlTable_2.4.3                
# [121] BiocNeighbors_2.2.0             xml2_1.3.8                      DelayedArray_0.34.1            
# [124] checkmate_2.3.2                 scales_1.4.0                    psych_2.5.6                    
# [127] quadprog_1.5-8                  digest_0.6.37                   rmarkdown_2.29                 
# [130] XVector_0.48.0                  jpeg_0.1-11                     htmltools_0.5.8.1              
# [133] pkgconfig_2.0.3                 base64enc_0.1-3                 sparseMatrixStats_1.20.0       
# [136] orca_1.1-3                      MatrixGenerics_1.20.0           fastmap_1.2.0                  
# [139] rlang_1.1.6                     htmlwidgets_1.6.4               UCSC.utils_1.4.0               
# [142] DelayedMatrixStats_1.30.0       farver_2.1.2                    zoo_1.8-14                     
# [145] jsonlite_2.0.0                  BiocParallel_1.42.1             BiocSingular_1.24.0            
# [148] magrittr_2.0.4                  Formula_1.2-5                   scuttle_1.18.0                 
# [151] GenomeInfoDbData_1.2.14         patchwork_1.3.1                 Rhdf5lib_1.30.0                
# [154] Rcpp_1.1.0                      ape_5.8-1                       ggnewscale_0.5.2               
# [157] viridis_0.6.5                   stringi_1.8.7                   rootSolve_1.8.2.4              
# [160] MASS_7.3-65                     plyr_1.8.9                      parallel_4.5.0                 
# [163] ggrepel_0.9.6                   doSNOW_1.0.20                   Biostrings_2.76.0              
# [166] splines_4.5.0                   gridtext_0.1.5                  multtest_2.64.0                
# [169] hms_1.1.3                       igraph_2.1.4                    fastcluster_1.3.0              
# [172] reshape2_1.4.4                  stats4_4.5.0                    ScaledMatrix_1.16.0            
# [175] evaluate_1.0.4                  tzdb_0.5.0                      foreach_1.5.2                  
# [178] BiocBaseUtils_1.10.0            rsvd_1.0.5                      xtable_1.8-4                   
# [181] tidytree_0.4.6                  glasso_1.11                     viridisLite_0.4.2              
# [184] ragg_1.4.0                      snow_0.4-4                      memoise_2.0.1                  
# [187] beeswarm_0.4.0                  AnnotationDbi_1.70.0            IRanges_2.42.0                 
# [190] cluster_2.1.8.1                 corrplot_0.95                   TreeSummarizedExperiment_2.16.1
# [193] timechange_0.3.0                mia_1.16.0

