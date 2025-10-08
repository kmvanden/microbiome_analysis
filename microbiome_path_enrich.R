# Workflow for Pathway Enrichment Analysis

# load libraries
library(ggplot2)
library(tidyverse)
library(readxl)
library(impute)
library(limma)
library(clusterProfiler)
library(KEGGREST)
library(gage)
library(pathview)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)

# subset meta to samples present in metabolomics data
sample_key <- read_excel("metabo_batch.xlsx", sheet = "Sample Meta Data") # sample name key
sample_key <- sample_key %>% mutate(CLIENT_SAMPLE_ID = sub("KMV_", "SMS_", CLIENT_SAMPLE_ID)) # rename sample ids
all(sample_key$CLIENT_SAMPLE_ID %in% meta$sample_id)
meta <- meta %>%
  left_join(sample_key %>% dplyr::select(CLIENT_SAMPLE_ID, PARENT_SAMPLE_NAME),
            by = c("sample_id" = "CLIENT_SAMPLE_ID")) %>% # add PARENT_SAMPLE_NAME to meta
  filter(!is.na(PARENT_SAMPLE_NAME)) # subset

group <- meta$condition # condition as a factor


### feature table
metabo <- read_excel("metabo_batch.xlsx", sheet = "Batch-normalized Data") # batch-normalized feature table
all(metabo$PARENT_SAMPLE_NAME %in% meta$PARENT_SAMPLE_NAME)

# replace names in PARENT_SAMPLE_NAME column in metabo with names from sample_name (from meta)
metabo <- metabo %>%
  left_join(meta %>% dplyr::select(PARENT_SAMPLE_NAME, sample_name), by = "PARENT_SAMPLE_NAME") %>%
  mutate(PARENT_SAMPLE_NAME = sample_name) %>%
  dplyr::select(-sample_name) %>%
  slice(match(meta$sample_name, PARENT_SAMPLE_NAME)) %>% # match order of samples in meta
  column_to_rownames(var = "PARENT_SAMPLE_NAME") # move sample id column to rownames
all(rownames(metabo) == meta$sample_name)


# replace CHEM_ID (colnames of metabo) with names from KEGG (from metabo_key)
metabo_key <- read_excel("metabo_batch.xlsx", sheet = "Chemical Annotation")
all(colnames(metabo) == metabo_key$CHEM_ID)
name_map <- metabo_key %>% dplyr::select(CHEM_ID, KEGG) %>% deframe() 
colnames(metabo) <- name_map[colnames(metabo)] %>% 
  str_split(",") %>% map_chr(1) # keep first KEGG id only

# remove metabolites without a KEGG id
metabo <- metabo[, !is.na(colnames(metabo))]

# assess missingness (zeros listed as NAs in data.frame)
missing_per_feature <- colMeans(is.na(metabo)) * 100 # missingness per feature
summary(missing_per_feature)
hist(missing_per_feature, breaks = 50, main = "Percent missing per metabolite", xlab = "Percent missing")

# remove features with greater than 30% missing values
metabo <- metabo[, missing_per_feature <= 30]

# impute remaining NAs with kNN
metabo_imputed <- impute.knn(as.matrix(t(metabo)))$data # impute.knn expects features as rows
metabo_imputed <- t(metabo_imputed)

# log2 transformation (with pseudocount)
metabo_log <- log2(metabo_imputed + 1e-6)


############################################################
########   DIFFERENTIAL ABUNDANCE ANALYSIS - LIMMA   #######
############################################################

# design matrix
# group <- factor(metadata$Group)
design <- model.matrix(~0 + group) # "groupdisease" "grouphealthy"

# transpose feature table
metabo_t <- t(metabo_log)

# fit linear model
limma_fit <- lmFit(metabo_t, design)

# define contrast
contrast <- makeContrasts(disease_vs_healthy = groupdisease - grouphealthy, levels = design)
fit <- contrasts.fit(limma_fit, contrast) # apply contrast to linear model
fit_bayes <- eBayes(fit) # empirical Bayes smoothing

# get results (logFC, p-value, padj)
results <- topTable(fit_bayes, coef = "disease_vs_healthy", adjust.method = "BH", number = Inf)
results$KEGG <- rownames(results)


########################################################
########   GSEA-LIKE PATHWAY ENRICHMENT - GAGE   #######
########################################################

# named numeric vector of log2 fold changes
logfc_vector <- results$logFC
names(logfc_vector) <- results$KEGG

# create pathway list
map <- keggLink("pathway", "compound")
map_table <- data.frame(TERM = unname(map),
                        GENE = sub("^cpd:", "", names(map)),
                        stringsAsFactors = FALSE)

kegg_path_sets <- split(map_table$GENE, map_table$TERM)

# run gage
gage_res <- gage(logfc_vector, gsets = kegg_path_sets)

# enriched pathways
head(as.data.frame(gage_res$greater), 10) # enriched in disease
head(as.data.frame(gage_res$less), 10) # enriched in healthy

# visualize with pathview
pathview(cpd.data = logfc_vector,
         pathway.id = "05230",
         species = "hsa",
         cpd.idtype = "kegg",
         kegg.native = TRUE,
         kegg.dir = ".")

pathview(cpd.data = logfc_vector,
         pathway.id = "00250",
         species = "hsa",
         cpd.idtype = "kegg",
         kegg.native = TRUE,
         kegg.dir = ".")

pathview(cpd.data = logfc_vector,
         pathway.id = "00310",
         species = "hsa",
         cpd.idtype = "kegg",
         kegg.native = TRUE,
         kegg.dir = ".")

pathview(cpd.data = logfc_vector,
         pathway.id = "04976",
         species = "hsa",
         cpd.idtype = "kegg",
         kegg.native = TRUE,
         kegg.dir = ".")


# data.frame for disease
gage_inc <- as.data.frame(gage_res$greater)
gage_inc$KEGG_id <- sub("^path:map", "", rownames(gage_inc))
gage_inc <- gage_inc %>% arrange(KEGG_id)

# data.frame for healthy
gage_dec <- as.data.frame(gage_res$less)
gage_dec$KEGG_id <- sub("^path:map", "", rownames(gage_dec))
gage_dec <- gage_dec %>% arrange(KEGG_id)

# pathways in same order
all(gage_inc$KEGG_id == gage_inc$KEGG_id)

# data.frame of KEGG ids and pathway names
pathways <- keggList("pathway", "hsa")
path_df <- data.frame(KEGG_id = sub("^hsa", "", names(pathways)),
                      path_name = sub(" - Homo sapiens \\(human\\)$", "", as.character(pathways)),
                      stringsAsFactors = FALSE)

# add pathway names to data.frames
gage_inc_name <- left_join(gage_inc, path_df, by = "KEGG_id") %>% filter(!is.na(path_name)) %>% arrange(desc(stat.mean))
gage_dec_name <- left_join(gage_dec, path_df, by = "KEGG_id") %>% filter(!is.na(path_name)) %>% arrange(stat.mean)

# plot stat.mean (mean of log2 fold changes for the metabolites in each KEGG pathway set)
ggplot(gage_inc_name[1:10, ], aes(x = reorder(path_name, stat.mean), y = stat.mean)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Top enriched pathways (disease)", y = "Stat.mean", x = "Pathway") +
  theme_minimal(base_size = 12) + theme(axis.text.y = element_text(size = 10))

ggplot(gage_dec_name[1:10, ], aes(x = reorder(path_name, -stat.mean), y = stat.mean)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Top enriched pathways (healthy)", y = "Stat.mean", x = "Pathway") +
  theme_minimal(base_size = 12) + theme(axis.text.y = element_text(size = 10))


##################################################################
########   OVERREPRESENTATION ANALYSIS - CLUSTERPROFILER   #######
##################################################################

# significant metabolites
sig_metabo_disease <- results %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  pull(KEGG) # metabolites upregualted in disease

sig_metabo_healthy <- results %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  pull(KEGG) # metabolites upregualtes in healthy

# background metabolites
back_metabo <- results$KEGG

# mapping table
map <- keggLink("pathway", "compound")
map_table <- data.frame(TERM = unname(map),
                        GENE = sub("^cpd:", "", names(map)),
                        stringsAsFactors = FALSE)

# enrichment analysis for disease
enrich_metabo_disease <- enricher(gene = sig_metabo_disease,
                                  universe = back_metabo,
                                  TERM2GENE = map_table,
                                  pvalueCutoff = 1)

# enrichment analysis for healthy
enrich_metabo_healthy <- enricher(gene = sig_metabo_healthy,
                                  universe = back_metabo,
                                  TERM2GENE = map_table,
                                  pvalueCutoff = 1)


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
#   [1] impute_1.82.0          pathview_1.48.0        gage_2.58.0            KEGGREST_1.48.1       
# [5] clusterProfiler_4.16.0 limma_3.64.3           readxl_1.4.5           lubridate_1.9.4       
# [9] forcats_1.0.1          stringr_1.5.2          dplyr_1.1.4            purrr_1.1.0           
# [13] readr_2.1.5            tidyr_1.3.1            tibble_3.3.0           tidyverse_2.0.0       
# [17] ggplot2_4.0.0         
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-9            DBI_1.2.3               gson_0.1.0              rlang_1.1.6            
# [5] magrittr_2.0.4          DOSE_4.2.0              compiler_4.5.0          RSQLite_2.4.3          
# [9] png_0.1-8               vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3        
# [13] crayon_1.5.3            fastmap_1.2.0           XVector_0.48.0          labeling_0.4.3         
# [17] tzdb_0.5.0              enrichplot_1.28.4       KEGGgraph_1.68.0        graph_1.86.0           
# [21] UCSC.utils_1.4.0        bit_4.6.0               cachem_1.1.0            aplot_0.2.9            
# [25] GenomeInfoDb_1.44.3     jsonlite_2.0.0          blob_1.2.4              BiocParallel_1.42.2    
# [29] parallel_4.5.0          R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3     
# [33] cellranger_1.1.0        GOSemSim_2.34.0         Rcpp_1.1.0              ggtangle_0.0.7         
# [37] R.utils_2.13.0          IRanges_2.42.0          Matrix_1.7-4            splines_4.5.0          
# [41] igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1        qvalue_2.40.0          
# [45] rstudioapi_0.17.1       dichromat_2.0-0.1       codetools_0.2-20        curl_7.0.0             
# [49] lattice_0.22-7          plyr_1.8.9              treeio_1.32.0           Biobase_2.68.0         
# [53] withr_3.0.2             S7_0.2.0                gridGraphics_0.5-1      Biostrings_2.76.0      
# [57] pillar_1.11.1           ggtree_3.16.3           stats4_4.5.0            ggfun_0.2.0            
# [61] generics_0.1.4          RCurl_1.98-1.17         S4Vectors_0.46.0        hms_1.1.3              
# [65] scales_1.4.0            tidytree_0.4.6          glue_1.8.0              lazyeval_0.2.2         
# [69] tools_4.5.0             data.table_1.17.8       fgsea_1.34.2            XML_3.99-0.19          
# [73] fs_1.6.6                fastmatch_1.1-6         cowplot_1.2.0           grid_4.5.0             
# [77] ape_5.8-1               AnnotationDbi_1.70.0    colorspace_2.1-2        nlme_3.1-168           
# [81] GenomeInfoDbData_1.2.14 patchwork_1.3.2         cli_3.6.5               rappdirs_0.3.3         
# [85] Rgraphviz_2.52.0        gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.2.1      
# [89] digest_0.6.37           BiocGenerics_0.54.0     ggrepel_0.9.6           ggplotify_0.1.3        
# [93] org.Hs.eg.db_3.21.0     farver_2.1.2            memoise_2.0.1           R.oo_1.27.1            
# [97] lifecycle_1.0.4         httr_1.4.7              GO.db_3.21.0            statmod_1.5.0          
# [101] bit64_4.6.0-1 

