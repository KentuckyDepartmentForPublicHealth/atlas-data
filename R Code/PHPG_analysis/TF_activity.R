# ================================
# 1. Load Required Libraries
# ================================
library(dorothea)
library(decoupleR)
library(dplyr)
library(tidyr)
library(tibble)

# ================================
# 2. Load and Prepare DoRothEA Regulon (A/B/C confidence)
# ================================
regulon <- get_dorothea(organism = "human", levels = c("A", "B", "C")) %>%
  mutate(target = toupper(target))

# ================================
# 3. Prepare Expression Matrix (VST-normalized TCGA)
# ================================
# Assumes you have:
# vst_data       = matrix of ENSEMBL IDs (genes x samples)
# exp_data       = SummarizedExperiment with gene symbols in rowData()
# annotation_tcga = sample annotations with Cluster column

# Convert Ensembl IDs to SYMBOL
ensembl_ids <- gsub("\\..*", "", rownames(vst_data))
gene_symbols <- rowData(exp_data)$gene_name

symbol_map <- data.frame(
  Ensembl = ensembl_ids,
  Symbol = toupper(gene_symbols),  # enforce consistent casing
  stringsAsFactors = FALSE
)
symbol_map <- symbol_map[!duplicated(symbol_map$Ensembl) & !is.na(symbol_map$Symbol), ]

# Subset vst_data to SYMBOL matrix
vst_symbol <- vst_data[symbol_map$Ensembl, ]
rownames(vst_symbol) <- symbol_map$Symbol
vst_symbol <- vst_symbol[!duplicated(rownames(vst_symbol)), ]

# ================================
# 4. Filter Regulon for Targets in vst_symbol
# ================================
regulon_filt <- regulon %>%
  filter(target %in% rownames(vst_symbol)) %>%
  group_by(source) %>%
  filter(n() >= 5) %>%
  ungroup()

# ================================
# 5. Subset vst_symbol to Valid Regulon Targets
# ================================
expr_mat <- vst_symbol[rownames(vst_symbol) %in% regulon_filt$target, ]

# ================================
# 6. Run TF Activity Estimation via Weighted Mean
# ================================
tf_activity <- run_wmean(
  mat = expr_mat,
  net = regulon_filt,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# ================================
# 7. Pivot Activity Output: TF x Sample
# ================================
# Collapse duplicates using mean
activity_matrix <- tf_activity %>%
  select(source, condition, score) %>%
  group_by(source, condition) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source")



# ================================
# 8. Merge with Cluster Annotations
# ================================
activity_df <- as.data.frame(t(activity_matrix))  # Samples x TFs
activity_df$Sample <- rownames(activity_df)

activity_df <- merge(activity_df, annotation_tcga, by = "Sample")

# ================================
# 9. Identify Top 5 TFs per Cluster
# ================================
top_tfs <- activity_df %>%
  pivot_longer(-c(Sample, Cluster), names_to = "TF", values_to = "Activity") %>%
  group_by(Cluster, TF) %>%
  summarise(mean_activity = mean(Activity), .groups = "drop") %>%
  group_by(Cluster) %>%
  slice_max(mean_activity, n = 25)

# ================================
# 10. Output
# ================================
print(top_tfs)




# ================================
# 3. Prepare Expression Matrix (Microarray: exp)
# ================================
# Assumes:
# - `exp` is samples x genes, with gene symbols as column names
# - Final 3 columns of exp are metadata and must be excluded
# - `annotation` contains Cluster assignments with rownames as Sample IDs

# Remove metadata columns
expr_micro <- exp[, 1:(ncol(exp) - 3)]
expr_micro <- t(expr_micro)  # Transpose to gene x sample

# Ensure gene symbols are uppercase and set as rownames
rownames(expr_micro) <- toupper(rownames(expr_micro))

# ================================
# 4. Filter Regulon for Targets in Microarray
# ================================
regulon_filt <- regulon %>%
  filter(target %in% rownames(expr_micro)) %>%
  group_by(source) %>%
  filter(n() >= 5) %>%
  ungroup()

# ================================
# 5. Subset Expression Matrix
# ================================
expr_mat <- expr_micro[rownames(expr_micro) %in% regulon_filt$target, ]

# ================================
# 6. Run TF Activity Estimation via Weighted Mean
# ================================
tf_activity <- run_wmean(
  mat = expr_mat,
  net = regulon_filt,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# ================================
# 7. Pivot Activity Output: TF x Sample
# ================================
activity_matrix <- tf_activity %>%
  select(source, condition, score) %>%
  group_by(source, condition) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source")

# ================================
# 8. Merge with Cluster Annotations (Cluster only)
# ================================
annotation_clean <- annotation %>%
  select(Cluster) %>%
  rownames_to_column("Sample")

activity_df <- as.data.frame(t(activity_matrix))  # Samples x TFs
activity_df$Sample <- rownames(activity_df)

activity_df_atlas <- merge(activity_df, annotation_clean, by = "Sample")

# ================================
# 9. Identify Top 10 TFs per Cluster
# ================================
top_tfs_atlas <- activity_df_atlas %>%
  pivot_longer(cols = -c(Sample, Cluster), names_to = "TF", values_to = "Activity") %>%
  group_by(Cluster, TF) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE), .groups = "drop") %>%
  group_by(Cluster) %>%
  slice_max(mean_activity, n = 25)


# ================================
# 10. Output
# ================================
print(top_tfs_atlas)






library(msigdbr)
library(clusterProfiler)
library(dplyr)
library(tibble)

msig_tf_targets <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%
  dplyr::select(gs_name, gene_symbol)


tcga_genes <- rownames(vst_symbol)
# Ensure consistent casing
tcga_genes <- toupper(tcga_genes)

# Filter cluster genes
cluster_genes <- lapply(specific_genes, function(gene_vec) {
  gene_vec <- toupper(gene_vec)
  intersect(gene_vec, tcga_genes)
})

tf_enrichment_results <- lapply(names(cluster_genes), function(cluster_name) {
  gene_list <- cluster_genes[[cluster_name]]
  
  enrichment <- enricher(
    gene = gene_list,
    TERM2GENE = msig_tf_targets
  )
  
  # Add cluster info for labeling
  if (!is.null(enrichment)) {
    enrichment@result$Cluster <- cluster_name
    return(enrichment@result)
  } else {
    return(NULL)
  }
})

# Combine all clusters into one data.frame
tf_enrichment_all <- do.call(rbind, tf_enrichment_results)

library(ggplot2)

top_terms <- tf_enrichment_all %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  slice_max(Count, n = 10)

ggplot(top_terms, aes(x = reorder(Description, Count), y = Count, fill = Cluster)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Cluster, scales = "free_y") +
  coord_flip() +
  labs(
    x = "TF Target Set",
    y = "Overlap Count",
    title = "TF Target Enrichment (MSigDB C3:TFT:GTRD)"
  ) +
  theme_minimal(base_size = 13)







ranked_tcga <- top_tfs %>%
  group_by(Cluster) %>%
  arrange(desc(mean_activity)) %>%
  mutate(Rank_TCGA = row_number()) %>%
  ungroup()

ranked_atlas <- top_tfs_atlas %>%
  group_by(Cluster) %>%
  arrange(desc(mean_activity)) %>%
  mutate(Rank_Atlas = row_number()) %>%
  ungroup()

merged_ranks <- inner_join(
  ranked_tcga %>% select(Cluster, TF, Rank_TCGA),
  ranked_atlas %>% select(Cluster, TF, Rank_Atlas),
  by = c("Cluster", "TF")
)
library(purrr)

cor_by_cluster <- merged_ranks %>%
  group_by(Cluster) %>%
  summarise(
    Spearman_rho = cor(Rank_TCGA, Rank_Atlas, method = "spearman"),
    Shared_TFs = n()
  )

top_tcga <- ranked_tcga %>%
  group_by(Cluster) %>%
  slice_min(Rank_TCGA, n = 10)

top_atlas <- ranked_atlas %>%
  group_by(Cluster) %>%
  slice_min(Rank_Atlas, n = 10)

shared_top <- inner_join(top_tcga, top_atlas, by = c("Cluster", "TF")) %>%
  arrange(Cluster, Rank_TCGA)

# Count shared top TFs per cluster
shared_tf_counts <- shared_top %>%
  group_by(Cluster) %>%
  summarise(Shared_Top10 = n())

# Plot as heatmap
ggplot(shared_tf_counts, aes(x = "Shared_Top10", y = Cluster, fill = Shared_Top10)) +
  geom_tile() +
  geom_text(aes(label = Shared_Top10), color = "white", size = 4) +
  scale_fill_gradient(low = "gray90", high = "darkgreen") +
  labs(title = "Overlap of Top 10 TFs Between TCGA and Atlas", x = "", y = "Cluster", fill = "Shared TFs") +
  theme_minimal(base_size = 13)






library(WGCNA)
options(stringsAsFactors = FALSE)

# Transpose to samples × genes
datExpr <- as.data.frame(t(vst_symbol))

# Keep top 5000 most variable genes
rm(var)
top_var_genes <- names(sort(apply(datExpr, 2, var), decreasing = TRUE))[1:20000]
datExpr_filtered <- datExpr[, top_var_genes]

# Remove genes with low variance
gsg <- goodSamplesGenes(datExpr_filtered, verbose = 3)
datExpr_filtered <- datExpr_filtered[gsg$goodSamples, gsg$goodGenes]

# Match cluster labels
annotation_tcga$Sample<-rownames(annotation_tcga)
datTraits <- annotation_tcga %>%
  filter(Sample %in% rownames(datExpr_filtered)) %>%
  arrange(match(Sample, rownames(datExpr_filtered)))

rownames(datTraits) <- datTraits$Sample
datTraits$Cluster <- as.factor(datTraits$Cluster)

powers <- c(seq(1, 10, 1), seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr_filtered, powerVector = powers, verbose = 5)

# Plot scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     type = "n", xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.85, col = "blue")

softPower <- 5  # (update this based on your plot)
adjacency <- adjacency(datExpr_filtered, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
moduleColors <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = 30)

# Convert colors to labels
moduleLabels <- labels2colors(moduleColors)

MEs <- moduleEigengenes(datExpr_filtered, colors = moduleLabels)$eigengenes
MEs <- orderMEs(MEs)
moduleTraitCor <- cor(MEs, model.matrix(~ Cluster - 1, data = datTraits), use = "p")

# Convert to tidy format
module_trait_df <- as.data.frame(as.table(moduleTraitCor))
colnames(module_trait_df) <- c("Module", "Cluster", "Correlation")

# TF list from DoRothEA or AnimalTFDB
library(dorothea)
regulon <- get_dorothea(organism = "human", levels = c("A", "B", "C")) %>%
  mutate(target = toupper(target))

# Extract unique TFs
human_tfs <- unique(toupper(regulon$source))

# For a module of interest (e.g., "blue"):
target_module <- "blue"
in_module <- moduleLabels == target_module
module_genes <- names(datExpr_filtered)[in_module]

# Compute kME (module membership)
kME <- signedKME(datExpr_filtered, MEs)
hub_genes <- data.frame(Gene = colnames(datExpr_filtered), kME = kME[, paste0("kME", target_module)])
hub_genes <- hub_genes %>% filter(Gene %in% module_genes & Gene %in% human_tfs) %>%
  arrange(desc(kME))


# Make design matrix (one-hot encoding of clusters)
cluster_matrix <- model.matrix(~ 0 + Cluster, data = datTraits)
colnames(cluster_matrix) <- gsub("Cluster", "", colnames(cluster_matrix))

# Correlation between modules and clusters
moduleTraitCor <- cor(MEs, cluster_matrix, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr_filtered))

# Format into tidy dataframe
library(tidyr)
library(dplyr)

mod_trait_df <- as.data.frame(as.table(moduleTraitCor)) %>%
  rename(Module = Var1, Cluster = Var2, Correlation = Freq)

mod_trait_df <- mod_trait_df %>%
  mutate(Module = gsub("^ME", "", Module))

mod_trait_df$Pvalue <- as.vector(moduleTraitP)

# Example: modules positively correlated with 'Mesenchymal' (r > 0.5, p < 0.05)
sig_mods_mesenchymal <- mod_trait_df %>%
  filter(Cluster == "Neuronal", Correlation > 0.5, Pvalue < 0.05) %>%
  pull(Module)

print(sig_mods_mesenchymal)

target_module <- sig_mods_mesenchymal[1]

# Get gene names from that module
mod_genes <- names(datExpr_filtered)[moduleLabels == target_module]

# Filter for TFs only
module_tfs <- intersect(mod_genes, human_tfs)

# Optional: rank by module membership (kME)
kME <- signedKME(datExpr_filtered, MEs)
tf_membership <- data.frame(
  TF = colnames(datExpr_filtered),
  kME = kME[, paste0("kME", target_module)]
) %>%
  filter(TF %in% module_tfs) %>%
  arrange(desc(kME))

# Wrapper function
get_tf_module <- function(cluster_name, min_corr = 0.5, max_p = 0.05) {
  mod_names <- mod_trait_df %>%
    filter(Cluster == cluster_name, Correlation > min_corr, Pvalue < max_p) %>%
    pull(Module)
  
  result <- list()
  
  for (mod in mod_names) {
    genes_in_mod <- names(datExpr_filtered)[moduleLabels == mod]
    tf_in_mod <- intersect(genes_in_mod, human_tfs)
    
    if (length(tf_in_mod) > 0) {
      kme_col <- paste0("kME", mod)
      if (!kme_col %in% colnames(kME)) {
        warning(paste("Module", mod, "not found in kME matrix — skipping"))
        next
      }
      tf_kme <- data.frame(
        TF = tf_in_mod,
        kME = kME[tf_in_mod, kme_col, drop = TRUE]
      )
      result[[mod]] <- tf_kme[order(-tf_kme$kME), ]
    }
  }
  
  return(result)
}


# Run for all clusters
tf_modules_by_cluster <- lapply(unique(datTraits$Cluster), get_tf_module)
names(tf_modules_by_cluster) <- unique(datTraits$Cluster)

# Flatten the nested list into a single data frame
tf_df_all <- purrr::imap_dfr(tf_modules_by_cluster, function(cluster_list, cluster_name) {
  purrr::imap_dfr(cluster_list, function(tf_df, module_name) {
    tf_df %>%
      dplyr::mutate(Module = module_name) %>%
      dplyr::arrange(desc(kME)) %>%
      dplyr::mutate(Cluster = cluster_name,
                    Rank = dplyr::row_number())
  })
})

top_tf_summary <- tf_df_all %>%
  dplyr::filter(Rank <= 30) %>%
  dplyr::select(Cluster, Module, Rank, TF, kME)

top_tf_summary <- top_tf_summary %>%
  arrange(Cluster, Module, desc(kME)) %>%
  mutate(TF_ordered = factor(TF, levels = unique(TF)))  # preserves new order

# Step 2: Plot
ggplot(top_tf_summary, aes(x = Cluster, y = TF_ordered)) +
  geom_point(aes(size = kME, color = Module)) +
  scale_size_continuous(range = c(2, 6)) +
  labs(title = "Top TFs per Cluster",
       y = "TF (grouped by Cluster → Module → kME)",
       x = "Cluster",
       size = "kME") +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(size = 9))






library(WGCNA)
library(dplyr)
library(tidyr)
library(dorothea)
library(decoupleR)
library(ggplot2)
options(stringsAsFactors = FALSE)

# 1. Remove last 3 metadata columns if present
datExpr_atlas <- exp[, 1:(ncol(exp) - 3)]

# 2. Subset top 20,000 variable genes
rm(var)
top_var_genes <- names(sort(apply(datExpr_atlas, 2, var), decreasing = TRUE))
datExpr_filtered_atlas <- datExpr_atlas[, top_var_genes]

# Ensure good samples and genes
gsg_atlas <- goodSamplesGenes(datExpr_filtered_atlas, verbose = 3)
datExpr_filtered_atlas <- datExpr_filtered_atlas[gsg_atlas$goodSamples, gsg_atlas$goodGenes]

# Match sample annotations
datTraits_atlas <- annotation[rownames(datExpr_filtered_atlas), , drop = FALSE]
datTraits_atlas$Cluster <- as.factor(datTraits_atlas$Cluster)

# Pick soft threshold
powers <- c(seq(1, 10, 1), seq(12, 20, 2))
sft_atlas <- pickSoftThreshold(datExpr_filtered_atlas, powerVector = powers, verbose = 5)

# Choose soft power manually based on plot
softPower_atlas <- 4
adjacency_atlas <- adjacency(datExpr_filtered_atlas, power = softPower_atlas)
TOM_atlas <- TOMsimilarity(adjacency_atlas)
dissTOM_atlas <- 1 - TOM_atlas

geneTree_atlas <- hclust(as.dist(dissTOM_atlas), method = "average")
moduleColors_atlas <- cutreeDynamic(dendro = geneTree_atlas, distM = dissTOM_atlas,
                                    deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
moduleLabels_atlas <- labels2colors(moduleColors_atlas)

MEs_atlas <- moduleEigengenes(datExpr_filtered_atlas, colors = moduleLabels_atlas)$eigengenes
MEs_atlas <- orderMEs(MEs_atlas)

# One-hot cluster design matrix
cluster_matrix_atlas <- model.matrix(~ 0 + Cluster, data = datTraits_atlas)
colnames(cluster_matrix_atlas) <- gsub("Cluster", "", colnames(cluster_matrix_atlas))

# Correlate modules with traits
moduleTraitCor_atlas <- cor(MEs_atlas, cluster_matrix_atlas, use = "p")
moduleTraitP_atlas <- corPvalueStudent(moduleTraitCor_atlas, nrow(datExpr_filtered_atlas))

# Tidy summary
mod_trait_df_atlas <- as.data.frame(as.table(moduleTraitCor_atlas)) %>%
  rename(Module = Var1, Cluster = Var2, Correlation = Freq) %>%
  mutate(Module = gsub("^ME", "", Module),
         Pvalue = as.vector(moduleTraitP_atlas))

# Load TFs
regulon_atlas <- get_dorothea(organism = "human", levels = c("A", "B", "C")) %>%
  mutate(target = toupper(target))
human_tfs_atlas <- unique(toupper(regulon_atlas$source))

# Compute kME
kME_atlas <- signedKME(datExpr_filtered_atlas, MEs_atlas)

# Wrapper to extract TFs for top modules per cluster
get_tf_module_atlas <- function(cluster_name, min_corr = 0.5, max_p = 0.05) {
  mod_names <- mod_trait_df_atlas %>%
    filter(Cluster == cluster_name, Correlation > min_corr, Pvalue < max_p) %>%
    pull(Module)
  
  result <- list()
  
  for (mod in mod_names) {
    genes_in_mod <- names(datExpr_filtered_atlas)[moduleLabels_atlas == mod]
    tf_in_mod <- intersect(genes_in_mod, human_tfs_atlas)
    
    if (length(tf_in_mod) > 0) {
      kme_col <- paste0("kME", mod)
      if (!kme_col %in% colnames(kME_atlas)) next
      
      tf_kme <- data.frame(
        TF = tf_in_mod,
        kME = kME_atlas[tf_in_mod, kme_col, drop = TRUE]
      )
      result[[mod]] <- tf_kme[order(-tf_kme$kME), ]
    }
  }
  
  return(result)
}

# Run per cluster
tf_modules_by_cluster_atlas <- lapply(unique(datTraits_atlas$Cluster), get_tf_module_atlas)
names(tf_modules_by_cluster_atlas) <- unique(datTraits_atlas$Cluster)

# Flatten to summary
tf_df_all_atlas <- purrr::imap_dfr(tf_modules_by_cluster_atlas, function(cluster_list, cluster_name) {
  purrr::imap_dfr(cluster_list, function(tf_df, module_name) {
    tf_df %>%
      mutate(Module = module_name) %>%
      arrange(desc(kME)) %>%
      mutate(Cluster = cluster_name,
             Rank = row_number())
  })
})

top_tf_summary_atlas <- tf_df_all_atlas %>%
  filter(Rank <= 30) %>%
  select(Cluster, Module, Rank, TF, kME) %>%
  arrange(Cluster, Module, desc(kME)) %>%
  mutate(TF_ordered = factor(TF, levels = unique(TF)))

# Plot
ggplot(top_tf_summary_atlas, aes(x = Cluster, y = TF_ordered)) +
  geom_point(aes(size = kME, color = Module)) +
  scale_size_continuous(range = c(2, 6)) +
  labs(title = "Top TFs per Cluster (Atlas)",
       y = "TF (Cluster → Module → kME)", x = "Cluster", size = "kME") +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(size = 9))


top_tf_summary$Dataset <- "TCGA"
top_tf_summary_atlas$Dataset <- "Current Atlas"

combined_tf <- bind_rows(top_tf_summary, top_tf_summary_atlas)

# Identify TFs that appear in top 20 for both datasets
common_tfs <- combined_tf %>%
  group_by(Cluster, TF) %>%
  summarise(n = n_distinct(Dataset), .groups = "drop") %>%
  filter(n == 2) %>%
  pull(TF)

# Subset combined data to common TFs
combined_common <- combined_tf %>%
  filter(TF %in% common_tfs) %>%
  arrange(Cluster, Module, desc(kME)) %>%
  mutate(TF_ordered = factor(TF, levels = unique(TF)))


library(ggh4x)

# Step 1: Ensure 'Dataset' is a factor with Atlas first
combined_common <- combined_common %>%
  mutate(Dataset = factor(Dataset, levels = c("Current Atlas", "TCGA")))

# Step 2: Define TF order based on Atlas kME values within Cluster
tf_order <- combined_common %>%
  filter(Dataset == "Current Atlas") %>%
  arrange(Cluster, kME) %>%
  pull(TF) %>%
  unique()

# Step 3: Apply TF order across entire dataframe
combined_common <- combined_common %>%
  mutate(TF_ordered = factor(TF, levels = tf_order))

# Step 4: Plot
hubtfs<-ggplot(combined_common, aes(x = Cluster, y = TF_ordered)) +
  geom_point(aes(size = abs(kME), color = kME)) +
  scale_size_continuous(breaks = c(0.34, 0.62, 0.91)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,limits=c(-max(abs(kME)), max(abs(kME)))) +
  labs(y = "Hub Transcription Factor",
       x = "Cluster",
       size = "Abs. Value (kME)",
       color = "kME") +
  theme_minimal() +
  theme(legend.position = "top",
        text = element_text(color = "black"),         # all general text
        axis.text.x = element_text(color = "black"),  # x-axis title specifically
        axis.text.y = element_text(color = "black"),  # y-axis title specifically
    strip.text = element_text(size = 13),
    panel.spacing = unit(1, "lines"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),                         # Remove y-axis line
    panel.grid.major.y = element_line(color = "grey90"),# Add horizontal grid
    panel.grid.major.x = element_line(color = "grey90"),# Add vertical grid
    panel.grid.minor = element_blank(),                    # Remove minor grids
    panel.border = element_blank()                         # Remove border
  ) +
  facet_wrap(~Dataset, scales = "free_y", ncol = 1, strip.position = "top") +
  #scale_x_discrete(guide = guide_axis(angle = 90))+
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_discrete(labels = function(x) x),  # Show labels for Atlas
      scale_y_discrete(labels = NULL)            # Hide labels for TCGA
    )
  )+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

