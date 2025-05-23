# ===================== Setup =====================
setwd("~/Library/CloudStorage/OneDrive-uoflhealth/Atlas/PHPG")

# Load required libraries
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(SummarizedExperiment)
library(ggplot2)

# ===================== Load Somatic Mutation Data =====================
maf_object <- read.maf("TCGA_PPGL_mutations.maf")

final_df <- maf_object@data %>%
  select(Sample = Tumor_Sample_Barcode, Gene = Hugo_Symbol) %>%
  bind_rows(data.frame(
    Sample = setdiff(
      unique(maf_object@clinical.data$Tumor_Sample_Barcode),
      unique(maf_object@data$Tumor_Sample_Barcode)
    ),
    Gene = "WT"
  )) %>%
  arrange(Sample) %>%
  mutate(Location = "Somatic")

# ===================== Load Expression and Germline Data =====================
query <- GDCquery(
  project = "TCGA-PCPG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")
GDCdownload(query)
exp_data <- GDCprepare(query)

germ2 <- read.csv("~/Downloads/mmc2.csv", header = TRUE, stringsAsFactors = TRUE) %>%
  select(-c(2:126, 5:8, 10, 12:48)) %>%
  mutate(
    gender = toupper(as.character(gender)),
    race = as.character(race),
    birth_days_to = as.integer(birth_days_to)
  )

clinical_sample_map <- as.data.frame(colData(exp_data)) %>%
  transmute(
    sample,
    gender = toupper(as.character(gender)),
    birth_days_to = as.integer(days_to_birth)
  )

germline <- left_join(germ2, clinical_sample_map, by = c("birth_days_to", "gender")) %>%
  select(sample, HUGO_Symbol) %>%
  left_join(final_df %>%
              transmute(sample_short = substr(Sample, 1, 16), Sample),
            by = c("sample" = "sample_short")) %>%
  distinct() %>%
  transmute(Sample, Gene = as.character(HUGO_Symbol), Location = "Germline")

# ===================== Combine Somatic and Germline Mutations =====================
muts_annotated <- bind_rows(final_df, germline) %>%
  distinct(Sample, Gene, .keep_all = TRUE) %>%
  mutate(sample_short = substr(Sample, 1, 16)) %>%
  left_join(annotation_tcga %>%
              rownames_to_column("Sample") %>%
              mutate(sample_short = substr(Sample, 1, 16)) %>%
              select(sample_short, Cluster),
            by = "sample_short")

# ===================== Manual Sample Annotations =====================
manual_samples <- tibble(
  sample_short = c(
    "TCGA-P8-A5KC-01A", "TCGA-PR-A5PF-01A", "TCGA-PR-A5PH-01A", "TCGA-QR-A70H-01A",
    "TCGA-RW-A686-01A", "TCGA-RW-A686-06A", "TCGA-RW-A8AZ-01A", "TCGA-WB-A80K-01A",
    "TCGA-WB-A80Q-01A", "TCGA-WB-A80Y-01A",  # MAML3 (Cluster 6)
    "TCGA-QR-A707-01A", "TCGA-QT-A5XN-01A", "TCGA-RT-A6YC-01A"  # Kinase (Cluster 1)
  ),
  Gene = c(rep("MAML3", 10), "NF1", "NGFR", "BRAF"),
  Theme = c(rep("Neurodevelopmental_Signaling", 10), rep("Kinase_Signaling", 3)),
  Cluster = c(rep(6, 10), rep(1, 3)),
  Present = 1
)

extra_df <- manual_samples %>%
  left_join(annotation_tcga %>%
              rownames_to_column("Sample") %>%
              mutate(sample_short = substr(Sample, 1, 16)) %>%
              select(Sample, sample_short),
            by = "sample_short") %>%
  mutate(
    Sample = coalesce(Sample, paste0(sample_short, "-FIXED")),
    Location = NA_character_,
    Cluster = factor(Cluster, levels = levels(muts_annotated$Cluster))
  ) %>%
  select(Sample, Gene, Location, sample_short, Cluster, Theme, Present)

# ===================== Append Manual Entries BEFORE Reassignment =====================
muts_annotated <- bind_rows(
  muts_annotated,
  extra_df %>% select(Sample, Gene, Location, sample_short, Cluster)
)

# ===================== Reassign Clusters Based on Gene =====================
muts_annotated <- muts_annotated %>%
  mutate(Cluster = case_when(
    Gene %in% c("NF1", "HRAS", "RET") ~ 1L,
    Gene == "EPAS1" ~ 2L,
    Gene == "MAML3" ~ 6L,
    Gene == "SDHD"  ~ 3L,
    TRUE ~ as.integer(as.character(Cluster))
  )) %>%
  mutate(Cluster = factor(Cluster))  # Ensure final format

# ===================== Prioritized Cluster Assignment =====================
# Step 1: Prioritize marker genes for determining cluster
marker_cluster_priority <- muts_annotated %>%
  filter(Gene %in% c("NF1", "HRAS", "RET", "EPAS1", "MAML3", "SDHD")) %>%
  select(sample_short, Cluster) %>%
  distinct()

# Step 2: Fallback to most common cluster otherwise
fallback_clusters <- muts_annotated %>%
  select(sample_short, Cluster) %>%
  distinct() %>%
  group_by(sample_short) %>%
  summarise(
    Updated_Cluster = as.integer(names(sort(table(Cluster), decreasing = TRUE)[1])),
    .groups = "drop"
  )

# Step 3: Use marker-based cluster if available, otherwise fallback
cluster_map_updated <- fallback_clusters %>%
  left_join(marker_cluster_priority, by = "sample_short") %>%
  mutate(
    Final_Cluster = coalesce(as.integer(as.character(Cluster)), Updated_Cluster)
  ) %>%
  select(sample_short, Updated_Cluster = Final_Cluster)

# ===================== Update annotation_tcga =====================
annotation_tcga <- annotation_tcga %>%
  rownames_to_column("Sample") %>%
  mutate(
    sample_short = substr(Sample, 1, 16),
    Original_Cluster = as.integer(as.character(Cluster))
  ) %>%
  left_join(cluster_map_updated, by = "sample_short") %>%
  mutate(
    Cluster = coalesce(Updated_Cluster, Original_Cluster)
  ) %>%
  select(Sample, Cluster) %>%
  column_to_rownames("Sample")



# ===================== Theme Assignment =====================
theme_list <- list(
  "DNA_Repair" = c( "RAD50", "RECQL", "PALB2", "PMS2", "BLM", "MLH3", "MDC1", "RAD54B"),
  "Hypoxia_Endothelial" = c("VHL", "EPAS1", "EGLN1"),
  "Kinase_Signaling" = c("NF1", "RET", "HRAS", "BRAF", "FGFR1","RYR1", "NGFR"),
  "Metabolic_TCA" = c("ATRX","BRCA1", "BRCA2", "SDHB", "SDHD", "IDH1", "DHTKD1", "TPI1", "FBP2","TGDS", "TRHDE", "FKBP9"),
  "Neurodevelopmental_Signaling" = c("MAML3","MAX","DLG5", "SOX5", "ZEB2", "GABRA2", "NPY1R", "TENM3", "NRXN3", "OLIG3", "BARHL2", "ZFHX3"),
  "Unassigned" = c("ABCA13", "CSDE1", "TTN", "ABCA12", "HUWE1", "KCNH5", "SETD2", "UNC79", "MUC16", "ALK", "TFAP2D", "IRX6")
)

theme_df <- enframe(theme_list, name = "Theme", value = "Gene") %>% unnest(Gene)

muts_annotated_with_theme <- muts_annotated %>%
  left_join(theme_df, by = "Gene") %>%
  mutate(Theme = ifelse(is.na(Theme), NA, Theme)) %>%
  bind_rows(extra_df) # ensure manual entries with theme + cluster are included

# ===================== Add Manual Sample-Gene-Theme Annotations =====================
manual_samples <- tibble(
  sample_short = c(
    "TCGA-P8-A5KC-01A", "TCGA-PR-A5PF-01A", "TCGA-PR-A5PH-01A", "TCGA-QR-A70H-01A",
    "TCGA-RW-A686-01A", "TCGA-RW-A686-06A", "TCGA-RW-A8AZ-01A", "TCGA-WB-A80K-01A",
    "TCGA-WB-A80Q-01A", "TCGA-WB-A80Y-01A",  # 10 MAML3
    "TCGA-QR-A707-01A", "TCGA-QT-A5XN-01A", "TCGA-RT-A6YC-01A"  # 3 Kinase samples
  ),
  Gene = c(
    rep("MAML3", 10),
    "NF1", "NGFR", "BRAF"
  ),
  Theme = c(
    rep("Neurodevelopmental_Signaling", 10),
    rep("Kinase_Signaling", 3)
  )
)

extra_df <- manual_samples %>%
  left_join(annotation_tcga %>%
              rownames_to_column("Sample") %>%
              mutate(sample_short = substr(Sample, 1, 16),
                     Cluster = as.character(Cluster)),
            by = "sample_short") %>%
  mutate(Present = 1)

muts_annotated_with_theme <- bind_rows(muts_annotated_with_theme, extra_df)

genes<-c("NF1","HRAS","RET","FGFR1","BRAF","NGRF",
         "VHL","EPAS1",
         "SDHB","SDHD","ATRX","BRCA1","BRCA2","TGDS",
         "MAML3",
         "TTN","CSDE1")
# Step 1: Subset only the relevant genes
tcga_muts <- subset(muts_annotated_with_theme, Gene %in% genes)

# Step 2: Keep only unique sample-gene combinations
tcga_muts <- distinct(tcga_muts[, c("Sample", "Gene")])

# Step 3: Merge gene labels
tcga_muts <- tcga_muts %>%
  mutate(Gene = case_when(
    Gene %in% c("BRCA1", "BRCA2") ~ "BRCA1/2",
    Gene %in% c("SDHB", "SDHD") ~ "SDHx",
    TRUE ~ Gene
  ))

# Step 4: Define structural vs non-structural genes
structural_genes <- c("ATRX", "TTN", "BRCA1/2")

tcga_muts <- tcga_muts %>%
  mutate(
    Gene_non_structural = if_else(!(Gene %in% structural_genes), Gene, NA_character_),
    Gene_structural = if_else(Gene %in% structural_genes, Gene, NA_character_)
  )

# Step 5: Collapse per sample
tcga_muts_collapsed <- tcga_muts %>%
  group_by(Sample) %>%
  summarise(
    Gene_non_structural = paste(na.omit(unique(Gene_non_structural)), collapse = ", "),
    Gene_structural = paste(na.omit(unique(Gene_structural)), collapse = ", "),
    .groups = "drop"
  )

colnames(tcga_muts_collapsed)[2:3]<-c("Primary_Mutations","Other")
rownames(tcga_muts_collapsed)<-tcga_muts_collapsed$Sample


# ===================== Enrichment: Theme × Cluster =====================
theme_matrix <- muts_annotated_with_theme %>%
  filter(!is.na(Theme)) %>%
  distinct(sample_short, Cluster, Theme) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Theme, values_from = value, values_fill = 0)

get_fisher_results_theme <- function(theme) {
  counts <- theme_matrix %>%
    select(Cluster, all_of(theme)) %>%
    group_by(Cluster) %>%
    summarise(Positive = sum(.data[[theme]] == 1),
              Negative = sum(.data[[theme]] == 0), .groups = "drop")
  map_dfr(1:6, function(k) {
    in_k <- counts %>% filter(Cluster == as.character(k))
    not_k <- counts %>% filter(Cluster != as.character(k)) %>%
      summarise(across(Positive:Negative, sum))
    mat <- matrix(c(in_k$Positive, in_k$Negative, not_k$Positive, not_k$Negative), nrow = 2, byrow = TRUE)
    tibble(Theme = theme, Cluster = k, P = fisher.test(mat)$p.value)
  })
}

all_theme_results <- map_dfr(setdiff(names(theme_matrix), c("sample_short", "Cluster")),
                             get_fisher_results_theme) %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

# ======= PLOT: Enriched Themes by Cluster =======
#top_themes_by_cluster <- all_theme_results %>%
#  group_by(Cluster) %>%
#  slice_min(FDR, n = 3)

#ggplot(top_themes_by_cluster %>%
#         mutate(logFDR = -log10(FDR),
#                Cluster = as.factor(Cluster)),
#       aes(x = Cluster, y = Theme)) +
#  geom_point(aes(size = logFDR, color = FDR)) +
#  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
#  scale_size_continuous(name = "-log10(FDR)") +
#  labs(
#   title = "Enriched Themes by Cluster",
#    x = "Cluster",
#    y = "Theme",
#    color = "FDR"
#  ) +
#  theme_minimal(base_size = 14)

# ===================== Reconstruct Plot Matrix =====================
annotation_tcga_fixed <- annotation_tcga %>%
  rownames_to_column("Sample") %>%
  mutate(sample_short = substr(Sample, 1, 16),
         Cluster = as.integer(Cluster)) %>%
  arrange(Cluster)

full_grid2 <- expand.grid(
  sample_short = annotation_tcga_fixed$sample_short,
  Gene = unique(theme_df$Gene)
) %>%
  left_join(theme_df, by = "Gene") %>%
  left_join(
    muts_annotated_with_theme %>%
      filter(Gene %in% theme_df$Gene) %>%
      mutate(Present = 1) %>%
      select(sample_short, Gene, Present),
    by = c("sample_short", "Gene")
  ) %>%
  mutate(Present = ifelse(is.na(Present), 0, Present)) %>%
  left_join(annotation_tcga_fixed %>% select(sample_short, Cluster), by = "sample_short") %>%
  distinct(sample_short, Gene, .keep_all = TRUE)

# ===================== Sample Ordering by Mutation Score =====================
mutation_matrix <- full_grid2 %>%
  select(sample_short, Cluster, Gene, Present) %>%
  pivot_wider(names_from = Gene, values_from = Present, values_fill = 0)

gene_cols <- setdiff(colnames(mutation_matrix), c("sample_short", "Cluster"))
gene_freqs <- colSums(mutation_matrix[, gene_cols])
genes_sorted <- names(sort(gene_freqs, decreasing = TRUE))

weights <- 2 ^ rev(seq_along(genes_sorted))
mutation_matrix$mutation_score <- as.numeric(as.matrix(mutation_matrix[, genes_sorted]) %*% weights)

sample_levels <- mutation_matrix %>%
  arrange(Cluster, desc(mutation_score)) %>%
  pull(sample_short)

# ===================== Gene Ordering within Theme =====================
gene_mutation_counts <- full_grid2 %>%
  group_by(Gene) %>%
  summarise(n_mutated = sum(Present), .groups = "drop")

gene_levels <- theme_df %>%
  left_join(gene_mutation_counts, by = "Gene") %>%
  arrange(Theme, n_mutated, Gene) %>%
  pull(Gene)

# ===================== Theme and Cluster Boundaries =====================
theme_boundaries <- tibble(Gene = gene_levels) %>%
  left_join(theme_df, by = "Gene") %>%
  mutate(row = row_number()) %>%
  group_by(Theme) %>%
  summarise(y_line = max(row) + 0.5, .groups = "drop")

cluster_boundaries <- annotation_tcga_fixed %>%
  mutate(col = match(sample_short, sample_levels)) %>%
  group_by(Cluster) %>%
  summarise(x_line = max(col) + 0.5, .groups = "drop")

# ======= PLOT: Final Heatmap =======
ggplot(full_grid2, aes(
  x = factor(sample_short, levels = sample_levels),
  y = factor(Gene, levels = gene_levels),
  fill = factor(Present)
)) +
  geom_tile(color = "white") +
  geom_hline(data = theme_boundaries, aes(yintercept = y_line), color = "white", size = 1) +
  geom_vline(data = cluster_boundaries, aes(xintercept = x_line), color = "white", size = 1) +
  scale_fill_manual(values = c("0" = "lightgrey", "1" = "steelblue"), name = "Mutation") +
  labs(
    title = "Theme-Associated Gene Mutations by Sample",
    x = "Sample (by Cluster, sorted by mutation score)",
    y = "Gene (by Theme, increasing mutation frequency)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )


# finalize - atlas
#- cluster atlas PG/PH
#- then heatmap tcga
#- then find genetic mutations in tcga 
#- re-arrange tcga heatmap based on samples arranged by common mutations / signature 
#(remember 4 samples don't have genetic information) "TCGA-P8-A5KD-11A" "TCGA-SQ-A6I4-11A" "TCGA-P8-A5KC-11A" "TCGA-S7-A7WU-01A"
# find DEG based on re-arranged TCGA -- that are commmon with atlas 
# finalize heatmap for both atlas and TCGA -- picking labeled genes and setting top annotation bars with mutation mutations for both heatmaps
# do GO / KEGG / Reactome based on commmon final DEGs to finalize naming of the clusters together with mutational sig
# make TCGA t-SNE based on less number of genes that match clustering -- consider 50% most MAD / variable

clin<-read.csv("~/Downloads/pcpg_clin.csv",header = T,stringsAsFactors = T)
# First, ensure Sample.ID and sample_short are character for matching
clin$Sample.ID <- as.character(clin$Sample.ID)

# Step 3: Create and apply annotation
annotation_tcga <- data.frame(Cluster = factor(sample_clusters))
rownames(annotation_tcga) <- names(sample_clusters)

# Step 1: Order samples by cluster assignment
ordered_samples <- rownames(annotation_tcga)[order(annotation_tcga$Cluster)]

# Step 2: Reorder expression matrix and annotation
expr_ordered <- expr_scaled[, ordered_samples]
annotation_tcga <- annotation_tcga[ordered_samples, , drop = FALSE]
# Format Cluster factor and assign Set2 colors
annotation_tcga$Cluster <- factor(annotation_tcga$Cluster, 
                                  labels = c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental"))

cluster_levels <- levels(annotation_tcga$Cluster)
cluster_colors <- setNames(brewer.pal(length(cluster_levels), "Set2"), cluster_levels)

# Prepare cluster mapping from annotation_tcga
cluster_map <- annotation_tcga %>%
  rownames_to_column("Sample") %>%
  mutate(sample_short = substr(Sample, 1, 16)) %>%
  select(sample_short, Cluster)

# Join the cluster column into the clinical data
clin <- clin %>%
  left_join(cluster_map, by = c("Sample.ID" = "sample_short"))

# Required packages
library(dplyr)
library(stats)

# Function for Kruskal-Wallis test (for continuous variable)
test_kruskal <- function(var) {
  kruskal.test(as.formula(paste(var, "~ Cluster")), data = clin)
}

# Function for Chi-square test (with fallback to Fisher's if needed)
test_chisq <- function(var) {
  tbl <- table(clin[[var]], clin$Cluster)
  test_chi <- suppressWarnings(chisq.test(tbl))
  
  # Check if any expected values are too low
  if (any(test_chi$expected < 5)) {
    test_fisher <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)
    return(list(p.value = test_fisher$p.value, method = "Fisher's exact test (simulated)"))
  } else {
    return(list(p.value = test_chi$p.value, method = "Chi-squared test"))
  }
}

# Collect results
results <- list()

# Categorical variables
cat_vars <- c("Gender", "Race", "Vital.Status",
              "Tumor", "Aggressive.or.Metastatic")

for (var in cat_vars) {
  tmp <- test_chisq(var)
  results[[var]] <- list(p.value = tmp$p.value, method = tmp$method)
}

# Continuous variable (non-parametric)
kw_test <- test_kruskal("Age.at.Initial.Pathologic.Diagnosis")
results[["Age.at.Initial.Pathologic.Diagnosis"]] <- list(p.value = kw_test$p.value, method = "Kruskal-Wallis test")

# Format results into a dataframe
tibble::tibble(
  Variable = names(results),
  Test = sapply(results, function(x) x$method),
  P.value = sapply(results, function(x) x$p.value)
) %>%
  arrange(P.value)

# Helper: proportion table for a categorical variable
cluster_proportions <- function(var) {
  prop_table <- table(clin[[var]], clin$Cluster)
  round(100 * prop.table(prop_table, margin = 2), 1)  # % within cluster
}

# 1. Tumor Type (Pheo vs PG)
cat("\n=== Pheochromocytoma or Paraganglioma ===\n")
print(cluster_proportions("Tumor"))

# 3. Age
cat("\n=== Age at Initial Diagnosis (median per cluster) ===\n")
print(clin %>%
        group_by(Cluster) %>%
        summarise(Median_Age = median(Age.at.Initial.Pathologic.Diagnosis, na.rm = TRUE)))

# 4. Clinically Aggressive or Metastatic
cat("\n=== Clinically Aggressive or Metastatic ===\n")
print(cluster_proportions("Aggressive.or.Metastatic"))

cat("\n=== Vital Status ===\n")
print(cluster_proportions("Vital.Status"))


library(dplyr)
library(tidyr)
library(ggpubr)

# Function to make a stacked bar plot ordered by % of "Yes" or category 2
make_stacked_bar <- function(data, var, yes_label = NULL) {
  # Ensure Cluster is a factor with full levels
  all_clusters <- sort(unique(na.omit(data$Cluster)))
  data <- data %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(Cluster = factor(Cluster, levels = all_clusters))
  
  # Count occurrences
  df <- data %>%
    count(Cluster, !!sym(var)) %>%
    complete(Cluster, !!sym(var), fill = list(n = 0)) %>%  # fill missing combos
    group_by(Cluster) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()
  
  # Reference level (to sort clusters)
  ref_level <- if (!is.null(yes_label)) {
    yes_label
  } else {
    levels(data[[var]])[2]
  }
  
  # Cluster order based on proportion of 'yes' category
  cluster_order <- df %>%
    filter(.data[[var]] == ref_level) %>%
    arrange(desc(pct)) %>%
    pull(Cluster)
  
  # Final plot
  ggplot(df, aes(x = factor(Cluster, levels = cluster_order), y = pct, fill = .data[[var]])) +
    geom_bar(stat = "identity", position = "fill", color = "black") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Cluster", y = "Proportion", fill = var,
         title = paste("Proportion of", var, "by Cluster")) +
    theme_pubr()
}

# Boxplot ordered by decreasing median
make_cluster_boxplot <- function(data, var) {
  medians <- data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(Cluster) %>%
    summarise(median_val = median(.data[[var]], na.rm = TRUE)) %>%
    arrange(desc(median_val))
  
  cluster_order <- medians$Cluster
  
  ggboxplot(data,x = "Cluster",y = var,order = cluster_order,
            color = "Cluster",add = "jitter",palette = "Dark2") +
    labs(x = "Cluster",y = var,
         title = paste(var, "by Cluster (ordered by median)")) +
    theme_pubr()
}


# Ensure Cluster is factor
clin$Cluster <- as.factor(clin$Cluster)

# Plot 1: Tumor type (PG/Pheo)
plot1 <- make_stacked_bar(clin, "Tumor", yes_label = "Paraganglioma")+scale_y_reverse()+
  scale_fill_manual(values = c("Paraganglioma" = "darkmagenta", "Pheochromocytoma" = "lightgreen"),
                    labels = c("Paraganglioma"="PG", "Pheochromocytoma"="PC"))+labs(title = NULL)+
  rremove("y.axis")+rremove("y.ticks")+rremove("y.text")+
  annotate("text", x = 5, y = -0.05, label = expression(chi^2 ~ "p < 0.0001"))+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "darkgrey")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot2 <- make_stacked_bar(clin, "Aggressive.or.Metastatic", yes_label = "Yes")+
  scale_fill_manual(values = c("Yes" = "darkred", "No" = "skyblue"))+
  labs(title = NULL, fill="Aggr./Met?")+rremove("y.axis")+rremove("y.ticks")+
  rremove("y.text")+
  annotate("text", x = 5, y = 1.05, label = expression(chi^2 ~ "p < 0.0007"))+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "darkgrey")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot3 <- make_cluster_boxplot(clin, "Age.at.Initial.Pathologic.Diagnosis")+
  labs(y="Age (Years)", title=NULL)+theme(legend.position = "none")+stat_compare_means(label.x=2,label.y=87)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# OR: arrange in grid (optional)
library(patchwork)
plot1+plot2+plot3
