# Load libraries
library(limma)
library(org.Hs.eg.db)
library(pheatmap)
library(FNN)
library(dbscan)
library(ggpubr)
library(patchwork)
library(grid)
library(gridExtra)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(gridGraphics)

# Load data
expdata <- readRDS("~/Library/CloudStorage/OneDrive-Personal/Research/Tumor - Basic/ExpData_BC.RDS")
d <- read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv", 
              header = TRUE, row.names = 1, stringsAsFactors = TRUE)
z <- d[d$Final.diagnosis_bc == "PH/PG", ]
p <- expdata[, rownames(z)]
exp <- t(p)
exp <- as.data.frame(exp)

# Map ENTREZ IDs to gene symbols
entrez_ids <- colnames(exp)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Clean mapping
valid <- !is.na(gene_symbols)
exp <- exp[, valid]
gene_symbols <- gene_symbols[valid]
colnames(exp) <- gene_symbols

# Add tSNE (assumed in z) and run DBSCAN
exp$t1 <- z$tsne1_bc
exp$t2 <- z$tsne2_bc

res <- dbscan(exp[, c("t1", "t2")], eps = 1, minPts = 7)
exp$class <- as.factor(res$cluster)

# Reassign noise points to nearest cluster
noise_idx <- which(res$cluster == 0)
if (length(noise_idx) > 0) {
  cluster_points <- exp[res$cluster != 0, c("t1", "t2")]
  cluster_labels <- res$cluster[res$cluster != 0]
  nearest <- get.knnx(cluster_points, exp[noise_idx, c("t1", "t2")], k = 1)
  reassigned_clusters <- cluster_labels[nearest$nn.index]
  res$cluster[noise_idx] <- reassigned_clusters
  exp$class <- as.factor(res$cluster)
}

#convert to class 5
exp$class[c(which(rownames(exp)=="AGZ_025_U133_2.CEL"),
which(rownames(exp)=="AGZ_188_U133_2.CEL"),
which(rownames(exp)=="AGZ_043_U133_2.CEL"),
which(rownames(exp)=="AGZ_170_U133_2.CEL"),
which(rownames(exp)=="GSM1638009_Path_FVN_030707_F121.CEL.gz"))]<-"5"

#convert to class 4
exp$class[c(which(rownames(exp)=="AGZ_090_U133_2.CEL"),
which(rownames(exp)=="AGZ_021_U133_2.CEL"),
which(rownames(exp)=="AGZ_121_U133_2.CEL"),
which(rownames(exp)=="AGZ_185_U133_2.CEL"),
which(rownames(exp)=="AGZ_023_U133_2.CEL"),
which(rownames(exp)=="AGZ_100_U133_2.CEL"))]<-"4"

#convert to class 3
exp$class[c(which(rownames(exp)=="GSM1638035_Path_FVN_071106_F170chip8.CEL.gz"),
            which(rownames(exp)=="AGZ_050_U133_2.CEL"),
            which(rownames(exp)=="AGZ_106b_U133_2.CEL"))]<-"3"

#convert to class 6
exp$class[c(which(rownames(exp)=="AGZ_073b_U133_2.CEL"),
            which(rownames(exp)=="GSM1638040_Path_FVN_120706_F303.CEL.gz"))]<-"6"

# ----- Differential Expression & Heatmap -----

# Prepare expression matrix with gene symbols as rownames
exprs_raw <- exp[, !(colnames(exp) %in% c("t1", "t2", "class"))]
exprs_mat <- t(as.matrix(exprs_raw))
rownames(exprs_mat) <- colnames(exprs_raw)

# Design matrix
group <- factor(exp$class)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(exprs_mat, design)

# All pairwise contrasts
contrast_names <- combn(levels(group), 2, FUN = function(x) paste(x[1], "vs", x[2], sep = ""), simplify = TRUE)
contrast_matrix <- sapply(contrast_names, function(comp) {
  g <- unlist(strsplit(comp, "vs"))
  v <- rep(0, length(levels(group)))
  names(v) <- levels(group)
  v[g[1]] <- 1
  v[g[2]] <- -1
  v
})
contrast_matrix <- as.matrix(contrast_matrix)
colnames(contrast_matrix) <- contrast_names

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Function to find genes upregulated in one cluster vs all others
get_up_genes_one_vs_rest <- function(cluster_label, exprs_mat, group, lfc_thresh = 1, fdr_thresh = 0.05) {
  # Create binary label: target cluster vs rest
  group_bin <- ifelse(group == cluster_label, paste0("C", cluster_label), "Other")
  group_bin <- factor(group_bin, levels = c("Other", paste0("C", cluster_label)))
  
  # Design matrix
  design <- model.matrix(~ 0 + group_bin)
  colnames(design) <- levels(group_bin)
  
  # Fit linear model
  fit <- lmFit(exprs_mat, design)
  contrast <- makeContrasts(contrasts = paste0("C", cluster_label, "-Other"), levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract upregulated genes
  tt <- topTable(fit2, coef = 1, number = Inf, adjust = "BH")
  rownames(tt)[tt$adj.P.Val < fdr_thresh & tt$logFC > lfc_thresh]
}

# For each cluster, get the top 10 upregulated genes
specific_genes <- lapply(levels(group), function(g) {
  top <- get_up_genes_one_vs_rest(g, exprs_mat, group)
  head(top, 30) #change this for G0 and KEGG
})
names(specific_genes) <- levels(group)

specific_genes[["6"]]<-c(specific_genes[["6"]], specific_genes[["4"]][c(21,27)])
specific_genes[["4"]]<-specific_genes[["4"]][-c(21,27)]

specific_genes[["4"]]<-(specific_genes[["4"]][-5]) #"ICA1-AS1" not in TCGA
specific_genes[["5"]]<-(specific_genes[["5"]][-4]) #"CIMIP2B" not in TCGA


#Neuronal
c1<-c("SHANK2", "NRGN", "FEV", "PNMT", "PPP1R1B", "NCAM2", "ATP8A2", "SALL4","RPH3A", "GABRG2")
c2<-c("TFAP2C", "EGLN3", "HK2", "NDRG1", "EGFL7","AQP1","DDIT4L","DEPP1","MIR210HG","PTPRB")
#Mesenchyml
#c2<-c("TFAP2C", "EGLN3", "HK2", "WNT3", "HEY2", "GLI3", "NDRG1", "ETS1", "RHOJ")
#Mesenchym/hypoxic
#c2<-c("EGLN3", "HK2", "NDRG1", "TFAP2C", "WNT3", "DDIT4L", "DEPP1", "GLI3", "RHOJ")
#Metabolic
c3<-c("SLITRK1", "NXPH4", "CNIH3", "SLC6A15", "STXBP6", "SLC7A11", "MTHFD2", "SHMT2", "ARG2", "CPS1")
#Indeterminate
c4<-c()
#Steroidogenic
c5<-c("CYP11A1", "FDXR", "FDX1", "ALAS1", "SULT2A1", "SOAT1", "MC2R", "NR0B1", "NR1H4", "SLC51A")
#Developmental
c6<-c("IRX4", "GLI2", "NTRK3", "MDGA1", "NOG", "WNT4", "GHR", "SST", "EPHB1","AR","MAMLD1")

genes_to_show<-c(c1,c2,c3,c4,c5,c6)

# Combine all unique genes
top_genes_all <- unique(unlist(specific_genes))

# Generate heatmap data
heat_data <- exp[, top_genes_all]
heat_data <- heat_data[order(group), ]               # sort samples by cluster
group_sorted <- group[order(group)]

# Scale expression by gene
scaled <- scale(as.matrix(heat_data))  # samples x genes, scale by gene (column)

# Transpose to get genes x samples
scaled <- t(scaled)

# Clip to improve color contrast
scaled[scaled > 2] <- 2
scaled[scaled < -2] <- -2

# Ensure sample annotation matches column order
group_sorted <- group[order(group)]
annotation <- data.frame(Cluster = group_sorted)
rownames(annotation) <- rownames(heat_data)  # original sample names

# Set column names of scaled matrix to sample names
colnames(scaled) <- rownames(annotation)

# Define row gaps (after these row indices)
row_gaps_indices <- c(30, 60, 90, 118, 148)
row_gaps <- rep(unit(2, "mm"), length(row_gaps_indices))

# Determine which rows to highlight
highlight_genes <- intersect(genes_to_show, rownames(scaled))
highlight_indices <- which(rownames(scaled) %in% highlight_genes)

# Create the right annotation with gene labels
row_anno <- rowAnnotation(
  Highlight = anno_mark(
    at = highlight_indices,
    labels = highlight_genes
  )
)

# Define color scale for heatmap expression values
col_fun <- colorRamp2(c(min(scaled), 0, max(scaled)),
                      c("navy", "white", "firebrick3"))

# Define row splits based on gene groups
n <- nrow(scaled)
split_labels <- character(n)
split_labels[1:30] <- "Neuronal"
split_labels[31:60] <- "Vascular"
split_labels[61:90] <- "Metabolic"
split_labels[91:115] <- "Indeterminate"
split_labels[116:144] <- "Steroidal"
split_labels[145:n] <- "Developmental"
split_labels <- factor(split_labels, levels = c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental"))

# Ensure annotation$Cluster is properly formatted
annotation$Cluster <- factor(annotation$Cluster, 
                             labels = c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental"))

# Define colors for top annotation using Set2
cluster_levels <- levels(annotation$Cluster)
cluster_colors <- setNames(brewer.pal(length(cluster_levels), "Set2"), cluster_levels)

annotation$Mutation <- z[rownames(annotation), 31]
annotation$Mutation <- factor(annotation$Mutation, levels = 
                                c("NF1", "RET", "pseudo RET", "VHL", "EPAS1","pseudo VHL", "SDHx",
                                  "pseudo SDHx", "TMEM127"))

mutation_colors <- c("NF1" = "darkgreen","RET" = "green","pseudo RET" = "olivedrab1",
  "VHL" = "orangered","EPAS1" = "orange","pseudo VHL" = "lightsalmon3",
  "SDHx" = "blue","pseudo SDHx" = "slateblue4")

# Top annotation with custom colors
top_annot <- HeatmapAnnotation(
  Cluster = annotation$Cluster,
  Mutation = annotation$Mutation,
  col = list(Cluster = cluster_colors,
             Mutation = mutation_colors),
  show_annotation_name = F
)

# Column split based on sample clusters
column_split <- annotation[colnames(scaled), "Cluster"]

# Generate the ComplexHeatmap
heatmap_obj<-Heatmap(scaled,
        name = "Expression",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = F,
        top_annotation = top_annot,
        #right_annotation = row_anno,
        row_split = split_labels,
        row_title_gp = gpar(fontsize = 10),
        row_gap = unit(1, "mm"),
        column_split = column_split,
        column_gap = unit(1, "mm"),
        show_heatmap_legend = F,
        column_title = "Current Atlas - PH/PG Samples (Microarray; N=240)")

heatmap_obj<-draw(heatmap_obj, show_annotation_legend = F)
heatmap_obj

# Plot tSNE
tsne_plot <- ggplot(exp, aes(t1, t2, color = class)) +
  geom_point() + xlim(70, 85) + scale_color_brewer(palette = "Set2") + ylim(min((exp$t2)[-232]), max((exp$t2)[-232]))+
  theme_pubr() + theme(legend.position = "none", plot.title = element_text(size = 13)) + 
  xlab("t-SNE1") + ylab("t-SNE2") + labs(title = "Current Atlas\nPH/PG Samples\nN=240")

ggplotly(tsne_plot)

# Combine plots
#grid.arrange(ggplotGrob(tsne_plot), heatmap_obj$gtable, ncol = 2)

################################PATHWAY ANALYSIS###################################

#change specific genes to 500 for GO/KEGG/Reactome analysis below
specific_genes <- lapply(levels(group), function(g) {
  top <- get_up_genes_one_vs_rest(g, exprs_mat, group)
  head(top, 500) #change this for G0 and KEGG
})
names(specific_genes) <- levels(group)

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ReactomePA)

# Step 1: Convert gene symbols to Entrez IDs for each cluster
gene_lists_entrez <- lapply(specific_genes, function(glist) {
  bitr(glist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
})



# Step 2: Run KEGG enrichment for each cluster
kegg_results <- lapply(gene_lists_entrez, function(entrez_ids) {
  enrichKEGG(gene = entrez_ids,
             organism = "hsa",
             pvalueCutoff = 0.05)
})

# Step 3: Extract top 3 enriched pathways per cluster
top_kegg_per_cluster <- lapply(kegg_results, function(res) {
  if (is.null(res) || nrow(res) == 0) return(data.frame())
  as.data.frame(res)[1:min(10, nrow(res)), c("ID", "Description", "p.adjust")]
})

# Step 4: Name the results by cluster
names(top_kegg_per_cluster) <- paste0("Cluster_", names(specific_genes))

# View example output
top_kegg_per_cluster[["Cluster_1"]]

#-----------------

# Step 2: Run GO:BP enrichment analysis
go_results <- lapply(gene_lists_entrez, function(entrez_ids) {
  enrichGO(gene = entrez_ids,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05)
})

# Step 3: Extract top 3 GO:BP terms per cluster
top_go_per_cluster <- lapply(go_results, function(res) {
  if (is.null(res) || nrow(res) == 0) return(data.frame())
  as.data.frame(res)[1:min(25, nrow(res)), c("ID", "Description", "p.adjust")]
})

# Step 4: Name results by cluster
names(top_go_per_cluster) <- paste0("Cluster_", names(specific_genes))

# View example output
top_go_per_cluster[["Cluster_1"]]

#------------
# Step 2: Run Reactome pathway enrichment per cluster
reactome_results <- lapply(gene_lists_entrez, function(entrez_ids) {
  enrichPathway(gene = entrez_ids,
                organism = "human",
                pvalueCutoff = 0.05,
                readable = TRUE)
})

# Step 3: Extract top 3 pathways per cluster
top_reactome_per_cluster <- lapply(reactome_results, function(res) {
  if (is.null(res) || nrow(res) == 0) return(data.frame())
  as.data.frame(res)[1:min(10, nrow(res)), c("ID", "Description", "p.adjust")]
})

# Step 4: Name results by cluster
names(top_reactome_per_cluster) <- paste0("Cluster_", names(specific_genes))

# Example output
top_reactome_per_cluster[["Cluster_1"]]

#-----------------------

# Load required libraries
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Step 1: Get MSigDB Hallmark gene sets as a data frame
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")[, c("gs_name", "gene_symbol")]

# Step 2: Run enrichment for each cluster’s gene list
msigdb_results <- lapply(specific_genes, function(glist) {
  enricher(gene = glist,
           TERM2GENE = hallmark_df,
           pvalueCutoff = 0.05)
})

# Step 3: Extract top 3 enriched Hallmark pathways per cluster
top_msigdb_per_cluster <- lapply(msigdb_results, function(res) {
  if (is.null(res) || nrow(res) == 0) return(data.frame())
  as.data.frame(res)[1:min(10, nrow(res)), c("ID", "Description", "p.adjust")]
})

# Step 4: Name results by cluster
names(top_msigdb_per_cluster) <- paste0("Cluster_", names(specific_genes))

# Example: View Cluster 1 results
top_msigdb_per_cluster[["Cluster_1"]]


top_kegg_per_cluster[["Cluster_6"]]
top_go_per_cluster[["Cluster_6"]]
top_reactome_per_cluster[["Cluster_6"]]
top_msigdb_per_cluster[["Cluster_6"]]


##############################. TCGA VALIDATION  #############################

# ────────────────────────────────────────────────────────────────
# 📦 Load packages
# ────────────────────────────────────────────────────────────────

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)

# ────────────────────────────────────────────────────────────────
# 📥 Step 1: Download and prepare TCGA-PCPG RNA-seq STAR-Counts
# ────────────────────────────────────────────────────────────────
setwd("~/Library/CloudStorage/OneDrive-uoflhealth/Atlas/PHPG")
query <- GDCquery(
  project = "TCGA-PCPG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
exp_data <- GDCprepare(query)  # SummarizedExperiment object

# ────────────────────────────────────────────────────────────────
# 🧬 Step 2: Extract full raw count matrix (Ensembl IDs)
# ────────────────────────────────────────────────────────────────
expr_matrix <- assay(exp_data)                   # genes x samples
rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))  # strip Ensembl version

# ────────────────────────────────────────────────────────────────
# 🔁 Step 3: Normalize full matrix using DESeq2::vst()
# ────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = expr_matrix,
                              colData = colData(exp_data),
                              design = ~1)

vst_data <- assay(vst(dds, blind = TRUE))    # genes x samples
rownames(vst_data) <- gsub("\\..*", "", rownames(vst_data))

# ────────────────────────────────────────────────────────────────
# 🧪 Step 4: Subset vst_data to cluster-specific genes
# ────────────────────────────────────────────────────────────────
#specific_genes is from phpg.R; 40 works well
# Assuming you have specific_genes list from earlier
top_genes_all <- unique(unlist(specific_genes))

# Get gene annotations for mapping Ensembl → SYMBOL
gene_annotations <- rowData(exp_data)
ensembl_to_symbol <- data.frame(
  Ensembl = gsub("\\..*", "", rownames(expr_matrix)),
  Symbol = gene_annotations$gene_name
)

# Filter and rename
gene_map <- ensembl_to_symbol[ensembl_to_symbol$Symbol %in% top_genes_all, ]
vst_filtered <- vst_data[gene_map$Ensembl, ]
rownames(vst_filtered) <- gene_map$Symbol

# ────────────────────────────────────────────────────────────────
# 🔥 Step 5: Generate scaled heatmap
# ────────────────────────────────────────────────────────────────
# Ensure the order of genes matches top_genes_all
top_genes_ordered <- intersect(top_genes_all, rownames(vst_filtered))
vst_filtered_ordered <- vst_filtered[top_genes_ordered, ]

# After scaling
expr_scaled <- scale(t(vst_filtered_ordered))
expr_scaled <- t(expr_scaled)

# ✅ Clip extreme values to improve color contrast
expr_scaled[expr_scaled > 2] <- 2
expr_scaled[expr_scaled < -2] <- -2

heat_out<-pheatmap::pheatmap(expr_scaled,
                   cluster_cols = TRUE,
                   cluster_rows = F,
                   clustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   show_colnames = FALSE,
                   fontsize_row = 8,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   main = "TCGA-PCPG Heatmap (Cluster-Specific Genes)",
                   silent = TRUE)

# Step 2: Cut column dendrogram into 6 clusters
col_dend <- heat_out$tree_col
sample_clusters <- cutree(col_dend, k = 6)

sample_clusters[sample_clusters == 3] <- 99     # temp label
sample_clusters[sample_clusters == 4] <- 3      # 4 → 3
sample_clusters[sample_clusters == 99] <- 4     # temp (was 3) → 4
sample_clusters["TCGA-RW-A68F-01A-11R-A35K-07"] <- 3

# Step 3: Create and apply annotation
annotation_tcga <- data.frame(Cluster = factor(sample_clusters))
rownames(annotation_tcga) <- names(sample_clusters)

# Step 1: Order samples by cluster assignment
ordered_samples <- rownames(annotation_tcga)[order(annotation_tcga$Cluster)]

# Step 2: Reorder expression matrix and annotation
expr_ordered <- expr_scaled[, ordered_samples]
annotation_tcga <- annotation_tcga[ordered_samples, , drop = FALSE]

# Step 3: Plot heatmap
tcga<-pheatmap::pheatmap(expr_ordered,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               annotation_col = annotation_tcga,
               show_colnames = FALSE,
               fontsize_row = 8,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
               main = "TCGA PH/PG Samples (RNA-seq)")


#T-sne of TCGA
set.seed(1)
source("~/Desktop/FIt-SNE-master/fast_tsne.R")
setwd("~/Library/CloudStorage/OneDrive-uoflhealth/Atlas/PHPG")

# Precompute MAD and order genes
gene_mad <- apply(vst_data, 1, function(x) mad(x, constant = 1))
gene_order <- names(sort(gene_mad, decreasing = TRUE))

# Define gene proportions and perplexities
gene_props <- c(0.25,0.33)
perplexities <- c(20,25)

# Create empty list to collect tSNE results
tsne_results <- list()

# Loop through combinations
for (p in perplexities) {
  for (gprop in gene_props) {
    ngenes <- round(length(gene_mad) * gprop)
    top_genes <- gene_order[1:ngenes]
    dat <- vst_data[top_genes, ]
    
    tsne_out <- fftRtsne(t(dat), perplexity = p)
    tsne_df <- as.data.frame(tsne_out)
    colnames(tsne_df) <- c("tSNE1", "tSNE2")
    tsne_df$Sample <- colnames(dat)
    tsne_df$Perplexity <- paste("Perplexity =", p)
    tsne_df$GeneProp <- paste0("Top ", gprop * 100, "% genes")
    
    tsne_results[[paste(p, gprop, sep = "_")]] <- tsne_df
  }
}

# Combine all results
tsne_all <- bind_rows(tsne_results)
annotation_tcga$Sample<-rownames(annotation_tcga)
tsne_all <- merge(tsne_all, annotation_tcga, by = "Sample")
tsne_all$sample_short <- substr(tsne_all$Sample, 1, 16)

# Factorize for plotting
tsne_all$Perplexity <- factor(tsne_all$Perplexity, levels = paste("Perplexity =", perplexities))
tsne_all$GeneProp <- factor(tsne_all$GeneProp, levels = paste0("Top ", gene_props * 100, "% genes"))

#tsne_all$Cluster[tsne_all$Sample == "TCGA-RW-A68F-01A-11R-A35K-07"] <- "Metabolic"

# Plot with facet grid
tcgatsne <- ggplot(tsne_all[tsne_all$Perplexity=="Perplexity = 20" & tsne_all$GeneProp=="Top 33% genes",], aes(tSNE1, tSNE2, color = factor(Cluster))) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Set2") + scale_x_reverse()+scale_y_reverse()+
  theme_pubr() + theme(legend.position = "none", plot.title = element_text(size = 13)) + 
  #facet_grid(Perplexity ~ GeneProp) +
  xlab("t-SNE1") + ylab("t-SNE2") + labs(title = "TCGA Atlas\nPH/PG Samples\nN=187")


# Define the row gaps
row_gaps_indices <- c(30, 60, 90, 118, 148)  # 90+28=118, 90+28+30=148
row_gaps <- rep(unit(2, "mm"), length(row_gaps_indices))

# Determine which rows to highlight
highlight_genes <- intersect(genes_to_show, rownames(expr_ordered))
highlight_indices <- which(rownames(expr_ordered) %in% highlight_genes)

# Create the right annotation with gene marks
row_anno <- rowAnnotation(
  Highlight = anno_mark(
    at = highlight_indices,
    labels = highlight_genes,
    labels_gp = gpar(fontsize = 8)
  )
)

# Define color scale for expression
col_fun <- colorRamp2(c(min(expr_ordered), 0, max(expr_ordered)),
                      c("navy", "white", "firebrick3"))

# Define row splits by gene subtype
n <- nrow(expr_ordered)
split_labels <- character(n)
split_labels[1:30] <- "Neuronal"
split_labels[31:60] <- "Vascular"
split_labels[61:90] <- "Metabolic"
split_labels[91:115] <- "Indeterminate"
split_labels[116:144] <- "Steroidal"
split_labels[145:n] <- "Developmental"
split_labels <- factor(split_labels, levels = c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental"))

# Format Cluster factor and assign Set2 colors
annotation_tcga$Cluster <- factor(annotation_tcga$Cluster, 
                                  labels = c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental"))

cluster_levels <- levels(annotation_tcga$Cluster)
cluster_colors <- setNames(brewer.pal(length(cluster_levels), "Set2"), cluster_levels)


#tcga_muts_collapsed from pcpg_muts_reassign.R
# Step 1: Add Sample column from rownames
annotation_tcga <- annotation_tcga %>%
  rownames_to_column(var = "Sample")

# Step 2: Create a short sample ID in both dataframes
annotation_tcga <- annotation_tcga %>%
  mutate(sample_short = substr(Sample, 1, 16))

tcga_muts_collapsed <- tcga_muts_collapsed %>%
  mutate(sample_short = substr(Sample, 1, 16))

# Step 3: Perform the join by short sample ID
annotation_tcga <- annotation_tcga %>%
  left_join(tcga_muts_collapsed %>% select(-Sample), by = "sample_short")

annotation_tcga<-annotation_tcga[!duplicated(annotation_tcga$sample_short),]

annotation_tcga$Primary_Mutations<-factor(annotation_tcga$Primary_Mutations, levels = c("NF1","RET","HRAS","BRAF","FGFR1","pseudo RET",
                                          "VHL","EPAS1","CSDE1", "pseudo VHL",
                                          "SDHx","TGDS","pseudo SDHx",
                                          "MAML3"),
                                          labels = c("NF1","RET","HRAS","BRAF","FGFR1","pseudo RET",
                                                     "VHL","EPAS1","CSDE1", "pseudo VHL",
                                                     "SDHx","TGDS","pseudo SDHx",
                                                     "MAML3 Fus."))

primcols <- c("NF1" = "darkgreen","RET" = "green","HRAS" = "olivedrab4","BRAF" = "olivedrab3","FGFR1" = "olivedrab2", "pseudo RET" = "olivedrab1",
              "VHL" = "orangered","EPAS1" = "orange","CSDE1" = "lightsalmon", "pseudo VHL" = "lightsalmon3",
              "SDHx" = "blue","TGDS" = "slateblue", "pseudo SDHx" = "slateblue4",
              "MAML3 Fus."="yellow")

# Create colored top annotation
top_annot <- HeatmapAnnotation(
  Cluster = annotation_tcga$Cluster,
  Mutation = annotation_tcga$Primary_Mutations,
  col = list(Cluster = cluster_colors,
             Mutation = primcols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10),  # 👈 controls "Cluster" label
  annotation_legend_param = list(
    Cluster = list(
      title_gp = gpar(fontsize = 10),       # legend title
      labels_gp = gpar(fontsize = 10)       # legend labels
    )
  )
)

# Match column order for column_split
# First, explicitly drop rownames
rownames(annotation_tcga) <- NULL
annotation_tcga_named <- annotation_tcga %>%
  column_to_rownames("Sample")
column_split_tcga <- annotation_tcga_named[colnames(expr_ordered), "Cluster"]

# Generate the final ComplexHeatmap
tcga<-Heatmap(expr_ordered,
              name = "Expression",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = F,
              show_column_names = FALSE,
              top_annotation = top_annot,
              right_annotation = row_anno,
              row_split = split_labels,
              row_title = NULL,
              row_gap = unit(1, "mm"),
              column_split = column_split_tcga,
              column_gap = unit(1, "mm"),
              column_title = "TCGA - PH/PG Samples (RNA-seq; N=187)",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10),   # heatmap legend title
                labels_gp = gpar(fontsize = 10)))

tcga<-draw(tcga,
           heatmap_legend_side = "right",           # place heatmap legend at the bottom
           annotation_legend_side = "right",        # place annotation legend at the bottom
           merge_legend = TRUE)          # ← THIS arranges them side by side



#gene expression in Affymetrix Microaaray Atlas normalized with RMA
ghr_expr <- data.frame(
  Sample = rownames(exp),
  GHR = exp[,"SST"],
  Cluster = annotation[rownames(exp), "Cluster"])
a<-ggplot(ghr_expr, aes(x = Cluster, y = GHR, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, aes(color = Cluster)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Atlas",
       x = "Cluster",
       y = "Expression") +
  stat_compare_means()+
  theme_pubr()+
  theme(legend.position = "none")

#gene expression in TCGA RNA-seq atlas normalized with DESeq2::vst()
# Transpose vst_filtered to have samples as rows
vst_df <- as.data.frame(t(vst_filtered))
vst_df$Sample <- rownames(vst_df)
plot_df <- merge(vst_df, annotation_tcga, by = "Sample")
b<-ggplot(plot_df, aes(x = Cluster, y = SST, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.6, aes(color = Cluster)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "TCGA",
       x = "Cluster",
       y = "Expression (VST)") +
  stat_compare_means()+
  theme_pubr() +
  theme(legend.position = "none")
a+b


library(cowplot)
library(grid)
library(ComplexHeatmap)

pdf_height <- 15.5
pdf_width <- 15.5

# tcga panel occupies 4/9 of height (from nested stacking math)
tcga_height <- pdf_height * (4 / 9)

# tcga panel width is 2.6 out of 5.6 columns
tcga_width <- pdf_width * (2.6 / 5.6)

# Set device height for first part
h <- 6

# Start PDF device (increase total height to fit everything)
pdf("PGPH_final.pdf", width = pdf_width, height = pdf_height, useDingbats = FALSE)

# Capture the heatmaps as grobs
heatmap_grob <- grid.grabExpr(draw(heatmap_obj), height = h)
tcga_grob <- grid.grabExpr(draw(tcga), width = tcga_width, height = tcga_height)

# Convert grobs to cowplot objects
heatmap_gg <- ggdraw() + draw_grob(heatmap_grob)
tcga_gg <- ggdraw() + draw_grob(tcga_grob)

# Combine tsne plots vertically (A and D)
tsne_combo <- plot_grid(
  tsne_plot, tcgatsne,
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("A", "D")
)

# Combine tsne_combo + heatmap + tcga heatmap horizontally (A–D)
top_four_panels <- plot_grid(
  tsne_combo, heatmap_gg, tcga_gg,
  labels = c("", "B", "C"),
  nrow = 1,
  rel_widths = c(1, 2, 2.6)
)

# Add wc underneath (Panel E)
top_five_panels <- plot_grid(
  top_four_panels,
  wc,
  labels = c("", "E"),
  ncol = 1,
  rel_heights = c(1, 0.5)  # adjust height of wc relative to above
)

# Now create the bottom row: plot1 + plot2 + plot3 + hubtfs (Panels F–I)
bottom_row <- plot_grid(
  hubtfs,plot1, plot2, plot3,
  labels = c("F", "G", "H", "I"),
  nrow = 1,
  rel_widths = c(2.5, 1, 1, 1)
)

# Final assembly: top_five_panels + bottom_row
full_final_plot <- plot_grid(
  top_five_panels,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 0.5)  # adjust bottom row size if needed
)

# Print to PDF
print(full_final_plot)

# Close device
dev.off()



### MUTATION P values

# Example: Perform Fisher's test for each gene
unique_mutations <- na.omit(unique(annotation_tcga$Primary_Mutations))

# Empty vector to store p-values
pvals <- c()

for (gene in unique_mutations) {
  # Build a 2xN table: Gene vs All Others
  gene_vs_rest <- table(annotation_tcga$Primary_Mutations == gene, annotation_tcga$Cluster)
  res <- fisher.test(gene_vs_rest)
  pvals <- c(pvals, res$p.value)
}
# Suppose pvals = your 11 raw p-values
# Add (97 - 11) = 86 dummy p-values of 1
pvals_extended <- c(pvals, rep(1, 97 - length(pvals)))

# Now adjust
pvals_fdr_extended <- p.adjust(pvals_extended, method = "fdr")

# Extract only the first 11 adjusted p-values
pvals_fdr_final <- pvals_fdr_extended[1:length(pvals)]

# Combine into a nice data frame
result_table <- data.frame(
  Gene = unique_mutations,
  Raw_p = pvals,
  FDR_p = pvals_fdr_final
)

# Round the p-values and avoid scientific notation
result_table$Raw_p <- format(round(result_table$Raw_p, 6), scientific = FALSE)
result_table$FDR_p <- format(round(result_table$FDR_p, 6), scientific = FALSE)

# Print
print(result_table, row.names = FALSE)



# Step 1: Define the Kinase genes
kinase_genes <- c("NF1", "HRAS", "BRAF", "RET", "FGFR1")

# Step 2: Create a new variable: does the sample have any kinase mutation?
annotation_tcga$Kinase_mut <- ifelse(annotation_tcga$Primary_Mutations %in% kinase_genes, "Yes", "No")

# Step 3: Build the 2xN contingency table
kinase_vs_cluster <- table(annotation_tcga$Kinase_mut, annotation_tcga$Cluster)

# Step 4: Perform Fisher's exact test
kinase_fisher_result <- fisher.test(kinase_vs_cluster, workspace = 2e7)  # 20 million workspace

# Step 5: Get the raw p-value
kinase_raw_p <- kinase_fisher_result$p.value

# Step 6: Pad with 96 dummy p-values of 1 (so total 97 tests)
kinase_pvals_extended <- c(kinase_raw_p, rep(1, 97 - 1))

# Step 7: Perform FDR adjustment
kinase_pvals_fdr_extended <- p.adjust(kinase_pvals_extended, method = "fdr")

# Step 8: Extract the FDR-adjusted p-value for Kinase
kinase_fdr_p <- kinase_pvals_fdr_extended[1]

# Step 9: Print both raw and FDR-adjusted p-values
cat("Kinase Raw p-value:", kinase_raw_p, "\n")
cat("Kinase FDR-adjusted p-value (97 tests):", kinase_fdr_p, "\n")

# Create a Mutation_Group variable
annotation_tcga$Mutation_Group[is.na(annotation_tcga$Mutation_Group)] <- "None"

# Calculate proportions
mutation_proportions <- annotation_tcga %>%
  group_by(Cluster, Mutation_Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Proportion = round(100 * n / sum(n), 1)) %>%
  arrange(Cluster, desc(Proportion))

# View the final table
print(mutation_proportions[mutation_proportions$Mutation_Group!="None",])
