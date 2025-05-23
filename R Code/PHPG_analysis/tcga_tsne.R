library(ggplot2)
library(ggpubr)
library(dplyr)

set.seed(1)
source("~/Desktop/FIt-SNE-master/fast_tsne.R")
setwd("~/Library/CloudStorage/OneDrive-uoflhealth/Atlas/PHPG")

# Precompute MAD and order genes
gene_mad <- apply(vst_data, 1, function(x) mad(x, constant = 1))
gene_order <- names(sort(gene_mad, decreasing = TRUE))

# Define gene proportions and perplexities
gene_props <- c(0.1,0.2,0.25,0.33)
perplexities <- c(10,15,20,25)

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

# Plot with facet grid
p <- ggplot(tsne_all[tsne_all$Perplexity=="Perplexity = 20" & tsne_all$GeneProp=="Top 33% genes",], aes(tSNE1, tSNE2, color = factor(Cluster))) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Set2") + scale_x_reverse()+
  theme_pubr() + theme(legend.position = "none", plot.title = element_text(size = 13)) + 
  #facet_grid(Perplexity ~ GeneProp) +
  xlab("t-SNE1") + ylab("t-SNE2") + labs(title = "TCGA Atlas\nPH/PG Samples\nN=187")


