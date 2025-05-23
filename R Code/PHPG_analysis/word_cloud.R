library(ggwordcloud)
library(dplyr)
library(tidyr)
library(textstem)
library(stopwords)
library(ggplot2)

#GO: using top 25

names(top_go_per_cluster)<-c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental")

# Extract all raw words from all clusters (no lemmatization or synonym mapping)
raw_word_list <- purrr::imap_dfr(top_go_per_cluster, function(df, cluster_id) {
  df %>%
    mutate(Cluster = as.character(cluster_id)) %>%
    select(Cluster, Description) %>%
    mutate(Description = tolower(Description)) %>%
    mutate(words = strsplit(Description, " ")) %>%
    unnest(words) %>%
    filter(!words %in% stopwords::stopwords("en"))
})

# Optional: Get overall word frequencies across all clusters
word_counts <- raw_word_list %>%
  count(words, sort = TRUE)

excluded_words <- c("regulation","process","response","acid","negative","positive","via",
                    "structure","system","compound","molecules","junction","levels","projection",
                    "activity","long-term","b","decreased","external","family","hydroxy","secondary",
                    "organization","structure","assembly")

word_counts_filtered <- word_counts %>%
  filter(!words %in% excluded_words)

#curate
synonym_map <- c(
  # Regulation terms
  "positive" = "regulation",
  "negative" = "regulation",
  
  # Development / differentiation / morphogenesis
  "morphogenesis" = "development",
  "differentiation" = "development",
  
  # Metabolism & biosynthesis
  "metabolic" = "metabolism",
  "biosynthetic" = "metabolism",
  "acid" = "amino_acid",
  "amino" = "amino_acid",
  "alpha-amino" = "amino_acid",
  "glutamine" = "amino_acid",
  "non-proteinogenic" = "amino_acid",
  "one-carbon" = "metabolism",
  "folic" = "folate",
  "pteridine-containing" = "folate",
  "alcohol" = "metabolism",
  "sterol" = "steroid",
  "cholesterol" = "steroid",
  "lipid" = "steroid",
  "secondary" = "metabolism",
  "organic" = "metabolism",
  "inorganic" = "metabolism",
  
  # Synapse-related
  "synaptic" = "synapse",
  "presynapse" = "synapse",
  "postsynaptic" = "synapse",
  "plasticity" = "synapse",
  "potentiation" = "synapse",
  "axon" = "synapse",
  "axonogenesis" = "synapse",
  "transmission" = "synapse",
  "signal" = "signaling",
  "transduction" = "signaling",
  
  # Signaling / pathway
  "signaling" = "signaling",
  "pathway" = "signaling",
  "jak-stat" = "signaling",
  "stat" = "signaling",
  "3-kinase/protein" = "signaling",
  "phosphatidylinositol" = "signaling",
  
  # Structural organization
  "matrix" = "ecm",
  "collagen" = "ecm",
  "basement" = "ecm",
  "extracellular" = "ecm",
  "encapsulating" = "ecm",
  
  # Adhesion / junctions
  "adhesion" = "adhesion",
  "junction" = "adhesion",
  "cell-cell" = "cell",
  "membrane" = "membrane",
  "plasma" = "membrane",
  "plasma-membrane" = "membrane",
  
  # Developmental anatomy
  "nephron" = "kidney",
  "glomerulus" = "kidney",
  "renal" = "kidney",
  "tubule" = "kidney",
  "tube" = "kidney",
  "endothelial" = "vascular",
  "endothelium" = "vascular",
  "epithelial" = "epithelium",
  "epithelium" = "epithelium",
  "heart" = "vascular",
  "retina" = "eye",
  "camera-type" = "eye",
  "artery" = "vascular",
  "vasculature" = "vascular",
  "vasculogenesis" = "vascular",
  "angiogenesis" = "vascular",
  
  # Hormones & peptides
  "neuropeptide" = "peptide",
  "peptide" = "peptide",
  "hormone" = "hormone",
  "insulin" = "hormone",
  "pituitary" = "hormone",
  
  # Nervous system
  "neuron" = "nervous_system",
  "nervous" = "nervous_system",
  "cognition" = "brain_function",
  "learning" = "brain_function",
  "memory" = "brain_function",
  "behavior" = "brain_function",
  
  # Others / specific functions
  "stress" = "stress",
  "hypoxia" = "oxygen",
  "decreased" = "response",
  "stimulus" = "signaling",
  "reproductive" = "reproduction",
  "ossification" = "bone",
  "bone" = "bone",
  "mineralization" = "bone",
  "zinc" = "metal",
  "copper" = "metal",
  "metal" = "metal",
  "toxic" = "xenobiotic",
  "xenobiotic" = "xenobiotic",
  "detoxification" = "xenobiotic",
  "fatty" = "steroid",
  "activity" = "function",
  "molecules" = "molecule",
  "levels" = "level",
  "projection" = "projection",
  "guidance" = "guidance",
  "transport" = "transport",
  "protein" = "protein",
  "compound" = "compound",
  "b" = "misc",
  "eye" = "eye",
  "system" = "system",
  "homophilic" = "adhesion",
  "external" = "external",
  "substance" = "compound",
  "family" = "group",
  "hydroxy" = "chemical",
  "humoral" = "immune",
  "antifungal" = "immune",
  "inhibitory" = "regulation"
)

# Re-count words within each cluster after excluding terms
word_counts_clustered <- raw_word_list %>%
  filter(!words %in% excluded_words) %>%
  mutate(mapped_word = recode(words, !!!synonym_map)) %>%
  count(Cluster, mapped_word, sort = TRUE) %>%
  rename(words = mapped_word)

highlight_words <- list(
  "1" = c("synapse", "brain_function"),
  "2" = c("vascular","oxygen"),
  "3" = c("metabolism", "amino_acid"),
  "4" = c(""),
  "5" = c("steroid", "metabolism"),
  "6" = c("development", "nervous_system"))

  library(RColorBrewer)
  dark2_colors <- brewer.pal(6, "Set2")
  names(dark2_colors) <- as.character(1:6)  # Match cluster IDs
  
  names(highlight_words) <- c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental")
  names(dark2_colors) <- names(highlight_words)  # Match exactly  
  
  # Mark words for highlight
  word_counts_clustered <- word_counts_clustered %>%
    mutate(highlight = mapply(function(word, Cluster) {
      word %in% highlight_words[[as.character(Cluster)]]
    }, words, Cluster),
    highlight_color = ifelse(highlight, dark2_colors[Cluster], "gray80"))

word_counts_clustered$Cluster<-factor(word_counts_clustered$Cluster, levels = c("Neuronal", "Vascular", "Metabolic", "Indeterminate", "Steroidal", "Developmental"))
  
ggplot(word_counts_clustered, aes(label = words, size = n, color = highlight_color)) +
  geom_text_wordcloud(area_corr_power = 1) +
  facet_wrap(~ Cluster, nrow =1) +
  scale_size_area(max_size = 10) +
  scale_color_identity() +  # Use the actual hex values as is
  theme_pubr()

library(ggh4x)
wc<-ggplot(word_counts_clustered, aes(label = words, size = n, color = highlight_color)) +
  geom_text_wordcloud() +
  facet_wrap2(
    ~ Cluster,
    nrow = 1,
    strip = strip_themed(
      background_x = lapply(unname(dark2_colors), function(col) element_rect(fill = col, color = "black")),
      text_x = lapply(seq_along(dark2_colors), function(i) element_text(color = "black", face = "bold"))
    )
  ) +
  scale_size_area(max_size = 10) +
  scale_color_identity() +
  theme_pubr()+
  ggtitle("    GO Enrichment Analysis: Biological Processes")
