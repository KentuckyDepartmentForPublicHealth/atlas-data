library(dplyr)
library(tidyr)

# Load data
d <- read.csv("inconst.csv", header = TRUE)

# Generate all unique unordered pairs per row
pair_list <- lapply(1:nrow(d), function(i) {
  row_values <- as.character(d[i, ])
  pairs <- combn(row_values, 2, simplify = FALSE)
  # Sort each pair to ensure (A,B) and (B,A) are treated the same
  lapply(pairs, function(p) sort(p))
})

# Flatten the list and convert to data frame
all_pairs <- do.call(rbind, lapply(pair_list, function(pairs) {
  do.call(rbind, lapply(pairs, function(p) data.frame(V1 = p[1], V2 = p[2])))
}))

# Keep only inconsistent pairs (V1 ≠ V2)
incon <- all_pairs %>% filter(V1 != V2)

# Count occurrences of each unique unordered pair
pair_counts <- incon %>%
  group_by(V1, V2) %>%
  summarise(COUNT = n(), .groups = "drop")

# Expand to one row per occurrence
final <- pair_counts %>% uncount(COUNT)

library(dplyr)
library(ggpubr)

# Step 1: Read training annotation and filter valid diagnoses
anno <- read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv",
                 header = TRUE, row.names = 1, stringsAsFactors = FALSE)

anno <- anno[anno$Train.diagnosis_bc != "#N/A", ]
training_freq <- table(anno$Train.diagnosis_bc)
training_prop <- prop.table(training_freq)

# Step 2: Count diagnosis occurrences in inconsistent pairs
incon_counts <- table(c(final$V1, final$V2))
incon_prop <- prop.table(incon_counts)

# Step 3: Normalize by training frequency
# Match and divide, only for common diagnoses
common_diagnoses <- intersect(names(incon_prop), names(training_prop))

# Step 3: Normalize and build the final dataframe
normalized_freq <- incon_prop[common_diagnoses] / training_prop[common_diagnoses]
z <- data.frame(
  Diagnosis = names(normalized_freq),
  NormalizedFreq = as.numeric(normalized_freq)
)


# Step 4: Dot plot
A<-ggdotchart(z,"Diagnosis", "NormalizedFreq",add="segments", sorting = "descending", 
           ylab = "Inconsistent Predictions Involving the Diagnosis\nNormalized By the Number of Training Samples")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#1000x550

#Generate a chord diagram
# Transform input data in a adjacency matrix
adjacencyData <- with(final, table(V1, V2)) #this is the same as f3, which can be entered into ChordDiagram

# Charge the circlize library
library(circlize)
chordDiagram(adjacencyData, transparency = 0.5)

#chord for only PA, GG, DIG, PXA - top 4 from above
pa<-final[final$V1=="DIG" | final$V2=="DIG" | 
          final$V1=="PXA" | final$V2=="PXA" |
          final$V1=="Ganglioglioma" | final$V2=="Ganglioglioma",]


# Define the 10 target diagnoses
target_dx <- c("DIG", "PXA", "Ganglioglioma")

# Reorder pairs so that the target diagnosis is in V1 when only one of them is in target
pa_fixed <- pa %>%
  rowwise() %>%
  mutate(
    # Count how many of the pair are in the target list
    in_target = sum(c(V1, V2) %in% target_dx),
    
    # If exactly one is in the list, put it in V1
    V1_new = if (in_target == 1) ifelse(V1 %in% target_dx, V1, V2) else V1,
    V2_new = if (in_target == 1) ifelse(V1 %in% target_dx, V2, V1) else V2
  ) %>%
  ungroup() %>%
  select(V1 = V1_new, V2 = V2_new)

# Check result
head(pa_fixed)

#pa<-read.csv("pa.csv", header = T, stringsAsFactors = T) #the excel file from above rearranged to have the above diagnosis in one column

adjacencyData <- with(pa_fixed, table(V1,V2))

#chordDiagram(adjacencyData, transparency = 0.5,  annotationTrack = c("name", "grid"),
#             annotationTrackHeight = c(0.03, 0.05), big.gap = 30)

chordDiagram(adjacencyData, transparency = 0.5, annotationTrack = "grid", big.gap = 30,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacencyData))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

#1200x1200
