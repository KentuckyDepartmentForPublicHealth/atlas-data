library(ggpubr)
library(dplyr)
library(patchwork)

d<-read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv", header=T, row.names = 1, stringsAsFactors = T)
d<-d[d$Train.diagnosis_bc!="#N/A",]

t<-d %>% 
  group_by(Train.diagnosis_bc) %>% 
  summarise(n_dataset = n_distinct(Dataset.ID)) %>%
  arrange(desc(n_dataset))

ba <- function(diagnosis) {
  # Subset data for the specific diagnosis
  z <- d[d$Train.diagnosis_bc == diagnosis,]
  
  # Calculate min and max values for x and y axes across both plots
  bxmin <- min(z$tsne1)
  bxmax <- max(z$tsne1)
  axmin <- min(z$tsne1_bc)
  axmax <- max(z$tsne1_bc)
  
  bymin <- min(z$tsne2)
  bymax <- max(z$tsne2)
  aymin <- min(z$tsne2_bc)
  aymax <- max(z$tsne2_bc)
  
  # Calculate the range of both axes
  bx_range <- bxmax - bxmin  # x-range for plot b
  by_range <- bymax - bymin  # y-range for plot b
  ax_range <- axmax - axmin  # x-range for plot a
  ay_range <- aymax - aymin  # y-range for plot a
  
  # Determine the largest x and y range across both plots
  max_x_range <- max(bx_range, ax_range)
  max_y_range <- max(by_range, ay_range)
  
  # Determine if plot 'b' or 'a' already has the max range for x or y axes
  if (bx_range < max_x_range) {
    # Plot 'b' needs expansion on the x-axis
    xlim_bmin <- (bxmin + bxmax) / 2 - max_x_range / 2
    xlim_bmax <- (bxmin + bxmax) / 2 + max_x_range / 2
  } else {
    # Plot 'b' already has the max x range, keep as is
    xlim_bmin <- bxmin
    xlim_bmax <- bxmax
  }
  
  if (ax_range < max_x_range) {
    # Plot 'a' needs expansion on the x-axis
    xlim_amin <- (axmin + axmax) / 2 - max_x_range / 2
    xlim_amax <- (axmin + axmax) / 2 + max_x_range / 2
  } else {
    # Plot 'a' already has the max x range, keep as is
    xlim_amin <- axmin
    xlim_amax <- axmax
  }
  
  if (by_range < max_y_range) {
    # Plot 'b' needs expansion on the y-axis
    ylim_bmin <- (bymin + bymax) / 2 - max_y_range / 2
    ylim_bmax <- (bymin + bymax) / 2 + max_y_range / 2
  } else {
    # Plot 'b' already has the max y range, keep as is
    ylim_bmin <- bymin
    ylim_bmax <- bymax
  }
  
  if (ay_range < max_y_range) {
    # Plot 'a' needs expansion on the y-axis
    ylim_amin <- (aymin + aymax) / 2 - max_y_range / 2
    ylim_amax <- (aymin + aymax) / 2 + max_y_range / 2
  } else {
    # Plot 'a' already has the max y range, keep as is
    ylim_amin <- aymin
    ylim_amax <- aymax
  }
  
  # Plot before SVA adjustment (plot b)
  b <- ggplot() +
    geom_point(data=z, aes(tsne1, tsne2, color=Dataset.ID), size=0.75) +
    theme_pubr() +
    theme(legend.position = "none") +  # No legend in this plot
    rremove("axis") + rremove("axis.text") + rremove("ticks") + rremove("xylab") +
    theme(panel.border = element_rect(colour = "black", fill=NA)) +
    ggtitle(diagnosis, subtitle = "Before ComBat") +
    xlim(xlim_bmin, xlim_bmax) +  # Set x-axis range symmetrically if needed
    ylim(ylim_bmin, ylim_bmax)    # Set y-axis range symmetrically if needed
  
  # Plot after SVA adjustment (plot a)
  a <- ggplot() +
    geom_point(data=z, aes(tsne1_bc, tsne2_bc, color=Dataset.ID), size=0.75) +
    theme_pubr() +
    theme(legend.position = "none",  # Legend position on the right
          legend.justification = c(0, 0.5)) +  # Justify the legend to the left
    rremove("axis") + rremove("axis.text") + rremove("ticks") + rremove("xylab") +
    theme(panel.border = element_rect(colour = "black", fill=NA)) +
    guides(color = guide_legend(nrow = 8)) +  # Force 8 rows in the legend
    ggtitle("", subtitle = "After ComBat") +
    xlim(xlim_amin, xlim_amax) +  # Set x-axis range symmetrically if needed
    ylim(ylim_amin, ylim_amax)    # Set y-axis range symmetrically if needed
  
  # Combine plots vertically with patchwork
  b + a
}


(ba("Diffuse Glioma-WT") / ba("PA") / ba("Neuroblastoma")  / ba("PFA") / ba("Meningioma") / ba("RELA") | 
  
  plot_spacer() / plot_spacer() / plot_spacer() / plot_spacer() / plot_spacer() / plot_spacer() |

ba("PCNSL") / ba("PitNET") / ba("Midline") / ba("MPNST") / ba("Lung") /ba("Cerebellum")) + plot_layout(widths = c(10,1,10))

#850 x 1400
