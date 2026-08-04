# Load necessary libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)
library(here)
library(plotrix)
library(scales)
library(car)
library(lattice)
library(multcompView)
library(FSA)
library(multcompView)
library(RColorBrewer)
library(cluster)
library(effectsize)
library(MASS)
library(vegan)
library(adegenet)
library(factoextra)
library(tidyverse)
library(ggridges)

#### Functions ####
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} # error bar function

generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY$label[,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$labels=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$label) , ]
  return(Tukey.labels)
} # function to generate Tukey groups for plotting

#### Read in habitat measurement dataset that has been filtered for populated with the handful of gap_fill measurements that were redone ####
mad_ann_seq_core_gfixed2 <- read.csv(here("3D_models/output_for_R/hab_measures_filtered_newgapfill.csv"), header = TRUE) #299 x 30

# Set color schemes for plotting
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#b0e0e6", "cornflowerblue")
colors_depth <- c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
colors_site <- c("#74D055FF", "#238A8DFF", "#453781FF")

# Perform PCA analysis: correlation matrix
pca_mad <- prcomp(mad_ann_seq_core_gfixed2[, c("vector_disp", "gap_fill_new", "elevation")], center = TRUE, scale. = TRUE) # correlation matrix pca

summ <- summary(pca_mad) # calculates variance
summ$importance

png(file=here("3D_models/figures/microhabitat_pca.png"), width=8, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=16, # point size, which is the font size
  mfrow = c(2,2)
)

plot(pca_mad$x, col = add.alpha(colors_mpha[factor(mad_ann_seq_core_gfixed2$Species_mpha)], 0.7), pch = 19, xlab = 'PC 1 (44.5%)', ylab = 'PC 2 (31.9%)', main = 'Taxon', xlim = c(-4, 4.8), las = 1)
legend('topright', legend = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA1", "MPHA2", "MSEN1", "MSEN2"), col = add.alpha(colors_mpha, 0.7), pch = 19, cex = 0.8)

plot(pca_mad$x, col = add.alpha(colors_depth[factor(mad_ann_seq_core_gfixed2$depth)], 0.7), pch = 19, xlab = 'PC 1 (44.5%)', ylab = 'PC 2 (31.9%)', main = 'Depth', xlim = c(-4, 4.8), las = 1)
legend('topright', legend = c('5', '10', '20', '40', '60'), col = add.alpha(colors_depth, 0.7), pch = 19, cex = 1)

plot(pca_mad$x, col = add.alpha(colors_site[factor(mad_ann_seq_core_gfixed2$site)], 0.7), pch = 19, xlab = 'PC 1 (44.5%)', ylab = 'PC 2 (31.9%)', main = 'Site', xlim = c(-4, 4.8), las = 1)
legend('topright', legend = c("est", "kal", "sna"), col = add.alpha(colors_site, 0.7), pch = 19,cex = 1)

dev.off()

