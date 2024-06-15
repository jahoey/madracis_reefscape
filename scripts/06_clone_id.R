## load packages
library(here)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(ape)
library(ggtree)
library(MASS)
library(adegenet)
library(dplyr)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### Read in data ####
m <- as.matrix(read.table(here("angsd_results/madracis672/all_snps.ibsMat"), header = FALSE)) # 28M SNPs after mapping to 14 somatic chromosomes & mt chromosome

# Load IDs file to assign species
ids <- read.csv(here("angsd_results/bam_lists/ids_recoded672.csv"), header = TRUE) # 672 colonies

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata$sample_id <- paste(metadata$collector, metadata$ID, sep = "_")

tip_names <- paste(metadata$species, metadata$collector, metadata$ID, sep = "_")

# Load IDs of known sample or technical replicates
reps <- read.csv(here("angsd_results/known_sample_and_tech_reps.csv"), header = TRUE)

#### Visualization ####
#heat map
heatmap(m, labRow = tip_names, labCol = tip_names, cexRow = 0.1, cexCol = 0.1)

#### Plot cluster dendrogram based on genetic distance (1-IBS)
# By species, PC1 vs PC2
png(file=here("angsd_results/madracis672_cluster_dendrogram.png"), width=40, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot(hclust(as.dist(m), "complete"), main = "Cluster Dendrogram: Madracis 672", labels = tip_names, cex = 0.3, ylab = 'Genetic distance (1-IBS)')

abline(h=0.05775805, col = '#ff9a81') #MDEC
abline(h=0.06650635, col = '#6b8e23') #MPHA
abline(h=0.0549302 , col = '#b0e0e6') #MSEN
abline(h=0.04816355, col = "#E6AB02") #MMIR

dev.off()

hist(as.numeric(m))

#### Subset IBS matrix to known sample or technical replicates by species by first getting extraction indices from metadata ####
# Get sample/technical rep indices from metadata for each species
# MDEC
which(reps$species == 'MDEC')
mdec.r <- match(reps$extraction_id1[2:5], metadata$extract)
mdec.c <- match(reps$extraction_id2[2:5], metadata$extract)

# Use mapply to extract the values from the matrix at the specified positions for MDEC
mdec.ibs <- mapply(function(r, c) m[r, c], mdec.r, mdec.c)

m[mdec.r[1],mdec.c[1]]
hist((mdec.ibs))
quantile(mdec.ibs, 0.95)

# MPHA
which(reps$species == 'MPHA')
mpha.r <- match(reps$extraction_id1[6:15], metadata$extract)
mpha.c <- match(reps$extraction_id2[6:15], metadata$extract)
mpha.ibs <- mapply(function(r, c) m[r, c], mpha.r, mpha.c)
hist((mpha.ibs))
quantile(mpha.ibs, 0.95)

# MSEN
which(reps$species == 'MSEN')
msen.r <- match(reps$extraction_id1[16:20], metadata$extract)
msen.c <- match(reps$extraction_id2[16:20], metadata$extract)
msen.ibs <- mapply(function(r, c) m[r, c], msen.r, msen.c)
hist((msen.ibs))
quantile(msen.ibs, 0.95)

# MMIR
which(reps$species == 'MMIR')
mmir.r <- match(reps$extraction_id1[21:22], metadata$extract)
mmir.c <- match(reps$extraction_id2[21:22], metadata$extract)
mmir.ibs <- mapply(function(r, c) m[r, c], mmir.r, mmir.c)
hist((mmir.ibs))
quantile(mmir.ibs, 0.95)

#### Figure out which additional colonies are likely to be clones based on sample/technical IBS distribution ####
# Find indices of matrix values below the cutoff
indices_below_cutoff <- which(m < 0.04816355, arr.ind = TRUE)

# Map indices to names in metadata
row_names_below_cutoff <- paste(metadata$species[indices_below_cutoff[, 1]], metadata$sample_id[indices_below_cutoff[, 1]], metadata$extract[indices_below_cutoff[, 1]], sep = '_')
col_names_below_cutoff <- paste(metadata$species[indices_below_cutoff[, 2]], metadata$sample_id[indices_below_cutoff[, 2]], metadata$extract[indices_below_cutoff[, 2]], sep = '_')

# Combine row and column names into a data frame
potential.clones <- data.frame(
  row_name = row_names_below_cutoff,
  col_name = col_names_below_cutoff
)

potential.clones.unique <- potential.clones[which(!potential.clones$row_name == potential.clones$col_name),] # potential clones without the self-self pairs

# Add a new column with sorted pairs
potential.clones.unique <- potential.clones.unique %>%
  mutate(sorted_pair = apply(potential.clones.unique, 1, function(x) paste(sort(x), collapse = "_")))

# Remove duplicate pairs
filtered_clones <- potential.clones.unique %>%
  distinct(sorted_pair, .keep_all = TRUE) %>%
  select(-sorted_pair)

# Display the filtered dataframe
filtered_clones

