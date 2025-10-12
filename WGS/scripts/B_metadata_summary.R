library(here)
library(stringr)
library(scales)
library(RColorBrewer)

ids <- read.csv(here("WGS/angsd_results/bam_lists/ids_recoded606_Oct2024.csv"), header = TRUE) # 606 Madracis, no clones

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata $id_original <- paste(metadata[,4], metadata[,5], sep = "_")
metadata$depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
metadata$sample_full <- ids$new_ID

# Get cluster assignments for coloring the points
data_cluster <- read.csv(here("madracis_master_merged_Oct2024.csv"), header = TRUE)

metadata$species_mpha <- data_cluster[match(metadata$id_original, data_cluster$ID), 7]

# Replace NAs with matching in species column
species_indices <- which(is.na(metadata$species_mpha) == TRUE) # indicies that have NAs
metadata$species_mpha[is.na(metadata$species_mpha)] <- metadata$species[species_indices] # replace the NAs with species in associated index

# Split out MPHA by depth to assign MPHA.ln3
metadata$species_mpha[which(metadata$species == 'MPHA' & metadata$depth >60)] <- "MPHA.ln3" # MPHA at depths >60

# Deal with the Pim MPHAs, which are mostly MSEN-like MPHAs
metadata[which(metadata$id_original == 'PB_4583'),9] <- 'MPHA.ln2'
metadata[which(metadata$id_original == 'PB_4578'),9] <- 'MPHA.ln4'
metadata[which(metadata$id_original == 'PB_4577'),9] <- 'MPHA.ln4'
metadata[which(metadata$id_original == 'PB_4576'),9] <- 'MPHA.ln4'
metadata[which(metadata$id_original == 'PB_4579'),9] <- 'MPHA.ln4'

metadata$species_mpha <- factor(metadata$species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MPHA.ln4', 'MSEN', 'UNKN')) #reorder the species factors

#### Set color scheme ####
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

colors_depth <- brewer.pal(8, 'Dark2')
show_col(colors_depth)

#### Start figuring out counts/percentages for species for a given depth ####
# What percent of species occur at a given depth
speciesbydepth <- as.data.frame.matrix(table(metadata$species_mpha, metadata$depth)) # table of species vs depths
count_by_depth <- as.vector(colSums(speciesbydepth)) # number of corals at each depth

speciesbydepth_percents <- data.frame(matrix(NA, nrow = 10, ncol = 8))
colnames(speciesbydepth_percents) <- colnames(speciesbydepth)
rownames(speciesbydepth_percents) <- rownames(speciesbydepth)

for (i in 1:ncol(speciesbydepth)) {
  speciesbydepth_percents[,i] <- speciesbydepth[,i]/count_by_depth[i]
}

#### Plot ####
# Barchart by proportion
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/species_prop_by_depth_bar.png", width=11, height=8.5, res=300, units="in")

par(
  mar=c(4, 4, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(3,3))

  for (j in 1:ncol(speciesbydepth_percents)) {
    barplot(height=as.matrix(t(speciesbydepth_percents))[j,], names=rownames(speciesbydepth_percents), col = colors_mpha, xlab = 'Species', ylab = 'Proportion', main = paste(colnames(speciesbydepth_percents)[j], 'meters', sep = ' '), ylim = c(0,1))
  }

dev.off()

# Pie chart
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/species_prop_by_depth_pie.png", width=9, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(3,3))

for (j in 1:ncol(speciesbydepth_percents)) {
  pie(speciesbydepth_percents[,j], labels = rownames(speciesbydepth_percents), col = colors_mpha, main = paste(colnames(speciesbydepth_percents)[j], 'meters', sep = ' '))
}

dev.off()

# Barchart by absolute count
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/species_count_by_depth_bar.png", width=4, height=18, res=300, units="in")

par(
  #mar=c(3, 4, 2, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=16, # point size, which is the font size
  bty = 'n',
  mfrow = c(8,1))

# Define a layout with equal panel heights but adjust space for the last plot
layout_matrix <- matrix(1:ncol(speciesbydepth), ncol = 1)
layout(layout_matrix, heights = c(rep(1, ncol(speciesbydepth) - 1), 1.2)) # Larger height for the last plot

for (j in 1:ncol(speciesbydepth)) {
  
  # Adjust margins for the last plot
  if (j == ncol(speciesbydepth)) {
    par(mar = c(7, 4.5, 2, 1)) # Increased bottom margin
  } else {
    par(mar = c(2, 4.5, 2, 1)) # Default margin for other plots
  }
  
  # Extract data for current depth
  counts <- as.matrix(t(speciesbydepth))[j,]
  
  # Create barplot and store midpoints
  bar_mids <- barplot(height = counts,
                      names = rownames(speciesbydepth),
                      col = colors_mpha,
                      xlab = '', ylab = 'Count',
                      ylim = c(0, 85), xaxt = 'n', las = 2)
  
  # Add depth label
  mtext(paste(colnames(speciesbydepth)[j], 'meters', sep = ' '), side = 3, adj = 0, line = 0.65)
  
  # Add sample size labels above bars (only if > 0)
  text(x = bar_mids[counts > 0],
       y = counts[counts > 0] + 4,  # Offset above the bar
       labels = counts[counts > 0],
       cex = 0.9)
}

# Add x-axis labels and axis title for the combined plot
text(x = c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9, 9.1, 10.3, 11.5),
     y = par("usr")[3] - 10,
     labels = c("MCAR", "MDEC", "MFOR", "MMIR","MPHA.ln1","MPHA.ln2", "MPHA.ln3", "MPHA.ln4","MSEN", "UNKN"),
     xpd = NA, srt = 35, adj = 0.965)

mtext('Genetic lineage', side = 1, line = 5.75)

dev.off()

#### Figure out counts/percentages for depths for a given species
# For a species, what is the depths at which they occur?
count_by_species <- as.vector(rowSums(speciesbydepth)) # number of corals of each species

depthbyspecies <- data.frame(matrix(NA, nrow = 10, ncol = 8))
colnames(depthbyspecies) <- colnames(speciesbydepth)
rownames(depthbyspecies) <- rownames(speciesbydepth)

for (i in 1:nrow(speciesbydepth)) {
  depthbyspecies[i,] <- speciesbydepth[i,]/count_by_species[i]
}

#### Plot ####
# Barplot
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/depth_prop_by_species_bar.png", width=12, height=16, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=17, # point size, which is the font size
  bty = 'n',
  mfrow = c(5,2))

for (j in 1:nrow(depthbyspecies)) {
  barplot(height=as.matrix(depthbyspecies[j,]), names=colnames(depthbyspecies), xlab = 'Depth', ylab = 'Proportion', main = paste(rownames(depthbyspecies)[j], ' (n = ', rowSums(speciesbydepth)[j], ')', sep = ''), col = colors_mpha[j], ylim = c(0,1), las = 1, cex.main = 1.5)
}

dev.off()

# Pie chart
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/depth_prop_by_species_pie.png", width=9, height=7, res=300, units="in")

par(
  mar=c(2, 1, 2, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(4,3))

for (j in 1:nrow(depthbyspecies)) {
  pie(as.numeric(depthbyspecies[j,]), labels = colnames(depthbyspecies), col = colors_depth, main = paste(rownames(depthbyspecies)[j], sep = ''))
}

dev.off()
