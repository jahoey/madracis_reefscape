# Author: Jennifer Hoey
# Copyright (C) 2024 Jennifer Hoey
# License: GPL

## Dependencies ==============================================================
library(here)
library(ade4)
library(adegenet)
library(hierfstat)
library(poppr)
library(ape)
library(ggtree)
library(ggplot2)
library(viridis)
library(stringr)
library(RColorBrewer)
library(scales)
library(ggstance)
library(reshape2)
library(dplyr)
library(grid)
library(gtable)
library(ggfun)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### Read in datasets ####
m <- as.matrix(read.table(here("WGS/angsd_results/madracis606/unlinked_snps.ibsMat"), header = FALSE))

madracis_tree <- nj(m)

#### Load IDs directly ####
ids <- read.csv(here("WGS/angsd_results/bam_lists/ids_recoded606.csv"), header = TRUE) # 606 Madracis, no clones

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata $id_original <- paste(metadata[,4], metadata[,5], sep = "_")
metadata$depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
metadata$sample_full <- ids$new_ID

# Get cluster assignments for coloring the points
data_cluster <- read.csv(here("madracis_master_merged_May2024.csv"), header = TRUE)

metadata$species_mpha <- data_cluster[match(metadata$id_original, data_cluster$ID), 7]

# Replace NAs with matching in species column
species_indices <- which(is.na(metadata$species_mpha) == TRUE) # indicies that have NAs
metadata$species_mpha[is.na(metadata$species_mpha)] <- metadata$species[species_indices] # replace the NAs with species in associated index

# Split out MPHA by depth to assign MPHA.ln3
metadata$species_mpha[which(metadata$species == 'MPHA' & metadata$depth >60)] <- "MPHA.ln3" # MPHA at depths >60

# Deal with the MPHAs, which are mostly MSEN-like MPHAs
metadata[which(metadata$id_original == 'PB_4583'),9] <- 'MPHA.ln2'
metadata[which(metadata$id_original == 'PB_4578'),9] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4577'),9] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4576'),9] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4579'),9] <- 'MSEN-like MPHA'

metadata$species_mpha <- factor(metadata$species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MSEN-like MPHA', 'MSEN', 'UNKN')) #reorder the species factors

#### Add labels ####
madracis_tree$tip.label <- ids$new_ID # when ids are directly read in

#### Data manipulation ####
# Reroot M. senaria?
is.rooted (madracis_tree)
outgroup <- c("MSEN_est_20_JH_0336_E0061", "MSEN_hul_20_JH_0684_E0228", "MSEN_dir_20_JH_0635_E0597", "MSEN_est_10_JH_0285_E0448",
              "MSEN_est_10_JH_0289_E0311", "MSEN_est_10_JH_0383_E0679", "MSEN_est_10_JH_0398_E0375", "MSEN_est_10_JH_0503_E0133",
              "MSEN_est_10_JH_0514_E0480", "MSEN_est_20_JH_0489_E0232", "MSEN_est_5_JH_0344_E0608", "MSEN_hul_20_JH_0668_E0525",
              "MSEN_hul_20_JH_0685_E0116", "MSEN_kal_10_JH_0134_E0542", "MSEN_sna_10_JH_0001_E0610", "MSEN_sna_10_JH_0061_E0358",
              "MSEN_sna_10_JH_0070_E0611", "MSEN_sna_10_JH_0072_E0544", "MSEN_sna_20_JH_0697_E0005", "MSEN_wat_20_JH_0608_E0710",
              "MSEN_wat_20_JH_0630_E0084", "MSEN_wat_5_JH_0626_E0186", "MSEN_dir_20_JH_0638_E0699","MSEN_dir_20_JH_0641_E0258",
              "MSEN_dir_20_JH_0650_E0230", "MSEN_est_10_JH_0266_E0579", "MSEN_est_10_JH_0271_E0021", "MSEN_est_10_JH_0274_E0437",
              "MSEN_est_10_JH_0283_E0075", "MSEN_est_10_JH_0290_E0143", "MSEN_est_10_JH_0382_E0039", "MSEN_est_10_JH_0385_E0240",
              "MSEN_est_10_JH_0392_E0490", "MSEN_est_10_JH_0396_E0408", "MSEN_est_10_JH_0426_E0463", "MSEN_est_10_JH_0511_E0661",
              "MSEN_est_20_JH_0335_E0602", "MSEN_est_20_JH_0473_E0256", "MSEN_est_20_JH_0496_E0288", "MSEN_est_5_JH_0414_E0181" ,
              "MSEN_hul_20_JH_0691_E0343", "MSEN_kal_10_JH_0090_E0578", "MSEN_kal_10_JH_0125_E0445", "MSEN_kal_10_JH_0128_E0381",
              "MSEN_kal_10_JH_0190_E0328", "MSEN_kal_10_JH_0548_E0550", "MSEN_kal_20_JH_0159_E0334", "MSEN_kal_20_JH_0161_E0647",
              "MSEN_kal_20_JH_0165_E0656", "MSEN_kal_20_JH_0204_E0460", "MSEN_kal_20_JH_0208_E0008",
              "MSEN_kal_20_JH_0226_E0019", "MSEN_kal_20_JH_0227_E0567", "MSEN_kal_20_JH_0231_E0177", "MSEN_kal_20_JH_0232_E0675",
              "MSEN_kal_20_JH_0240_E0193", "MSEN_kal_20_JH_0245_E0663", "MSEN_kal_20_JH_0538_E0569", "MSEN_kal_5_JH_0564_E0669" ,
              "MSEN_sea_10_JH_0530_E0384", "MSEN_sea_10_JH_0536_E0371", "MSEN_sea_20_JH_0518_E0120", "MSEN_sea_20_JH_0519_E0022",
              "MSEN_sea_20_JH_0523_E0497", "MSEN_sea_20_JH_0524_E0600", "MSEN_sea_20_JH_0525_E0051", "MSEN_sea_20_JH_0528_E0549",
              "MSEN_sna_10_JH_0024_E0652", "MSEN_sna_10_JH_0030_E0461", "MSEN_sna_10_JH_0590_E0053", "MSEN_sna_20_JH_0040_E0503",
              "MSEN_sna_20_JH_0043_E0326", "MSEN_sna_20_JH_0045_E0199", "MSEN_sna_20_JH_0580_E0531", "MSEN_sna_20_JH_0586_E0391",
              "MSEN_sna_20_JH_0712_E0372", "MSEN_sna_5_JH_0601_E0329", "MSEN_wat_20_JH_0607_E0368","MSEN_wat_20_JH_0622_E0303",
              "MSEN_wat_20_JH_0623_E0103", "MSEN_wat_20_JH_0631_E0162", "MSEN_wat_5_JH_0627_E0351", "MSEN_est_10_JH_0311_E0564",
              "MSEN_est_10_JH_0314_E0547", "MSEN_est_10_JH_0316_E0218", "MSEN_est_10_JH_0327_E0624", "MSEN_est_20_JH_0472_E0670",
              "MSEN_est_40_JH_0439_E0654", "MSEN_kal_10_JH_0250_E0481", "MSEN_kal_10_JH_0259_E0118", "MSEN_kal_10_JH_0263_E0207",
              "MSEN_kal_20_JH_0155_E0449", "MSEN_kal_20_JH_0207_E0691", "MSEN_kal_20_JH_0233_E0213", "MSEN_sna_20_JH_0719_E0093",
              "MSEN_sna_20_JH_0720_E0320", "MSEN_est_20_JH_0307_E0090", "MSEN_est_20_JH_0462_E0399", "MSEN_kal_10_JH_0124_E0066",
              "MSEN_kal_20_JH_0213_E0017", "MSEN_sna_10_JH_0011_E0361", "MSEN_sna_10_JH_0012_E0470", "MSEN_sna_10_JH_0031_E0621",
              "MSEN_sna_20_JH_0709_E0100", "MSEN_est_20_JH_0356_E0616", "MSEN_est_20_JH_0363_E0468", "MSEN_kal_10_JH_0085_E0667",
              "MSEN_kal_10_JH_0550_E0349", "MSEN_sna_20_JH_0041_E0584", "MSEN_sna_20_JH_0052_E0357", "MPHA_kal_10_JH_0255_E0038",
              "MPHA_kal_20_JH_0196_E0115", "MPHA_hul_20_JH_0690_E0140", "MPHA_kal_20_JH_0224_E0119", "MDEC_sna_10_JH_0068_E0030",
              "MPHA_kal_10_JH_0123_E0006", "MPHA_kal_10_JH_0176_E0284") # MSEN_kal_20_JH_0170_E0206 grouping with MDEC

madracis_tree_reroot <- root(madracis_tree, outgroup, resolve.root = TRUE)
is.rooted(madracis_tree_reroot)

# Prepare tree-specific labels
species.table <- data.frame(sample_full = madracis_tree_reroot$tip.label, species_mpha = metadata$species_mpha, depth = metadata$depth, site = metadata$site)

#### Set color scheme ####
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70")
show_col(colors_mpha)

#### Produce individual plots ####
njtree <- ggtree(madracis_tree_reroot, layout="rectangular")

tree_reroot_names <- get_taxa_name(njtree) # get exact tip order of rerooted tree

species.table2 <- species.table[match(tree_reroot_names, species.table$sample_full),] # order data by new exact tip order

p1 <- njtree %<+% species.table2 +
  geom_tiplab(aes(color = factor(species_mpha)), align=TRUE, size=2, linesize=0.5, hjust = -0.05, show.legend = FALSE) + 
  geom_tippoint(aes(color = factor(species_mpha), shape = factor(species_mpha)), align=TRUE, size=2)  +
  scale_color_manual(values = colors_mpha, limits = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MSEN-like MPHA', 'MSEN', 'UNKN')) + scale_shape_manual(values = rep(19,10)) +  theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, 'cm')) + xlim_tree(0.3)

print(p1)
