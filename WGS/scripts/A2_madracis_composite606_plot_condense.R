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
library(ggnewscale)

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

# Cluster assignments (from nsgadmix K=6)
madracis_clusters_k6 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K6_5_out.qopt"), header = FALSE, sep = "", col.names = c("cluster3", "cluster1", "cluster5", "cluster6", "cluster4", "cluster2"))

# Cluster assignments (from nsgadmix K=7)
# madracis_clusters_k7 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K7_9_out.qopt"), header = FALSE, sep = "", col.names = c("cluster7", "cluster5", "cluster1", "cluster6", "cluster2", "cluster4", "cluster3"))
madracis_clusters_k7 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K7_2_out.qopt"), header = FALSE, sep = "", col.names = c("cluster7", "cluster4", "cluster1", "cluster3", "cluster6", "cluster5", "cluster2"))

# Reorder K7 columns to match K6
madracis_clusters_k7_ordered <- madracis_clusters_k7[,c(1,3:4,6,5,7,2)]

# Cluster assignments (from nsgadmix K=8)
madracis_clusters_k8 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K8_3_out.qopt"), header = FALSE, sep = "", col.names = c("cluster2", "cluster7", "cluster5", "cluster6", "cluster1", "cluster3", "cluster8", "cluster4"))

# Reorder K8 columns to match K6-7
madracis_clusters_k8_ordered <- madracis_clusters_k8[,c(5,2:4,6:8,1)]

# DAPC results
madracis_dapc_pc5_post <- read.delim(here('WGS/angsd_results/madracis606/madracis606_dapc_5pcs_posteriors.txt'), header = TRUE, sep = '\t')

#### Option 1: Load bam names to assign species ####
# bam_paths <- read.table(here("angsd_results/bam_lists/angsd_bam_list_653.txt"), header = FALSE, sep = '\t', col.names = c('path'))
# 
# # Split file paths to extract metadata
# names <- str_split_fixed(bam_paths$path, '.dedup.overlapclipped.bam', 2) # Get rid of suffix
# names2 <- str_split_fixed(names[,1], '/', 8)
# names3 <- sub('_tech_', '_', names2[,8]) # remove 'tech' prefix of technical replicates so that str_split_fixed will work well
# metadata <- as.data.frame(str_split_fixed(names3, '_', 7))
# colnames(metadata) <- c('merged', 'species', 'site', 'depth', 'collector', 'ID', 'extract')
# 
# species <- metadata[,2]
# site <- metadata[,3]
# depth <- as.numeric(metadata[,4]) # Must be numeric or the subsetting gets very weird
# id_original <- paste(metadata[,5], metadata[,6], sep = "_")
# extract <- metadata[,7]
# 
# names4 <- paste(species, site, depth, id_original, extract, sep = "_")

#### Option 2: Load IDs directly ####
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

# Deal with the Pim MPHAs, which are mostly MSEN-like MPHAs (MPHA.ln4)
metadata[which(metadata$id_original == 'PB_4583'),9] <- 'MPHA.ln2'
metadata[which(metadata$id_original == 'PB_4578'),9] <- 'MPHA.ln4'
metadata[which(metadata$id_original == 'PB_4577'),9] <- 'MPHA.ln4'
metadata[which(metadata$id_original == 'PB_4576'),9] <- 'MPHA.ln4'
metadata[which(metadata$id_original == 'PB_4579'),9] <- 'MPHA.ln4'

metadata$species_mpha <- factor(metadata$species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MPHA.ln4', 'MSEN', 'UNKN')) #reorder the species factors

#### Add labels ####
madracis_tree$tip.label <- ids$new_ID # when ids are directly read in

# madracis_clusters_k5 <- cbind(names4, madracis_clusters_k5)
# madracis_clusters_k6 <- cbind(names4, madracis_clusters_k6)
# madracis_clusters_k7 <- cbind(names4, madracis_clusters_k7)

madracis_clusters_k6 <- cbind(ids$new_ID, madracis_clusters_k6)
madracis_clusters_k7 <- cbind(ids$new_ID, madracis_clusters_k7_ordered)
madracis_clusters_k8 <- cbind(ids$new_ID, madracis_clusters_k8_ordered)

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
              "MSEN_kal_10_JH_0550_E0349", "MSEN_sna_20_JH_0041_E0584", "MSEN_sna_20_JH_0052_E0357", "MSEN_kal_10_JH_0255_E0038",
              "MSEN_kal_20_JH_0196_E0115", "MSEN_hul_20_JH_0690_E0140", "MSEN_kal_20_JH_0224_E0119", "MSEN_sna_10_JH_0068_E0030",
              "MSEN_kal_10_JH_0123_E0006", "MSEN_kal_10_JH_0176_E0284")

madracis_tree_reroot <- root(madracis_tree, outgroup, resolve.root = TRUE)
is.rooted(madracis_tree_reroot)

# Prepare tree-specific labels
species.table <- data.frame(sample_full = madracis_tree_reroot$tip.label, species_mpha = metadata$species_mpha, depth = metadata$depth, site = metadata$site)

# Format ngsadmix data for ggplot
madracis_k6_melt <- melt(madracis_clusters_k6)
madracis_k7_melt <- melt(madracis_clusters_k7)
madracis_k8_melt <- melt(madracis_clusters_k8)

# Format DAPC data for ggplot
madracis_dapc_pc5_post_melt <- melt(madracis_dapc_pc5_post)

# # Prepare site categories
# species.table$site[which(species.table$site == 'buy')] <- 1 #buoy1
# species.table$site[which(species.table$site == 'dir')] <- 2 #directors bay
# species.table$site[which(species.table$site == 'est')] <- 3 #coral estate
# species.table$site[which(species.table$site == 'hul')] <- 4 #hulu
# species.table$site[which(species.table$site == 'kal')] <- 5 #kalki
# species.table$site[which(species.table$site == 'sea')] <- 6 #seaquarium
# species.table$site[which(species.table$site == 'sna')] <- 7 #snake bay
# species.table$site[which(species.table$site == 'wat')] <- 8 #water factor
# 
# # Prepare depth categories
# species.table$depth[which(species.table$depth <=10)] <- 1 # depths less than 10
# species.table$depth[which(species.table$depth >11 & species.table$depth <=20)] <- 2 # depths between 11-20
# species.table$depth[which(species.table$depth >21 & species.table$depth <=40)] <- 3 # depths between 21-40
# species.table$depth[which(species.table$depth >41 & species.table$depth <=60)] <- 4 # depths between 41-60
# species.table$depth[which(species.table$depth >61)] <- 5 # depths greater than 61

#### Set color scheme ####
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

# For STRUCTURE plots
colors_struc <- c("#E6AB02", "#6b8e23", "#ff9a81", "gray70", "#b0e0e6", "#b6c590", "#79d2a6", "#E7298A", "#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70")
show_col(colors_struc)

colors_heatmap <- c("#DCE318FF", "#B8DE29FF", "#74D055FF", "#20A386FF", "#238A8DFF", "#2D718EFF", "#453781FF", "#481568FF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
show_col(colors_heatmap)

colors_all <- c(colors_mpha, colors_heatmap)

#### Produce individual plots ####
njtree <- ggtree(madracis_tree_reroot, layout="rectangular")

tree_reroot_names <- get_taxa_name(njtree) # get exact tip order of rerooted tree

species.table2 <- species.table[match(tree_reroot_names, species.table$sample_full),] # order data by new exact tip order

# cluster.table <- madracis653_assigns_k8[match(tree_reroot_names, madracis653_assigns_k8$ID),] # order species table by new exact tip order
# cluster.table$max.cluster.num <- gsub("cluster", "", cluster.table$max.cluster)# remove 'cluster' from max.cluster

# Tree for composite figure
# p1 <- njtree %<+% species.table2 +
#       geom_tiplab(aes(color = factor(species_mpha)), align=TRUE, size=2, linesize=0.5, hjust = -0.05, show.legend = FALSE) + 
#       geom_tippoint(aes(color = factor(species_mpha), shape = factor(species_mpha)), align=TRUE, size=2)  +
#       scale_color_manual(values = colors_mpha, limits = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MPHA.ln4', 'MSEN', 'UNKN')) + scale_shape_manual(values = rep(19,10)) +  theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, 'cm')) + xlim_tree(0.3)

# njtree + geom_nodelab(aes(label = node), hjust = -0.5) # to view node labels

p1 <- njtree %<+% species.table2 + 
  geom_tippoint(aes(color = factor(species_mpha), shape = factor(species_mpha)), align=TRUE, size=1)  +
  scale_color_manual(values = colors_mpha, limits = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MPHA.ln4', 'MSEN', 'UNKN')) + scale_shape_manual(values = rep(19,10)) +  theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, 'cm'), strip.text.x = element_blank())

p1.1 <- rotate(p1, 926)
p1.2 <- rotate(p1.1, 1047)

# Because I reordered the tree again, I need to find the new tip order
tree_reroot_names2 <- get_taxa_name(p1.2) # get exact tip order of reordered tree

species.table3 <- species.table2[match(tree_reroot_names2, species.table2$sample_full),] # order data by new exact tip order

species.table3$depth <- species.table3$depth + 10 # turning sites into numbers is causing a conflict between 5 as depth and 5 as a site number, so I have to be clever about how I deal with this

heatmap_data <- cbind.data.frame(species.table3$sample_full, species.table3$site, species.table3$depth)
colnames(heatmap_data) <- c('sample_full', 'site', 'depth')
rownames(heatmap_data) <- species.table3$sample_full
# heatmap_data$site <- as.numeric(factor(heatmap_data$site)) # for Method 1 only
heatmap_data$depth <- as.factor(heatmap_data$depth) # for Method 2 only

p2 <- p1.2 + new_scale_fill()

#### Two different ways of visualizing the metadata. Method 1 consists of a column for every possible category. Method 2 is a single column with categories represented by different colors ####
# Method 1
p2.a <-  p2 + 
  geom_facet(
    panel = "site", 
    data = heatmap_data, 
    geom = geom_point, 
    mapping = aes(x = rev(heatmap_data$site), color = rev(factor(heatmap_data$site))),
    stat='identity'
  ) + 
  scale_color_manual(values=c(colors_all)) 

p2.b <-  p2.a + 
  geom_facet(
    panel = "depth", 
    data = heatmap_data, 
    geom = geom_point, 
    mapping = aes(x = rev(heatmap_data$depth), color = rev(factor(heatmap_data$depth))),
    stat='identity'
  ) + 
  scale_color_manual(values=c(colors_all)) 

# Method 2
p2.1 <- gheatmap(p2, heatmap_data[,-1], offset=0, width=0.5, colnames = FALSE) +
        scale_fill_manual(values = c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B", "#DCE318FF", "#B8DE29FF", "#74D055FF", "#20A386FF", "#238A8DFF", "#2D718EFF", "#453781FF", "#481568FF")) # cannot have the samples names, for whichever reason

#### Plot the rest of the genetic data ####
p3 <- p2.1 + new_scale_fill() # depending on metadata visualization method, either p2.1 or p2.b
p3.1 <- facet_plot(p3, panel= "z.structure a6", data=madracis_k6_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity') +
  scale_fill_manual(values = colors_struc) + theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.5, 'cm'), strip.text.x = element_blank()) # for whatever reason, the facet plots are being plotted in alphabetical order regardless of the order they are called in, so need to add letters to get all subsequent facets to plot afterwards and then relabel at the end

# # p2.2 <- p2 + facet_set(label = c('Tree'= '1.tree', 'z.structure k5'= '2.structure k5')) # thought I needed this to get the ordering of the facets right, but I don't think I do

p4 <- facet_plot(p3.1, panel="z.structure b7", data=madracis_k7_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p5 <- facet_plot(p4, panel="z.structure c8", data=madracis_k8_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p6 <- facet_plot(p5, panel="zz.dapc", data=madracis_dapc_pc5_post_melt,
                geom=geom_barh, width = 1.0,
                mapping = aes(x = value, fill = as.factor(variable)),
                stat='identity')

# p6 <- facet_plot(p5, panel="zz.site", data=species.table2,
#                  geom = geom_point, aes(x = rev(as.numeric(species.table2$site)), color = rep(NA, 606)), shape = 15, size=2, show.legend = F) # use dataframe reordered by rerooted tree
# 
# p7 <- facet_plot(p6, panel="zzz.<10, 11-20, 21-40, 41-60, 60+", data=species.table2,
#            geom = geom_point, aes(x = rev(as.numeric(species.table2$depth)), color = rep(NA,606)), shape = 7, size=2, show.legend = F) # use dataframe reordered by rerooted tree

# p7 <- facet_labeller(p6, label = c('Tree' = 'tree', 'z.structure a6'= 'k6', 'z.structure b7'= 'k7', 'z.structure c8'= 'k8', 'zz.dapc' = 'DAPC'))

# Change the width of the individual facet_plot panels
gt <- ggplot_gtable(ggplot_build(p6))
gt$widths[7] <- 0.3*gt$widths[7] #tree (0.5 vs 0.3)
gt$widths[9] <- 0.12*gt$widths[9]
gt$widths[11] <- 0.12*gt$widths[11]
gt$widths[13] <- 0.1*gt$widths[13] 
gt$widths[15] <- 0.1*gt$widths[15]
# With method 1 for metadata
gt$widths[17] <- 0.1*gt$widths[17]
gt$widths[19] <- 0.1*gt$widths[19]

grid.draw(gt)

ggsave(gt,file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/figures/madracis606_composite_condensed_alt.png", width = 23.5, height = 28, units = "cm", limitsize = FALSE)

