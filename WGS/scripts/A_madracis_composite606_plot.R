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

# Cluster assignments (from nsgadmix K=5)
madracis_clusters_k5 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K5_1_out.qopt"), header = FALSE, sep = "", col.names = c("cluster5",  "cluster1", "cluster3", "cluster2", "cluster4"))

# Cluster assignments (from nsgadmix K=6)
madracis_clusters_k6 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K6_5_out.qopt"), header = FALSE, sep = "", col.names = c("cluster3", "cluster1", "cluster5", "cluster6", "cluster4", "cluster2"))

# Cluster assignments (from nsgadmix K=7)
madracis_clusters_k7 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K7_9_out.qopt"), header = FALSE, sep = "", col.names = c("cluster7", "cluster5", "cluster1", "cluster6", "cluster2", "cluster4", "cluster3"))

# Cluster assignments (from nsgadmix K=8)
madracis_clusters_k8 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K8_3_out.qopt"), header = FALSE, sep = "", col.names = c("cluster2", "cluster8", "cluster5", "cluster6", "cluster1", "cluster3", "cluster7", "cluster4"))

# # Cluster assignments (from nsgadmix K=9)
madracis_clusters_k9 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K9_9_out.qopt"), header = FALSE, sep = "", col.names = c("cluster2", "cluster5", "cluster6", "cluster8", "cluster4", "cluster7", "cluster9", "cluster3", "cluster1"))

# Cluster assignments (from nsgadmix K=10)
madracis_clusters_k10 <- read.delim(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K10_6_out.qopt"), header = FALSE, sep = "", col.names = c("cluster7", "cluster6", "cluster8", "cluster5", "cluster3", "cluster1", "cluster4", "cluster9", "cluster10", "cluster2"))

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

# madracis_clusters_k5 <- cbind(names4, madracis_clusters_k5)
# madracis_clusters_k6 <- cbind(names4, madracis_clusters_k6)
# madracis_clusters_k7 <- cbind(names4, madracis_clusters_k7)

madracis_clusters_k5 <- cbind(ids$new_ID, madracis_clusters_k5)
madracis_clusters_k6 <- cbind(ids$new_ID, madracis_clusters_k6)
madracis_clusters_k7 <- cbind(ids$new_ID, madracis_clusters_k7)
madracis_clusters_k8 <- cbind(ids$new_ID, madracis_clusters_k8)
madracis_clusters_k9 <- cbind(ids$new_ID, madracis_clusters_k9)
madracis_clusters_k10 <- cbind(ids$new_ID, madracis_clusters_k10)

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

# Root by deep MPHA 
# outgroup <- c("MPHA_sea_60_PB_4577_E0775", "MPHA_sea_60_PB_4576_E0774", "MPHA_sea_60_PB_4579_E0777", "MPHA_sna_40_JH_0034_E0138", "UNKN_est_40_JH_0431_E0168", "MPHA_est_20_JH_0332_E0083", "MPHA_sea_60_PB_4578_E0776", "UNKN_est_60_JH_0304_E0593")

madracis_tree_reroot <- root(madracis_tree, outgroup, resolve.root = TRUE)
is.rooted(madracis_tree_reroot)

# Prepare tree-specific labels
species.table <- data.frame(sample_full = madracis_tree_reroot$tip.label, species_mpha = metadata$species_mpha, depth = metadata$depth, site = metadata$site)

# Format ngsadmix data for ggplot
madracis_k5_melt <- melt(madracis_clusters_k5)
madracis_k6_melt <- melt(madracis_clusters_k6)
madracis_k7_melt <- melt(madracis_clusters_k7)
madracis_k8_melt <- melt(madracis_clusters_k8)
madracis_k9_melt <- melt(madracis_clusters_k9)
madracis_k10_melt <- melt(madracis_clusters_k10)

# Format DAPC data for ggplot
madracis_dapc_pc5_post_melt <- melt(madracis_dapc_pc5_post)

# Prepare site categories
species.table$site[which(species.table$site == 'buy')] <- 1 #buoy1
species.table$site[which(species.table$site == 'dir')] <- 2 #directors bay
species.table$site[which(species.table$site == 'est')] <- 3 #coral estate
species.table$site[which(species.table$site == 'hul')] <- 4 #hulu
species.table$site[which(species.table$site == 'kal')] <- 5 #kalki
species.table$site[which(species.table$site == 'sea')] <- 6 #seaquarium
species.table$site[which(species.table$site == 'sna')] <- 7 #snake bay
species.table$site[which(species.table$site == 'wat')] <- 8 #water factor

# Prepare depth categories
species.table$depth[which(species.table$depth <=10)] <- 1 # depths less than 10
species.table$depth[which(species.table$depth >11 & species.table$depth <=20)] <- 2 # depths between 11-20
species.table$depth[which(species.table$depth >21 & species.table$depth <=40)] <- 3 # depths between 21-40
species.table$depth[which(species.table$depth >41 & species.table$depth <=60)] <- 4 # depths between 41-60
species.table$depth[which(species.table$depth >61)] <- 5 # depths greater than 61

#### Set color scheme ####
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

# For STRUCTURE plots
# colors_struc <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#0D0887FF", "#5402A3FF", "lightblue")
colors_struc <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#0D0887FF", "#5402A3FF", "lightblue", "#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70")
# c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", plasma(8))
show_col(colors_struc)

#### Produce individual plots ####
njtree <- ggtree(madracis_tree_reroot, layout="rectangular")
# njtree <- ggtree(madracis_tree_reroot, layout="rectangular", branch.length = 'none')

tree_reroot_names <- get_taxa_name(njtree) # get exact tip order of rerooted tree

species.table2 <- species.table[match(tree_reroot_names, species.table$sample_full),] # order data by new exact tip order
# cluster.table <- madracis653_assigns_k8[match(tree_reroot_names, madracis653_assigns_k8$ID),] # order species table by new exact tip order
# cluster.table$max.cluster.num <- gsub("cluster", "", cluster.table$max.cluster)# remove 'cluster' from max.cluster

p1 <- njtree %<+% species.table2 +
      geom_tiplab(aes(color = factor(species_mpha)), align=TRUE, size=2, linesize=0.5, hjust = -0.05, show.legend = FALSE) + 
      geom_tippoint(aes(color = factor(species_mpha), shape = factor(species_mpha)), align=TRUE, size=2)  +
      scale_color_manual(values = colors_mpha, limits = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MSEN-like MPHA', 'MSEN', 'UNKN')) + scale_shape_manual(values = rep(19,10)) +  theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, 'cm')) + xlim_tree(0.3)

p2 <- facet_plot(p1, panel= "z.structure a5", data=madracis_k5_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity') +
  scale_fill_manual(values = colors_struc) + theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.5, 'cm')) # for whatever reason, the facet plots are being plotted in alphabetical order regardless of the order they are called in, so need to add letters to get all subsequent facets to plot afterwards and then relabel at the end

# # p2.2 <- p2 + facet_set(label = c('Tree'= '1.tree', 'z.structure k5'= '2.structure k5')) # thought I needed this to get the ordering of the facets right, but I don't think I do

p3 <- facet_plot(p2, panel="z.structure b6", data=madracis_k6_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p4 <- facet_plot(p3, panel="z.structure c7", data=madracis_k7_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p5 <- facet_plot(p4, panel="z.structure d8", data=madracis_k8_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p6 <- facet_plot(p5, panel="z.structure e9", data=madracis_k9_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p7 <- facet_plot(p6, panel="z.structure f10", data=madracis_k10_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p8 <- facet_plot(p7, panel="zz.dapc", data=madracis_dapc_pc5_post_melt,
                geom=geom_barh, width = 1.0,
                mapping = aes(x = value, fill = as.factor(variable)),
                stat='identity')

p9 <- facet_plot(p8, panel="zz.site", data=species.table2,
                 geom = geom_point, aes(x = rev(as.numeric(species.table2$site)), color = rep(NA, 606)), shape = 15, size=2, show.legend = F) # use dataframe reordered by rerooted tree

p10 <- facet_plot(p9, panel="zzz.<10, 11-20, 21-40, 41-60, 60+", data=species.table2,
           geom = geom_point, aes(x = rev(as.numeric(species.table2$depth)), color = rep(NA,606)),shape = 7, size=2, show.legend = F) # use dataframe reordered by rerooted tree

p11 <- facet_labeller(p10, label = c('Tree' = 'tree', 'z.structure a5'= 'k5', 'z.structure b6'= 'k6', 'z.structure c7'= 'k7', 'z.structure d8'= 'k8', 'z.structure e9'= 'k9', 'z.structure f10'= 'k10', 'zz.dapc' = 'DAPC', 'zz.site' = 'Site', 'zzz.<10, 11-20, 21-40, 41-60, 60+'= 'Depth'))

# Change the width of the individual facet_plot panels
gt <- ggplot_gtable(ggplot_build(p11))
# gt$widths[5] <- 0.3*gt$widths[5] #tree
gt$widths[7] <- 0.3*gt$widths[7] #k5
gt$widths[9] <- 0.08*gt$widths[9]
gt$widths[11] <- 0.08*gt$widths[11]
gt$widths[13] <- 0.08*gt$widths[13] 
gt$widths[15] <- 0.08*gt$widths[15]
gt$widths[17] <- 0.08*gt$widths[17] #k10
gt$widths[19] <- 0.08*gt$widths[19] #dapc?
gt$widths[21] <- 0.08*gt$widths[21] #site/dapc?
gt$widths[23] <- 0.2*gt$widths[23] #site/depth?
gt$widths[25] <- 0.2*gt$widths[25] #depth?
grid.draw(gt)

ggsave(gt,file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/figures/madracis606_composite.png", width = 38, height = 130, units = "cm", limitsize = FALSE)












# p2 <- p1 %>% collapse(node=1066) + 
#   geom_point2(aes(subset=(node==1066)), shape=21, size=3, fill='#374618') # very deep MPHA
# p3 <- collapse(p2, node=1164) + 
#   geom_point2(aes(subset=(node==1164)), shape=21, size=3, fill='#ff6347') # MCAR
# p4 <- collapse(p3, node=658) + 
#   geom_point2(aes(subset=(node==658)), shape=21, size=3, fill='#b6c590') # shallow MPHA
# p5 <- collapse(p4, node=1005) + 
#   geom_point2(aes(subset=(node==1005)), shape=21, size=3, fill='#6b8e23') # deep MPHA
# p6 <- collapse(p5, node=1075) + 
#   geom_point2(aes(subset=(node==1075)), shape=21, size=3, fill='#E7298A') # MFOR
# p7 <- collapse(p6, node=821) + 
#   geom_point2(aes(subset=(node==821)), shape=21, size=3, fill='#E6AB02') # MMIR
# p8 <- collapse(p7, node=746) + 
#   geom_point2(aes(subset=(node==746)), shape=21, size=3, fill='#ff9a81') # MDEC
# p9 <- collapse(p8, node=1169) + 
#   geom_point2(aes(subset=(node==1169)), shape=21, size=3, fill='#b0e0e6') # MSEN
# p10 <- collapse(p9, node=1176) + 
#   geom_point2(aes(subset=(node==1176)), shape=21, size=3, fill='#6b8e23') # MSEN like MPHA







p7 <- facet_plot(p6, panel="z.structure f11", data=madracis_k11_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p8 <- facet_plot(p7, panel="z.structure g12", data=madracis_k12_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p9 <- facet_plot(p8, panel="z.structure h13", data=madracis_k13_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

# p10 <- facet_plot(p9, panel="z.structure i14", data=madracis_k14_melt,
#                  geom=geom_barh, width = 1.0,
#                  mapping = aes(x = value, fill = as.factor(variable)),
#                  stat='identity')
# 
# p11 <- facet_plot(p10, panel="z.structure j15", data=madracis_k15_melt,
#                   geom=geom_barh, width = 1.0,
#                   mapping = aes(x = value, fill = as.factor(variable)),
#                   stat='identity')

# p6 <- facet_plot(p5, panel="zz.assignment k7", data=madracis_dapc_k7_melt,
#                 geom=geom_barh, width = 1.0,
#                 mapping = aes(x = V1, fill = as.factor(V3)),
#                 stat='identity')

# #p7 <- facet_plot(p6, panel="zz.dapc k7", data=madracis_dapc_k7_melt, 
# #                 geom=geom_barh, width = 1.0,
# #                 mapping = aes(x = value, fill = as.factor(Var2)), 
# #                 stat='identity')

# p6 <- facet_plot(p5, panel="zz.assign k8", data=cluster.table,
#                  geom = geom_point, aes(x = rev(as.numeric(cluster.table$max.cluster.num)), color = rep(NA, 653)), shape = 15, size=2, show.legend = F)

# p8 <- facet_plot(p5, panel="zzz.site", data=species.table.merged2,
#                  geom = geom_point, aes(x = rev(as.numeric(species.table.merged2$site)), color = rep(NA, 87)), shape = 15, size=2, show.legend = F) # use dataframe reordered by rerooted tree
# 
# p9 <- facet_plot(p8, panel="zzzz.<10, 11-20, 21-40, 41-60, 60+", data=species.table.merged2,
#            geom = geom_point, aes(x = rev(as.numeric(species.table.merged2$depth)), color = rep(NA,87)),shape = 15, size=2, show.legend = F) # use dataframe reordered by rerooted tree
# 
#  
p10 <- facet_labeller(p9, label = c('Tree' = 'tree', 'z.structure a6'= 'k6', 'z.structure b7'= 'k7', 'z.structure c8'= 'k8', 'z.structure d9'= 'k9', 'z.structure e10'= 'k10', 'z.structure f11'= 'k11', 'z.structure g12'= 'k12', 'z.structure h13'= 'k13'))
                                    
# 'zzz.site' = 'site', 'zzzz.<10, 11-20, 21-40, 41-60, 60+' = 'depth')) #relabels tree facets once the facets are in the desired order



