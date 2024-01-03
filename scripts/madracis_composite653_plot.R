# Author: Jennifer Hoey
# Copyright (C) 2023 Jennifer Hoey
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
m <- as.matrix(read.table(here("angsd_results/coral653_snps.ibsMat"), header = FALSE))

madracis_tree <- nj(m)

# # Cluster assignments (from nsgadmix K=5)
# madracis_clusters_k5 <- read.delim(here("admixture_results/qfiles/K5/madracis87_ngsAdmix_K5_1_out.qopt"), header = FALSE, sep = " ", col.names = c("cluster1",  "cluster2", "cluster3", "cluster4", "cluster5", "extra"))
# madracis_clusters_k5 <- madracis_clusters_k5[,-6]
# 
# Cluster assignments (from nsgadmix K=6)
madracis_clusters_k6 <- read.delim(here("admixture_results/madracis653/qfiles/K6/madracis653_ngsAdmix_K6_2_out.qopt"), header = FALSE, sep = " ", col.names = c("cluster2", "cluster1", "cluster6", "cluster3", "cluster5", "cluster4", "extra"))
madracis_clusters_k6 <- madracis_clusters_k6[,-7]

# Cluster assignments (from nsgadmix K=7)
madracis_clusters_k7 <- read.delim(here("admixture_results/madracis653/qfiles/K7/madracis653_ngsAdmix_K7_2_out.qopt"), header = FALSE, sep = " ", col.names = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster7", "cluster6", "extra"))
madracis_clusters_k7 <- madracis_clusters_k7[,-8]

# Cluster assignments (from nsgadmix K=8)
madracis_clusters_k8_big <- read.delim(here("admixture_results/madracis653/qfiles/K8/madracis653_ngsAdmix_K8_1_out.qopt"), header = FALSE, sep = " ", col.names = c("cluster4", "cluster8", "cluster6", "cluster5", "cluster7", "cluster2", "cluster1", "cluster3", "extra"))
madracis_clusters_k8_big <- madracis_clusters_k8_big[,-9]

madracis_clusters_k8 <- read.delim(here("admixture_results/madracis653/qfiles/K8/madracis653_ngsAdmix_K8_2_out.qopt"), header = FALSE, sep = " ", col.names = c("cluster8", "cluster1", "cluster2", "cluster5", "cluster6", "cluster7", "cluster4", "cluster3", "extra"))
madracis_clusters_k8 <- madracis_clusters_k8[,-9]

# Cluster assignments (from nsgadmix K=9)
madracis_clusters_k9 <- read.delim(here("admixture_results/madracis653/qfiles/K9/madracis653_ngsAdmix_K9_2_out.qopt"), header = FALSE, sep = " ", col.names = c("cluster8", "cluster6", "cluster9", "cluster5", "cluster1", "cluster7", "cluster3", "cluster4", "cluster2", "extra"))
madracis_clusters_k9 <- madracis_clusters_k9[,-10]

# DAPC results
# madracis_dapc <- read.csv(here("data/dapc_results.csv"), header = FALSE, col.names = c("sample_full", "sample", "sample2", "species", "dapc_unsup"))
# load(here("results/dapc_k5.RData")) # K5
# load(here("results/dapc_k6.RData")) # K6
# load(here("results/dapc_k7.RData")) # K7

# Modify column names so that colors match later
# colnames(dapc1$posterior) <- c("1", "2", "3", "4", "5")
# colnames(dapc2$posterior) <- c("4", "2", "6", "1", "3", "5")
# colnames(dapc3$posterior) <- c("3", "2", "4", "7", "6", "5", "1")

# Symbiont data


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
ids <- read.csv(here("angsd_results/bam_lists/ids_recoded.csv"), header = TRUE)

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')

species <- metadata[,1]
site <- metadata[,2]
depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
id_original <- paste(metadata[,4], metadata[,5], sep = "_")
extract <- metadata[,6]

#### Organize MPHA ####                        
# Split out MPHA by depth
species_mpha <- species
species_mpha[which(species == 'MPHA' & depth <=20)] <- "MPHA-shallow" # MPHA at depths of 20m or shallower
species_mpha[which(species == 'MPHA' & depth >20 & depth <=60)] <- "MPHA-deep" 
species_mpha[which(species == 'MPHA' & depth >60)] <- "MPHA-very deep" # MPHA at depths >60

species_mpha <- factor(species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'SINT', 'UNKN')) #reorder the species factors

#### Add labels ####
# madracis_tree$tip.label <- names4 # when bamlist is used
madracis_tree$tip.label <- ids$new_ID # when ids are directly read in

# madracis_clusters_k5 <- cbind(names4, madracis_clusters_k5)
# madracis_clusters_k6 <- cbind(names4, madracis_clusters_k6)
# madracis_clusters_k7 <- cbind(names4, madracis_clusters_k7)

madracis_clusters_k6 <- cbind(ids$new_ID, madracis_clusters_k6)
madracis_clusters_k7 <- cbind(ids$new_ID, madracis_clusters_k7)
madracis_clusters_k8_big <- cbind(ids$new_ID, madracis_clusters_k8_big)
madracis_clusters_k8 <- cbind(ids$new_ID, madracis_clusters_k8)
madracis_clusters_k9 <- cbind(ids$new_ID, madracis_clusters_k9)

#### Data manipulation ####
# Reroot M. senaria
is.rooted (madracis_tree)
# outgroup <- c("MCAR_est_20_JH_0336_E0061", "MDEC_hul_20_JH_0684_E0228", "MDEC_kal_10_JH_0089_E0485", "MFOR_sna_10_JH_0068_E0030", "MPHA_dir_20_JH_0635_E0597", "MPHA_est_10_JH_0285_E0448", "MPHA_est_10_JH_0289_E0311", "MPHA_est_10_JH_0383_E0679", "MPHA_est_10_JH_0398_E0375", "MPHA_est_10_JH_0503_E0133", "MPHA_est_10_JH_0514_E0480", "MPHA_est_20_JH_0462_E0399", "MPHA_est_20_JH_0489_E0232", "MPHA_est_5_JH_0344_E0608",  "MPHA_hul_20_JH_0668_E0525", "MPHA_hul_20_JH_0685_E0116", "MPHA_hul_20_JH_0690_E0140", "MPHA_kal_10_JH_0134_E0542", "MPHA_kal_20_JH_0196_E0115", "MPHA_kal_20_JH_0224_E0119", "MPHA_sna_10_JH_0001_E0610", "MPHA_sna_10_JH_0011_E0361", "MPHA_sna_10_JH_0012_E0470", "MPHA_sna_10_JH_0061_E0358", "MPHA_sna_10_JH_0070_E0611", "MPHA_sna_10_JH_0071_E0089", "MPHA_sna_10_JH_0072_E0544", "MPHA_sna_20_JH_0697_E0005", "MPHA_sna_20_JH_0707_E0201", "MPHA_wat_20_JH_0608_E0710", "MPHA_wat_20_JH_0630_E0084", "MPHA_wat_5_JH_0626_E0186",  "MPHA_est_20_JH_0356_E0616", "MSEN_dir_20_JH_0638_E0699", "MSEN_dir_20_JH_0641_E0258", "MSEN_dir_20_JH_0650_E0230", "MSEN_est_10_JH_0266_E0579", "MSEN_est_10_JH_0271_E0021", "MSEN_est_10_JH_0274_E0437", "MSEN_est_10_JH_0283_E0075", "MSEN_est_10_JH_0290_E0143", "MSEN_est_10_JH_0382_E0039", "MSEN_est_10_JH_0385_E0240", "MSEN_est_10_JH_0392_E0490", "MSEN_est_10_JH_0396_E0408", "MSEN_est_10_JH_0426_E0463", "MSEN_est_10_JH_0511_E0661", "MSEN_est_20_JH_0335_E0602", "MSEN_est_20_JH_0473_E0256", "MSEN_est_20_JH_0496_E0288", "MSEN_est_5_JH_0414_E0181", "MSEN_hul_20_JH_0691_E0343", "MSEN_kal_10_JH_0090_E0578", "MSEN_kal_10_JH_0122_E0515", "MSEN_kal_10_JH_0124_E0066", "MSEN_kal_10_JH_0125_E0445", "MSEN_kal_10_JH_0128_E0381", "MSEN_kal_10_JH_0190_E0328", "MSEN_kal_10_JH_0548_E0550", "MSEN_kal_10_JH_0550_E0349", "MSEN_kal_20_JH_0159_E0334", "MSEN_kal_20_JH_0161_E0647", "MSEN_kal_20_JH_0165_E0656", "MSEN_kal_20_JH_0204_E0460", "MSEN_kal_20_JH_0208_E0008", "MSEN_kal_20_JH_0226_E0019", "MSEN_kal_20_JH_0227_E0567", "MSEN_kal_20_JH_0231_E0177", "MSEN_kal_20_JH_0232_E0675", "MSEN_kal_20_JH_0240_E0193", "MSEN_kal_20_JH_0245_E0663", "MSEN_kal_20_JH_0538_E0569", "MSEN_kal_5_JH_0564_E0669",  "MSEN_sea_10_JH_0530_E0384", "MSEN_sea_10_JH_0536_E0371", "MSEN_sea_20_JH_0518_E0120", "MSEN_sea_20_JH_0519_E0022", "MSEN_sea_20_JH_0523_E0497", "MSEN_sea_20_JH_0524_E0600", "MSEN_sea_20_JH_0525_E0051", "MSEN_sea_20_JH_0528_E0549", "MSEN_sna_10_JH_0024_E0652", "MSEN_sna_10_JH_0030_E0461", "MSEN_sna_10_JH_0031_E0621", "MSEN_sna_10_JH_0590_E0053", "MSEN_sna_20_JH_0040_E0503", "MSEN_sna_20_JH_0043_E0326", "MSEN_sna_20_JH_0045_E0199", "MSEN_sna_20_JH_0580_E0531", "MSEN_sna_20_JH_0586_E0391", "MSEN_sna_20_JH_0709_E0100", "MSEN_sna_20_JH_0712_E0372", "MSEN_sna_5_JH_0601_E0329",  "MSEN_wat_20_JH_0607_E0368", "MSEN_wat_20_JH_0622_E0303", "MSEN_wat_20_JH_0623_E0103", "MSEN_wat_20_JH_0631_E0162", "MSEN_wat_5_JH_0627_E0351", "MSEN_est_20_JH_0363_E0468", "MSEN_kal_10_JH_0085_E0667", "MSEN_sna_20_JH_0041_E0584", "MSEN_sna_20_JH_0052_E0357", "SINT_sna_20_JH_0013_E0713", "UNKN_est_10_JH_0311_E0564", "UNKN_est_10_JH_0314_E0547", "UNKN_est_10_JH_0316_E0218", "UNKN_est_10_JH_0327_E0624", "UNKN_est_20_JH_0306_E0261", "UNKN_est_20_JH_0307_E0090", "UNKN_est_20_JH_0463_E0322", "UNKN_est_20_JH_0464_E0572", "UNKN_est_20_JH_0472_E0670", "UNKN_est_40_JH_0439_E0654", "UNKN_kal_10_JH_0250_E0481", "UNKN_kal_10_JH_0255_E0038", "UNKN_kal_10_JH_0259_E0118", "UNKN_kal_10_JH_0263_E0207", "UNKN_kal_20_JH_0155_E0449", "UNKN_kal_20_JH_0207_E0691", "UNKN_kal_20_JH_0212_E0533", "UNKN_kal_20_JH_0213_E0017", "UNKN_kal_20_JH_0233_E0213", "UNKN_sna_20_JH_0719_E0093", "UNKN_sna_20_JH_0720_E0320")

outgroup <- c("SINT_kal_40_JH_0105_E0657", "SINT_kal_40_JH_0147_E0236", "SINT_sna_60_JH_0075_E0215") # seem like SINT

madracis_tree_reroot <- root(madracis_tree, outgroup, resolve.root = TRUE)
is.rooted(madracis_tree_reroot)

# Prepare tree-specific labels
species.table <- data.frame(sample_full = madracis_tree$tip.label, species = species, species_mpha = species_mpha, depth = depth, site = site)

# Format ngsadmix data for ggplot
# madracis_k5_melt <- melt(madracis_clusters_k5)
madracis_k6_melt <- melt(madracis_clusters_k6)
madracis_k7_melt <- melt(madracis_clusters_k7)
madracis_k8_melt_big <- melt(madracis_clusters_k8_big)
madracis_k8_melt <- melt(madracis_clusters_k8)
madracis_k9_melt <- melt(madracis_clusters_k9)

# Format DAPC data for ggplot
# madracis_dapc_k5_melt <- melt(dapc1$posterior)
# madracis_dapc_k6_melt <- melt(dapc2$posterior)

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
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#b0e0e6", "#2664d8", "gray70") # new colors
show_col(colors_mpha)

# For STRUCTURE plots
colors7 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", plasma(7))
show_col(colors7)

#### Produce individual plots ####
njtree <- ggtree(madracis_tree_reroot, layout="rectangular")
# njtree <- ggtree(madracis_tree_reroot, layout="rectangular", branch.length = 'none')

tree_reroot_names <- get_taxa_name(njtree) # get exact tip order of rerooted tree

species.table <- species.table[match(tree_reroot_names, species.table$sample_full),] # order species table by new exact tip order
cluster.table <- madracis653_assigns_k8[match(tree_reroot_names, madracis653_assigns_k8$ID),] # order species table by new exact tip order
cluster.table$max.cluster.num <- gsub("cluster", "", cluster.table$max.cluster)# remove 'cluster' from max.cluster

p1 <- njtree %<+% species.table +
      geom_tiplab(aes(color = factor(species_mpha)), align=TRUE, size=2, linesize=0.5, hjust = -0.05, show.legend = FALSE) + # without 'show.legend = FALSE' random 'a' was in the species legend
      geom_tippoint(aes(color = factor(species_mpha), shape = factor(species_mpha)), align=TRUE, size=2)  +
      scale_color_manual(values = colors_mpha, limits = c('MCAR', 'MDEC','MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'SINT', 'UNKN')) + scale_shape_manual(values = rep(19,10)) +  theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, 'cm')) + xlim_tree(0.3) # other ways to add more space for long tip labels are 1) scale_x_continuous(expand=expand_scale(0.7)) and 2) xlim_expand(c(0,2)

p2 <- facet_plot(p1, panel= "z.structure k6", data=madracis_k6_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity') +
                 scale_fill_manual(values = colors7) + theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.5, 'cm')) # for whatever reason, the facet plots are being plotted in alphabetical order regardless of the order they are called in, so need to add letters to get all subsequent facets to plot afterwards and then relabel at the end

# # p2.2 <- p2 + facet_set(label = c('Tree'= '1.tree', 'z.structure k5'= '2.structure k5')) # thought I needed this to get the ordering of the facets right, but I don't think I do
# 
p3 <- facet_plot(p2, panel="z.structure k7", data=madracis_k7_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p4 <- facet_plot(p3, panel="z.structure k8", data=madracis_k8_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

p5 <- facet_plot(p4, panel="z.structure k8 big", data=madracis_k8_melt_big,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')


p6 <- facet_plot(p5, panel="z.structure k9", data=madracis_k9_melt,
                 geom=geom_barh, width = 1.0,
                 mapping = aes(x = value, fill = as.factor(variable)),
                 stat='identity')

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
p10 <- facet_labeller(p6, label = c('Tree' = 'tree', 'z.structure k6'= 'k6 (12454 snps)', 'z.structure k7'= 'k7 (12454 snps)', 'z.structure k8'= 'k8 (12454 snps)', 'z.structure k8 big'= 'k8 (7588346 snps)', 'z.structure k9'= 'k9 (12454 snps)')) #, 'zz.assign k8' = 'assignment'))
                                    
# 'zzz.site' = 'site', 'zzzz.<10, 11-20, 21-40, 41-60, 60+' = 'depth')) #relabels tree facets once the facets are in the desired order

# Change the width of the individual facet_plot panels
# gt$layout$l[grep('panel-4', gt$layout$name)] # obtain id for panel-4
gt <- ggplot_gtable(ggplot_build(p10))
gt$widths[5] <- 0.15*gt$widths[5] #tree
gt$widths[7] <- 0.05*gt$widths[7]
gt$widths[9] <- 0.05*gt$widths[9]
gt$widths[11] <- 0.05*gt$widths[11]
gt$widths[13] <- 0.05*gt$widths[13] 
gt$widths[15] <- 0.05*gt$widths[15]

# gt$widths[15] <- 0.05*gt$widths[15] # site
# gt$widths[17] <- 0.05*gt$widths[17] # depth
grid.draw(gt)

ggsave(p10,file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/figures/madracis653_composite.png", width = 50, height = 150, units = "cm", limitsize = FALSE)


# Selecting the best K
#### Load NGSadmixture log files to determine best K ####
ngsadmix_loglist <- list.files(here("admixture_results/madracis653/logs/"), pattern = ".log", full.names = T)

#use lapply to read in all our log files at once
loglines <- lapply(1:4, FUN = function(i) readLines(ngsadmix_loglist[i]))

#this will pull out the line that starts with "b" from each file and return it as a list
foundset <- sapply(1:4, FUN= function(x) loglines[[x]][which(str_sub(loglines[[x]], 1, 1) == 'b')])

as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ) # this should extract the log likelihood

# store the log likelihoods in a dataframe
# make a dataframe with an index 6:9, this corresponds to our K values
logs <- data.frame(K = rep(6:9, each=1)) # 'each' modifies number of replicates for each K

# add to it our likelihood values
logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

# Calculate our delta K and probability
tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

