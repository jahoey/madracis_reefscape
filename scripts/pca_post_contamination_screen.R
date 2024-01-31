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

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

##########################################
#### Initial analysis using 87 corals ####
##########################################

#### Read in necessary files ####
#Load the covariance matrix for 87 corals with low, middle and high number of reads
cov_pcangsd <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis87_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd

# Load bam names to assign species
bam_paths <- read.table(here("angsd_results/bam_lists/angsd_bam_list.txt"), header = FALSE, sep = '\t', col.names = c('path'))

# Split file paths to extract metadata
names <- str_split_fixed(bam_paths$path, '.dedup.overlapclipped.bam', 2) # Get rid of suffix
names2 <- str_split_fixed(names[,1], '/', 8)
names3 <- sub('_tech_', '_', names2[,8]) # remove 'tech' prefix of technical replicates so that str_split_fixed will work well
metadata <- str_split_fixed(names3, '_', 7)

species <- metadata[,2]
site <- metadata[,3]
depth <- as.numeric(metadata[,4]) # Must be numeric or the subsetting gets very weird
id_original <- paste(metadata[,5], metadata[,6], sep = "_")
extract <- metadata[,7]

names4 <- paste(species, site, depth, id_original, extract, sep = "_")

# Split out MPHA by depth
species_mpha <- species
species_mpha[which(species == 'MPHA' & depth <=20)] <- "MPHA-shallow" # MPHA at depths of 20m or shallower
species_mpha[which(species == 'MPHA' & depth >20 & depth <=60)] <- "MPHA-deep" # MPHA at depths between 20-60
species_mpha[which(species == 'MPHA' & depth >60)] <- "MPHA-very deep" # MPHA at depths >60

species_mpha <- factor(species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'SINT', 'UNKN')) #reorder the species factors

# Combine back into a dataframe
data_count <- cbind.data.frame(species, species_mpha, site, depth, id_original, extract)

# Import 'low', 'medium' and 'high' categories for plotting
# read_cat <- read.table(here("angsd_results/input/madracis90.full.txt"), header = F, sep = '\t')

# Set color scheme
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#b0e0e6", "#2664d8", "gray70") # new colors
show_col(colors_mpha)

# This is covariance matrix created by pcangsd
madracis.pca <- eigen(cov_pcangsd) #perform the pca using the eigen function. 

eigenvectors <- madracis.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(data_count$species_mpha, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(madracis.pca$values) #sum of eigenvalues

varPC1 <- (madracis.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (madracis.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (madracis.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (madracis.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = data_count$species_mpha)) + geom_point()
plot(pca)

#### Base plot for pcangsd covariance matrix
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis87_coralpca_pcangsd.png", width=8, height=8, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), xlim = c(-0.05, 0.25), ylim = c(-0.25,0.1), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.05,0.25, by=0.05), labels=seq(-0.05,0.25, by= 0.05), line = 0)
axis(2, at=seq(-0.25,0.1, by = 0.1), labels=seq(-0.25,0.1, by= 0.1), line = -3, las = 2)
mtext("PC1 (36.25%)", side = 1, line = 2.5)
mtext("PC2 (6.31%)", side = 2, line = 1)

legend(0.15,-0.1, 
       legend = levels(factor(data_count$species_mpha)),
       pch = 19,
       col = add.alpha(colors_mpha, 0.7),
       cex = 0.9)

dev.off()

# PC1-PC5
plot(pca.vectors[,2:6], col = alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), pch = 16, lower.panel = NULL)

#### Tree ####
## heatmap / clustering / trees from angsd wiki
m <- as.matrix(read.table(here("angsd_results/madracis87_coral_snps.ibsMat"), header = FALSE))

#heat map
heatmap(m)

#neighbour joining
plot(ape::nj(m))
plot(hclust(dist(m), "ave"), main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.6)

## Using pcangsd tree method ##
# tree file
# madracis_tree <- ape::read.tree(here("angsd_results/PCAngsd_madracis87_covmat.tree"))
# plot(madracis_tree, main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.6)
# plot(hclust(madracis_tree, "ave"), main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.6)

# create tree from distance matrix
m <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis87_covmat.tree.cov"), header = FALSE))
heatmap(m)

#neighbour joining
plot(ape::nj(m))
plot(hclust(dist(m), "ave"), main = "", labels = paste(species, "_", site, "_", depth, sep = ''), cex = 0.6)


##########################################
madracis_tree3 <- ape::read.tree(here("angsd_results/PCAngsd_madracis87_covmat.tree"))
tip.order <- as.numeric(madracis_tree3$tip.label)

names4.ordered <- names4[tip.order]
madracis_tree3$tip.label <- names4.ordered
plot(madracis_tree3)

m <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis87_covmat.tree.cov"), header = FALSE)) # distribution of edge lengths is really weird
madracis_tree <- nj(m)
madracis_tree$tip.label <- names4

# Reroot M. senaria
is.rooted (madracis_tree3)
outgroup <- c("MSEN_est_10_JH_0392_E0490", "MSEN_kal_10_JH_0555_E0295", "MSEN_kal_20_JH_0240_E0193", "MSEN_sea_10_JH_0530_E0384", "MSEN_sna_10_JH_0030_E0461", "MSEN_sna_10_JH_0590_E0053", "MSEN_wat_20_JH_0607_E0368", "MSEN_wat_20_JH_0622_E0303", "SINT_sna_20_JH_0013_E0713", "MDEC_kal_10_JH_0089_E0485",
              "MPHA_kal_10_JH_0176_E0284", "MPHA_est_10_JH_0383_E0679",
              "MPHA_est_5_JH_0344_E0608", "MPHA_est_10_JH_0503_E0133",
              "MPHA_est_20_JH_0332_E0083","UNKN_est_40_JH_0431_E0168",
              "UNKN_est_10_JH_0314_E0547") # UNKN_est_60_JH_0304_E0593 is new from first iteration
madracis_tree_reroot <- root(madracis_tree3, outgroup, resolve.root = TRUE)
# madracis_tree_reroot <- root(madracis_tree, outgroup, interactive = TRUE, resolve.root = TRUE)
is.rooted(madracis_tree_reroot)

tree_reroot_names <- get_taxa_name(madracis_tree)

njtree <- ggtree(madracis_tree_reroot, layout="rectangular")

tree_reroot_names <- get_taxa_name(njtree) # get exact tip order of rerooted tree

species.table2 <- species.table[match(tree_reroot_names, species.table$sample_full),] # order species table by new exact tip order

p1 <- njtree %<+% species.table2 +
  geom_tiplab(aes(color = factor(species_mpha)), align=TRUE, size=2, linesize=0.5, hjust = -0.05, show.legend = FALSE) + # without 'show.legend = FALSE' random 'a' was in the species legend
  geom_tippoint(aes(color = factor(species_mpha), shape = factor(species_mpha)), align=TRUE, size=2)  +
  scale_color_manual(values = colors_mpha, limits = c('MCAR', 'MDEC','MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'SINT', 'UNKN')) + scale_shape_manual(values = rep(19,10)) +  theme(legend.text = element_text(size=15), legend.title = element_blank(), legend.key.size = unit(0.8, 'cm')) + xlim_tree(0.3) # other ways to add more space for long tip labels are 1) scale_x_continuous(expand=expand_scale(0.7)) and 2) xlim_expand(c(0,2)









#####################################
#### Full analysis of 653 corals ####
#####################################

#### Read in necessary files ####
#Load the covariance matrix for 87 corals with low, middle and high number of reads
# angsd_coral653 <- as.matrix(read.table(here("angsd_results/coral653_snps.covMat"), header = F)) # This is the covariance matrix created by angsd
# pcangsd_coral653 <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis653_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd using minMaf 0.05 default
pcangsd_coral653 <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis653_covmat_0.01.cov"), header = F)) # This is the covariance matrix created by pcangsd using minMaf 0.01 

# Load id names to assign metadata
ids <- read.csv(here("angsd_results/bam_lists/ids_recoded.csv"), header = TRUE)

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')

species <- metadata[,1]
site <- metadata[,2]
depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
id_original <- paste(metadata[,4], metadata[,5], sep = "_")
extract <- metadata[,6]

# Split out MPHA by depth
species_mpha <- species
species_mpha[which(species == 'MPHA' & depth <=20)] <- "MPHA-shallow" # MPHA at depths of 20m or shallower
species_mpha[which(species == 'MPHA' & depth >20 & depth <=60)] <- "MPHA-deep" # MPHA at depths between 20-60
species_mpha[which(species == 'MPHA' & depth >60)] <- "MPHA-very deep" # MPHA at depths >60

species_mpha <- factor(species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'SINT', 'UNKN')) #reorder the species factors

# Combine back into a dataframe
data_count <- cbind.data.frame(species, species_mpha, site, depth, id_original, extract)

# Set color scheme
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#b0e0e6", "#2664d8", "gray70") # new colors
show_col(colors_mpha)

# Decide which covariance matrix to use 
# madracis.pca <- eigen(angsd_coral653) #perform the pca using the eigen function. This is the angsd covariance matrix
madracis.pca <- eigen(pcangsd_coral653)# this is the pcangsd covariance matrix

eigenvectors <- madracis.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(data_count$species_mpha, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(madracis.pca$values) #sum of eigenvalues

varPC1 <- (madracis.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (madracis.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (madracis.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (madracis.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = data_count$species_mpha)) + geom_point()
plot(pca)

#### Base plot for angsd covariance matrix ####
# By species, PC1 vs PC2
# png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis653_coralpca_angsd.png", width=8, height=8, res=300, units="in")
# 
# par(
#   mar=c(4, 3, 3, 1), # panel margin size in "line number" units
#   mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#   tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#   cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#   ps=14, # point size, which is the font size
#   bty = 'n'
# )
# 
# plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), xlim = c(-0.03, 0.1), ylim = c(-0.6,0.1), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
# axis(1, at=seq(-0.03,0.1, by=0.01), labels=seq(-0.03,0.1, by= 0.01), line = 1)
# axis(2, at=seq(-0.6,0.1, by = 0.1), labels=seq(-0.6,0.1, by= 0.1), line = -1, las = 2)
# mtext("PC1 (7.86%)", side = 1, line = 3)
# mtext("PC2 (2.12%)", side = 2, line = 1.5)
# 
# legend(0.05,-0.1, 
#        legend = levels(factor(data_count$species_mpha)),
#        pch = 19,
#        col = add.alpha(colors_mpha, 0.7),
#        cex = 0.9)
# 
# dev.off()

# PC1-PC10
# plot(pca.vectors[,2:11], col = alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), pch = 16, lower.panel = NULL)

####################################################################################################################################################
#### Base plot for pcangsd covariance matrix ####
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis653_coralpca_pcangsd_0.01.png", width=9, height=9, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), xlim = c(-0.1, 0.025), ylim = c(-0.1,0.05), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.1,0.025, by=0.01), labels=seq(-0.1,0.025, by= 0.01), line = -1)
axis(2, at=seq(-0.1,0.05, by = 0.05), labels=seq(-0.1,0.05, by= 0.05), line = -2, las = 2)
mtext("PC1 (20.86%)", side = 1, line = 1.5) # modify variance explained to match dataset
mtext("PC2 (4.32%)", side = 2, line = 1) # modify variance explained to match dataset

legend(-0.09,-0.04, 
       legend = levels(factor(data_count$species_mpha)),
       pch = 19,
       col = add.alpha(colors_mpha, 0.7),
       cex = 1)

dev.off()

# PC1-PC10
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis653_coralpca_pcangsd_0.01_pcs1-10.png", width=14, height=9, res=300, units="in")

plot(pca.vectors[,2:11], col = alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), pch = 16, lower.panel = NULL)

dev.off()

#### Tree (exploratory) ####
## heatmap / clustering / trees from angsd wiki
m <- as.matrix(read.table(here("angsd_results/coral653_snps.ibsMat"), header = FALSE))

#heat map
heatmap(m)

#neighbour joining
plot(ape::nj(m))
plot(hclust(dist(m), "ave"), main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.3)

## Using pcangsd tree method ##
# tree file
madracis_tree <- ape::read.tree(here("angsd_results/PCAngsd_madracis653_covmat_0.01.tree"))
plot(madracis_tree, main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.6)
plot(hclust(madracis_tree, "ave"), main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.6)

# create tree from distance matrix
m <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis653_covmat_0.01.tree.cov"), header = FALSE))
colnames(m) <- ids$new_ID
rownames(m) <- ids$new_ID
heatmap(m, cexRow = 0.07, cexCol = 0.07)

#neighbour joining
plot(ape::nj(m))
pdf(file = "~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/tree653_angsd_cov.pdf", width = 50, height = 8)
plot(hclust(dist(m), "ave"), main = "", labels = paste(species, "_", site, "_", depth, sep = ''), cex = 0.5)
dev.off()

# from pcangsd tree file
madracis_tree3 <- ape::read.tree(here("angsd_results/PCAngsd_madracis653_covmat_0.01.tree"))
tip.order <- as.numeric(madracis_tree3$tip.label)

names4.ordered <- ids$new_ID[tip.order]
madracis_tree3$tip.label <- names4.ordered
plot(madracis_tree3)


#####################################################
#### Full analysis of 650 corals, 3 SINT removed ####
#####################################################

#### Read in necessary files ####
#Load the pcangsdcovariance matrix
pcangsd_coral650 <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis650_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd using minMaf 0.05 default

# Load id names to assign metadata
ids <- read.csv(here("angsd_results/bam_lists/ids_recoded.csv"), header = TRUE) # This is 653 colonies
ids2<- ids[-c(626,633,648),] #remove the 3 SINT colonies

metadata <- as.data.frame(str_split_fixed(ids2$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')

species <- metadata[,1]
site <- metadata[,2]
depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
id_original <- paste(metadata[,4], metadata[,5], sep = "_")
extract <- metadata[,6]

# Split out MPHA by depth
species_mpha <- species
species_mpha[which(species == 'MPHA' & depth <=20)] <- "MPHA-shallow" # MPHA at depths of 20m or shallower
species_mpha[which(species == 'MPHA' & depth >20 & depth <=60)] <- "MPHA-deep" # MPHA at depths between 20-60
species_mpha[which(species == 'MPHA' & depth >60)] <- "MPHA-very deep" # MPHA at depths >60

species_mpha <- factor(species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'UNKN')) #reorder the species factors

# Combine back into a dataframe
data_count <- cbind.data.frame(species, species_mpha, site, depth, id_original, extract)

# Set color scheme
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

# Decide which covariance matrix to use 
madracis.pca <- eigen(pcangsd_coral650)# this is the pcangsd covariance matrix

eigenvectors <- madracis.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(data_count$species_mpha, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(madracis.pca$values) #sum of eigenvalues

varPC1 <- (madracis.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (madracis.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (madracis.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (madracis.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = data_count$species_mpha)) + geom_point()
plot(pca)

####################################################################################################################################################
#### Base plot for pcangsd covariance matrix ####
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis650_coralpca_pcangsd.png", width=9, height=9, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), xlim = c(-0.1, 0.025), ylim = c(-0.1,0.05), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.1,0.025, by=0.01), labels=seq(-0.1,0.025, by= 0.01), line = -1)
axis(2, at=seq(-0.1,0.05, by = 0.05), labels=seq(-0.1,0.05, by= 0.05), line = -2, las = 2)
mtext("PC1 (31.40%)", side = 1, line = 1.5) # modify variance explained to match dataset
mtext("PC2 (5.90%)", side = 2, line = 1) # modify variance explained to match dataset

legend(-0.09,-0.04, 
       legend = levels(factor(data_count$species_mpha)),
       pch = 19,
       col = add.alpha(colors_mpha, 0.7),
       cex = 1)

dev.off()

# PC1-PC10
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis650_coralpca_pcangsd_pcs1-10.png", width=14, height=9, res=300, units="in")

plot(pca.vectors[,2:11], col = alpha(colors_mpha[factor(data_count$species_mpha)], 0.7), pch = 16, lower.panel = NULL)

dev.off()

#### Tree (exploratory) ####
## heatmap / clustering / trees from angsd wiki
#heat map
heatmap(pcangsd_coral650)

#neighbour joining
plot(ape::nj(pcangsd_coral650))
plot(hclust(dist(pcangsd_coral650), "ave"), main = "", labels = paste(species, "_", depth, sep = ''), cex = 0.3)
