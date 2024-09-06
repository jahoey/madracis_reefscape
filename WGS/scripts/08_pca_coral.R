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

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

###############################################################################
#### Full analysis of 606 Madracis, clones removed, SNPs are not LD pruned ####
###############################################################################

#### Read in necessary files ####
#Load the pcangsdcovariance matrix
pcangsd_coral606_linked <- as.matrix(read.table(here("WGS/angsd_results/madracis606/PCAngsd_madracis606_allsnps_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd using minMaf 0.05 default

# Load id names to assign metadata
ids <- read.csv(here("WGS/angsd_results/bam_lists/ids_recoded606.csv"), header = TRUE) # This is 606 colonies, clones removed

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata $id_original <- paste(metadata[,4], metadata[,5], sep = "_")
metadata$depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird

# Get cluster assignments for coloring the points
data_cluster <- read.csv(here("madracis_master_merged_May2024.csv"), header = TRUE)

metadata$species_mpha <- data_cluster[match(metadata$id_original, data_cluster$ID), 7]

# Replace NAs with matching in species column
species_indices <- which(is.na(metadata$species_mpha) == TRUE) # indicies that have NAs
metadata$species_mpha[is.na(metadata$species_mpha)] <- metadata$species[species_indices] # replace the NAs with species in associated index

# Split out MPHA by depth to assign MPHA.ln3
# species_mpha <- species
# species_mpha[which(species == 'MPHA' & depth <=20)] <- "MPHA-shallow" # MPHA at depths of 20m or shallower
# species_mpha[which(species == 'MPHA' & depth >20 & depth <=60)] <- "MPHA-deep" # MPHA at depths between 20-60
metadata$species_mpha[which(metadata$species == 'MPHA' & metadata$depth >60)] <- "MPHA.ln3" # MPHA at depths >60

# Deal with the MPHAs, which are mostly MSEN-like MPHAs
metadata[which(metadata$id_original == 'PB_4583'),8] <- 'MPHA.ln2'
metadata[which(metadata$id_original == 'PB_4578'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4577'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4576'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4579'),8] <- 'MSEN-like MPHA'

metadata$species_mpha <- factor(metadata$species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MSEN-like MPHA', 'MSEN', 'UNKN')) #reorder the species factors

# Set color scheme
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

# Decide which covariance matrix to use 
madracis.pca <- eigen(pcangsd_coral606_linked)# this is the pcangsd covariance matrix

eigenvectors <- madracis.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(metadata$species_mpha, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(madracis.pca$values) #sum of eigenvalues

varPC1 <- (madracis.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (madracis.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (madracis.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (madracis.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = metadata$species_mpha)) + geom_point()
plot(pca)

######################################################################################################################################
#### Base plot for pcangsd covariance matrix ####
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/madracis606_coralpca_pcangsd.png", width=9.5, height=9, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors_mpha[factor(metadata$species_mpha)], 0.7), xlim = c(-0.1, 0.025), ylim = c(-0.1,0.05), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.1,0.025, by=0.01), labels=seq(-0.1,0.025, by= 0.01), line = -1)
axis(2, at=seq(-0.1,0.05, by = 0.05), labels=seq(-0.1,0.05, by= 0.05), line = -2, las = 2)
mtext("PC1 (33.75%)", side = 1, line = 1.5) # modify variance explained to match dataset
mtext("PC2 (5.91%)", side = 2, line = 1) # modify variance explained to match dataset

legend(-0.09,-0.04, 
       legend = levels(factor(metadata$species_mpha)),
       pch = 19,
       col = add.alpha(colors_mpha, 0.7),
       cex = 1)

dev.off()

# PC1-PC10
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/madracis606_coralpca_pcangsd_pcs1-10.png", width=14, height=9, res=300, units="in")

plot(pca.vectors[,2:11], col = alpha(colors_mpha[factor(metadata$species_mpha)], 0.7), pch = 16, lower.panel = NULL)

dev.off()

#### Tree (exploratory) ####
## heatmap / clustering / trees from angsd wiki
#heat map
heatmap(pcangsd_coral606)

#neighbour joining
plot(hclust(dist(pcangsd_coral606), "ave"), main = "", labels = paste(metadata$species_mpha, "_", metadata$depth, sep = ''), cex = 0.3)

#######################################################################
#### Full analysis of 606 Madracis, clones removed, LD pruned SNPs ####
#######################################################################

#### Read in necessary files ####
#Load the pcangsdcovariance matrix
pcangsd_coral606_unlinked <- as.matrix(read.table(here("WGS/angsd_results/madracis606/PCAngsd_madracis606_unlinkedsnps_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd using minMaf 0.05 default

# Load id names to assign metadata
ids <- read.csv(here("WGS/angsd_results/bam_lists/ids_recoded606.csv"), header = TRUE) # This is 606 colonies, clones removed

metadata <- as.data.frame(str_split_fixed(ids$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata $id_original <- paste(metadata[,4], metadata[,5], sep = "_")
metadata$depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird

# Get cluster assignments for coloring the points
data_cluster <- read.csv(here("madracis_master_merged_May2024.csv"), header = TRUE)

metadata$species_mpha <- data_cluster[match(metadata$id_original, data_cluster$ID), 7]

# Replace NAs with matching in species column
species_indices <- which(is.na(metadata$species_mpha) == TRUE) # indicies that have NAs
metadata$species_mpha[is.na(metadata$species_mpha)] <- metadata$species[species_indices] # replace the NAs with species in associated index

# Split out MPHA by depth to assign MPHA.ln3
# species_mpha <- species
# species_mpha[which(species == 'MPHA' & depth <=20)] <- "MPHA-shallow" # MPHA at depths of 20m or shallower
# species_mpha[which(species == 'MPHA' & depth >20 & depth <=60)] <- "MPHA-deep" # MPHA at depths between 20-60
metadata$species_mpha[which(metadata$species == 'MPHA' & metadata$depth >60)] <- "MPHA.ln3" # MPHA at depths >60

# Deal with the MPHAs, which are mostly MSEN-like MPHAs
metadata[which(metadata$id_original == 'PB_4583'),8] <- 'MPHA.ln2'
metadata[which(metadata$id_original == 'PB_4578'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4577'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4576'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4579'),8] <- 'MSEN-like MPHA'

metadata$species_mpha <- factor(metadata$species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MSEN-like MPHA', 'MSEN', 'UNKN')) #reorder the species factors

# Set color scheme
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#374618", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

# Decide which covariance matrix to use 
madracis.pca <- eigen(pcangsd_coral606_unlinked)# this is the pcangsd covariance matrix

eigenvectors <- madracis.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(metadata$species_mpha, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(madracis.pca$values) #sum of eigenvalues

varPC1 <- (madracis.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (madracis.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (madracis.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (madracis.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = metadata$species_mpha)) + geom_point()
plot(pca)

######################################################################################################################################
#### Base plot for pcangsd covariance matrix ####
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/madracis606_coralpca_pcangsd_unlinked.png", width=7.5, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors_mpha[factor(metadata$species_mpha)], 0.7), xlim = c(-0.025, 0.1), ylim = c(-0.1,0.05), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.02,0.1, by=0.02), labels=seq(-0.02,0.1, by= 0.02), line = -1)
axis(2, at=seq(-0.1,0.05, by = 0.05), labels=seq(-0.1,0.05, by= 0.05), line = -2, las = 2)
mtext("PC1 (16.58%)", side = 1, line = 2) # modify variance explained to match dataset
mtext("PC2 (5.18%)", side = 2, line = 1) # modify variance explained to match dataset

legend(0.060,-0.03, 
       legend = levels(factor(metadata$species_mpha)),
       pch = 19,
       col = add.alpha(colors_mpha, 0.7),
       cex = 1,
       bty = 'n')

dev.off()

# PC1-PC10
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/madracis606_coralpca_pcangsd_unlinked_pcs1-10.png", width=14, height=9, res=300, units="in")

plot(pca.vectors[,2:11], col = alpha(colors_mpha[factor(metadata$species_mpha)], 0.7), pch = 16, lower.panel = NULL)

dev.off()






#### DAPC for 650 Madracis ####
# Prepare dataframe for DAPC
mad650_dapc_data <- cbind(ids2$new_ID, species_mpha, data.frame(madracis.pca$vectors))
colnames(mad650_dapc_data)[3:(dim(mad650_dapc_data)[2])] <- paste0("PC",1:(dim(mad650_dapc_data)[2]-2))
colnames(mad650_dapc_data)[1:2] <- c("individual", "population")

# Evaluate number of PCs to retain
plot(cumsum(eig_percent)[1:15], type = 'b')
plot(1:15, madracis.pca$values[1:15], xlab = 'principal component', ylab = 'eigenvalue') #screeplot

n <- 6 # the number of principal components to use for discriminant analysis

#### Do PCA-LDA ####
set.seed(175)
train.index <- sample(c(1:650), 449, replace=FALSE) #455/650 = 70%
train.index <- c(train.index,8, 125, 127, 128, 561, 564) # There are only two MCAR, so put 1 in training dataset; 5 UNKN so put 2 in training; 6 MFOR so put 3 in training

# Split into training-test datasets
train <- mad650_dapc_data[train.index,]
test <- mad650_dapc_data[-train.index,]

# Use 6 PCs to fit training data & then use that to predict the test datasets. Finally get training coordinates
pc6.fit <- lda(population ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data=train, na.action="na.omit", CV=F) # needed for predict step
pc6.pred <- predict(object = pc6.fit, newdata = test)

pc6.train <- predict(object = pc6.fit, newdata = train)

# Extract eigenvalues (similar to PCA)
eigenvalues_lda <- pc6.fit$svd^2

# Calculate the proportion of variance explained by each linear discriminant
variance_explained_lda <- eigenvalues_lda / sum(eigenvalues_lda)

# Assess accuracy of the prediction
# percent correct for each species
table(pc6.pred$class, test$population) # correct assignment rate of test dataset (rows are 1st entry, columns are 2nd)
ct1 <- table(pc6.pred$class, test$population)
diag(prop.table(ct1,1))

# total percent correct
sum(diag(prop.table(ct1)))

# Assess accuracy of entire training-test dataset
all.class <- as.factor(append(as.character(pc6.pred$class), as.character(pc6.train$class)))
all.class.ordered <- factor(all.class, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'UNKN'))
all.pop <- as.factor(append(as.character(test$population), as.character(train$population)))
all.pop.ordered <- factor(all.pop, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA-shallow', 'MPHA-deep', 'MPHA-very deep', 'MSEN', 'UNKN'))

ct2 <- table(all.class.ordered, all.pop.ordered)# correct assignment rate of training-test dataset
diag(prop.table(ct2,1))
sum(diag(prop.table(ct2)))

#### Plot LDA ####
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/madracis650_dapc.png", width=10, height=9, res=300, units="in")

par(
  mar=c(4, 4, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  par(mfrow = c(2,2))
)

plot(pc6.pred$x[,"LD1"], pc6.pred$x[,"LD2"], col = add.alpha(colors_mpha[factor(test$population)], 0.7), pch = 16, xlab = 'LD1 (77.3%)', ylab = 'LD2 (16.6%)', xlim = c(-100,10), ylim = c(-25,33))
points(pc6.train$x[,"LD1"], pc6.train$x[,"LD2"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

plot(pc6.pred$x[,"LD2"], pc6.pred$x[,"LD3"], col = add.alpha(colors_mpha[factor(test$population)], 0.7), pch = 16, xlab = 'LD2 (16.6%)', ylab = 'LD3 (4.1%)', xlim = c(-22,32), ylim = c(-4,7))
points(pc6.train$x[,"LD2"], pc6.train$x[,"LD3"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

plot(pc6.pred$x[,"LD1"], pc6.pred$x[,"LD3"], col = add.alpha(colors_mpha[factor(test$population)], 0.7), pch = 16, xlab = 'LD1 (77.3%)', ylab = 'LD3 (4.1%)', xlim = c(-100,20), ylim = c(-4,7))
points(pc6.train$x[,"LD1"], pc6.train$x[,"LD3"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

dev.off()

