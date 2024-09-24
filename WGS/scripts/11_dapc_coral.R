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

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot to confirm all looks good
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = metadata$species_mpha)) + geom_point()
plot(pca)

#### DAPC for 606 Madracis ####
# Prepare dataframe for DAPC
mad606_dapc_data <- cbind(ids$new_ID, metadata$species_mpha, data.frame(madracis.pca$vectors))
colnames(mad606_dapc_data)[3:(dim(mad606_dapc_data)[2])] <- paste0("PC",1:(dim(mad606_dapc_data)[2]-2)) #rename columns
colnames(mad606_dapc_data)[1:2] <- c("individual", "population")

# Evaluate number of PCs to retain
plot(cumsum(eig_percent)[1:15], type = 'b', ylab = 'Cumulative variance (%)', xlab = 'Number of retained PCs') #variance explained by PCA; look where curve starts to flatten out, maybe 4-5
plot(1:15, madracis.pca$values[1:15], xlab = 'principal component', ylab = 'eigenvalue', type = 'b') #screeplot; look for elbow, 4-5 PCs

n <- 5 # the number of principal components to use for discriminant analysis

#### Do PCA-LDA ####
# set.seed(175)
# train.index <- sample(c(1:608), 449, replace=FALSE) #455/606 = 70%
# train.index <- c(train.index,8, 125, 127, 128, 561, 564) # There are only two MCAR, so put 1 in training dataset; 5 UNKN so put 2 in training; 6 MFOR so put 3 in training

set.seed(175)

train.index <- c(586,112,114,442,301,302,303,294,296,145,116,258,239,474,328,233,147,1,272) #1 for MCAR,3 MFOR, 3 MPHA.ln3, 2 MSEN-like MPHA, 10 for UNKN = 19
train.index2 <- sample(setdiff(1:606,c(586,112,114,442,301,302,303,294,296,145,116,258,239,474,328,233,147,1,272)), 406, replace=FALSE) #425/606 = 70%: sample 417, excluding c(586,112,114,442,301,302,303,294,296,145,116,258,239,474,328,233,147,1,272)
train.indices <- c(train.index, train.index2)

# Split into training-test datasets
train <- mad606_dapc_data[train.indices,]
test <- mad606_dapc_data[-train.indices,]

# Use 5 PCs to fit training data & then use that to predict the test datasets. Finally get training coordinates
pc5.fit <- lda(population ~ PC1 + PC2 + PC3 + PC4 + PC5 , data=train, na.action="na.omit", CV=F) # needed for predict step
pc5.pred <- predict(object = pc5.fit, newdata = test)

pc5.train <- predict(object = pc5.fit, newdata = train)

# Extract eigenvalues (similar to PCA)
eigenvalues_lda <- pc5.fit$svd^2

# Calculate the proportion of variance explained by each linear discriminant
variance_explained_lda <- round((eigenvalues_lda / sum(eigenvalues_lda))*100,2)

# Assess accuracy of the prediction
# percent correct for each species
ct1 <- table(pc5.pred$class, test$population) # correct assignment rate of test dataset (rows are 1st entry, columns are 2nd)
diag(prop.table(ct1,1))

# total percent correct
sum(diag(prop.table(ct1)))

# Assess accuracy of entire training-test dataset
all.class <- as.factor(append(as.character(pc5.pred$class), as.character(pc5.train$class)))
all.class.ordered <- factor(all.class, levels = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA.ln1", "MPHA.ln2", "MPHA.ln3", "MSEN-like MPHA", "MSEN", "UNKN"))
all.pop <- as.factor(append(as.character(test$population), as.character(train$population)))
all.pop.ordered <- factor(all.pop, levels = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA.ln1", "MPHA.ln2", "MPHA.ln3", "MSEN-like MPHA", "MSEN", "UNKN"))

ct2 <- table(all.class.ordered, all.pop.ordered)# correct assignment rate of training-test dataset
diag(prop.table(ct2,1))
sum(diag(prop.table(ct2)))

#### Plot LDA as a scatterplot ####
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/madracis606_dapc.png", width=10, height=9, res=300, units="in")

par(
  mar=c(4, 4, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  par(mfrow = c(2,2))
)

plot(pc5.pred$x[,"LD1"], pc5.pred$x[,"LD2"], col = add.alpha(colors_mpha[-1][factor(test$population)], 0.7), pch = 16, xlab = 'LD1 (51.5%)', ylab = 'LD2 (38.7%)', xlim = c(-40,20), ylim = c(-20,35))
points(pc5.train$x[,"LD1"], pc5.train$x[,"LD2"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

plot(pc5.pred$x[,"LD2"], pc5.pred$x[,"LD3"], col = add.alpha(colors_mpha[-1][factor(test$population)], 0.7), pch = 16, xlab = 'LD2 (38.7%)', ylab = 'LD3 (5.5%)', xlim = c(-20,40), ylim = c(-35,10))
points(pc5.train$x[,"LD2"], pc5.train$x[,"LD3"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

plot(pc5.pred$x[,"LD1"], pc5.pred$x[,"LD3"], col = add.alpha(colors_mpha[-1][factor(test$population)], 0.7), pch = 16, xlab = 'LD1 (51.5%)', ylab = 'LD3 (5.5%)', xlim = c(-40,20), ylim = c(-35,10))
points(pc5.train$x[,"LD1"], pc5.train$x[,"LD3"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

legend(-39,-10, 
       legend = levels(factor(train$population)),
       pch = 16,
       col = add.alpha(colors_mpha, 0.7),
       cex = 0.75,
       bty = 'n')

# plot(pc5.pred$x[,"LD3"], pc5.pred$x[,"LD4"], col = add.alpha(colors_mpha[-1][factor(test$population)], 0.7), pch = 16, xlab = 'LD3 (5.5%)', ylab = 'LD4 (2.9%)', xlim = c(-40,10), ylim = c(-25,15))
# points(pc5.train$x[,"LD3"], pc5.train$x[,"LD4"], col = add.alpha(colors_mpha[factor(train$population)], 0.7), pch = 16)

dev.off()

#### Plot LDA as a barplot ####
posterior_names <- rbind(test[,1:2], train[,1:2])
posteriors <- rbind(pc5.pred$posterior, pc5.train$posterior) # join the group membership posteriors for the training and test datasets
# row.names(posteriors) <- posterior_names
posteriors2 <- as.data.frame(cbind(posterior_names, posteriors))
rownames(posteriors2) <- posterior_names$individual
write.table(posteriors2, here('WGS/angsd_results/madracis606/madracis606_dapc_5pcs_posteriors.txt'), quote = FALSE, sep = '\t', row.names = FALSE) # write text file that contains DAPC posterior values so that it can be used for the big composite plot

# Set plotting order
ord<-order(posteriors2$population)

png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/madracis606_dapc_barplot.png", width=40, height=5, res=300, units="in")

par(mar = c(0.5, 3.5, 0.5, 1),
    mgp = c(3, 1, 0),
    oma = c(10,3,0,1) +0.1,
    xpd=NA)

barplot(t(posteriors2)[3:12,ord], # dont use rows 1 and 2
        col=colors_mpha,
        cex.names = 0.3,
        las=2,
        space=0,
        border=NA,
        ylab="",
        yaxt = 'n',
        xaxt = 'n')

axis(1, at=seq(0.5,605.5, by=1), labels=paste(posteriors2$individual[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 7)

axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -6, las = 2)
mtext("membership probability", side = 2, line = -2)

dev.off()
