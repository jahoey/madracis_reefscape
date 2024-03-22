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

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### MDEC intraspecies population structure ####
# # Read in necessary files #
# pcangsd_coral_MDEC161 <- as.matrix(read.table(here("angsd_results/MDEC/setMinDepth85/PCAngsd_madracis_MDEC161_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd using default minMaf 0.05 & setMinDepth 85
# 
# # Load id names to assign metadata
# ids_MDEC <- read.csv(here("angsd_results/bam_lists/MDEC.csv"), header = TRUE)
# 
# metadata <- as.data.frame(str_split_fixed(ids_MDEC$new_ID, '_', 6))
# colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
# metadata$depth <- as.numeric(metadata$depth)
# 
# species <- metadata[,1]
# site <- metadata[,2]
# depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
# id_original <- paste(metadata[,4], metadata[,5], sep = "_")
# extract <- metadata[,6]
# 
# # Set color schemes
# colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666") #'Dark2' R color brewer palette
# show_col(colors)
# 
# # Decide which covariance matrix to use 
# mdec.pca <- eigen(pcangsd_coral_MDEC161)# this is the pcangsd covariance matrix
# 
# eigenvectors <- mdec.pca$vectors #extract eigenvectors 
# pca.vectors <- as_tibble(cbind(ids_MDEC$new_ID, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities
# 
# # Calculate variance explained by each axis
# pca.eigenval.sum <- sum(mdec.pca$values) #sum of eigenvalues
# 
# varPC1 <- (mdec.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
# varPC2 <- (mdec.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
# varPC3 <- (mdec.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
# varPC4 <- (mdec.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
# 
# eig_percent <- round((mdec.pca$values)/(pca.eigenval.sum)*100,2)
# eig_percent[1:10]
# 
# # Quick ggplot
# pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = metadata$depth)) + geom_point()
# plot(pca)
# 
# #### Base plot for pcangsd covariance matrix ####
# # By species, PC1 vs PC2
# png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/MDEC/setMinDepth85/mdec161_coral_site_pcangsd.png", width=9, height=8, res=300, units="in")
# 
# par(
#   mar=c(4, 3, 3, 1), # panel margin size in "line number" units
#   mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#   tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#   cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#   ps=14, # point size, which is the font size
#   bty = 'n'
#  )
# 
# plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors[factor(metadata$site)], 0.7), xlim = c(-0.1, 0.7), ylim = c(-0.75,0.25), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
# axis(1, at=seq(-0.1,0.65, by=0.25), labels=seq(-0.1,0.65, by= 0.25), line = -1)
# axis(2, at=seq(-0.75,0.25, by = 0.25), labels=seq(-0.75,0.25, by= 0.25), line = -1.5, las = 2)
# mtext("PC1 (2.01%)", side = 1, line = 1.5)
# mtext("PC2 (1.05%)", side = 2, line = 1.5)
# 
# legend(0.5,-0.4,
#        legend = levels(factor(metadata$site)),
#        pch = 19,
#        col = add.alpha(colors, 0.7),
#        cex = 0.9)
# 
# dev.off()
# 
# # PC1-PC10
# plot(pca.vectors[,2:11], col = alpha(colors[factor(metadata$site)], 0.7), pch = 16, lower.panel = NULL)
# 
# #### Remove all outliers ####
# pcangsd_coral_MDEC161_rmouts <- pcangsd_coral_MDEC161[-c(18,20,56,64,138,149,93,99,6:8), -c(18,20,56,64,138,149,93,99,6:8)]
# mdec.pca_rmouts <- eigen(pcangsd_coral_MDEC161_rmouts)
# 
# eigenvectors.rmouts <- mdec.pca_rmouts$vectors #extract eigenvectors 
# pca.vectors.rmouts <- as_tibble(cbind(ids_MDEC$new_ID[-c(18,20,56,64,138,149,93,99,6:8)], data.frame(eigenvectors.rmouts)))
# 
# # Calculate variance explained by each axis
# pca.eigenval.sum.rmouts <- sum(mdec.pca_rmouts$values) #sum of eigenvalues
# 
# eig_percent <- round((mdec.pca_rmouts$values)/(pca.eigenval.sum.rmouts)*100,2)
# eig_percent[1:10]
# 
# pca <- ggplot(data = pca.vectors.rmouts, aes(x=X1, y=X2, colour = metadata$depth[-c(18,20,56,64,138,149,93,99,6:8)])) + geom_point()
# plot(pca)
# 
# plot(pca.vectors.rmouts$X1, pca.vectors.rmouts$X2, col = add.alpha(colors[factor(metadata$site[-c(18,20,56,64,138,149,93,99,6:8)])], 0.7), xlim = c(-0.25, 0.25), ylim = c(-0.2,0.35), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
# axis(1, at=seq(-0.25,0.25, by=0.25), labels=seq(-0.25,0.25, by= 0.25), line = -1)
# axis(2, at=seq(-0.2,0.35, by = 0.25), labels=seq(-0.2,0.35, by= 0.25), line = -1.5, las = 2)
# mtext("PC1 (0.79%)", side = 1, line = 1.5)
# mtext("PC2 (0.75%)", side = 2, line = 1.5)
# 
# legend(0.2,0.25,
#        legend = levels(factor(metadata$site[-c(18,20,56,64,138,149,93,99,6:8)])),
#        pch = 19,
#        col = add.alpha(colors, 0.7),
#        cex = 0.9)
# 
# # PC1-PC10
# plot(pca.vectors.rmouts[,2:11], col = alpha(colors[factor(metadata$site[-c(64,149,93,99,6:8)])], 0.7), pch = 16, lower.panel = NULL)
# 
# #### Remove all outliers, more lax ####
# pcangsd_coral_MDEC161_rmouts_lax <- pcangsd_coral_MDEC161[-c(64,149,93,99,6:8), -c(64,149,93,99,6:8)]
# mdec.pca_rmouts_lax <- eigen(pcangsd_coral_MDEC161_rmouts_lax)
# 
# eigenvectors.rmouts.lax <- mdec.pca_rmouts_lax$vectors #extract eigenvectors 
# pca.vectors.rmouts.lax <- as_tibble(cbind(ids_MDEC$new_ID[-c(64,149,93,99,6:8)], data.frame(eigenvectors.rmouts.lax)))
# 
# # Calculate variance explained by each axis
# pca.eigenval.sum.rmouts.lax <- sum(mdec.pca_rmouts_lax$values) #sum of eigenvalues
# 
# eig_percent <- round((mdec.pca_rmouts_lax$values)/(pca.eigenval.sum.rmouts.lax)*100,2)
# eig_percent[1:10]
# 
# # PC1-PC10
# plot(pca.vectors.rmouts.lax[,2:11], col = alpha(colors[factor(metadata$site[-c(64,149,93,99,6:8)])], 0.7), pch = 16, lower.panel = NULL)

##############################################################################################################
# Read in necessary files #
pcangsd_coral_MDEC161 <- as.matrix(read.table(here("angsd_results/MDEC/setMinDepth5/PCAngsd_madracis_MDEC161_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd using default minMaf 0.05 & setMinDepth 85

# Load id names to assign metadata
ids_MDEC <- read.csv(here("angsd_results/bam_lists/MDEC.csv"), header = TRUE)

metadata <- as.data.frame(str_split_fixed(ids_MDEC$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata$depth <- as.numeric(metadata$depth)

species <- metadata[,1]
site <- metadata[,2]
depth <- as.numeric(metadata[,3]) # Must be numeric or the subsetting gets very weird
id_original <- paste(metadata[,4], metadata[,5], sep = "_")
extract <- metadata[,6]

# Set color schemes
colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666") #'Dark2' R color brewer palette
show_col(colors)

# Decide which covariance matrix to use 
mdec.pca <- eigen(pcangsd_coral_MDEC161)# this is the pcangsd covariance matrix

eigenvectors <- mdec.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(ids_MDEC$new_ID, data.frame(eigenvectors))) #combine with our fine-scale MPHA species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(mdec.pca$values) #sum of eigenvalues

varPC1 <- (mdec.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mdec.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mdec.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mdec.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((mdec.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = metadata$depth)) + geom_point()
plot(pca)

#### Base plot for pcangsd covariance matrix ####
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/MDEC/setMinDepth5/mdec161_coral_site_pcangsd.png", width=9, height=8, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors[factor(metadata$site)], 0.7), xlim = c(-0.1, 0.7), ylim = c(-0.75,0.25), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.1,0.65, by=0.25), labels=seq(-0.1,0.65, by= 0.25), line = -1)
axis(2, at=seq(-0.75,0.25, by = 0.25), labels=seq(-0.75,0.25, by= 0.25), line = -1.5, las = 2)
mtext("PC1 (2.01%)", side = 1, line = 1.5)
mtext("PC2 (1.05%)", side = 2, line = 1.5)

legend(0.5,-0.4,
       legend = levels(factor(metadata$site)),
       pch = 19,
       col = add.alpha(colors, 0.7),
       cex = 0.9)

dev.off()

# PC1-PC10
plot(pca.vectors[,2:11], col = alpha(colors[factor(metadata$site)], 0.7), pch = 16, lower.panel = NULL)

#### DAPC for 161 MDEC ####
# Prepare dataframe for DAPC
mdec161_dapc_data <- cbind(ids_MDEC$new_ID, site, data.frame(mdec.pca$vectors))
colnames(mdec161_dapc_data)[3:(dim(mdec161_dapc_data)[2])] <- paste0("PC",1:(dim(mdec161_dapc_data)[2]-2))
colnames(mdec161_dapc_data)[1:2] <- c("individual", "site")

# Evaluate number of PCs to retain
plot(cumsum(eig_percent)[1:15], type = 'b')
plot(1:15, mdec.pca$values[1:15], xlab = 'principal component', ylab = 'eigenvalue') #screeplot

n <- 6 # the number of principal components to use for discriminant analysis

#### Do PCA-LDA ####
set.seed(34)
train.index <- sample(c(1:161), 113, replace=FALSE) #113/161 = 70%

# Split into training-test datasets
train <- mdec161_dapc_data[train.index,]
test <- mdec161_dapc_data[-train.index,]

# Use 6 PCs to fit training data & then use that to predict the test datasets. Finally get training coordinates
pc6.fit <- lda(site ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data=train, na.action="na.omit", CV=F) # needed for predict step
pc6.pred <- predict(object = pc6.fit, newdata = test)

pc6.train <- predict(object = pc6.fit, newdata = train)

# Extract eigenvalues (similar to PCA)
eigenvalues_lda <- pc6.fit$svd^2

# Calculate the proportion of variance explained by each linear discriminant
variance_explained_lda <- eigenvalues_lda / sum(eigenvalues_lda)

# Assess accuracy of the prediction
# percent correct for each species
table(pc6.pred$class, test$site) # correct assignment rate of test dataset (rows are 1st entry, columns are 2nd)
ct1 <- table(pc6.pred$class, test$site)
diag(prop.table(ct1,1))

# total percent correct
sum(diag(prop.table(ct1)))

# Assess accuracy of entire training-test dataset
all.class <- factor(append(as.character(pc6.pred$class), as.character(pc6.train$class)), levels = c("dir", "est", "hul", "kal", "sea", "sna", "wat"))
all.pop <- factor(append(as.character(test$site), as.character(train$site)), levels = c("dir", "est", "hul", "kal", "sea", "sna", "wat"))

ct2 <- table(all.class, all.pop) # correct assignment rate of training-test dataset
ct2.1 <- diag(prop.table(ct2,1))
ct2.1[is.na(ct2.1)] <- 0
sum(diag(prop.table(ct2.1)))

#### Plot LDA ####
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/MDEC/mdec161_dapc.png", width=10, height=9, res=300, units="in")

par(
  mar=c(4, 4, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  par(mfrow = c(2,2))
)

plot(pc6.pred$x[,"LD1"], pc6.pred$x[,"LD2"], col = add.alpha(colors[factor(test$site)], 0.7), pch = 16, xlab = 'LD1 (83.5%)', ylab = 'LD2 (9.7%)', xlim = c(-10,0.8), ylim = c(-4,5))
points(pc6.train$x[,"LD1"], pc6.train$x[,"LD2"], col = add.alpha(colors[factor(train$site)], 0.7), pch = 16)

plot(pc6.pred$x[,"LD2"], pc6.pred$x[,"LD3"], col = add.alpha(colors[factor(test$site)], 0.7), pch = 16, xlab = 'LD2 (9.7%)', ylab = 'LD3 (5.2%)', xlim = c(-4,5), ylim = c(-4,11))
points(pc6.train$x[,"LD2"], pc6.train$x[,"LD3"], col = add.alpha(colors[factor(train$site)], 0.7), pch = 16)

legend(3.5,11,
       legend = levels(factor(test$site)),
       pch = 19,
       col = add.alpha(colors, 0.7),
       cex = 0.9)

plot(pc6.pred$x[,"LD1"], pc6.pred$x[,"LD3"], col = add.alpha(colors[factor(test$site)], 0.7), pch = 16, xlab = 'LD1 (83.5%)', ylab = 'LD3 (5.2%)', xlim = c(-9,2), ylim = c(-4,10))
points(pc6.train$x[,"LD1"], pc6.train$x[,"LD3"], col = add.alpha(colors[factor(train$site)], 0.7), pch = 16)

dev.off()



