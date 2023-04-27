## load packages
library(here)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(RcppCNPy)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### Read in necessary files ####
#Load the covariance matrix
cov_angsd <- as.matrix(read.table(here("angsd_results/madracis_lib1_angsd_pca_global_snps.covMat"), header = F)) # This is the covariance matrix created by angsd

cov_pcangsd <- as.matrix(read.table(here("angsd_results/PCAngsd_madracis_lib1_global_covmat.cov"), header = F)) # This is the covariance matrix created by pcangsd

# Load bam names to assign species
bam_paths <- read.table(here("angsd_results/bam_lists/library1_bamlist.txt"), header = FALSE, sep = '\t', col.names = c('path'))

# Split file paths to extract metadata
names <- str_split_fixed(bam_paths$path, '.fixmate.sorted.markdup.', 2) # Get rid of suffix
names2 <- str_split_fixed(names[,1], '/', 7)
metadata <- str_split_fixed(names2[,7], '_', 6)

species <- metadata[,1]
site <- metadata[,2]
depth <- metadata[,3]
id_original <- paste(metadata[,4], metadata[,5], sep = "_")
extract <- metadata[,6]

# Set colors
colors <- brewer.pal(n=6, name = 'Dark2')
show_col(colors)

#### Do PCA using different covariance matricies ####
# This is covariance matrix created by angsd
madracis.pca <- eigen(cov_angsd) #perform the pca using the eigen function. 

# This is covariance matrix created by pcangsd
madracis.pca <- eigen(cov_pcangsd) 

eigenvectors <- madracis.pca$vectors #extract eigenvectors 
pca.vectors <- as_tibble(cbind(species, data.frame(eigenvectors))) #combine with our species identities

# Calculate variance explained by each axis
pca.eigenval.sum <- sum(madracis.pca$values) #sum of eigenvalues

varPC1 <- (madracis.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (madracis.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (madracis.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (madracis.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

eig_percent <- round((madracis.pca$values)/(pca.eigenval.sum)*100,2)
eig_percent[1:10]

# Quick ggplot
pca <- ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = species)) + geom_point()
plot(pca)

#### Base plot for angsd covariance matrix
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/pca_angsd_library1_global_snps.png", width=9, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors[factor(species)], 0.8), xlim = c(-0.1, 0.2), ylim = c(-0.2,0.4), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.1,0.2, by=0.05), labels=seq(-0.1,0.2, by= 0.05), line = -1)
axis(2, at=seq(-0.2,0.4, by = 0.1), labels=seq(-0.2,0.4, by= 0.1), line = -1.5, las = 2)
mtext("PC1 (15.77%)", side = 1, line = 1.5)
mtext("PC2 (3.36%)", side = 2, line = 1.5)

legend('topright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors,
       cex = 0.9)

dev.off()

# PCA comparison for PC1-5
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/pca_library1_global_snps_pc1-5.png", width=11, height=9, res=300, units="in")
plot(pca.vectors[,c(2:6)], col = alpha(colors[factor(species)], 0.8), pch = 16)
dev.off()


#### Base plot for pcangsd covariance matrix
# By species, PC1 vs PC2
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/pca_pcangsd_library1_global_snps.png", width=9, height=7, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n'
)

plot(pca.vectors$X1, pca.vectors$X2, col = add.alpha(colors[factor(species)], 0.8), xlim = c(-0.1, 0.2), ylim = c(-0.2,0.4), pch = 16, xaxt = 'n', xlab = '', ylab = '', yaxt = 'n', frame.plot = F)
axis(1, at=seq(-0.1,0.2, by=0.05), labels=seq(-0.1,0.2, by= 0.05), line = -1)
axis(2, at=seq(-0.2,0.4, by = 0.1), labels=seq(-0.2,0.4, by= 0.1), line = -1.5, las = 2)
mtext("PC1 (39.36%)", side = 1, line = 1.5)
mtext("PC2 (6.81%)", side = 2, line = 1.5)

legend('topright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors,
       cex = 0.9)

dev.off()

# PCA comparison for PC1-5
png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_WGS/angsd_results/pca_pcangsd_library1_global_snps_pc1-5.png", width=11, height=9, res=300, units="in")
plot(pca.vectors[,c(2:6)], col = alpha(colors[factor(species)], 0.8), pch = 16)
dev.off()
