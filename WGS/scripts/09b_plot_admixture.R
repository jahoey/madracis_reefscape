## load packages
library(here)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### Load NGSadmixture log files to determine best K ####
ngsadmix_loglist <- list.files(here("WGS/angsd_results/madracis606/admixture/logs"), pattern = ".log", full.names = T)

#use lapply to read in all our log files at once
loglines <- lapply(1:120, FUN = function(i) readLines(ngsadmix_loglist[i]))

#this will pull out the line that starts with "b" from each file and return it as a list
foundset <- sapply(1:120, FUN= function(x) loglines[[x]][which(str_sub(loglines[[x]], 1, 1) == 'b')])

as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ) # this should extract the log likelihood

# store the log likelihoods in a dataframe
# Get K value from the filename, since R orders numbers weirdly
path <- as.data.frame(str_split_fixed(ngsadmix_loglist, '/', 13))
file <- as.data.frame(str_split_fixed(path$V13, '_', 5))
colnames(file) <- c('file1', 'file2', 'K', 'rep', 'file3')
logs <- file[,3:4]

# add to it our likelihood values
logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

# Generate a two column text file (K, likelihood) for input into CLUMPAK
for_clumpak <- sapply(str_split(logs$K, "K"), "[[",2)
for_clumpak2 <- cbind(for_clumpak, logs$like)
write.table(for_clumpak2, here("WGS/angsd_results/madracis606/admixture/likes_clumpak_input.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

# Calculate our delta K and probability
bestK <- data.frame(
  K = as.integer(unique(str_split_fixed(logs$K, 'K', 2))[,2]),  # Convert K to character
  LK = tapply(logs$like, logs$K, mean),
  sd = tapply(logs$like, logs$K, sd)
)

# Order the Ks in ascending order
bestK.ordered <- bestK[order(bestK[,1]),]

# Use likelihoods to calculate delta K
bestK.ordered <- bestK.ordered %>%
  mutate(LK_diff = LK - lag(LK),
         LK_diff2 = lead(LK_diff - lag(LK_diff))) %>%  # Shift LK_diff2 up by using lead()
  mutate(deltaK = abs(LK_diff2) / sd)                 # Calculate deltaK after shifting LK_diff2

# Plot the likelihoods & delta Ks to pick optimal K using Evanno method
plot(bestK.ordered$K, log(bestK.ordered$LK), type = 'b', xlab = 'K', ylab = 'Log Likelihoods')

plot(bestK.ordered$K, bestK.ordered$deltaK, type = 'b', xlab = 'K', ylab = 'deltaK') # highest value is the optimal K using Evanno method

#### Load NGSadmixture results for K=6-8 ####
# For when multiple replicates per K are read in. Array dimensions for each K must match number of replicates
# q6_list <- list.files(here("admixture_results/qfiles/K6/"), pattern = ".qopt", full.names = T)
# q6_files <- lapply(q6_list, FUN = function(x)read.table(x, header=F, fill = T))
# q6_files_array <- array(as.numeric(unlist(q6_files)), dim=c(606, 6, 1))

# For when the best run for a given K are read in
k6 <- read.table(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K6_5_out.qopt"), header = F)
k7 <- read.table(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K7_9_out.qopt"), header = F)
k8 <- read.table(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K8_3_out.qopt"), header = F)
k9 <- read.table(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K9_9_out.qopt"), header = F)
k10 <- read.table(here("WGS/angsd_results/madracis606/admixture/qfiles/madracis606_ngsAdmix_K10_6_out.qopt"), header = F)

#### Load bam names to assign species ####
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
metadata$species_mpha[which(metadata$species == 'MPHA' & metadata$depth >60)] <- "MPHA.ln3" # MPHA at depths >60

# Deal with the MPHAs, which are mostly MSEN-like MPHAs
metadata[which(metadata$id_original == 'PB_4583'),8] <- 'MPHA.ln2'
metadata[which(metadata$id_original == 'PB_4578'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4577'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4576'),8] <- 'MSEN-like MPHA'
metadata[which(metadata$id_original == 'PB_4579'),8] <- 'MSEN-like MPHA'

metadata$species_mpha <- factor(metadata$species_mpha, levels = c('MCAR','MDEC', 'MFOR', 'MMIR', 'MPHA.ln1', 'MPHA.ln2', 'MPHA.ln3', 'MSEN-like MPHA', 'MSEN', 'UNKN')) #reorder the species factors

# For STRUCTURE plots
colors10 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#0D0887FF", "#5402A3FF", "lightblue")
show_col(colors10)

#### Plot ####
# Set plotting order
ord<-order(metadata$species_mpha)

png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/admixture/k6-10.png", width=40, height=15, res=300, units="in")

par(mfrow = c(5,1),
    mar = c(0.5, 3.5, 0.5, 1),
    mgp = c(3, 1, 0),
    oma = c(10,3,0,1) +0.1,
    #omi=c(1,0,0,0), 
    xpd=NA)

# barplot(t(q6_files_array[,,1])[,ord],
#         col=colors_mpha,
#         #names=species[ord],
#         cex.names = 0.5,
#         las=2,
#         space=0,
#         border=NA,
#         ylab="K=6")

barplot(t(k6)[,ord],
        col=colors10,
        cex.names = 0.5,
        las=2,
        space=0,
        border=NA,
        ylab="",
        yaxt = 'n')
axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -10, las = 2)
mtext("K=6", side = 2, line = -6)

barplot(t(k7)[c(7,3,2,1,6,5,4),ord], # order the cluster proportions based on K=6 to match the colors
        col=colors10,
        cex.names = 0.5,
        las=2,
        space=0,
        border=NA,
        ylab="",
        yaxt = 'n')
axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -10, las = 2)
mtext("K=7", side = 2, line = -6)

barplot(t(k8)[c(6,5,3,7,8,1,4,2),ord],
        col=colors10,
        cex.names = 0.5,
        las=2,
        space=0,
        border=NA,
        ylab="",
        yaxt = 'n')
axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -10, las = 2)
mtext("K=8", side = 2, line = -6)

mtext("Ancestry proportions", side = 2, line = -1)

barplot(t(k9)[c(8,3,2,6,5,1,7,4,9),ord],
        col=colors10,
        cex.names = 0.5,
        las=2,
        space=0,
        border=NA,
        ylab="",
        yaxt = 'n')
axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -10, las = 2)
mtext("K=9", side = 2, line = -6)

barplot(t(k10)[c(5,2,4,1,7,10,9,3,6,8),ord],
        col=colors10,
        cex.names = 0.5,
        las=2,
        space=0,
        border=NA,
        ylab="",
        yaxt = 'n')
axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -10, las = 2)
mtext("K=10", side = 2, line = -6)

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.5)
mtext("Individual", side = 1, line = 9)

dev.off()
