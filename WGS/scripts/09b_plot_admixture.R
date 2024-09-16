## load packages
library(here)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(pophelper)

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
write.table(for_clumpak2, here("WGS/angsd_results/madracis606/admixture/likes_clumpak_input.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

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

##################################################################################
#### Load 'best' NGSadmixture results for K=6-10, for visualization across Ks ####
##################################################################################
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

##############################################################################################################################
#### Load NGSadmixture replicate results for K=5-16, for visualization within and across Ks ####
# For when multiple replicates per K are read in. Array dimensions for each K must match number of replicates
q5_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K5*.qopt"), full.names = T)
q5_path <- as.data.frame(str_split_fixed(q5_list, '/', 13))
q5_rep <- as.data.frame(str_split_fixed(q5_path$V13, '_', 5))
# q5_files <- lapply(q5_list, FUN = function(x)read.table(x, header=F, fill = T))
# q5_files_array <- array(as.numeric(unlist(q5_files)), dim=c(606, 5, 10))
qlist5 <- readQ(files = q5_list)

q6_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K6*.qopt"), full.names = T)
q6_path <- as.data.frame(str_split_fixed(q6_list, '/', 13))
q6_rep <- as.data.frame(str_split_fixed(q6_path$V13, '_', 5))
# q6_files <- lapply(q6_list, FUN = function(x)read.table(x, header=F, fill = T))
# q6_files_array <- array(as.numeric(unlist(q6_files)), dim=c(606, 6, 10))
qlist6 <- readQ(files = q6_list)

q7_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K7*.qopt"), full.names = T)
q7_path <- as.data.frame(str_split_fixed(q7_list, '/', 13))
q7_rep <- as.data.frame(str_split_fixed(q7_path$V13, '_', 5))
# q7_files <- lapply(q7_list, FUN = function(x)read.table(x, header=F, fill = T))
# q7_files_array <- array(as.numeric(unlist(q7_files)), dim=c(606, 7, 10))
qlist7 <- readQ(files = q7_list)

q8_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K8*.qopt"), full.names = T)
q8_path <- as.data.frame(str_split_fixed(q8_list, '/', 13))
q8_rep <- as.data.frame(str_split_fixed(q8_path$V13, '_', 5))
# q8_files <- lapply(q8_list, FUN = function(x)read.table(x, header=F, fill = T))
# q8_files_array <- array(as.numeric(unlist(q8_files)), dim=c(606, 8, 10))
qlist8 <- readQ(files = q8_list)

q9_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K9*.qopt"), full.names = T)
q9_path <- as.data.frame(str_split_fixed(q9_list, '/', 13))
q9_rep <- as.data.frame(str_split_fixed(q9_path$V13, '_', 5))
# q9_files <- lapply(q9_list, FUN = function(x)read.table(x, header=F, fill = T))
# q9_files_array <- array(as.numeric(unlist(q9_files)), dim=c(606, 9, 10))
qlist9 <- readQ(files = q9_list)

q10_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K10*.qopt"), full.names = T)
q10_path <- as.data.frame(str_split_fixed(q10_list, '/', 13))
q10_rep <- as.data.frame(str_split_fixed(q10_path$V13, '_', 5))
# q10_files <- lapply(q10_list, FUN = function(x)read.table(x, header=F, fill = T))
# q10_files_array <- array(as.numeric(unlist(q10_files)), dim=c(606, 10, 10))
qlist10 <- readQ(files = q10_list)

q11_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K11*.qopt"), full.names = T)
q11_path <- as.data.frame(str_split_fixed(q11_list, '/', 13))
q11_rep <- as.data.frame(str_split_fixed(q11_path$V13, '_', 5))
qlist11 <- readQ(files = q11_list)

q12_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K12*.qopt"), full.names = T)
q12_path <- as.data.frame(str_split_fixed(q12_list, '/', 13))
q12_rep <- as.data.frame(str_split_fixed(q12_path$V13, '_', 5))
qlist12 <- readQ(files = q12_list)

q13_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K13*.qopt"), full.names = T)
q13_path <- as.data.frame(str_split_fixed(q13_list, '/', 13))
q13_rep <- as.data.frame(str_split_fixed(q13_path$V13, '_', 5))
qlist13 <- readQ(files = q13_list)

q14_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K14*.qopt"), full.names = T)
q14_path <- as.data.frame(str_split_fixed(q14_list, '/', 13))
q14_rep <- as.data.frame(str_split_fixed(q14_path$V13, '_', 5))
qlist14 <- readQ(files = q14_list)

q15_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K15*.qopt"), full.names = T)
q15_path <- as.data.frame(str_split_fixed(q15_list, '/', 13))
q15_rep <- as.data.frame(str_split_fixed(q15_path$V13, '_', 5))
qlist15 <- readQ(files = q15_list)

q16_list <- list.files(here("WGS/angsd_results/madracis606/admixture/qfiles"), pattern = glob2rx("*K16*.qopt"), full.names = T)
q16_path <- as.data.frame(str_split_fixed(q16_list, '/', 13))
q16_rep <- as.data.frame(str_split_fixed(q16_path$V13, '_', 5))
qlist16 <- readQ(files = q16_list)

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
# colors10 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#0D0887FF", "#5402A3FF", "lightblue")
colors16 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#0D0887FF", "#5402A3FF", "lightblue", "gray30", "burlywood4", "darkcyan", "plum2", "pink1", "wheat")
show_col(colors16)

#### Manipulate the data to ensure consistent ordering across x and y axes ####
# Set plotting order for x axis for all Ks
ord<-order(metadata$species_mpha)

# Set y axis order for K5
# q5_files_array_ord <- q5_files_array
# q5_files_array_ord[,,2] <- q5_files_array[, c(1,4,2,3,5), 2]  # Reorder each array slice individually based on slice1
# q5_files_array_ord[,,3] <- q5_files_array[, c(5,2,1,4,3), 3]
# q5_files_array_ord[,,4] <- q5_files_array[, c(3,4,5,2,1), 4]
# q5_files_array_ord[,,5] <- q5_files_array[, c(3,1,2,4,5), 5]
# q5_files_array_ord[,,6] <- q5_files_array[, c(5,4,3,1,2), 6]
# q5_files_array_ord[,,7] <- q5_files_array[, c(3,5,4,1,2), 7]
# q5_files_array_ord[,,8] <- q5_files_array[, c(2,4,5,1,3), 8]
# q5_files_array_ord[,,9] <- q5_files_array[, c(3,2,1,5,4), 9]
# q5_files_array_ord[,,10] <- q5_files_array[, c(2,3,5,1,4), 10]

# Set y axis order for K6
# q6_files_array_ord <- q6_files_array
# q6_files_array_ord[,,2] <- q6_files_array[, c(6,3,5,2,4,1), 2]  # Reorder each array slice individually based on slice1
# q6_files_array_ord[,,3] <- q6_files_array[, c(2,1,5,3,4,6), 3]
# q6_files_array_ord[,,4] <- q6_files_array[, c(4,5,2,1,3,6), 4]
# q6_files_array_ord[,,5] <- q6_files_array[, c(1,2,3,4,5,6), 5]
# q6_files_array_ord[,,6] <- q6_files_array[, c(1,2,3,4,5,6), 6]
# q6_files_array_ord[,,7] <- q6_files_array[, c(1,2,3,4,5,6), 7]
# q6_files_array_ord[,,8] <- q6_files_array[, c(1,2,3,4,5,6), 8]
# q6_files_array_ord[,,9] <- q6_files_array[, c(5,1,6,2,3,4), 9]
# q6_files_array_ord[,,10] <- q6_files_array[, c(6,2,3,1,5,4), 10]

aligned_qlist5 <- alignK(qlist5)
q5_files_array_aligned <- array(unlist(aligned_qlist5), dim=c(606, 5, 10))

aligned_qlist6 <- alignK(qlist6)
q6_files_array_aligned <- array(unlist(aligned_qlist6), dim=c(606, 6, 10))

aligned_qlist7 <- alignK(qlist7)
q7_files_array_aligned <- array(unlist(aligned_qlist7), dim=c(606, 7, 10))

aligned_qlist8 <- alignK(qlist8)
q8_files_array_aligned <- array(unlist(aligned_qlist8), dim=c(606, 8, 10))

aligned_qlist9 <- alignK(qlist9)
q9_files_array_aligned <- array(unlist(aligned_qlist9), dim=c(606, 9, 10))

aligned_qlist10 <- alignK(qlist10)
q10_files_array_aligned <- array(unlist(aligned_qlist10), dim=c(606, 10, 10))

aligned_qlist11 <- alignK(qlist11)
q11_files_array_aligned <- array(unlist(aligned_qlist11), dim=c(606, 11, 10))

aligned_qlist12 <- alignK(qlist12)
q12_files_array_aligned <- array(unlist(aligned_qlist12), dim=c(606, 12, 10))

aligned_qlist13 <- alignK(qlist13)
q13_files_array_aligned <- array(unlist(aligned_qlist13), dim=c(606, 13, 10))

aligned_qlist14 <- alignK(qlist14)
q14_files_array_aligned <- array(unlist(aligned_qlist14), dim=c(606, 14, 10))

aligned_qlist15 <- alignK(qlist15)
q15_files_array_aligned <- array(unlist(aligned_qlist15), dim=c(606, 15, 10))

aligned_qlist16 <- alignK(qlist16)
q16_files_array_aligned <- array(unlist(aligned_qlist16), dim=c(606, 16, 10))

### Plot ####
pdf(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/admixture/k5-16_10reps.pdf", width=30, height=18)

par(mfrow = c(10,1),
    mar = c(0.5, 4.5, 0.5, 1),
    mgp = c(3, 1, 0),
    oma = c(10,3,0,1) +0.1,
    #omi=c(1,0,0,0), 
    xpd=NA)

# Loop over the 10 slices of the array for K5
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q5_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    #names=ids$new_ID[ord],             # Names of individuals
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q5_rep$V3[i], q5_rep$V4[i], sep = '_'), side = 2, line = -4)
}

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 6)

# Loop over the 10 slices of the array for K6
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q6_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q6_rep$V3[i], q6_rep$V4[i], sep = '_'), side = 2, line = -4)
}

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 6)

# Loop over the 10 slices of the array for K7
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q7_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q7_rep$V3[i], q7_rep$V4[i], sep = '_'), side = 2, line = -4)
}

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 6)

# Loop over the 10 slices of the array for K8
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q8_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q8_rep$V3[i], q8_rep$V4[i], sep = '_'), side = 2, line = -4)
}

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 6)

# Loop over the 10 slices of the array for K9
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q9_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q9_rep$V3[i], q9_rep$V4[i], sep = '_'), side = 2, line = -4)
}

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 6)

# Loop over the 10 slices of the array for K10
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q10_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q10_rep$V3[i], q10_rep$V4[i], sep = '_'), side = 2, line = -4)
}

# Loop over the 10 slices of the array for K11
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q11_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q11_rep$V3[i], q11_rep$V4[i], sep = '_'), side = 2, line = -4)
}

# Loop over the 10 slices of the array for K12
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q12_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q12_rep$V3[i], q12_rep$V4[i], sep = '_'), side = 2, line = -4)
}

# Loop over the 10 slices of the array for K13
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q13_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q13_rep$V3[i], q13_rep$V4[i], sep = '_'), side = 2, line = -4)
}

# Loop over the 10 slices of the array for K14
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q14_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q14_rep$V3[i], q14_rep$V4[i], sep = '_'), side = 2, line = -4)
}

# Loop over the 10 slices of the array for K15
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q15_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q15_rep$V3[i], q15_rep$V4[i], sep = '_'), side = 2, line = -4)
}

# Loop over the 10 slices of the array for K16
for (i in 1:10) {
  # Generate the barplot for the i-th slice
  barplot(
    t(q16_files_array_aligned[,,i])[, ord],  # Transpose and subset the array
    col = colors16,              # Colors for the stacked bars
    cex.names = 0.5,                # Font size for the names
    las = 2,                        # Rotate axis labels
    space = 0,                      # No space between bars
    border = NA,                    # No borders around bars
    ylab = "",
    yaxt = 'n'
  )
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), line = -7, las = 2)
  mtext(paste(q16_rep$V3[i], q16_rep$V4[i], sep = '_'), side = 2, line = -4)
}

axis(1, at=seq(0.5,605.5, by=1), labels=paste(ids$new_ID[ord]), line = 0.5, las = 2, cex.axis = 0.3)
mtext("Individual", side = 1, line = 6)

dev.off()



















barplot(t(q5_files_array[,,1])[,ord],
        col=colors_mpha,
        #names=species[ord],
        cex.names = 0.5,
        las=2,
        space=0,
        border=NA,
        ylab="K=6")




