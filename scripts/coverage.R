library(here)
library(tidyr)
library(viridis)
library(stringr)
library(RColorBrewer)
library(scales)

add.alpha <- function(col, alpha=1){
        if(missing(col))
                stop("Please provide a vector of colours.")
        apply(sapply(col, col2rgb)/255, 2, 
              function(x) 
                      rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Set color scheme
# col <- viridis(6)
colors <- brewer.pal(n=8, name = 'Dark2')
show_col(colors)

#### Coverage and quality data when reads for 672 Madracis mapped to ASG MMIR genome ####
asg_mmir_mapping_stats <- read.table(here("mapping/mapping2mmir_asg_June2024_stats_headers.txt"), header = TRUE, sep = '\t') # this is sample level mapping stats merged across sequencing runs

asg_mmir_mapping_stats$proportion_mapping <- asg_mmir_mapping_stats$number_of_mapped_reads/asg_mmir_mapping_stats$number_of_reads

# Extract metadata from filenames
names <- asg_mmir_mapping_stats$bam
names_split <- str_split_fixed(names, '_', 7)
names_split2 <- str_split_fixed(names_split[,7], '.dedup.overlapclipped.bam', 2)

species <- names_split[,2]
site <- names_split[,3]
depth <- names_split[,4]
id_original <- paste(names_split[,5], names_split[,6], sep = "_")
extract <- names_split2[,1]

#### Read in recoded ids for 672 Madracis ####
ids672 <- read.csv(here("angsd_results/bam_lists/ids_recoded672.csv"), header = TRUE) # This is 672 colonies

metadata <- as.data.frame(str_split_fixed(ids672$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')

# What is the average coverage for 672 samples?
mean(asg_mmir_mapping_stats$coverage_mean, na.rm = TRUE) # 2.09X
mean(asg_mmir_mapping_stats$coverage_sd, na.rm = TRUE) # 4.13X
mean(asg_mmir_mapping_stats$number_of_reads, na.rm = TRUE) # 8752428
mean(asg_mmir_mapping_stats$number_of_mapped_reads, na.rm = TRUE) # 8197947
mean(asg_mmir_mapping_stats$mean_mapping_quality, na.rm = TRUE) # 53.81X
mean(asg_mmir_mapping_stats$proportion_mapping, na.rm = TRUE) # 0.927

# Total reads vs coverage
# Total reads = reads remaining after trimming
plot(jitter(asg_mmir_mapping_stats$number_of_reads), jitter(asg_mmir_mapping_stats$coverage_mean), col = add.alpha(colors[factor(metadata$species)], 0.8), pch = 19, las = 1, xlab = "Number of reads after trimming", ylab = "Average coverage across ASG MMIR")
legend('topleft',
       legend = levels(factor(metadata$species)),
       pch = 19,
       col = colors)

# Mapped reads vs coverage
plot(asg_mmir_mapping_stats$number_of_mapped_reads, asg_mmir_mapping_stats$coverage_mean, col = add.alpha(colors[factor(metadata$species)], 0.8), pch = 19, las = 1, xlab = "Number of mapped reads", ylab = "Average coverage across ASG MMIR")
legend('bottomright',
       legend = levels(factor(metadata$species)),
       pch = 19,
       col = colors)

# Zoomed in a bit
plot(jitter(asg_mmir_mapping_stats$number_of_mapped_reads), jitter(asg_mmir_mapping_stats$coverage_mean), col = add.alpha(colors[factor(metadata$species)], 0.8), pch = 19, las = 1, xlab = "Number of mapped reads", ylab = "Average coverage across ASG MMIR", xlim = c(0,1000000), ylim = c(0,0.25))
legend('topleft',
       legend = levels(factor(metadata$species)),
       pch = 19,
       col = colors)

# Even more zoomed in
plot(jitter(asg_mmir_mapping_stats$number_of_mapped_reads), jitter(asg_mmir_mapping_stats$coverage_mean), col = add.alpha(colors[factor(metadata$species)], 0.8), pch = 19, las = 1, xlab = "Number of mapped reads", ylab = "Average coverage across ASG MMIR", xlim = c(0,450000), ylim = c(0,0.10))
legend('topleft',
       legend = levels(factor(metadata$species)),
       pch = 19,
       col = colors)

# Total vs mapped
plot(asg_mmir_mapping_stats$number_of_reads, asg_mmir_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(metadata$species)], 0.8), pch = 19, xlab = "Number of reads after trimming", ylab = "Number of mapped reads")
legend('topleft',
       legend = levels(factor(metadata$species)),
       pch = 19,
       col = colors)

# Proportion mapped
i <- aggregate(asg_mmir_mapping_stats$proportion_mapping, by = list(metadata$species), FUN = mean) # proportion mapped by species
barplot(i$x, col = colors, las = 1, names.arg = i$Group.1, xlab = 'Species ID', main = 'Proportion mapping to ASG MMIR')
