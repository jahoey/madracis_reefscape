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
colors <- brewer.pal(n=6, name = 'Dark2')
show_col(colors)

#### Coverage and quality data when reads mapped to M. senaria de novo genome ####
msen_denovo_mapping_stats <- read.table(here("genotyping/mapping2msen_denovo_stats_headers_lib1.txt"), header = TRUE, sep = '\t')

# Calculate proportion of reads mapping to M. senaria de novo
msen_denovo_mapping_stats$proportion_mapping <- msen_denovo_mapping_stats$number_of_mapped_reads/msen_denovo_mapping_stats$number_of_reads

# Clean up names & add to coverage data
names <- str_split_fixed(msen_denovo_mapping_stats$bam, '.fixmate.sorted.markdup.', 2)
names_split <- str_split_fixed(names[,1], '_', 6)

species <- names_split[,1]
site <- names_split[,2]
depth <- names_split[,3]
id_original <- paste(names_split[,4], names_split[,5], sep = "_")
extract <- names_split[,6]

# What is the average coverage across all 96 samples in library 1?
mean(msen_denovo_mapping_stats$coverage_mean) # 2.15X
mean(msen_denovo_mapping_stats$coverage_sd) # 13.7X
mean(msen_denovo_mapping_stats$number_of_reads) # 10687847
mean(msen_denovo_mapping_stats$number_of_mapped_reads) # 9622868
mean(msen_denovo_mapping_stats$mean_mapping_quality) # 20.05

# Total reads vs coverage
# Total reads = reads remaining after trimming
plot(jitter(msen_denovo_mapping_stats$number_of_reads), jitter(msen_denovo_mapping_stats$coverage_mean), col = add.alpha(col[factor(species)], 0.8), pch = 19, las = 1, xlim = c(407000, 68000000), ylim = c(0,15), xlab = "Number of reads after trimming", ylab = "Average coverage across M. senaria de novo")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = col)

mod1 <- lm(coverage_mean ~ number_of_reads, data = msen_denovo_mapping_stats) # y = (1.994e-07 * #of reads) + 2.214e-02
abline(mod1)
text(2.5e+07, 14, 'y = (1.994e-07 * # of reads) + 2.214e-02', cex = 0.7)

# Inset for low coverage
plot(jitter(msen_denovo_mapping_stats$number_of_reads), jitter(msen_denovo_mapping_stats$coverage_mean), col = add.alpha(col[factor(species)], 0.8), pch = 19, las = 1, xlim = c(407000, 22000000), ylim = c(0,4), xlab = "Number of reads after trimming", ylab = "Average coverage across M. senaria de novo")

# Mapped reads vs coverage
plot(msen_denovo_mapping_stats$number_of_mapped_reads, msen_denovo_mapping_stats$coverage_mean, col = add.alpha(col[factor(species)], 0.8), pch = 19, las = 1, ylim = c(0,15), xlab = "Number of mapped reads", ylab = "Average coverage across M. senaria de novo")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = col)

mod2 <- lm(coverage_mean ~ number_of_mapped_reads, data = msen_denovo_mapping_stats) # y = (2.218e-07 * #of mapped reads) + 1.877e-02
abline(mod2)
text(2.9e+07, 14, 'y = (2.218e-07 * #of mapped reads) + 1.877e-02', cex = 0.7)

# Inset for low coverage
plot(msen_denovo_mapping_stats$number_of_mapped_reads, msen_denovo_mapping_stats$coverage_mean, col = add.alpha(col[factor(species)], 0.8), pch = 19, las = 1, xlim = c(407000, 20000000), ylim = c(0,4), xlab = "Number of mapped reads", ylab = "Average coverage across M. senaria de novo")

# Total vs mapped
plot(msen_denovo_mapping_stats$number_of_reads, msen_denovo_mapping_stats$number_of_mapped_reads, col = add.alpha(col[factor(species)], 0.8), pch = 19, xlab = "Number of reads after trimming", ylab = "Number of mapped reads")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = col)

mod3 <- lm(number_of_mapped_reads ~ number_of_reads, data = msen_denovo_mapping_stats) # mapped reads = (0.8994 * # of total reads) + 9711.3684
abline(mod3)
text(2.7e+07, 6e+07, 'y = (0.8994 * # of total reads) + 9711.3684', cex = 0.7)

# Do the earlier sampled Madracis have worse mapping/coverage?
plot(names_split[,5], msen_denovo_mapping_stats$number_of_mapped_reads, col = add.alpha(col[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Number of mapped reads")
plot(names_split[,5], msen_denovo_mapping_stats$coverage_mean, col = add.alpha(col[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Average coverage across M. senaria de novo")
legend('topleft', 
       legend = levels(factor(species)),
       pch = 19,
       col = col)

# How do the mapping stats differ by species identity?
a <- aggregate(msen_denovo_mapping_stats$number_of_mapped_reads, by = list(species), FUN = mean) # number of mapped reads by species
b <- aggregate(msen_denovo_mapping_stats$mean_mapping_quality, by = list(species), FUN = mean) # mapping quality by species
c <- aggregate(msen_denovo_mapping_stats$coverage_mean, by = list(species), FUN = mean) # average coverage by species
d <- aggregate(msen_denovo_mapping_stats$proportion_mapping, by = list(species), FUN = mean) # proportion mapped by species

barplot(a$x, col = colors, las = 1, names.arg = a$Group.1, xlab = 'Species ID - Field', main = 'Average number of mapped reads')
barplot(b$x, col = colors, las = 1, names.arg = b$Group.1, xlab = 'Species ID - Field', main = 'Average mapping quality')
barplot(c$x, col = colors, las = 1, names.arg = c$Group.1, xlab = 'Species ID - Field', main = 'Average coverage')
barplot(d$x, col = colors, las = 1, names.arg = d$Group.1, xlab = 'Species ID - Field', main = 'Average proportion mapping')


temp <- data.frame(cbind(as.numeric(names_split[,5]), as.numeric(msen_denovo_mapping_stats$coverage_mean)))


#### Coverage and quality data when reads mapped to M. mirabilis mito genome ####
mmir_mito_mapping_stats <- read.table(here("genotyping/mapping2mmir_mito_stats_headers_lib1.txt"), header = TRUE, sep = '\t')

# Calculate proportion of reads mapping to M. senaria de novo
mmir_mito_mapping_stats$proportion_mapping <- mmir_mito_mapping_stats$number_of_mapped_reads/mmir_mito_mapping_stats$number_of_reads

# Clean up names & add to coverage data
names <- str_split_fixed(mmir_mito_mapping_stats$bam, '.fixmate.sorted.markdup.', 2)
names_split <- str_split_fixed(names[,1], '_', 6)

species <- names_split[,1]
site <- names_split[,2]
depth <- names_split[,3]
id_original <- paste(names_split[,4], names_split[,5], sep = "_")
extract <- names_split[,6]

# What is the average coverage across all 96 samples in library 1?
mean(mmir_mito_mapping_stats$coverage_mean) # 38X
mean(mmir_mito_mapping_stats$coverage_sd) # 25X
mean(mmir_mito_mapping_stats$number_of_reads) # 9916774
mean(mmir_mito_mapping_stats$number_of_mapped_reads) # 6220.344
mean(mmir_mito_mapping_stats$mean_mapping_quality) # 57X

# Total reads vs coverage
# Total reads = reads remaining after trimming
plot(jitter(mmir_mito_mapping_stats$number_of_reads), jitter(mmir_mito_mapping_stats$coverage_mean), col = add.alpha(colors[factor(species)], 0.8), pch = 19, las = 1, xlab = "Number of reads after trimming", ylab = "Average coverage across M. mirabilis mito")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Mapped reads vs coverage
plot(mmir_mito_mapping_stats$number_of_mapped_reads, mmir_mito_mapping_stats$coverage_mean, col = add.alpha(colors[factor(species)], 0.8), pch = 19, las = 1, xlab = "Number of mapped reads", ylab = "Average coverage across M. mirabilis mito")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Total vs mapped
plot(mmir_mito_mapping_stats$number_of_reads, mmir_mito_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Number of reads after trimming", ylab = "Number of mapped reads")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Do the earlier sampled Madracis have worse mapping/coverage?
plot(names_split[,5], mmir_mito_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Number of mapped reads")
plot(names_split[,5], mmir_mito_mapping_stats$coverage_mean, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Average coverage across M. mirabilis mito")
legend('topleft', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Proportion mapped
e <- aggregate(mmir_mito_mapping_stats$proportion_mapping, by = list(species), FUN = mean) # proportion mapped by species
barplot(e$x, col = colors, las = 1, names.arg = e$Group.1, xlab = 'Species ID - Field', main = 'Proportion mapping to M. mirabilis mito')

#### Coverage and quality data when reads mapped to P. verrucosa genome ####
pver_mapping_stats <- read.table(here("genotyping/mapping2pver_stats_headers_lib1.txt"), header = TRUE, sep = '\t')

# Calculate proportion of reads mapping to M. senaria de novo
pver_mapping_stats$proportion_mapping <- pver_mapping_stats$number_of_mapped_reads/pver_mapping_stats$number_of_reads

# Clean up names & add to coverage data
names <- str_split_fixed(pver_mapping_stats$bam, '.fixmate.sorted.markdup.', 2)
names_split <- str_split_fixed(names[,1], '_', 6)

species <- names_split[,1]
site <- names_split[,2]
depth <- names_split[,3]
id_original <- paste(names_split[,4], names_split[,5], sep = "_")
extract <- names_split[,6]

# What is the average coverage across all 96 samples in library 1?
mean(pver_mapping_stats$coverage_mean) # 0.36X
mean(pver_mapping_stats$coverage_sd) # 13X
mean(pver_mapping_stats$number_of_reads) # 10011968
mean(pver_mapping_stats$number_of_mapped_reads) # 2383365
mean(pver_mapping_stats$mean_mapping_quality) # 7.13X

# Total reads vs coverage
# Total reads = reads remaining after trimming
plot(jitter(pver_mapping_stats$number_of_reads), jitter(pver_mapping_stats$coverage_mean), col = add.alpha(colors[factor(species)], 0.8), pch = 19, las = 1, xlab = "Number of reads after trimming", ylab = "Average coverage across P. verrucosa")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Mapped reads vs coverage
plot(pver_mapping_stats$number_of_mapped_reads, pver_mapping_stats$coverage_mean, col = add.alpha(colors[factor(species)], 0.8), pch = 19, las = 1, xlab = "Number of mapped reads", ylab = "Average coverage across P. verrucosa")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Total vs mapped
plot(pver_mapping_stats$number_of_reads, pver_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Number of reads after trimming", ylab = "Number of mapped reads")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Do the earlier sampled Madracis have worse mapping/coverage?
plot(names_split[,5], pver_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Number of mapped reads")
plot(names_split[,5], pver_mapping_stats$coverage_mean, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Average coverage across P.verrucosa")
legend('topleft', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Proportion mapped
f <- aggregate(pver_mapping_stats$proportion_mapping, by = list(species), FUN = mean) # proportion mapped by species
barplot(f$x, col = colors, las = 1, names.arg = f$Group.1, xlab = 'Species ID - Field', main = 'Proportion mapping to P. verrucosa')

#### Coverage and quality data when reads mapped to Breviolum minutum mito genome ####
bmin_mito_mapping_stats <- read.table(here("genotyping/mapping2bmin_mito_stats_headers_lib1.txt"), header = TRUE, sep = '\t')

# Calculate proportion of reads mapping to Breviolum minutum mito
bmin_mito_mapping_stats$proportion_mapping <- bmin_mito_mapping_stats$number_of_mapped_reads/bmin_mito_mapping_stats$number_of_reads

# Clean up names & add to coverage data
names <- str_split_fixed(bmin_mito_mapping_stats$bam, '.fixmate.sorted.markdup.', 2)
names_split <- str_split_fixed(names[,1], '_', 6)

species <- names_split[,1]
site <- names_split[,2]
depth <- names_split[,3]
id_original <- paste(names_split[,4], names_split[,5], sep = "_")
extract <- names_split[,6]

# What is the average coverage across all 96 samples in library 1?
mean(bmin_mito_mapping_stats$coverage_mean) # 0.13X
mean(bmin_mito_mapping_stats$coverage_sd) # 12X
mean(bmin_mito_mapping_stats$number_of_reads) # 9965031
mean(bmin_mito_mapping_stats$number_of_mapped_reads) # 1440123
mean(bmin_mito_mapping_stats$mean_mapping_quality) # 16.16X

# Total reads vs coverage
# Total reads = reads remaining after trimming
plot(jitter(bmin_mito_mapping_stats$number_of_reads), jitter(bmin_mito_mapping_stats$coverage_mean), col = add.alpha(colors[factor(species)], 0.8), pch = 19, las = 1, xlab = "Number of reads after trimming", ylab = "Average coverage across Breviolum minutum mito")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Mapped reads vs coverage
plot(bmin_mito_mapping_stats$number_of_mapped_reads, bmin_mito_mapping_stats$coverage_mean, col = add.alpha(colors[factor(species)], 0.8), pch = 19, las = 1, xlab = "Number of mapped reads", ylab = "Average coverage across Breviolum minutum mito")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Total vs mapped
plot(bmin_mito_mapping_stats$number_of_reads, bmin_mito_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Number of reads after trimming", ylab = "Number of mapped reads")
legend('bottomright', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Do the earlier sampled Madracis have worse mapping/coverage?
plot(names_split[,5], bmin_mito_mapping_stats$number_of_mapped_reads, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Number of mapped reads")
plot(names_split[,5], bmin_mito_mapping_stats$coverage_mean, col = add.alpha(colors[factor(species)], 0.8), pch = 19, xlab = "Sample number", ylab = "Average coverage across Breviolum minutum mito")
legend('topleft', 
       legend = levels(factor(species)),
       pch = 19,
       col = colors)

# Proportion mapped
g <- aggregate(bmin_mito_mapping_stats$proportion_mapping, by = list(species), FUN = mean) # proportion mapped by species
barplot(g$x, col = colors, las = 1, names.arg = g$Group.1, xlab = 'Species ID - Field', main = 'Proportion mapping to Breviolum minutum mito')

