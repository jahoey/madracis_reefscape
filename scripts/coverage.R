library(here)
library(tidyr)
library(viridis)
library(stringr)

add.alpha <- function(col, alpha=1){
        if(missing(col))
                stop("Please provide a vector of colours.")
        apply(sapply(col, col2rgb)/255, 2, 
              function(x) 
                      rgb(x[1], x[2], x[3], alpha=alpha))  
}

#### Coverage and quality data when reads mapped to M. senaria de novo genome ####
msen_denovo_mapping_stats <- read.table(here("genotyping/mapping2msen_denovo_stats_headers_lib1.txt"), header = TRUE, sep = '\t')

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

# Set color scheme
col <- viridis(6)

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


temp <- data.frame(cbind(as.numeric(names_split[,5]), as.numeric(msen_denovo_mapping_stats$coverage_mean)))

