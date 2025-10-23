# Load necessary libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)
library(here)
library(plotrix)
library(scales)
library(car)
library(lattice)
library(multcompView)
library(FSA)
library(multcompView)
library(RColorBrewer)
library(cluster)
library(effectsize)
library(MASS)
library(vegan)
library(adegenet)
library(factoextra)
library(tidyverse)
library(ggridges)

#### Functions ####
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} # error bar function

generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY$label[,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$labels=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$label) , ]
  return(Tukey.labels)
} # function to generate Tukey groups for plotting

#### Read in metadata to check genotyping status ####
metadata <- read.csv(here("madracis_master_merged_Oct2025.csv"), header = TRUE)

#### Read in all habitat measurements & name according to their filename ####
# Set the directory containing the CSV files
directory <- here("3D_models/output/hab_measures")

# List all CSV files with the pattern *_measurements.csv; modify depending on desired radius
csv_files <- list.files(path = directory, pattern = "measurements_quality.csv", full.names = TRUE)

# Function to extract the desired parts of the file name
extract_parts <- function(file_path) {
  file_name <- basename(file_path)
  # Extract the part of the file name (e.g., "est_05" from "est_05m_example_measurements_jah_quality.csv" or "cur_est_05m_measurements.csv")
  #site <- str_extract(file_name, "^[a-z]+")  # Extract 'est' when it is first
  site <- str_match(file_name, "cur_(\\w+?)_")[,2] # When cur is first, extract the second set of letters (est)
  depth <- str_extract(file_name, "\\d+")     # Extract '05'
  return(list(site = site, depth = depth))
}

# Loop through the list of CSV files, read them, and assign to variables with new columns
for (file in csv_files) {
  parts <- extract_parts(file)
  site <- parts$site
  depth <- parts$depth
  object_name <- paste(site, depth, sep = "_")
  
  # Read the CSV file into a data frame
  data <- read.csv(file)
  
  # Add new columns to the data frame
  data$site <- site
  data$depth <- depth
  
  # Assign the data frame to a variable with the extracted name
  assign(object_name, data)
}

# Verify that the objects have been created and they have the correct dimensions; these are for radius 0.2
dim(dir_05) #5x25
dim(est_05) #9x25
dim(est_10) #66x25
dim(est_20) #54x25
dim(est_40) #32x25
dim(est_60) #14x25
dim(kal_05) #13x25
dim(kal_10) #55x25
dim(kal_20) #48x25
dim(kal_40) #19x25
dim(sea_10) #6x25
dim(sea_20) #12x25
dim(sna_05) #10x25
dim(sna_10) #15x25
dim(sna_20) #46x25
dim(sna_40) #5x25

################
#### Part 1 ####
################

#### Combine the separate objects into one, manipulate and filter ####
# Data set 1 for est, kal & sna 20m only
coral_ann_eks20 <- rbind(est_20, kal_20, sna_20) #148 x 25; only for est, kal & sna 20m supplement figure

# Add labels from metadata
# Merge metadata$Species_mpha into coral_ann by matching 'id' to 'ID'
coral_ann_eks20.2 <- merge(coral_ann_eks20, metadata[, c("ID", "Species_mpha")],
                    by.x = "id", by.y = "ID", all.x = TRUE)

# Update the 'label' column with Species_mpha
coral_ann_eks20.2$label <- coral_ann_eks20.2$Species_mpha

mad_ann_eks20 <- coral_ann_eks20.2[!coral_ann_eks20.2$label %in% c('PPOR', 'SINT'),] # 138x26 remove corals that aren't Madracis

mad_ann_seq_eks20 <- mad_ann_eks20[mad_ann_eks20$id %in% (metadata[which(metadata$Davis_library == TRUE),5]),] # work with only genotyped Madracis, 123x26

# Make sure any clones are removed
Madracis606 <- read.csv(here("WGS/angsd_results/bam_lists/madracis606_names.csv"), header = FALSE) # list of unique ids (no clones)
Madracis606_split <- as.data.frame(str_split_fixed(Madracis606$V1, '_', 6))
Madracis606_names <- paste(Madracis606_split$V4, Madracis606_split$V5, sep = "_")

mad_ann_seq_eks20$id %in% Madracis606_names # identify potential clones and carefully check 06a_merge_sample_bams.sh to see if one or both are in annotation dataset. If 2, remove 1.

mad_ann_seq_eks20[!mad_ann_seq_eks20$id %in% Madracis606_names,]$id # potential clones. Clone pair can be seen in 06a_merge_sample_bams.sh
mad_ann_seq_eks20[mad_ann_seq_eks20$id %in% Madracis606_names,]$id

clones <- c('JH_0089', 'JH_0212', 'JH_0221', 'JH_0306', 'JH_0381', 'JH_0387', 'JH_0475', 'JH_0483', 'JH_0694', 'JH_0696') # clones in est, kal & sna all depths

mad_ann_seq_eks20.2 <- mad_ann_seq_eks20[!mad_ann_seq_eks20$id %in% clones, ] # remove the clones

mad_ann_seq_eks20_core <- mad_ann_seq_eks20.2[!mad_ann_seq_eks20.2$id %in% c('JH_0306', 'JH_0307', 'JH_0165', 'JH_0233', 'JH_0013', 'JH_0052', 'JH_0586'),] # remove edge cases for exploration (some already removed in filters above), 110 x 26 for est, kal & sna 20m only

# See if the unknowns have a more updated ID. Should be empty.
mad_ann_seq_eks20_core[which(mad_ann_seq_eks20_core$label == 'UNKN'),1]
metadata[metadata$ID %in% mad_ann_seq_eks20_core[which(mad_ann_seq_eks20_core$label == 'UNKN'),1],]

# Change depth categories into integers
mad_ann_seq_eks20_core$depth <- as.integer(mad_ann_seq_eks20_core$depth)

# Look at means
aggregate(x=mad_ann_seq_eks20_core$gapF_fill,
          by=list(mad_ann_seq_eks20_core$label,mad_ann_seq_eks20_core$depth, mad_ann_seq_eks20_core$site),
          FUN=mean)

#### For est, kal and sna 20m only
plot_boxplot_facet_anova <- function(df, y_vars) {
  # Reshape data to long format
  df_long <- df %>%
    pivot_longer(cols = all_of(y_vars), names_to = "measure", values_to = "value")
  
  # Tukey test results
  tukey_results <- df_long %>%
    group_by(site, depth, measure) %>%
    group_split() %>%
    lapply(function(sub_df) {
      if (length(unique(sub_df$label)) < 2) return(NULL)
      anova_model <- aov(value ~ label, data = sub_df)
      if (length(coef(anova_model)) < 2) return(NULL)
      
      num_groups <- length(unique(sub_df$label))
      if (num_groups < 2) return(NULL)
      
      tukey_test <- tryCatch(TukeyHSD(anova_model), error = function(e) NULL)
      if (is.null(tukey_test) || is.null(tukey_test$label)) return(NULL)
      
      tukey_values <- tukey_test$label[, "p adj"]
      valid_names <- grep("-", names(tukey_values), value = TRUE)
      tukey_values <- tukey_values[valid_names]
      if (length(tukey_values) == 0) return(NULL)
      
      tukey_letters <- multcompLetters(tukey_values < 0.05)$Letters
      data.frame(label = names(tukey_letters),
                 Letters = tukey_letters,
                 site = unique(sub_df$site),
                 depth = unique(sub_df$depth),
                 measure = unique(sub_df$measure))
    }) %>%
    bind_rows()
  
  # Max y for letters
  y_max <- df_long %>%
    group_by(site, depth, label, measure) %>%
    summarize(max_y = max(value, na.rm = TRUE), .groups = "drop")
  
  tukey_results <- left_join(tukey_results, y_max, 
                             by = c("site", "depth", "label", "measure"))
  
  # Plot
  ggplot(df_long, aes(x = label, y = value, fill = label)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_grid(measure ~ site, scales = "free_y", switch = "y") +  # Removed depth facetting
    geom_text(data = tukey_results, 
              aes(x = label, y = max_y * 1.15, label = Letters),  # Lowered letters
              inherit.aes = FALSE, color = "black", size = 4, fontface = "bold", na.rm = TRUE) +
    theme_minimal() +
    labs(y = NULL, x = "Genetic lineage", fill = 'lineage') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(margin = margin(r = 10)), # space for y-axis title
      panel.grid = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_text(face = "bold"), # y-axis facet labels (measure)
      strip.placement = "outside", # moves measure labels to outside left
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    scale_fill_manual(values = c("#ff9a81", "#E6AB02", 
                                 "#b6c590", "#6b8e23", "#b0e0e6", "gray70"))
}

png(file=here("3D_models/figures/VD_GF_SA180_rad0.2_facet_anova_20m.png"), width=8, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_eks20_core, y_vars = c("gapF_fill", "elevation", "vector_disp"))

dev.off()

################
#### Part 2 ####
################
#### Aggregate data across sites by depth, with an option for ANOVA ####
# Data set 2 for est, kal & sna all depths
coral_ann <- rbind(est_05, est_10, est_20, est_40, est_60, kal_05, kal_10, kal_20, kal_40, sna_10, sna_20, sna_40) #376 x 25

# Add labels from metadata
# Merge metadata$Species_mpha into coral_ann by matching 'id' to 'ID'
coral_ann2 <- merge(coral_ann, metadata[, c("ID", "Species_mpha")],
                    by.x = "id", by.y = "ID", all.x = TRUE)

# Update the 'label' column with Species_mpha
coral_ann2$label <- coral_ann2$Species_mpha

mad_ann <- coral_ann2[!coral_ann2$label %in% c('PPOR', 'SINT'),] #363x26 remove corals that aren't Madracis, change input data set depending on elevation range

mad_ann_seq <- mad_ann[mad_ann$id %in% (metadata[which(metadata$Davis_library == TRUE),5]),] # work with only genotyped Madracis, 308x26

# Make sure any clones are removed
Madracis606 <- read.csv(here("WGS/angsd_results/bam_lists/madracis606_names.csv"), header = FALSE) # list of unique ids (no clones)
Madracis606_split <- as.data.frame(str_split_fixed(Madracis606$V1, '_', 6))
Madracis606_names <- paste(Madracis606_split$V4, Madracis606_split$V5, sep = "_")

mad_ann_seq$id %in% Madracis606_names # identify potential clones and carefully check 06a_merge_sample_bams.sh to see if one or both are in annotation dataset. If 2, remove 1.

mad_ann_seq[!mad_ann_seq$id %in% Madracis606_names,]$id # potential clones. Clone pair can be seen in 06a_merge_sample_bams.sh
mad_ann_seq[mad_ann_seq$id %in% Madracis606_names,]$id

clones <- c('JH_0089', 'JH_0212', 'JH_0221', 'JH_0306', 'JH_0381', 'JH_0387', 'JH_0475', 'JH_0483', 'JH_0694', 'JH_0696') # clones in est, kal & sna all depths

mad_ann_seq2 <- mad_ann_seq[!mad_ann_seq$id %in% clones, ] # remove the clones

mad_ann_seq_core <- mad_ann_seq2[!mad_ann_seq2$id %in% c('JH_0306', 'JH_0307', 'JH_0165', 'JH_0233', 'JH_0013', 'JH_0052', 'JH_0586'),] # remove edge cases for exploration (some already removed in filters above), 292x26

# See if the unknowns have a more updated ID. Should be empty.
mad_ann_seq_core[which(mad_ann_seq_core$label == 'UNKN'),1]
metadata[metadata$ID %in% mad_ann_seq_core[which(mad_ann_seq_core$label == 'UNKN'),1],]

# Change depth categories into integers
mad_ann_seq_core$depth <- as.integer(mad_ann_seq_core$depth)

# Reduce to only MDEC, MPHA1, MPHA2 & MSEN 
mad_ann_seq_core_reduced <- mad_ann_seq_core[mad_ann_seq_core$label %in% c('MDEC', 'MPHA.ln1', 'MPHA.ln2', 'MSEN.ln1'),] # 286 x 26
mad_ann_seq_core_reduced2 <- mad_ann_seq_core_reduced[-c(180,190),] # remove MSEN from 5 (row 180) and 40m (row 190) #284 x 26

# Plotting function for multiple y variables
plot_boxplot_anova_multi_pub <- function(df, y_vars) {
  
  # Reshape to long format
  df_long <- df %>%
    pivot_longer(cols = all_of(y_vars), names_to = "variable", values_to = "y_value")
  
  # Perform t-test or ANOVA + Tukey for each variable and depth, as appropriate
  tukey_results <- df_long %>%
    group_by(variable, depth) %>%
    group_split() %>%
    lapply(function(sub_df) {
      n_groups <- length(unique(sub_df$label))
      
      # Case 1: fewer than 2 groups → skip
      if (n_groups < 2) return(NULL)
      
      # Case 2: exactly 2 groups → run t-test
      if (n_groups == 2) {
        ttest <- t.test(y_value ~ label, data = sub_df)
        sig <- ttest$p.value < 0.05
        labs <- if (sig) c("a", "b") else c("a", "a")
        names(labs) <- unique(sub_df$label)
        return(data.frame(
          label = names(labs),
          Letters = labs,
          depth = unique(sub_df$depth),
          variable = unique(sub_df$variable)
        ))
      }
      
      # Case 3: 3+ groups → ANOVA + Tukey
      anova_model <- aov(y_value ~ label, data = sub_df)
      tukey_test <- tryCatch(
        TukeyHSD(anova_model),
        error = function(e) NULL
      )
      if (is.null(tukey_test)) return(NULL)
      
      tukey_values <- tukey_test$label[, "p adj"]
      valid_names <- grep("-", names(tukey_values), value = TRUE)
      tukey_values <- tukey_values[valid_names]
      if (length(tukey_values) == 0) return(NULL)
      
      tukey_letters <- multcompLetters(tukey_values < 0.05)$Letters
      
      data.frame(
        label = names(tukey_letters),
        Letters = tukey_letters,
        depth = unique(sub_df$depth),
        variable = unique(sub_df$variable)
      )
    }) %>%
    bind_rows()
  
  y_max <- df_long %>%
    group_by(variable, depth, label) %>%
    summarize(max_y = max(y_value, na.rm = TRUE) * 1.05, .groups = "drop")
  
  tukey_results <- left_join(tukey_results, y_max,
                             by = c("variable", "depth", "label"))
  
  ggplot(df_long, aes(x = label, y = y_value, fill = label)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_grid(variable ~ depth, scales = "free", space = "free_x", switch = "y") +
    geom_text(data = tukey_results, 
              aes(x = label, y = max_y * 1.05, label = Letters),
              color = "black", size = 4, fontface = "bold", na.rm = TRUE) +
    theme_minimal() +
    labs(y = NULL, x = "Genetic lineage", fill = 'lineage') + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # tick labels
      axis.text.y = element_text(size = 12),
      panel.grid = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold"),
      strip.placement = "outside",
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title.x = element_text(size = 14, face = "bold"),          # axis titles
      axis.title.y = element_text(size = 14, face = "bold"),
      strip.text.x = element_text(size = 14, face = "bold"),  # facet labels along x
      strip.text.y = element_text(size = 14, face = "bold")   # facet labels along y
    ) +
    scale_fill_manual(values = c("#ff9a81", "#b6c590", 
                                 "#6b8e23", "#b0e0e6", "gray70"))
}

# Nice plot summarizing habitat measures across depths and genetic lineages
png(file=here("3D_models/figures/VD_GF_SA180_bydepth_eks_anova_simple.png"), width=10, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova_multi_pub(mad_ann_seq_core_reduced2, c("vector_disp", "gapF_fill", "elevation"))

dev.off()
