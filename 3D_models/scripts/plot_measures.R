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
metadata <- read.csv(here("madracis_master_merged_Oct2024.csv"), header = TRUE)

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
# dim(dir_20)
dim(est_05) #9x25
dim(est_10) #66x25
dim(est_20) #54x25
dim(est_40) #32x25
dim(est_60) #14x25
#dim(hul_05)
#dim(hul_20) #21x25
dim(kal_05) #13x25
dim(kal_10) #55x25
dim(kal_20) #48x25
dim(kal_40) #19x25
dim(sea_10) #6x25
dim(sea_20) #12x25
# dim(sna_05) #10x25
dim(sna_10) #15x25
dim(sna_20) #46x25
dim(sna_40) #5x25
#dim(wat_05) #2x25
#dim(wat_20)

#### Combine the separate objects into one, manipulate and filter ####
# coral_ann <- rbind(dir_05, est_05, est_10, est_20, est_40, est_60, hul_20, kal_05, kal_10, kal_20, kal_40, sea_10, sea_20, sna_10, sna_20, sna_40, wat_05) #422x24

coral_ann <- rbind(dir_05, est_05, est_10, est_20, est_40, est_60, kal_05, kal_10, kal_20, kal_40, sea_10, sea_20, sna_10, sna_20, sna_40) #399 x 26

# Add labels from metadata
# Merge metadata$Species_mpha into coral_ann by matching 'id' to 'ID'
coral_ann2 <- merge(coral_ann, metadata[, c("ID", "Species_mpha")],
                    by.x = "id", by.y = "ID", all.x = TRUE)

# Update the 'label' column with Species_mpha
coral_ann2$label <- coral_ann2$Species_mpha

# # Manipulate elevation/slope angle by subtracting any angles greater than 90 degrees by 180. No longer needed as of April 2025, since we now want to know if coral is angled up or angled down.
# coral_ann3 <- coral_ann2 %>%
#   mutate(elevation = if_else(elevation > 90, 180 - elevation, elevation))

mad_ann <- coral_ann2[!coral_ann2$label %in% c('PPOR', 'SINT'),] #386x26, remove corals that aren't Madracis

mad_ann_seq <- mad_ann[mad_ann$id %in% (metadata[which(metadata$Davis_library == TRUE),5]),] # work with only genotyped Madracis, 331x26

mad_ann_seq_core <- mad_ann_seq[!mad_ann_seq$id %in% c('JH_0306', 'JH_0307', 'JH_0165', 'JH_0233', 'JH_0013', 'JH_0052', 'JH_0586'),] # remove edge cases for exploration, 324x26

# See if the unknowns have a more updated ID. Should be empty.
mad_ann_seq_core[which(mad_ann_seq_core$label == 'UNKN'),1]
metadata[metadata$ID %in% mad_ann_seq_core[which(mad_ann_seq_core$label == 'UNKN'),1],]

# Change depth categories into integers
mad_ann_seq_core$depth <- as.integer(mad_ann_seq_core$depth)

# Look at means
aggregate(x=mad_ann_seq_core$gapF_fill,
          by=list(mad_ann_seq_core$label,mad_ann_seq_core$depth, mad_ann_seq_core$site),
          FUN=mean)

###############################
#### Univariate approaches ####
###############################

#Quick facet plot for each unique site/depth combination. There is some variation in the pattern between sites
# plot_boxplot_facet <- function(df, y_var) {
#   ggplot(df, aes(x = label, y = .data[[y_var]], fill = label)) +
#     geom_boxplot(alpha = 1, outlier.shape = NA) +
#     geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
#     facet_grid(depth ~ site) +  # Facet by site and depth
#     theme_minimal() +
#     labs(y = y_var, title = paste(y_var), x = "Genetic lineage") +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid = element_blank(),  # Remove all grid lines
#       panel.background = element_blank(),  # Remove background color
#       strip.background = element_blank(),  # Remove facet background
#       strip.text = element_text(face = "bold"),  # Keep facet text bold
#       axis.line = element_line(color = "black"),  # Show x and y axes
#       axis.ticks = element_line(color = "black"), # Show axis ticks
#       panel.border = element_rect(color = "black", fill = NA) 
#     ) +
#     scale_fill_manual(values = c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", 
#                                  "#6b8e23", "#79d2a6", "#b0e0e6", "gray70"))
# }

# Facet plot with ANOVA for each unique site/depth
plot_boxplot_facet_anova <- function(df, y_var) {
  tukey_results <- df %>%
    group_by(site, depth) %>%
    group_split() %>%
    lapply(function(sub_df) {
      # Skip if there are fewer than two unique labels
      if (length(unique(sub_df$label)) < 2) {
        return(NULL)
      }
      
      # Perform ANOVA
      anova_model <- aov(sub_df[[y_var]] ~ sub_df$label)
      
      # Skip if ANOVA couldn't run
      if (length(coef(anova_model)) < 2) {
        return(NULL)
      }
      
      # Check if there are enough levels for Tukey test
      num_groups <- length(unique(sub_df$label))
      if (num_groups < 2) {
        return(NULL)
      }
      
      # Perform Tukey’s HSD test
      tukey_test <- tryCatch(
        TukeyHSD(anova_model),
        error = function(e) NULL
      )
      
      # Skip if Tukey test failed or if there are fewer than 2 groups
      if (is.null(tukey_test) || nrow(tukey_test$`sub_df$label`) < 1) return(NULL)
      
      # Extract p-values and perform multiple comparison letters
      tukey_values <- tukey_test$`sub_df$label`[, "p adj"]
      valid_names <- grep("-", names(tukey_values), value = TRUE)  # Ensure valid comparisons
      tukey_values <- tukey_values[valid_names]
      
      # Skip if no valid comparisons
      if (length(tukey_values) == 0) return(NULL)
      
      # Assign letters based on Tukey test significance
      tukey_letters <- multcompLetters(tukey_values < 0.05)$Letters
      
      # Return formatted results
      data.frame(label = names(tukey_letters), Letters = tukey_letters, 
                 site = unique(sub_df$site), depth = unique(sub_df$depth))
    }) %>%
    bind_rows()
  
  # Compute maximum y-values for each group to position Tukey letters properly
  y_max <- df %>%
    group_by(site, depth, label) %>%
    summarize(max_y = max(.data[[y_var]], na.rm = TRUE) * 1.05, .groups = "drop")
  
  # Merge Tukey results with max y-values for plotting
  tukey_results <- left_join(tukey_results, y_max, by = c("site", "depth", "label"))
  
  # Create boxplot with Tukey letters
  ggplot(df, aes(x = label, y = .data[[y_var]], fill = label)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_grid(depth ~ site) +  # Facet by site and depth
    geom_text(data = tukey_results, aes(x = label, y = max_y *1.2, label = Letters),
              color = "black", size = 4, fontface = "bold", na.rm = TRUE) +  # Tukey letters
    theme_minimal() +
    labs(y = y_var, title = paste(y_var), x = "Genetic lineage") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),  # Remove grid lines
      panel.background = element_blank(),  # Remove background color
      strip.background = element_blank(),  # Remove facet background
      strip.text = element_text(face = "bold"),  # Keep facet text bold
      axis.line = element_line(color = "black"),  # Show axes
      axis.ticks = element_line(color = "black"), # Show axis ticks
      panel.border = element_rect(color = "black", fill = NA)  # Add borders
    ) +
    scale_fill_manual(values = c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", 
                                 "#6b8e23", "#79d2a6", "#b0e0e6", "gray70"))
}

# Plot each of the measurements and save
png(file=here("3D_models/figures/dev_rug_rad0.2_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "dev_rug")

dev.off()

png(file=here("3D_models/figures/gapF_fill_rad0.2_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "gapF_fill")

dev.off()

png(file=here("3D_models/figures/gapF_raw_rad0.2_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "gapF_raw")

dev.off()

png(file=here("3D_models/figures/luminance_rad0.2_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "luminance")

dev.off()

png(file=here("3D_models/figures/vector_disp_rad0.2_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "vector_disp")

dev.off()

png(file=here("3D_models/figures/substrate_angle_rad0.2_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "elevation")

dev.off()

# Aggregate data across sites by depth, with an option for ANOVA
# plot_boxplot_by_depth <- function(df, y_var) {
#   ggplot(df, aes(x = label, y = .data[[y_var]], fill = label)) +
#     geom_boxplot(alpha = 1, outlier.shape = NA) +
#     geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
#     facet_wrap(~depth, nrow = 1) +  # Facet by depth only, scales = "free_y" lets y-axis vary
#     theme_minimal() +
#     labs(y = y_var, title = paste(y_var), x = "Genetic lineage") +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid = element_blank(),  # Remove all grid lines
#       panel.background = element_blank(),  # Remove background color
#       strip.background = element_blank(),  # Remove facet background
#       strip.text = element_text(face = "bold"),  # Keep facet text bold
#       axis.line = element_line(color = "black"),  # Show x and y axes
#       axis.ticks = element_line(color = "black"), # Show axis ticks
#       panel.border = element_rect(color = "black", fill = NA)  # Add borders
#     ) +
#     scale_fill_manual(values = c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", 
#                                  "#6b8e23", "#79d2a6", "#b0e0e6", "gray70"))
# }

plot_boxplot_anova <- function(df, y_var) {
  tukey_results <- df %>%
    group_by(depth) %>%
    group_split() %>%
    lapply(function(sub_df) {
      # Perform ANOVA only if there are at least 2 unique labels
      if (length(unique(sub_df$label)) < 2) {
        return(NULL)  # Skip ANOVA if only one group exists
      }
      
      anova_model <- aov(sub_df[[y_var]] ~ sub_df$label)
      
      # Perform Tukey’s HSD test
      tukey_test <- tryCatch(
        TukeyHSD(anova_model),
        error = function(e) NULL  # Handle error cases gracefully
      )
      
      # Skip if Tukey test failed
      if (is.null(tukey_test)) return(NULL)
      
      # Extract p-values safely
      tukey_values <- tukey_test$`sub_df$label`[, "p adj"]
      valid_names <- grep("-", names(tukey_values), value = TRUE)  # Ensure valid comparisons
      tukey_values <- tukey_values[valid_names]
      
      # Skip if no valid comparisons
      if (length(tukey_values) == 0) return(NULL)
      
      # Assign letters
      tukey_letters <- multcompLetters(tukey_values < 0.05)$Letters
      
      # Return formatted dataframe
      data.frame(label = names(tukey_letters), Letters = tukey_letters, depth = unique(sub_df$depth))
    }) %>%
    bind_rows()
  
  # Compute max y values per depth to position letters properly
  y_max <- df %>%
    group_by(depth, label) %>%
    summarize(max_y = max(.data[[y_var]], na.rm = TRUE) * 1.05, .groups = "drop")
  
  # Merge Tukey results with max y values
  tukey_results <- left_join(tukey_results, y_max, by = c("depth", "label"))
  
  # Plot boxplots with significance letters
  ggplot(df, aes(x = label, y = .data[[y_var]], fill = label)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_wrap(~depth, nrow = 1) +  # Facet by depth in a single row
    geom_text(data = tukey_results, aes(x = label, y = max_y * 1.05, label = Letters),
              color = "black", size = 4, fontface = "bold", na.rm = TRUE) +  # Tukey letters
    theme_minimal() +
    labs(y = y_var, x = "Genetic lineage") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),  # Remove all grid lines
      panel.background = element_blank(),  # Remove background color
      strip.background = element_blank(),  # Remove facet background
      strip.text = element_text(face = "bold"),  # Keep facet text bold
      axis.line = element_line(color = "black"),  # Show x and y axes
      axis.ticks = element_line(color = "black"), # Show axis ticks
      panel.border = element_rect(color = "black", fill = NA)  # Add borders
    ) +
    scale_fill_manual(values = c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", 
                                 "#6b8e23", "#79d2a6", "#b0e0e6", "gray70"))
}

# Plot each of the measurements and save
png(file=here("3D_models/figures/dev_rug_bydepth_anova.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "dev_rug")

dev.off()

png(file=here("3D_models/figures/gapF_fill_bydepth_anova.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "gapF_fill")

dev.off()

png(file=here("3D_models/figures/gapF_raw_bydepth_anova.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
plot_boxplot_anova(mad_ann_seq_core, "gapF_raw")

dev.off()

png(file=here("3D_models/figures/luminance_bydepth_anova.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "luminance")

dev.off()

png(file=here("3D_models/figures/vector_disp_bydepth_anova.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "vector_disp")

dev.off()

png(file=here("3D_models/figures/substrate_angle_bydepth_anova.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "elevation")

dev.off()

#################################
#### Multivariate approaches ####
#################################
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#79d2a6", "#b0e0e6")
show_col(colors_mpha)
colors_depth <- c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
show_col(colors_depth)
colors_site <- c("#B8DE29FF", "#74D055FF", "#20A386FF", "#238A8DFF", "#2D718EFF", "#453781FF", "#481568FF")
show_col(colors_site)

# MANOVA
dependent_vars <- cbind(mad_ann_seq_core$dev_rug, mad_ann_seq_core$elevation, mad_ann_seq_core$gapF_fill, mad_ann_seq_core$luminance, mad_ann_seq_core$vector_disp)
independent_vars <- cbind.data.frame(mad_ann_seq_core$label, mad_ann_seq_core$site, mad_ann_seq_core$depth)

manova_results <- manova(cbind(dev_rug, elevation, gapF_fill, luminance, vector_disp) ~ label + site + depth, data = mad_ann_seq_core)
summary(manova_results)
summary.aov(manova_results)
eta_squared(manova_results) #depth has the biggest effect size

# MANOVA post-hoc test (LDA) for each independent variable
png(file=here("3D_models/figures/manova_posthoc_lda.png"), width=8, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow=c(2,2)
)

# Depth
depth_lda <- lda(mad_ann_seq_core$depth ~ dependent_vars, CV = F)
depth_lda_df <- data.frame(depth = mad_ann_seq_core[, "depth"],
                     lda = predict(depth_lda)$x)

plot(depth_lda_df$lda.LD1, depth_lda_df$lda.LD2, col = add.alpha(colors_depth[factor(depth_lda_df$depth)], 0.7), pch = 19, xlab = 'LD1', ylab = 'LD2', main = 'Depth')
legend('topleft', legend = c('5', '10', '20', '40', '60'), col = add.alpha(colors_depth, 0.7), pch = 19)

# Species
species_lda <- lda(mad_ann_seq_core$label ~ dependent_vars, CV = F)
species_lda_df <- data.frame(label = mad_ann_seq_core[, "label"],
                           lda = predict(species_lda)$x)

plot(species_lda_df$lda.LD1, species_lda_df$lda.LD2, col = add.alpha(colors_mpha[factor(species_lda_df$label)], 0.7), pch = 19, xlab = 'LD1', ylab = 'LD2', main = 'Species', xlim = c(-5,4))
legend('topleft', legend = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA.ln1", "MPHA.ln2", "MPHA.ln4", "MSEN"), col = add.alpha(colors_mpha, 0.7), pch = 19, cex = 0.8)

# Site
site_lda <- lda(mad_ann_seq_core$site ~ dependent_vars, CV = F)
site_lda_df <- data.frame(site = mad_ann_seq_core[, "site"],
                             lda = predict(site_lda)$x)

plot(site_lda_df$lda.LD1, site_lda_df$lda.LD2, col = add.alpha(colors_site[factor(site_lda_df$site)], 0.7), pch = 19, xlab = 'LD1', ylab = 'LD2', main = 'Site')
legend('bottomright', legend = c("dir", "est", "hul", "kal", "sea", "sna", "wat"), col = add.alpha(colors_site, 0.7), pch = 19)

dev.off()

# MDS analysis
dist_matrix <- dist(mad_ann_seq_core[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")], method = "euclidean")

mds_mad <- cmdscale(dist_matrix, eig = TRUE)

mds_distances <- mds_mad$points
shepard_df <- data.frame(Original = as.vector(dist_matrix), MDS = as.vector(dist(mds_distances)))

ggplot(shepard_df, aes(x = Original, y = MDS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Shepard Plot: Original vs. MDS Distances") +
  theme_minimal() # checks linearity, linear correlation so MDS is appropriate

# Fit environmental variables to the MDS space
env_fit <- envfit(mds_mad$points, mad_ann_seq_core[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")], permutations = 999)

# Extract arrow scores
arrow_scores <- as.data.frame(scores(env_fit, display = "vectors"))  # Extract coordinates
arrow_scores$variable <- rownames(arrow_scores)  # Add variable names

# Filter for significant vetors
significant_arrows <- arrow_scores[env_fit$vectors$pvals < 0.05, ]

# All colored by species
# Interpretation: elevation is perfectly correlated with MDS1, gapF_fill is perfectly correlated with MDS2, rest of the variables have low correlation with both axes
png(file=here("3D_models/figures/mds_all_species.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 10), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  xpd = TRUE
)

plot(mds_mad$points[,1], mds_mad$points[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_mpha[factor(mad_ann_seq_core$label)], 0.7), pch = 19, ylim = c(-0.5, 1.2))
arrows(0, 0, significant_arrows$Dim1, significant_arrows$Dim2, length = 0.1, col = "black", lwd = 1.5)
text(significant_arrows$Dim1, significant_arrows$Dim2, labels = significant_arrows$variable, pos = 3, cex = 0.8)
legend('topright', inset = c(-0.3,0), legend = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA.ln1", "MPHA.ln2", "MPHA.ln4", "MSEN"), col = add.alpha(colors_mpha, 0.7), pch = 19, title = 'Species')

dev.off()

# All colored by depth
png(file=here("3D_models/figures/mds_all_depth.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 6), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  xpd = TRUE
)

plot(mds_mad$points[,1], mds_mad$points[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_depth[factor(mad_ann_seq_core$depth)], 0.7), pch = 19, ylim = c(-0.5, 1.2))
arrows(0, 0, significant_arrows$Dim1, significant_arrows$Dim2, length = 0.1, col = "black", lwd = 1.5)
text(significant_arrows$Dim1, significant_arrows$Dim2, labels = significant_arrows$variable, pos = 3, cex = 0.8)
legend('topright', inset = c(-0.15,0), legend = c('5', '10', '20', '40', '60'), col = add.alpha(colors_depth, 0.7), pch = 19, title = 'Depth')

dev.off()

# By depth
mad5 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 5),]
mad10 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 10),]
mad20 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 20),]
mad40 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 40),]
mad60 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 60),]

mds_mad5 <- cmdscale(dist(mad5[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")]), k = 2)
mds_mad10 <- cmdscale(dist(mad10[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")]), k = 2)
mds_mad20 <- cmdscale(dist(mad20[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")]), k = 2)
mds_mad40 <- cmdscale(dist(mad40[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")]), k = 2)
mds_mad60 <- cmdscale(dist(mad60[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")]), k = 2)

png(file=here("3D_models/figures/mds_by_depth.png"), width=8, height=12, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow=c(3,2)
)

plot(mds_mad5[,1], mds_mad5[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_mpha[-c(1,3,4,6,7)][factor(mad5$label)], 0.7), pch = 19, main = '5m')
plot(mds_mad10[,1], mds_mad10[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_mpha[-c(3,4,7)][factor(mad10$label)], 0.7), pch = 19, main = '10m')
plot(mds_mad20[,1], mds_mad20[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_mpha[-c(1,3,7)][factor(mad20$label)], 0.7), pch = 19, main = '20m')
plot(mds_mad40[,1], mds_mad40[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_mpha[-c(1,2,4,5)][factor(mad40$label)], 0.7), pch = 19, main = '40m')
plot(mds_mad60[,1], mds_mad60[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_mpha[-c(1:5)][factor(mad60$label)], 0.7), pch = 19, main = '60m')
plot.new()
legend('topleft', legend = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA.ln1", "MPHA.ln2", "MPHA.ln4", "MSEN"), col = add.alpha(colors_mpha, 0.7), pch = 19)

dev.off()

# NMDS
dist_nmds <- vegdist(mad_ann_seq_core[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")], method = "bray")
nmds_results <- metaMDS(dist_nmds, k = 2, trymax = 100)
nmds_results$stress # check stress value: < 0.05 is excellent
plot(nmds_results, type = "t")
nmds_scores <- as.data.frame(scores(nmds_results))
nmds_scores$depth <- mad_ann_seq_core$depth

ggplot(nmds_scores, aes(NMDS1, NMDS2, color = depth)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS Ordination by depth")

# PCA analysis: correlation matrix
pca_mad <- prcomp(mad_ann_seq_core[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")], center = TRUE, scale. = TRUE) # correlation matrix pca

summ <- summary(pca_mad) # calculates variance
summ$importance

png(file=here("3D_models/figures/pca.png"), width=8, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow = c(2,2)
)

plot(pca_mad$x, col = add.alpha(colors_depth[factor(mad_ann_seq_core$depth)], 0.7), pch = 19, xlab = 'PC 1 (39.9%)', ylab = 'PC 2 (22.2%)', main = 'Depth')
legend('topleft', legend = c('5', '10', '20', '40', '60'), col = add.alpha(colors_depth, 0.7), pch = 19, cex = 0.9)

plot(pca_mad$x, col = add.alpha(colors_mpha[factor(mad_ann_seq_core$label)], 0.7), pch = 19, xlab = 'PC 1 (39.9%)', ylab = 'PC 2 (22.2%)', main = 'Species')
legend('topleft', legend = c("MCAR", "MDEC", "MFOR", "MMIR", "MPHA.ln1", "MPHA.ln2", "MPHA.ln4", "MSEN"), col = add.alpha(colors_mpha, 0.7), pch = 19, cex = 0.6)

plot(pca_mad$x, col = add.alpha(colors_site[factor(mad_ann_seq_core$site)], 0.7), pch = 19, xlab = 'PC 1 (39.9%)', ylab = 'PC 2 (22.2%)', main = 'Sites')
legend('topleft', legend = c("dir", "est", "hul", "kal", "sea", "sna", "wat"), col = add.alpha(colors_site, 0.7), pch = 19, cex = 0.7)

dev.off()

# DAPC analysis by species
# Step 1: scale data & perform PCA
model3d.data.scaled <- scale(mad_ann_seq_core[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")])
model3d.dapc.species <- dapc(model3d.data.scaled, mad_ann_seq_core$label)# by species

# Choose number of PCs to retain
a_score <- optim.a.score(model3d.dapc.species) 

# Step 2: perform DAPC
model3d.dapc.species <- dapc(model3d.data.scaled, mad_ann_seq_core$label, n.pca = 3, n.da = 2) # slightly more separation when fewer pcs are retained

# Step 3: plot
scatter(model3d.dapc.species)        # Visualizes the DAPC space
compoplot(model3d.dapc.species)      # Shows posterior probabilities of group assignment

# DAPC analysis by depth
# Step 1: scale data & perform PCA
model3d.data.scaled <- scale(mad_ann_seq_core[, c("dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")])
model3d.dapc.depth <- dapc(model3d.data.scaled, mad_ann_seq_core$depth)# by depth

# Choose number of PCs to retain
a_score <- optim.a.score(model3d.dapc.depth) 

# Step 2: perform DAPC
model3d.dapc.depth <- dapc(model3d.data.scaled, mad_ann_seq_core$depth, n.pca = 3, n.da = 2)

# Step 3: plot
scatter(model3d.dapc.depth)        # Visualizes the DAPC space
compoplot(model3d.dapc.depth)      # Shows posterior probabilities of group assignment

# Plot a priori DAPCs for species and depth together
png(file=here("3D_models/figures/dapc.png"), width=9, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow = c(1,2)
)

scatter(model3d.dapc.species)  
scatter(model3d.dapc.depth)  

dev.off()

#######################################################
#### Multivariate for MPHA1 vs MPHA2 vs MPHA4 only ####
#######################################################
mpha_data <- mad_ann_seq_core[mad_ann_seq_core$label %in% c('MPHA.ln1', 'MPHA.ln2', 'MPHA.ln4'), c("label", "depth", "site", "dev_rug", "elevation", "gapF_fill", "luminance", "vector_disp")] # subsets data to only MPHA1, MPHA2 and MPHA4; 173 x 8

# MDS analysis
dist_matrix <- dist(mpha_data[,-c(1:3)], method = "euclidean")

mds_mpha <- cmdscale(dist_matrix, eig = TRUE)

mds_distances <- mds_mpha$points
shepard_df <- data.frame(Original = as.vector(dist_matrix), MDS = as.vector(dist(mds_distances)))

ggplot(shepard_df, aes(x = Original, y = MDS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Shepard Plot: Original vs. MDS Distances") +
  theme_minimal() # checks linearity, linear correlation so MDS is appropriate

# Fit environmental variables to the MDS space
env_fit <- envfit(mds_mpha$points, mpha_data[,-c(1:3)], permutations = 999)

# Extract arrow scores
arrow_scores <- as.data.frame(scores(env_fit, display = "vectors"))  # Extract coordinates
arrow_scores$variable <- rownames(arrow_scores)  # Add variable names

# Filter for significant vetors
significant_arrows <- arrow_scores[env_fit$vectors$pvals < 0.05, ]

# All colored by species
# Interpretation: elevation is perfectly correlated with MDS1, gapF_fill is perfectly correlated with MDS2, rest of the variables have low correlation with both axes
png(file=here("3D_models/figures/mds_mpha_species.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 10), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  xpd = TRUE
)

plot(mds_mpha$points[,1], mds_mpha$points[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(c("#b6c590", "#6b8e23", "#79d2a6")[factor(mpha_data$label)], 0.7), pch = 19, ylim = c(-1, 0.6))
arrows(0, 0, significant_arrows$Dim1, significant_arrows$Dim2, length = 0.1, col = "black", lwd = 1.5)
text(significant_arrows$Dim1, significant_arrows$Dim2, labels = significant_arrows$variable, pos = 3, cex = 0.8)
legend('topright', inset = c(-0.3,0), legend = c("MPHA.ln1", "MPHA.ln2", "MPHA.ln4"), col = add.alpha(c("#b6c590", "#6b8e23", "#79d2a6"), 0.7), pch = 19, title = 'Species')

dev.off()

# All colored by depth
png(file=here("3D_models/figures/mds_mpha_depth.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 6), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  xpd = TRUE
)

plot(mds_mpha$points[,1], mds_mpha$points[,2], xlab = 'MDS 1', ylab = 'MDS 2', col = add.alpha(colors_depth[factor(mpha_data$depth)], 0.7), pch = 19, ylim = c(-1, 0.6))
arrows(0, 0, significant_arrows$Dim1, significant_arrows$Dim2, length = 0.1, col = "black", lwd = 1.5)
text(significant_arrows$Dim1, significant_arrows$Dim2, labels = significant_arrows$variable, pos = 3, cex = 0.8)
legend('topright', inset = c(-0.15,0), legend = c('5', '10', '20', '40', '60'), col = add.alpha(colors_depth, 0.7), pch = 19, title = 'Depth')

dev.off()

# PCA analysis: correlation matrix
pca_mpha <- prcomp(mpha_data[,-c(1:3)], center = TRUE, scale. = TRUE) # correlation matrix pca

summ <- summary(pca_mpha) # calculates variance
summ$importance

png(file=here("3D_models/figures/pca_mpha.png"), width=8, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow = c(2,2)
)

plot(pca_mpha$x, col = add.alpha(colors_depth[factor(mpha_data$depth)], 0.7), pch = 19, xlab = 'PC 1 (43.7%)', ylab = 'PC 2 (20.2%)', main = 'Depth')
legend('bottomright', legend = c('5', '10', '20', '40', '60'), col = add.alpha(colors_depth, 0.7), pch = 19, cex = 0.8)

plot(pca_mpha$x, col = add.alpha(c("#b6c590", "#6b8e23", "#79d2a6")[factor(mpha_data$label)], 0.7), pch = 19, xlab = 'PC 1 (43.7%)', ylab = 'PC 2 (20.2%)', main = 'Species')
legend('bottomright', legend = c("MPHA.ln1", "MPHA.ln2", "MPHA.ln4"), col = add.alpha(c("#b6c590", "#6b8e23", "#79d2a6"), 0.7), pch = 19, cex = 0.9)

plot(pca_mpha$x, col = add.alpha(colors_site[factor(mpha_data$site)], 0.7), pch = 19, xlab = 'PC 1 (43.7%)', ylab = 'PC 2 (20.2%)', main = 'Sites')
legend('bottomright', legend = c("dir", "est", "hul", "kal", "sea", "sna", "wat"), col = add.alpha(colors_site, 0.7), pch = 19, cex = 0.65)

dev.off()

# DAPC analysis by species
# Step 1: scale data & perform PCA
model3d.data.mpha.scaled <- scale(mpha_data[, -c(1:3)])
model3d.dapc.mpha.species <- dapc(model3d.data.mpha.scaled, mpha_data$label)# by species

# Choose number of PCs to retain
a_score <- optim.a.score(model3d.dapc.mpha.species) 

# Step 2: perform DAPC
model3d.dapc.mpha.species <- dapc(model3d.data.mpha.scaled, mpha_data$label, n.pca = 4, n.da = 2)

# Step 3: plot
scatter(model3d.dapc.mpha.species)        # Visualizes the DAPC space
compoplot(model3d.dapc.mpha.species)      # Shows posterior probabilities of group assignment

# DAPC analysis by depth
# Step 1: scale data & perform PCA
model3d.data.mpha.scaled <- scale(mpha_data[, -c(1:3)])
model3d.dapc.mpha.depth <- dapc(model3d.data.mpha.scaled, mpha_data$depth)# by depth

# Choose number of PCs to retain
a_score <- optim.a.score(model3d.dapc.mpha.depth) 

# Step 2: perform DAPC
model3d.dapc.mpha.depth <- dapc(model3d.data.mpha.scaled, mpha_data$depth, n.pca = 4, n.da = 2)

# Step 3: plot
scatter(model3d.dapc.mpha.depth)        # Visualizes the DAPC space
compoplot(model3d.dapc.mpha.depth)      # Shows posterior probabilities of group assignment

# Plot a priori DAPCs for MPHA and depth together
png(file=here("3D_models/figures/dapc_mpha.png"), width=9, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow = c(1,2)
)

scatter(model3d.dapc.mpha.species)  
scatter(model3d.dapc.mpha.depth)  

dev.off()

# DAPC analysis by species AND depth
# Data manipulation: create a new column for unique species & depth
mpha_data$label_depth <- paste(mpha_data$label, mpha_data$depth, sep = '_')
mpha_data$label_depth <- factor(mpha_data$label_depth, levels = c("MPHA.ln1_5", "MPHA.ln1_10", "MPHA.ln1_20", "MPHA.ln2_10", "MPHA.ln2_20", "MPHA.ln2_40", "MPHA.ln2_60", "MPHA.ln4_40", "MPHA.ln4_60"))

# Step 1: scale data & perform PCA
model3d.data.mpha.scaled <- scale(mpha_data[, -c(1:3, 9)])
model3d.dapc.mpha.speciesdepth <- dapc(model3d.data.mpha.scaled, mpha_data$label_depth)# by species/depth

# Choose number of PCs to retain
a_score <- optim.a.score(model3d.dapc.mpha.speciesdepth) 

# cross validation, dataset too small for this method to work super well
set.seed(999)
datax <- xvalDapc(model3d.data.mpha.scaled, mpha_data$label_depth)
system.time(datax <- xvalDapc(model3d.data.mpha.scaled, mpha_data$label_depth,
                              n.pca = 1:5, n.rep = 1000,
                              parallel = "multicore", ncpus = 4L)) #4 PCs

# Step 2: perform DAPC
model3d.dapc.mpha.speciesdepth <- dapc(model3d.data.mpha.scaled, mpha_data$label_depth, n.pca = 4, n.da = 2)

# Step 3: plot
scatter(model3d.dapc.mpha.speciesdepth, clabel = FALSE)        # Visualizes the DAPC space
compoplot(model3d.dapc.mpha.speciesdepth)      # Shows posterior probabilities of group assignment

# Step 4: better plot
table(mpha_data$label, mpha_data$depth)
colors_mpha_depth <- c("#A6CEE3", "#1F78B4", "blue4", "greenyellow", "#B2DF8A", "#33A02C", "darkgreen", "#FDBF6F", "#FF7F00")
var_explained <- round(100 * model3d.dapc.mpha.speciesdepth$eig / sum(model3d.dapc.mpha.speciesdepth$eig), 1)

# Step 5: how correct is the DAPC assignment?
summary(model3d.dapc.mpha.speciesdepth)
table(model3d.dapc.mpha.speciesdepth$grp, model3d.dapc.mpha.speciesdepth$assign)
assignplot(model3d.dapc.mpha.speciesdepth) #blue crosses (prior) on red rectangles (posterior; red=1, white=0) means DAPC posterior classification is consistent with genetic species identify

# Plot a priori DAPCs for MPHA and depth together
png(file=here("3D_models/figures/dapc_mphabydepth.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 10), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  xpd = TRUE
)

plot(model3d.dapc.mpha.speciesdepth$ind.coord[,1], model3d.dapc.mpha.speciesdepth$ind.coord[,2], col = add.alpha(colors_mpha_depth[factor(mpha_data$label_depth)], 0.7), pch = 19, las = 1, xlab = 'LD1, variance explained 82.7%', ylab = 'LD2, variance explained 13.2%')

legend('topright', inset = c(-0.35,0), legend = c("MPHA.ln1_5", "MPHA.ln1_10", "MPHA.ln1_20",  "MPHA.ln2_10", "MPHA.ln2_20", "MPHA.ln2_40", "MPHA.ln2_60", "MPHA.ln4_40", "MPHA.ln4_60"), col = add.alpha(colors_mpha_depth, 0.7), pch = 19, title = 'Species by Depth')

dev.off()

# Separate DAPC by depth (doesn't work super well because only two groups, so therefore 1 DF)
# Data by depth
mpha5 <-  mpha_data[which(mpha_data$depth == 5),] # only MPHA.ln1
mpha10 <- mpha_data[which(mpha_data$depth == 10),]
mpha20 <- mpha_data[which(mpha_data$depth == 20),]
mpha40 <- mpha_data[which(mpha_data$depth == 40),]
mpha60 <- mpha_data[which(mpha_data$depth == 60),]

# Step 1: scale data & perform PCA
mpha10.scaled <- scale(mpha10[, -c(1:3,9)])
mpha10.dapc <- dapc(mpha10.scaled, mpha10$label)# by species

mpha20.scaled <- scale(mpha20[, -c(1:3,9)])
mpha20.dapc <- dapc(mpha20.scaled, mpha20$label)# by species

mpha40.scaled <- scale(mpha40[, -c(1:3,9)])
mpha40.dapc <- dapc(mpha40.scaled, mpha40$label)# by species

mpha60.scaled <- scale(mpha60[, -c(1:3,9)])
mpha60.dapc <- dapc(mpha60.scaled, mpha60$label)# by species

# Choose number of PCs to retain
a_score <- optim.a.score(mpha10.dapc) 
a_score2 <- optim.a.score(mpha20.dapc) 
a_score3 <- optim.a.score(mpha40.dapc) 
a_score4 <- optim.a.score(mpha60.dapc) 

# Step 2: perform DAPC
mpha10.dapc <- dapc(mpha10.scaled, mpha10$label, n.pca = 3, n.da = 1)
mpha20.dapc <- dapc(mpha20.scaled, mpha20$label, n.pca = 3, n.da = 1)
mpha40.dapc <- dapc(mpha40.scaled, mpha40$label, n.pca = 3, n.da = 1)
mpha60.dapc <- dapc(mpha60.scaled, mpha60$label, n.pca = 3, n.da = 1)

# Step 3: plot
colors_mpha10 <- c("#1F78B4", "greenyellow")
colors_mpha20 <- c("blue4", "#B2DF8A")
colors_mpha40 <- c("#33A02C", "#FDBF6F")
colors_mpha60 <- c("darkgreen", "#FF7F00")

scatter(mpha10.dapc, legend = TRUE, col = colors_mpha10)        # Visualizes the DAPC space; doesn't work because only 1 individual in a group
compoplot(mpha10.dapc, col = colors_mpha10)      # Shows posterior probabilities of group assignment

scatter(mpha20.dapc, legend = TRUE, col = colors_mpha20)        # Visualizes the DAPC space
compoplot(mpha20.dapc, col = colors_mpha20)      # Shows posterior probabilities of group assignment

scatter(mpha40.dapc, legend = TRUE, col = colors_mpha40)        # Visualizes the DAPC space
compoplot(mpha40.dapc, col = colors_mpha40)      # Shows posterior probabilities of group assignment

scatter(mpha60.dapc, legend = TRUE, col = colors_mpha60)        # Visualizes the DAPC space; doesn't work because only 1 individual in a group
compoplot(mpha60.dapc, col = colors_mpha60)      # Shows posterior probabilities of group assignment

# Step 4: how correct is the DAPC assignment?
summary(mpha20.dapc)
table(mpha20.dapc$grp, mpha20.dapc$assign)
assignplot(mpha20.dapc) #blue crosses (prior) on red rectangles (posterior; red=1, white=0) means DAPC posterior classification is consistent with genetic species identify

# Plot a priori DAPCs for MPHA and depth together
png(file=here("3D_models/figures/dapc_mphaeverydepth_density.png"), width=10, height=5, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow=c(1,2)
)

scatter(mpha20.dapc, legend = TRUE, col = colors_mpha20)
mtext("20m", side = 3, line = 1, cex = 1.5, font = 2)

scatter(mpha40.dapc, legend = TRUE, col = colors_mpha40)
mtext("40m", side = 3, line = 1, cex = 1.5, font = 2)

dev.off()

png(file=here("3D_models/figures/dapc_mphaeverydepth_compoplot.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(2, 5, 2.5, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  mfrow=c(2,2)
)

compoplot(mpha10.dapc, col = colors_mpha10)
mtext("10m", side = 3, line = 1, cex = 1, font = 2)

compoplot(mpha20.dapc, col = colors_mpha20)
mtext("20m", side = 3, line = 1, cex = 1, font = 2)

compoplot(mpha40.dapc, col = colors_mpha40)
mtext("20m", side = 3, line = 1, cex = 1, font = 2)

compoplot(mpha60.dapc, col = colors_mpha60)
mtext("60m", side = 3, line = 1, cex = 1, font = 2)

dev.off()










































###########################################
#### Validate ANOVA results for kal_20 ####
###########################################
mad_ann_seq_core_20 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == '20' & mad_ann_seq_core$site == 'kal'),]

# Visualization #
boxplot(gapF_fill ~ site + label, data = mad_ann_seq_core_20)

aggregate(gapF_fill ~ label + site, data = mad_ann_seq_core_20,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

table(mad_ann_seq_core_20$depth, mad_ann_seq_core_20$label)

#### ANOVA ####
gapF_fill_lm <- lm(gapF_fill ~ label, data = mad_ann_seq_core_20) #just 20m for kal
summary(gapF_fill_lm)

par(mfrow = c(1,2)) # combine plots
hist(gapF_fill_lm$residuals) # histogram
qqPlot(gapF_fill_lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(gapF_fill_lm$residuals) # formal Shapiro-Wilk normality test

boxplot(gapF_fill ~ label, data = mad_ann_seq_core_20)
dotplot(gapF_fill ~ label, data = mad_ann_seq_core_20)

boxplot(gapF_fill ~ site, data = mad_ann_seq_core_20)
dotplot(gapF_fill ~ site, data = mad_ann_seq_core_20)

leveneTest(gapF_fill ~ label, data = mad_ann_seq_core_20) # formal Levene's test for equality of variances

Anova(gapF_fill_lm)
TUKEY <- TukeyHSD(aov(gapF_fill_lm))

plot(TUKEY)

png(file=here("figures/mad_kal20_species_boxplots.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

# Basic boxplot
colors_mpha <- c("#ff9a81", "#b6c590", "#6b8e23", "#b0e0e6", "gray70")
# colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

a <- boxplot(mad_ann_seq_core_20$gapF_fill ~ mad_ann_seq_core_20$label, col = colors_mpha, ylim = c(0,1.1), xlab = 'Genetic lineage', ylab = 'Gap fraction', main = "20 meters", outline = FALSE) # don't draw outliers so that they aren't doubled when the points go on top

# Add data points
# Points
stripchart(mad_ann_seq_core_20$gapF_fill ~ mad_ann_seq_core_20$label, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "mad_ann_seq_core_20$label")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:4) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col= 'black' )

dev.off()

#####################################
#### ANOVA, look at main effects ####
#####################################
mad_ann_seq_core <- mad_ann_seq_core[which(mad_ann_seq_core$site %in% c('est', 'kal', 'sna')),] # looking at the main sites

# Gap Fraction Fill
gapF_fill_lm <- lm(gapF_fill ~ label + site + depth, data = mad_ann_seq_core) #all 20m corals
summary(gapF_fill_lm)

par(mfrow = c(1,2)) # combine plots
hist(gapF_fill_lm$residuals) # histogram
qqPlot(gapF_fill_lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(gapF_fill_lm$residuals) # formal Shapiro-Wilk normality test

boxplot(gapF_fill ~ label, data = mad_ann_seq_core)
dotplot(gapF_fill ~ label, data = mad_ann_seq_core)

boxplot(gapF_fill ~ site, data = mad_ann_seq_core)
dotplot(gapF_fill ~ site, data = mad_ann_seq_core)

leveneTest(gapF_fill ~ label, data = mad_ann_seq_core) # formal Levene's test for equality of variances

Anova(gapF_fill_lm) # everything has a main effect, but site just barely

# Dev rug
dev_rug_lm <- lm(dev_rug ~ label + site + depth, data = mad_ann_seq_core) #all 20m corals
summary(dev_rug_lm)

par(mfrow = c(1,2)) # combine plots
hist(dev_rug_lm$residuals) # histogram
qqPlot(dev_rug_lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(dev_rug_lm$residuals) # formal Shapiro-Wilk normality test

boxplot(dev_rug ~ label, data = mad_ann_seq_core)
dotplot(dev_rug ~ label, data = mad_ann_seq_core)

boxplot(dev_rug ~ site, data = mad_ann_seq_core)
dotplot(dev_rug ~ site, data = mad_ann_seq_core)

leveneTest(dev_rug ~ label, data = mad_ann_seq_core) # formal Levene's test for equality of variances

Anova(dev_rug_lm) # depth has a main effect

# Vector dispersion
vector_disp_lm <- lm(vector_disp ~ label + site + depth, data = mad_ann_seq_core) #all 20m corals
summary(vector_disp_lm)

par(mfrow = c(1,2)) # combine plots
hist(vector_disp_lm$residuals) # histogram
qqPlot(vector_disp_lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(vector_disp_lm$residuals) # formal Shapiro-Wilk normality test

boxplot(vector_disp ~ label, data = mad_ann_seq_core)
dotplot(vector_disp ~ label, data = mad_ann_seq_core)

boxplot(vector_disp ~ site, data = mad_ann_seq_core)
dotplot(vector_disp ~ site, data = mad_ann_seq_core)

leveneTest(vector_disp ~ label, data = mad_ann_seq_core) # formal Levene's test for equality of variances

Anova(vector_disp_lm) # site and depth has a main effect, but site just barely

# Luminance
luminance_lm <- lm(luminance ~ label + site + depth, data = mad_ann_seq_core) #all 20m corals
summary(luminance_lm)

par(mfrow = c(1,2)) # combine plots
hist(luminance_lm$residuals) # histogram
qqPlot(luminance_lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(luminance_lm$residuals) # formal Shapiro-Wilk normality test

boxplot(luminance ~ label, data = mad_ann_seq_core)
dotplot(luminance ~ label, data = mad_ann_seq_core)

boxplot(luminance ~ site, data = mad_ann_seq_core)
dotplot(luminance ~ site, data = mad_ann_seq_core)

leveneTest(luminance ~ label, data = mad_ann_seq_core) # formal Levene's test for equality of variances

Anova(luminance_lm) # everything has a main effect

#### Scatterplots ####
cols <- brewer.pal(n = 3, name = 'Dark2')
plot(mad_ann_seq_core$dev_rug, mad_ann_seq_core$vector_disp, pch = 19, col = cols[factor(mad_ann_seq_core$site)], xlab = 'dev_rug', ylab = 'vector_disp')
legend('bottomright', legend = c('est', 'kal', 'sna'), pch = 19, col = cols)

plot(mad_ann_seq_core$dev_rug, mad_ann_seq_core$gapF_fill, pch = 19, col = cols[factor(mad_ann_seq_core$site)], xlab = 'dev_rug', ylab = 'gapF_fill')
legend('topright', legend = c('est', 'kal', 'sna'), pch = 19, col = cols)

#############################################################
#### Try to recreate ECRS results for 20m est, kal & sna ####
#############################################################
# Additional samples for all sites are drawing down gap_frac for MPHA.ln1, in particular
mad_ann_seq_core_20 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == '20'),] # try to recreate ECRS results that just used 20m kal, est, sna

mad_ann_seq_core_20_ecrs <- mad_ann_seq_core_20[which(mad_ann_seq_core_20$label %in% c('MDEC', 'MPHA.ln1', 'MPHA.ln2', 'MSEN')),] # 144 corals

gapF_fill_ecrs <- lm(gapF_fill ~ label, data = mad_ann_seq_core_20_ecrs) #just 20m for kal, est and sna
summary(gapF_fill_ecrs)

par(mfrow = c(1,2)) # combine plots
hist(gapF_fill_ecrs$residuals) # histogram
qqPlot(gapF_fill_ecrs$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(gapF_fill_ecrs$residuals) # formal Shapiro-Wilk normality test

boxplot(gapF_fill ~ label, data = mad_ann_seq_core_20_ecrs)
dotplot(gapF_fill ~ label, data = mad_ann_seq_core_20_ecrs)

boxplot(gapF_fill ~ site, data = mad_ann_seq_core_20_ecrs)
dotplot(gapF_fill ~ site, data = mad_ann_seq_core_20_ecrs)

leveneTest(gapF_fill ~ label, data = mad_ann_seq_core_20_ecrs) # formal Levene's test for equality of variances

Anova(gapF_fill_ecrs)
TUKEY <- TukeyHSD(aov(gapF_fill_ecrs))

plot(TUKEY)

# Basic boxplot
colors_mpha <- c("#ff9a81", "#b6c590", "#6b8e23", "#b0e0e6")
show_col(colors_mpha)

a <- boxplot(mad_ann_seq_core_20_ecrs$gapF_fill ~ mad_ann_seq_core_20_ecrs$label, col = colors_mpha, ylim = c(0,1.1), xlab = 'Genetic lineage', ylab = 'Gap fraction', main = "20 meters", outline = FALSE) # don't draw outliers so that they aren't doubled when the points go on top

# Add data points
# Points
stripchart(mad_ann_seq_core_20_ecrs$gapF_fill ~ mad_ann_seq_core_20_ecrs$label, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "mad_ann_seq_core_20_ecrs$label")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:4) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col= 'black' )






























#### By individual site and depth combinations ####
a <- boxplot(kal_10$gapF_fill ~ kal_10$label, col = colors_mpha, ylim = c(0,1.1), xlab = 'Genetic lineage', ylab = 'Gap fraction', main = "", outline = FALSE) # don't draw outliers so that they aren't doubled when the points go on top

# Add data points
# Points
stripchart(mad_ann_seq_20$gapF_fill ~ mad_ann_seq_20$label, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over





























# # Read in Madracis master datasheet
# mad_master <- read.csv("~/Documents/CalAcademy_postdoc/Annotation/madracis_master_merged_May2024.csv", header = TRUE)
# 
# # Read in site-specific habitat measurements
# kal_20 <- read.csv(here("output/hab_measures/kal_20m_example_measurements_jah.csv"), header = TRUE)
# 
# kal_20_merged <- merge(kal_20, mad_master[, c(5,21,24)], by.x = 'id', by.y = 'ID', all.x = TRUE)
# 
# kal_20_merged[which(kal_20_merged$gapF_fill < 0.1),]
# 
# # Annotation quality
# good <- kal_20_merged[which(kal_20_merged$Annotation_point_on_actual_coral == 'yes'),]
# bad <- kal_20_merged[which(kal_20_merged$Annotation_point_on_actual_coral == 'no'),]
# 
# hist(good$gapF_fill)
# hist(bad$gapF_fill)
# 
# bad[which(bad$gapF_fill < 0.1),]
# 
# boxplot(good$gapF_fill ~ good$label, xlab = 'Species')


kal_20 <- read.csv(here("output/hab_measures/kal_20m_example_measurements_jah_quality.csv"), header = TRUE)


hist(kal_20[which(kal_20$label == 'MPHA.ln1'), 'gapF_fill'], col = 'blue')
hist(kal_20[which(kal_20$label == 'MPHA.ln2'), 'gapF_fill'], add = TRUE, col = 'tomato')
lines(density(kal_20[which(kal_20$label == 'MPHA.ln1'), 'gapF_fill']), col = 'blue', lwd = 2)
lines(density(kal_20[which(kal_20$label == 'MPHA.ln2'), 'gapF_fill']), col = 'tomato', lwd = 2)

kal_20_means <- aggregate(kal_20[9:15], list(kal_20$label), FUN = function(x) c(mean = mean(x), se = std.error(x)))
# gapF_fill <- barplot(kal_20_means$gapF_fill[,1] ~ kal_20_means$Group.1, xlab = 'Species', ylab = 'Average gapF_fill')
# error.bar(gapF_fill, kal_20_means$gapF_fill[,1], kal_20_means$gapF_fill[,2])




# Loop through each metric
for (i in 1:length(colnames(kal_20_means)[-1])) {
  # Get group mean and standard error
  metric_mean <- kal_20_means[,i+1][,1]
  metric_se <- kal_20_means[,i+1][,2]
  
  # Plot the bar
  bar_x <- barplot(metric_mean ~ kal_20_means$Group.1, xlab = 'Species', ylab = paste('Average ', colnames(kal_20_means)[i+1]), ylim = c(0, max(metric_mean)+0.1))
  
  # Add error bars
  error.bar(bar_x, metric_mean, metric_se)
}


