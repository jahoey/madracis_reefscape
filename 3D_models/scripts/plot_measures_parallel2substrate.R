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
library(multcomp)
library(RColorBrewer)
library(cluster)
library(effectsize)
library(MASS)
library(vegan)
library(adegenet)
library(factoextra)
library(lme4)

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
directory <- here("3D_models/output/hab_measures_parallel2substrate/")

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
dim(dir_05) #5x24
dim(est_05) #9x24
dim(est_10) #66x24
dim(est_20) #54x24
dim(est_40) #32x24
dim(est_60) #14x24
dim(hul_20) #21x24
dim(kal_05) #13x24
dim(kal_10) #55x24
dim(kal_20) #48x24
dim(kal_40) #19x24
dim(sea_10) #6x24
dim(sea_20) #12x24
dim(sna_10) #15x24
dim(sna_20) #46x24
dim(sna_40) #5x24
dim(wat_05) #2x24

#### Combine the separate objects into one, manipulate and filter ####
coral_ann <- rbind(dir_05, est_05, est_10, est_20, est_40, est_60, hul_20, kal_05, kal_10, kal_20, kal_40, sea_10, sea_20, sna_10, sna_20, sna_40, wat_05) #422x24

# Add labels from metadata
# Merge metadata$Species_mpha into coral_ann by matching 'id' to 'ID'
coral_ann2 <- merge(coral_ann, metadata[, c("ID", "Species_mpha")],
                     by.x = "id", by.y = "ID", all.x = TRUE)

# Update the 'label' column with Species_mpha
coral_ann2$label <- coral_ann$Species_mpha

# Manipulate elevation/slope angle by subtracting any angles greater than 90 degrees by 180
coral_ann3 <- coral_ann2 %>%
  mutate(elevation = if_else(elevation > 90, 180 - elevation, elevation))

mad_ann <- coral_ann3[!coral_ann2$Species_mpha %in% c('PPOR', 'SINT'),] #409x24, remove corals that aren't Madracis

mad_ann_seq <- mad_ann[mad_ann$id %in% (metadata[which(metadata$Davis_library == TRUE),5]),] # work with only genotyped Madracis, 349x24

mad_ann_seq_core <- mad_ann_seq[!mad_ann_seq$id %in% c('JH_0306', 'JH_0307', 'JH_0165', 'JH_0233', 'JH_0013', 'JH_0052', 'JH_0586'),] # remove edge cases for exploration, 342x24

# See if the unknowns have a more updated ID
mad_ann_seq_core[which(mad_ann_seq_core$Species_mpha == 'UNKN'),1]
metadata[metadata$ID %in% mad_ann_seq_core[which(mad_ann_seq_core$Species_mpha == 'UNKN'),1],]

# Change depth categories into integers
mad_ann_seq_core$depth <- as.integer(mad_ann_seq_core$depth)

# Look at means
aggregate(x=mad_ann_seq_core$gapF_fill,
          by=list(mad_ann_seq_core$Species_mpha,mad_ann_seq_core$depth, mad_ann_seq_core$site),
          FUN=mean)

###############################
#### Univariate approaches ####
###############################
#Quick facet plot for each unique site/depth combination. There is some variation in the pattern between sites
# Facet plot with ANOVA for each unique site/depth
plot_boxplot_facet_anova <- function(df, y_var) {
  tukey_results <- df %>%
    group_by(site, depth) %>%
    group_split() %>%
    lapply(function(sub_df) {
      # Skip if there are fewer than two unique Species_mphas
      if (length(unique(sub_df$Species_mpha)) < 2) {
        return(NULL)
      }
      
      # Perform ANOVA
      anova_model <- aov(sub_df[[y_var]] ~ sub_df$Species_mpha)
      
      # Skip if ANOVA couldn't run
      if (length(coef(anova_model)) < 2) {
        return(NULL)
      }
      
      # Check if there are enough levels for Tukey test
      num_groups <- length(unique(sub_df$Species_mpha))
      if (num_groups < 2) {
        return(NULL)
      }
      
      # Perform Tukey’s HSD test
      tukey_test <- tryCatch(
        TukeyHSD(anova_model),
        error = function(e) NULL
      )
      
      # Skip if Tukey test failed or if there are fewer than 2 groups
      if (is.null(tukey_test) || nrow(tukey_test$`sub_df$Species_mpha`) < 1) return(NULL)
      
      # Extract p-values and perform multiple comparison letters
      tukey_values <- tukey_test$`sub_df$Species_mpha`[, "p adj"]
      valid_names <- grep("-", names(tukey_values), value = TRUE)  # Ensure valid comparisons
      tukey_values <- tukey_values[valid_names]
      
      # Skip if no valid comparisons
      if (length(tukey_values) == 0) return(NULL)
      
      # Assign letters based on Tukey test significance
      tukey_letters <- multcompLetters(tukey_values < 0.05)$Letters
      
      # Return formatted results
      data.frame(Species_mpha = names(tukey_letters), Letters = tukey_letters, 
                 site = unique(sub_df$site), depth = unique(sub_df$depth))
    }) %>%
    bind_rows()
  
  # Compute maximum y-values for each group to position Tukey letters properly
  y_max <- df %>%
    group_by(site, depth, Species_mpha) %>%
    summarize(max_y = max(.data[[y_var]], na.rm = TRUE) * 1.05, .groups = "drop")
  
  # Merge Tukey results with max y-values for plotting
  tukey_results <- left_join(tukey_results, y_max, by = c("site", "depth", "Species_mpha"))
  
  # Create boxplot with Tukey letters
  ggplot(df, aes(x = Species_mpha, y = .data[[y_var]], fill = Species_mpha)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_grid(depth ~ site) +  # Facet by site and depth
    geom_text(data = tukey_results, aes(x = Species_mpha, y = max_y *1.2, label = Letters),
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
png(file=here("3D_models/figures/dev_rug_rad0.2_facet_anova_parallel2substrate.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "dev_rug")

dev.off()

png(file=here("3D_models/figures/gapF_fill_rad0.2_facet_anova_parallel2substrate.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
plot_boxplot_facet_anova(mad_ann_seq_core, "gapF_fill")

dev.off()

png(file=here("3D_models/figures/gapF_raw_rad0.2_facet_anova_parallel2substrate.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
plot_boxplot_facet_anova(mad_ann_seq_core, "gapF_raw")

dev.off()

png(file=here("3D_models/figures/luminance_rad0.2_facet_anova_parallel2substrate.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "luminance")

dev.off()

png(file=here("3D_models/figures/vector_disp_rad0.2_facet_anova_parallel2substrate.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "vector_disp")

dev.off()

png(file=here("3D_models/figures/substrate_angle_rad0.2_facet_anova_parallel2substrate.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "elevation")

dev.off()

# Aggregate data across sites by depth, with an option for ANOVA
plot_boxplot_anova <- function(df, y_var) {
  tukey_results <- df %>%
    group_by(depth) %>%
    group_split() %>%
    lapply(function(sub_df) {
      # Perform ANOVA only if there are at least 2 unique Species_mphas
      if (length(unique(sub_df$Species_mpha)) < 2) {
        return(NULL)  # Skip ANOVA if only one group exists
      }
      
      anova_model <- aov(sub_df[[y_var]] ~ sub_df$Species_mpha)
      
      # Perform Tukey’s HSD test
      tukey_test <- tryCatch(
        TukeyHSD(anova_model),
        error = function(e) NULL  # Handle error cases gracefully
      )
      
      # Skip if Tukey test failed
      if (is.null(tukey_test)) return(NULL)
      
      # Extract p-values safely
      tukey_values <- tukey_test$`sub_df$Species_mpha`[, "p adj"]
      valid_names <- grep("-", names(tukey_values), value = TRUE)  # Ensure valid comparisons
      tukey_values <- tukey_values[valid_names]
      
      # Skip if no valid comparisons
      if (length(tukey_values) == 0) return(NULL)
      
      # Assign letters
      tukey_letters <- multcompLetters(tukey_values < 0.05)$Letters
      
      # Return formatted dataframe
      data.frame(Species_mpha = names(tukey_letters), Letters = tukey_letters, depth = unique(sub_df$depth))
    }) %>%
    bind_rows()
  
  # Compute max y values per depth to position letters properly
  y_max <- df %>%
    group_by(depth, Species_mpha) %>%
    summarize(max_y = max(.data[[y_var]], na.rm = TRUE) * 1.05, .groups = "drop")
  
  # Merge Tukey results with max y values
  tukey_results <- left_join(tukey_results, y_max, by = c("depth", "Species_mpha"))
  
  # Plot boxplots with significance letters
  ggplot(df, aes(x = Species_mpha, y = .data[[y_var]], fill = Species_mpha)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_wrap(~depth, nrow = 1) +  # Facet by depth in a single row
    geom_text(data = tukey_results, aes(x = Species_mpha, y = max_y * 1.05, label = Letters),
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
png(file=here("3D_models/figures/dev_rug_bydepth_anova_parallel2substrate.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "dev_rug")

dev.off()

png(file=here("3D_models/figures/gapF_fill_bydepth_anova_parallel2substrate.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "gapF_fill")

dev.off()

png(file=here("3D_models/figures/gapF_raw_bydepth_anova_parallel2substrate.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
plot_boxplot_anova(mad_ann_seq_core, "gapF_raw")

dev.off()

png(file=here("3D_models/figures/luminance_bydepth_anova_parallel2substrate.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "luminance")

dev.off()

png(file=here("3D_models/figures/substrate_angle_bydepth_anova_parallel2substrate.png"), width=16, height=4, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis Species_mpha, tick Species_mpha, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_anova(mad_ann_seq_core, "elevation")

dev.off()

#### Test main effects for each variable ####
# Convert variables to factors
mad_ann_seq_core$Species_mpha <- as.factor(mad_ann_seq_core$Species_mpha)
mad_ann_seq_core$site <- as.factor(mad_ann_seq_core$site)

# Subset data by depth
mad_ann_seq_core_05 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 5),]
mad_ann_seq_core_10 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 10),]
mad_ann_seq_core_20 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 20),]
mad_ann_seq_core_40 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 40),]
mad_ann_seq_core_60 <- mad_ann_seq_core[which(mad_ann_seq_core$depth == 60),]

# Gap fill
gapfill_lm <- lm(gapF_fill ~ Species_mpha + depth + site, data = mad_ann_seq_core)
summary(gapfill_lm) # site (est & sea), depth and lineage
Anova(gapfill_lm)
gapfill_species_tukey <- glht(gapfill_lm, linfct = mcp(Species_mpha = "Tukey"))
summary(gapfill_species_tukey)
gapfill_site_tukey <- glht(gapfill_lm, linfct = mcp(site = "Tukey"))
summary(gapfill_site_tukey)

gapfill_lm20 <- lm(gapF_fill ~ Species_mpha + site, data = mad_ann_seq_core_20)
summary(gapfill_lm20)
Anova(gapfill_lm20)
gapfill_species_tukey20 <- glht(gapfill_lm20, linfct = mcp(Species_mpha = "Tukey"))
summary(gapfill_species_tukey20)

gapfill_lmer <- lmer(gapF_fill ~ Species_mpha * depth + (1 | site), data = mad_ann_seq_core)
summary(gapfill_lmer)

gapfill_lmer20 <- lmer(gapF_fill ~ Species_mpha + (1 | site), data = mad_ann_seq_core_20)
summary(gapfill_lmer20)

# Elevation
elevation_lm <- lm(mad_ann_seq_core$elevation ~ mad_ann_seq_core$Species_mpha + mad_ann_seq_core$depth + mad_ann_seq_core$site) # only depth
summary(elevation_lm)
Anova(elevation_lm)

# Rugosity
devrug_lm <- lm(mad_ann_seq_core$dev_rug ~ mad_ann_seq_core$Species_mpha + mad_ann_seq_core$depth + mad_ann_seq_core$site) # only depth
summary(devrug_lm)
Anova(devrug_lm)

# Vector dispersion
vectordist_lm <- lm(mad_ann_seq_core$vector_disp ~ mad_ann_seq_core$Species_mpha + mad_ann_seq_core$depth + mad_ann_seq_core$site) # site & depth (all sites except wat)
summary(vectordist_lm)
Anova(vectordist_lm)

# Luminance
luminance_lm <- lm(mad_ann_seq_core$luminance ~ mad_ann_seq_core$Species_mpha + mad_ann_seq_core$depth + mad_ann_seq_core$site) # site (overall, but individual sites don't contribute a lot), depth & lineage
summary(luminance_lm)
Anova(luminance_lm)

