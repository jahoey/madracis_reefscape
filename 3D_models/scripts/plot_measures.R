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

#### Functions ####
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

# List all CSV files with the pattern *_measurements.csv
csv_files <- list.files(path = directory, pattern = "*_measurements.csv", full.names = TRUE)

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

# Verify that the objects have been created and they have the correct dimensions
dim(est_05) #19x17
dim(est_10) #104x17
dim(est_20) #68x17
dim(est_40) #33x17
dim(est_60) #14x17
dim(hul_20) #21x17
dim(kal_05) #19x17
dim(kal_10) #77x17
dim(kal_20) #66x17
dim(kal_40) #24x17
dim(sea_10) #8x17
dim(sea_20) #13x17
# dim(sna_05) #10x20
dim(sna_10) #23x17
dim(sna_20) #51x17
dim(sna_40) #5x17

#### Combine the separate objects into one and filter ####
coral_ann <- rbind(est_05, est_10, est_20, est_40, est_60, hul_20, kal_05, kal_10, kal_20, kal_40, sna_10, sna_20, sna_40) #524x17

mad_ann <- coral_ann[!coral_ann$label %in% c( 'SINT'),] #512x17, remove SINT

mad_ann_seq <- mad_ann[mad_ann$id %in% (metadata[which(metadata$Davis_library == TRUE),5]),] # work with only genotyped Madracis, 431x17

mad_ann_seq_core <- mad_ann_seq[!mad_ann_seq$id %in% c('JH_0306', 'JH_0307', 'JH_0165', 'JH_0233', 'JH_0013', 'JH_0052', 'JH_0586'),] # remove edge cases for exploration

# See if the unknowns have a more updated ID
mad_ann_seq_core[which(mad_ann_seq_core$label == 'UNKN'),1]
metadata[metadata$ID %in% mad_ann_seq_core[which(mad_ann_seq_core$label == 'UNKN'),1],]

# Change depth categories into integers
mad_ann_seq_core$depth <- as.integer(mad_ann_seq_core$depth)

#### Quick facet plot for each unique site/depth combination. There is some variation in the pattern between sites ####
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
png(file=here("3D_models/figures/dev_rug_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "dev_rug")

dev.off()

png(file=here("3D_models/figures/gapF_fill_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
plot_boxplot_facet_anova(mad_ann_seq_core, "gapF_fill")

dev.off()

png(file=here("3D_models/figures/gapF_raw_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)
plot_boxplot_facet_anova(mad_ann_seq_core, "gapF_raw")

dev.off()

png(file=here("3D_models/figures/luminance_facet_anova.png"), width=12, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_facet_anova(mad_ann_seq_core, "luminance")

dev.off()

#### Aggregate data across sites by depth, with an option for ANOVA ####
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
























mad_ann_seq_20 <- mad_ann_seq[which(mad_ann_seq$depth == '20'),]

mad_ann_seq_20 <- mad_ann_seq_20[!mad_ann_seq_20$id %in% c('JH_0306', 'JH_0307', 'JH_0165', 'JH_0233', 'JH_0013', 'JH_0052', 'JH_0586'),]# remove edge cases

#### Visualization ####
boxplot(gapF_fill ~ site + label, data = mad_ann_seq_20)

aggregate(gapF_fill ~ label + site, data = mad_ann_seq_20,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

table(mad_ann_seq$depth, mad_ann_seq$label)


#### ANOVA ####
gapF_fill_lm <- lm(gapF_fill ~ label, data = mad_ann_seq_20) #just 20m for kal, est and sna
summary(gapF_fill_lm)

par(mfrow = c(1,2)) # combine plots
hist(gapF_fill_lm$residuals) # histogram
qqPlot(gapF_fill_lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(gapF_fill_lm$residuals) # formal Shapiro-Wilk normality test

boxplot(gapF_fill ~ label, data = mad_ann_seq_20)
dotplot(gapF_fill ~ label, data = mad_ann_seq_20)

boxplot(gapF_fill ~ site, data = mad_ann_seq_20)
dotplot(gapF_fill ~ site, data = mad_ann_seq_20)

leveneTest(gapF_fill ~ label, data = mad_ann_seq_20) # formal Levene's test for equality of variances

Anova(gapF_fill_lm)
TUKEY <- TukeyHSD(aov(gapF_fill_lm))

plot(TUKEY)

png(file=here("figures/mad_20m_species_boxplots.png"), width=8, height=6, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

# Basic boxplot
# colors_mpha <- c("#ff9a81", "#b6c590", "#6b8e23", "#b0e0e6")
colors_mpha <- c("#ff6347", "#ff9a81", "#E7298A", "#E6AB02", "#b6c590", "#6b8e23", "#79d2a6", "#b0e0e6", "gray70") # new colors, "#2664d8"
show_col(colors_mpha)

a <- boxplot(mad_ann_seq_20$gapF_fill ~ mad_ann_seq_20$label, col = colors_mpha, ylim = c(0,1.1), xlab = 'Genetic lineage', ylab = 'Gap fraction', main = "20 meters", outline = FALSE) # don't draw outliers so that they aren't doubled when the points go on top

# Add data points
# Points
stripchart(mad_ann_seq_20$gapF_fill ~ mad_ann_seq_20$label, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "mad_ann_seq_20$label")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:4) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col= 'black' )

dev.off()

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


