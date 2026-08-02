library(here)
library(ggplot2)
library(plotrix)
library(scales)
library(car)
library(lattice)
library(multcompView)
library(FSA)
library(visreg)
library(stringr)
library(tidyverse)

#### Functions ####
# Function to draw error bar later
error.bar.y <- function(x, y, upper, lower=upper, length=0.01,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} # error bar function for y axis

error.bar.x <- function(x, y, right, left=right, length=0.01,...){
  arrows(x+right,y, x-left, y, angle=90, code=3, length=length, ...)
} # error bar function for x axis

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Function to generate Tukey labels
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  #Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.levels <- TUKEY$Species_mpha[,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$Species_mpha=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Species_mpha) , ]
  return(Tukey.labels)
}

#### Read in Nelly's & genotyped data ####
data <- read.delim(here("morphology/data/Nelly_data_May2024_updatedspeciesMPHA_rescaled_2025.csv"), header = TRUE, sep = ',') #691 x 44

jah <- read.delim("~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/bam_lists/ids_recoded672.csv", header = TRUE, sep = ',') #672 x 3

# Extract metadata from filenames
names <- jah$new_ID
names_split <- str_split_fixed(names, '_', 6)

id_original <- paste(names_split[,4], names_split[,5], sep = "_")

#### Subset the data ####
nelly_geno <- data[data$ID %in% id_original,]

# Decide if all of Nelly's data will be used or just genotyped photos
# madracis_data <- data[!data$Species_mpha %in% c('ORB', 'PPOR', 'SINT'),] #686 x 44, includes non-genotyped colonies

madracis_data <- nelly_geno[!nelly_geno$Species_mpha %in% c('ORB', 'PPOR', 'SINT'),] #592 x 44, only genotyped colonies

madracis_data$Yellow.Ruler[which(madracis_data$Yellow.Ruler == "")] <- 'no'

# Make sure any clones are removed
Madracis606 <- read.csv(here("WGS/angsd_results/bam_lists/madracis606_names.csv"), header = FALSE) # list of unique ids (no clones)
Madracis606_split <- as.data.frame(str_split_fixed(Madracis606$V1, '_', 6))
Madracis606_names <- paste(Madracis606_split$V4, Madracis606_split$V5, sep = "_")

madracis_data$ID %in% Madracis606_names # identify potential clones and carefully check 06a_merge_sample_bams.sh to see if one or both are in morphology dataset. If 2, remove 1.

madracis_data[!madracis_data$ID %in% Madracis606_names,]$ID # potential clones. Clone pair can be seen in 06a_merge_sample_bams.sh. Check to see if pair is on the list below. If yes, remove one
madracis_data[madracis_data$ID %in% Madracis606_names,]$ID # were genotyped, but could have a clone

clones <- c('JH_0033', 'JH_0062', 'JH_0071', 'JH_0079', 'JH_0089', 'JH_0122', 'JH_0143', 'JH_0192', 'JH_0212', 'JH_0221', 'JH_0287','JH_0306','JH_0381','JH_0387','JH_0404','JH_0405','JH_0409','JH_0420', 'JH_0463','JH_0464','JH_0466', 'JH_0475', 'JH_0476','JH_0483','JH_0517','JH_0539','JH_0540','JH_0544','JH_0545','JH_0547','JH_0555','JH_0593','JH_0617','JH_0618','JH_0651','JH_0657', 'JH_0661','JH_0663','JH_0677', 'JH_0694','JH_0696','JH_0707', 'JH_0384') # clones; MCAR clones are JH_0136 and JH_0143 - keep JH_0136 b/c measured using white square; remove JH_0384 because JH_0385 actually got measured in this photo

madracis_data2 <- madracis_data[!madracis_data$ID %in% clones, ] # remove the clones; 549 x 44

#### Divide the dataframe into 2 separate dataframes: one for septa and one for IC distance ####
# Intercorallite distance
ic.dist <- madracis_data2[,-c(10:24)]

# Let's figure out how much missing data there are in the data columns. 
which(rowSums(is.na(ic.dist[10:24])) > 0) # this returns the rows of IC distance containing NAs

# Remove the rows containing any NAs for IC distance
ic.dist.trimmed <- ic.dist[-which(rowSums(is.na(ic.dist[10:24])) > 0),] #555 x 29

# Septa
septa <- madracis_data2[,-c(25:39)]

# Let's figure out how much missing data there are in the data columns. 
which(rowSums(is.na(septa[10:24])) > 5) # this returns the rows of septa containing more than 5 NAs

# Remove the rows containing more than 5 NAs for septa 
septa.trimmed <- septa[-which(rowSums(is.na(septa[10:24])) > 5),] #525 x 30

#### Calculate the mean of columns 9 through 23 and store it in a new column ####
septa.trimmed$Mean_Septa <- rowMeans(septa.trimmed[, 10:24], na.rm = TRUE)

ic.dist.trimmed$Mean_IC <- rowMeans(ic.dist.trimmed[, 10:24], na.rm = TRUE)

#### Calculate the SD of columns 9 through 23 and store it in a new column ####
septa.trimmed$SD_Septa <- apply(septa.trimmed[, 10:24], 1, sd, na.rm = TRUE)

ic.dist.trimmed$SD_IC <- apply(ic.dist.trimmed[, 10:24], 1, sd, na.rm = TRUE)

#### Subset datasets ####
# Both septa and IC distance
septa.and.icdist.trimmed <- merge(septa.trimmed, ic.dist.trimmed, by = c("ID", "Species", "Species_mpha", "Depth", "Date", "Dive_Num", "Site", "Original_macrophoto_file_name", "Diver", "Notes", "Yellow.Ruler", "Angled_Ruler", "Blurry_Square", "Yellow_Ruler_Backgroud")) # these are the colonies that have both IC dist and septa measurements: 517 x 46

# IC distance (white & yellow), without blurry white squares
ic.dist.trimmed.quality <- ic.dist.trimmed[!ic.dist.trimmed$Blurry_Square %in% 'Blurry',] #remove blurry square

# IC distance with only blurry white squares
ic.dist.trimmed.blurry <- ic.dist.trimmed[ic.dist.trimmed$Blurry_Square %in% 'Blurry',]

# IC distance only, by method
ic.dist.trimmed.white <- ic.dist.trimmed[which(ic.dist.trimmed$Yellow.Ruler == 'no'),]
ic.dist.trimmed.yellow <- ic.dist.trimmed[which(ic.dist.trimmed$Yellow.Ruler == 'yes'),]

# High quality white squares
ic.dist.trimmed.white.quality <- ic.dist.trimmed.white[!ic.dist.trimmed.white$Blurry_Square %in% 'Blurry',] #remove blurry square

# High quality yellow ruler
ic.dist.trimmed.yellow.quality <- ic.dist.trimmed.yellow[!ic.dist.trimmed.yellow$Angled_Ruler %in% 'Angled',] #remove angled

#### Measurements by species ####
septa.by.species <- aggregate(Mean_Septa ~ Species_mpha, data = septa.trimmed, 
                              FUN = function(x) c(mean = mean(x), se = std.error(x))) # septa only

ic.by.species <- aggregate(Mean_IC ~ Species_mpha, data = ic.dist.trimmed, 
                              FUN = function(x) c(mean = mean(x), se = std.error(x))) # ic dist only

septa.and.icdist.by.species <- aggregate(list(septa = septa.and.icdist.trimmed$Mean_Septa, ic = septa.and.icdist.trimmed$Mean_IC), by = list(septa.and.icdist.trimmed$Species_mpha), FUN = function(x) c(mean = mean(x), se = std.error(x))) # both septa and ic dist

#### Visualizing data: Scatterplots & boxplots ####
# colors_mpha <- c("#ff6347", "#ff9a81", "#32a852", "#b6c590", "#6b8e23", "#b0e0e6", "gray70") #use when all of Nelly's data are used for extra MPHA group
colors_mpha <- c("#ff6347", "#ff9a81", "#b6c590", "#6b8e23", "#b0e0e6", "cornflowerblue", "gray70") #for genotyped data only
show_col(colors_mpha)

# Yellow ruler vs white square
cols <- c('tomato', 'blue')
plot(septa.and.icdist.trimmed$Mean_IC, septa.and.icdist.trimmed$Mean_Septa, col = add.alpha(cols[factor(septa.and.icdist.trimmed$Yellow.Ruler)], 0.7))

boxplot(Mean_IC ~ Yellow.Ruler + Species_mpha, data = ic.dist.trimmed, col = cols)

ggplot(ic.dist.trimmed) +
  aes(x = Species_mpha, y = Mean_IC, fill = Yellow.Ruler) +
  geom_boxplot()

# Average species IC dist vs septa
plot(septa.and.icdist.by.species$ic[,1], septa.and.icdist.by.species$septa[,1], col = alpha('white', 0.8), pch = 19, xlab = 'Average IC dist (mm)', ylab = 'Average septa count', xlim = c(1.1,2), ylim = c(6,10))
error.bar.y(septa.and.icdist.by.species$ic[,1], septa.and.icdist.by.species$septa[,1], septa.and.icdist.by.species$septa[,2], col = colors_mpha, lwd = 2) # septa error bars
error.bar.x(septa.and.icdist.by.species$ic[,1], septa.and.icdist.by.species$septa[,1], septa.and.icdist.by.species$ic[,2], col = colors_mpha, lwd = 2) # ic dist error bars

# Scatterplot of IC dist vs species, depth, site & method
ggplot(ic.dist.trimmed) +
  aes(x = Species_mpha, y = Mean_IC, colour = Yellow.Ruler) +
  geom_point() +
  scale_color_continuous() +
  labs(
    y = "Intercorallite Distance",
    x = "Species",
    color = "Yellow Ruler"
  ) +
  theme_minimal()

# IC distance: Yellow vs quality white squares
ggplot(c) +
  aes(x = Species_mpha, y = Mean_IC, colour = Yellow.Ruler) +
  geom_point() +
  scale_color_discrete() +
  labs(
    y = "Intercorallite Distance",
    x = "Species",
    color = "Yellow Ruler"
  ) +
  theme_minimal()

boxplot(Mean_IC ~ Yellow.Ruler + Species_mpha, data = ic.dist.trimmed.quality, col = cols)

ggplot(ic.dist.trimmed.quality) +
  aes(x = Species_mpha, y = Mean_IC, fill = Yellow.Ruler) +
  geom_boxplot()

# White square photo quality
ggplot(ic.dist.trimmed.white) +
  aes(x = Species_mpha, y = Mean_IC, fill = Blurry_Square) +
  geom_boxplot()

ggplot(ic.dist.trimmed.white) +
  aes(x = Species_mpha, y = Mean_IC, colour = Blurry_Square) +
  geom_point() +
  scale_color_discrete() +
  labs(
    y = "Intercorallite Distance",
    x = "Species",
    color = "Blurry_Square"
  ) +
  theme_minimal()

# Yellow quality
ggplot(ic.dist.trimmed.yellow) +
  aes(x = Species_mpha, y = Mean_IC, fill = Yellow_Ruler_Backgroud) +
  geom_boxplot()

ggplot(ic.dist.trimmed.yellow) +
  aes(x = Species_mpha, y = Mean_IC, fill = Angled_Ruler) +
  geom_boxplot()

#### Statistics ####
## Inter corallite distance for white squares & yellow ruler ##
# ANOVA with interactions
ic.lm <- lm(Mean_IC ~ Species_mpha + Yellow.Ruler, data = ic.dist.trimmed) # + Blurry_Square + Angled_Ruler + Yellow_Ruler_Backgroud * Depth * Site
# First tried: Mean_IC ~ Species_mpha * Yellow.Ruler * Depth * Site, data = ic.dist.trimmed. Yellow.Ruler:Site, Yellow.Ruler:Depth:Site and Species_mpha:Yellow.Ruler:Depth:Site were significant. Adjusted R-squared:  0.4082
# Next: Mean_IC ~ Species_mpha * Yellow.Ruler, data = ic.dist.trimmed. Interaction not sig, so removed. Adjusted R-squared w/o interaction:  0.3845

summary(ic.lm)
Anova(ic.lm, type = 3)

visreg(ic.lm, "Depth", by = "Yellow.Ruler")
visreg(ic.lm, "Species_mpha", by = "Yellow.Ruler")
visreg(ic.lm, "Site", by = "Yellow.Ruler")

## Inter corallite distance for white squares ##
ic.white.lm <- lm(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white) # first included depth, then removed b/c not significant

summary(ic.white.lm)
Anova(ic.white.lm, type = 3)

visreg(ic.white.lm, "Species_mpha", by = "Blurry_Square")

par(mfrow = c(1,2)) # combine plots
hist(ic.white.lm$residuals) # histogram
qqPlot(ic.white.lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(ic.white.lm$residuals) # formal Shapiro-Wilk normality test

boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white)
dotplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white)

boxplot(Mean_IC ~ Site, data = ic.dist.trimmed.white)
dotplot(Mean_IC ~ Site, data = ic.dist.trimmed.white)

leveneTest(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white) # formal Levene's test for equality of variances
leveneTest(Mean_IC ~ Site, data = ic.dist.trimmed.white) # formal Levene's test for equality of variances

aggregate(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

summary(ic.white.lm)

Anova(ic.white.lm)
TUKEY <- TukeyHSD(aov(ic.white.lm))

plot(TUKEY)

png(file=here("morphology/figures/ic_dist_white_geno.png"), width=8, height=6, res=300, units="in")

# Basic boxplot
a <- boxplot(ic.dist.trimmed.white$Mean_IC ~ ic.dist.trimmed.white$Species_mpha, col = colors_mpha, ylim = c(0.4,2.5), xlab = 'Genetic lineage', ylab = 'Average Intercorallite Distance (mm)', outline = FALSE)

# Add data points
# Points
stripchart(ic.dist.trimmed.white$Mean_IC ~ ic.dist.trimmed.white$Species_mpha,              # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  #Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.levels <- TUKEY$Species_mpha[,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$Species_mpha=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Species_mpha) , ]
  return(Tukey.labels)
  }

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "ic.dist.trimmed.white$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:7) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col= 'black' )

dev.off()

## Inter corallite distance for yellow ruler ##
# ANOVA
ic.yellow.lm <- lm(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.yellow.quality) # first included depth, then removed b/c not significant
# First: Mean_IC ~ Species_mpha * Depth * Site, data = ic.dist.trimmed.yellow. Adjusted R-squared:  0.02988 
# Second: Mean_IC ~ Species_mpha * Angled_Ruler * Yellow_Ruler_Backgroud, data = ic.dist.trimmed.yellow. Angled_Ruler is sig.
# Third: Remove Angled ruler from the dataset

summary(ic.yellow.lm)
Anova(ic.yellow.lm, type = 3)

visreg(ic.yellow.lm, "Species_mpha", by = "Depth")
visreg(ic.yellow.lm, "Species_mpha", by = "Site")
visreg(ic.yellow.lm, "Site", by = "Yellow.Ruler")

par(mfrow = c(1,2)) # combine plots
hist(ic.yellow.lm$residuals) # histogram
qqPlot(ic.yellow.lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(ic.yellow.lm$residuals) # formal Shapiro-Wilk normality test

boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.yellow.quality)
dotplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.yellow.quality)

boxplot(Mean_IC ~ Site, data = ic.dist.trimmed.yellow.quality)
dotplot(Mean_IC ~ Site, data = ic.dist.trimmed.yellow.quality)

leveneTest(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.yellow.quality) # formal Levene's test for equality of variances

aggregate(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.yellow.quality,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

summary(ic.yellow.lm)

Anova(ic.yellow.lm)
TUKEY2 <- TukeyHSD(aov(ic.yellow.lm))

plot(TUKEY2)

png(file=here("morphology/figures/ic_dist_yellow_quality_geno.png"), width=8, height=6, res=300, units="in")

# Basic boxplot
b <- boxplot(ic.dist.trimmed.yellow.quality$Mean_IC ~ ic.dist.trimmed.yellow.quality$Species_mpha, col = colors_mpha[-1], ylim = c(0.8,3.5), xlab = 'Genetic lineage', ylab = 'Average Intercorallite Distance (mm)', outline = FALSE)

# Add data points
# Points
stripchart(ic.dist.trimmed.yellow.quality$Mean_IC ~ ic.dist.trimmed.yellow.quality$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS2 <- generate_label_df(TUKEY2 , "ic.dist.trimmed.yellow.quality$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.15*max( b$stats[nrow(b$stats),] )

#Add the labels
text( c(1:7) , b$stats[nrow(b$stats),]+over , LABELS2[,1]  , col= 'black' )

dev.off()

## Intercorallite distance using only high quality white squares ##
# ANOVA
ic.white.quality.lm <- lm(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality)
# First: Mean_IC ~ Species_mpha*Depth*Site, data = ic.dist.trimmed.white.quality. Species_mpha and Species_mpha:Depth are sig. Adjusted R-squared:  0.1735 
# Second: Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality. Species_mpha is sig. Adjusted R-squared:  0.1384

summary(ic.white.quality.lm)
Anova(ic.white.quality.lm, type = 3)

visreg(ic.white.quality.lm, "Species_mpha", by = "Depth")
visreg(ic.white.quality.lm, "Species_mpha", by = "Site")

par(mfrow = c(1,2)) # combine plots
hist(ic.white.quality.lm$residuals) # histogram
qqPlot(ic.white.quality.lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(ic.white.quality.lm$residuals) # formal Shapiro-Wilk normality test

boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality)
dotplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality)

leveneTest(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality) # formal Levene's test for equality of variances

aggregate(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

summary(ic.white.quality.lm)

Anova(ic.white.quality.lm, type = '2')
TUKEY3 <- TukeyHSD(aov(ic.white.quality.lm))

plot(TUKEY3)

png(file=here("morphology/figures/ic_dist_white_quality_geno.png"), width=8, height=6, res=300, units="in")

# Basic boxplot
c <- boxplot(ic.dist.trimmed.white.quality$Mean_IC ~ ic.dist.trimmed.white.quality$Species_mpha, col = colors_mpha, xlab = 'Genetic lineage', ylab = 'Average Intercorallite Distance (mm)', ylim = c(0.8, 2.8), outline = FALSE)

# Add data points
# Points
stripchart(ic.dist.trimmed.white.quality$Mean_IC ~ ic.dist.trimmed.white.quality$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS3 <- generate_label_df(TUKEY3 , "ic.dist.trimmed.white.quality$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.2*max(c$stats[nrow(c$stats),] )

#Add the labels
text( c(1:7) , c$stats[nrow(c$stats),]+over , LABELS3[,1]  , col= 'black' )

dev.off()

# Plot combination plot of IC distance split by quality yellow ruler and quality white squares
png(file=here("morphology/figures/ic_dist_yellow vs_quality_white_geno.png"), width=9, height=5, res=300, units="in")

par(
  mar=c(5, 4, 2, 0.5), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(1,2)
)

# Basic boxplot for yellow ruler
b <- boxplot(ic.dist.trimmed.yellow.quality$Mean_IC ~ ic.dist.trimmed.yellow.quality$Species_mpha, col = colors_mpha[-1], ylim = c(1,3), xlab = '', ylab = '', outline = FALSE, main = 'Quality yellow ruler (n = 105)', las = 1, xaxt = 'n', yaxt = 'n')
axis(1, labels = FALSE)
mtext("Genetic lineage", side = 1, line = 4)
axis(2, las = 2)
mtext("Average Intercorallite Distance (mm)", side = 2, line = 3)

text(x = 1:5,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.1,
     ## Use names from the data list.
     labels = c("MDEC","MPHA.ln1","MPHA.ln2","MSEN.ln1", "MSEN.ln2"), # "MPHA"
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

# Add data points
# Points
stripchart(ic.dist.trimmed.yellow.quality$Mean_IC ~ ic.dist.trimmed.yellow.quality$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS2 <- generate_label_df(TUKEY2 , "ic.dist.trimmed.yellow.quality$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.23*max( b$stats[nrow(b$stats),] )

#Add the labels
text( c(1:5) , b$stats[nrow(b$stats),]+over , LABELS2[,1]  , col= 'black' )

# Basic boxplot for quailty white squares
c <- boxplot(ic.dist.trimmed.white.quality$Mean_IC ~ ic.dist.trimmed.white.quality$Species_mpha, col = colors_mpha, xlab = '', ylab = '', ylim = c(1, 3), outline = FALSE, main = 'Quality white square (n = 252)', xaxt = 'n', yaxt = 'n')
axis(1, labels = FALSE)
mtext("Genetic lineage", side = 1, line = 4)
axis(2, las = 2)
mtext("Average Intercorallite Distance (mm)", side = 2, line = 3)

text(x = 1:6,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.1,
     ## Use names from the data list.
     labels = c("MCAR","MDEC","MPHA.ln1","MPHA.ln2","MSEN.ln1", "MSEN.ln2"), #"MPHA"
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

# Add data points
# Points
stripchart(ic.dist.trimmed.white.quality$Mean_IC ~ ic.dist.trimmed.white.quality$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS3 <- generate_label_df(TUKEY3 , "ic.dist.trimmed.white.quality$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.25*max(c$stats[nrow(c$stats),] )

#Add the labels
text( c(1:6) , c$stats[nrow(c$stats),]+over , LABELS3[,1]  , col= 'black' )

dev.off()

## Septa ##
# ANOVA
septa.lm <- lm(Mean_Septa ~ Species_mpha, data = septa.trimmed) # + Depth + Site + Blurry_Square
summary(septa.lm)

par(mfrow = c(1,2)) # combine plots
hist(septa.lm$residuals) # histogram
qqPlot(septa.lm$residuals, id = TRUE) # qq plot: id = FALSE to remove point identification

shapiro.test(septa.lm$residuals) # formal Shapiro-Wilk normality test

boxplot(Mean_Septa ~ Species_mpha, data = septa.trimmed)
dotplot(Mean_Septa ~ Species_mpha, data = septa.trimmed)

boxplot(Mean_Septa ~ Site, data = septa.trimmed)
dotplot(Mean_Septa ~ Site, data = septa.trimmed)

leveneTest(Mean_Septa ~ Species_mpha, data = septa.trimmed) # formal Levene's test for equality of variances
leveneTest(Mean_Septa ~ Site, data = septa.trimmed) # formal Levene's test for equality of variances

aggregate(Mean_Septa ~ Species_mpha, data = septa.trimmed,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

summary(septa.lm)

Anova(septa.lm)
TUKEY4 <- TukeyHSD(aov(septa.lm))

plot(TUKEY4)

png(file=here(here("morphology/figures/septa_geno.png")), width=8, height=6, res=300, units="in")

par(
  mar=c(5, 4, 2, 0.5), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(1,1)
)

# Basic boxplot
d <- boxplot(septa.trimmed$Mean_Septa ~ septa.trimmed$Species_mpha, col = colors_mpha, ylim = c(5,12.5), xlab = '', ylab = '', outline = FALSE, main = 'Septa (n = 490)', xaxt = 'n', yaxt = 'n')
axis(1, labels = FALSE)
mtext("Genetic lineage", side = 1, line = 4)
axis(2, las = 2)
mtext("Average Septa Count", side = 2, line = 3)

text(x = 1:7,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.4,
     ## Use names from the data list.
     labels = c("MCAR", "MDEC","MPHA.ln1","MPHA.ln2","MSEN.ln1", "MSEN.ln2","UNKN"),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

# Add data points
# Points
stripchart(septa.trimmed$Mean_Septa ~ septa.trimmed$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS4 <- generate_label_df(TUKEY4 , "septa.trimmed$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.11*max( d$stats[nrow(d$stats),] )

#Add the labels
text( c(1:7) , d$stats[nrow(d$stats),]+over , LABELS4[,1]  , col= 'black' )

dev.off()

# Plot combination plot of IC distance for quality white squares and septa count
png(file=here("morphology/figures/septa_and_ic_dist_quality_white_geno.png"), width=9, height=5, res=300, units="in")

par(
  mar=c(5, 4, 2, 0.5), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(1,2)
)

# Basic boxplot for septa count
d <- boxplot(septa.trimmed$Mean_Septa ~ septa.trimmed$Species_mpha, col = colors_mpha, ylim = c(5,12.5), xlab = '', ylab = '', outline = FALSE, main = 'Septa (n = 490)', xaxt = 'n', yaxt = 'n')
axis(1, labels = FALSE)
mtext("Genetic taxon", side = 1, line = 4)
axis(2, las = 2)
mtext("Average Septa Count", side = 2, line = 3)

text(x = 1:7,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.4,
     ## Use names from the data list.
     labels = c("MCAR", "MDEC","MPHA1","MPHA2","MSEN1", "MSEN2","UNKN"),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

# Add data points
# Points
stripchart(septa.trimmed$Mean_Septa ~ septa.trimmed$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS4 <- generate_label_df(TUKEY4 , "septa.trimmed$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.11*max( d$stats[nrow(d$stats),] )

#Add the labels
text( c(1:7) , d$stats[nrow(d$stats),]+over , LABELS4[,1]  , col= 'black' )

# Basic boxplot for quailty white squares
c <- boxplot(ic.dist.trimmed.white.quality$Mean_IC ~ ic.dist.trimmed.white.quality$Species_mpha, col = colors_mpha, xlab = '', ylab = '', ylim = c(0.75, 2.8), outline = FALSE, main = 'Inter-columella  Distance (n = 252)', xaxt = 'n', yaxt = 'n')
axis(1, labels = FALSE)
mtext("Genetic taxon", side = 1, line = 4)
axis(2, las = 2)
mtext("Average Inter-columella  Distance (mm)", side = 2, line = 3)

text(x = 1:6,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.1,
     ## Use names from the data list.
     labels = c("MCAR","MDEC","MPHA1","MPHA2","MSEN1", "MSEN2"), #"MPHA"
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

# Add data points
# Points
stripchart(ic.dist.trimmed.white.quality$Mean_IC ~ ic.dist.trimmed.white.quality$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS3 <- generate_label_df(TUKEY3 , "ic.dist.trimmed.white.quality$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.14*max(c$stats[nrow(c$stats),] )

#Add the labels
text( c(1:6) , c$stats[nrow(c$stats),]+over , LABELS3[,1]  , col= 'black' )

dev.off()

## Calculate summary stats for a table
summary_stats <- function(x) {
  c(
    Mean   = mean(x, na.rm = TRUE),
    SD     = sd(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Min    = min(x, na.rm = TRUE),
    Max    = max(x, na.rm = TRUE)
  )
}

# IC distance
aggregate(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality, FUN = summary_stats)

# Septa
aggregate(Mean_Septa ~ Species_mpha, data = septa.trimmed, FUN = summary_stats)






#### Measurements broken out by depth, regardless of scale used (both yellow & white) ####
# septa.and.icdist.trimmed # colonies with both data for septa & IC dist
# ic.dist.trimmed # colonies with data for IC dist
ic.dist.trimmed.white.quality # colonies with quality data for IC dist

# All IC dist data
# ic.dist.trimmed$Species_mpha <- factor(ic.dist.trimmed$Species_mpha, levels = c("MCAR", "MDEC", "MPHA.ln1", "MPHA.ln2", "MPHA.ln4", "MSEN", "UNKN"))
# 
# ic.dist.trimmed.5 <- ic.dist.trimmed[which(ic.dist.trimmed$Depth == 5),]
# ic.dist.trimmed.10 <- ic.dist.trimmed[which(ic.dist.trimmed$Depth == 10),]
# ic.dist.trimmed.20 <- ic.dist.trimmed[which(ic.dist.trimmed$Depth == 20),]
# ic.dist.trimmed.40 <- ic.dist.trimmed[which(ic.dist.trimmed$Depth == 40),]
# ic.dist.trimmed.60 <- ic.dist.trimmed[which(ic.dist.trimmed$Depth == 60),]
# 
# png(file=here(here("morphology/figures/ic_dist_all_by_depth_geno.png")), width=4, height=11.25, res=300, units="in")
# 
# par(
#   mar=c(2, 4, 2, 1), # panel margin size in "line number" units
#   mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#   tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#   cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#   ps=14, # point size, which is the font size
#   bty = 'n',
#   mfrow = c(5,1)
# )
# 
# # Define a layout with equal panel heights but adjust space for the last plot
# layout_matrix <- matrix(1:5, ncol = 1)
# layout(layout_matrix, heights = c(rep(1, 5 - 1), 1.2)) # Larger height for the last plot
# 
# e1 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.5, col = colors_mpha, ylim = c(0,4), outline = FALSE, xlab = '', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
# axis(1, labels = FALSE)
# mtext("5 meters", side=3, adj = 0, line = 0.7)
# 
# # Points
# stripchart(ic.dist.trimmed.5$Mean_IC ~ ic.dist.trimmed.5$Species_mpha, # Data
#            method = "jitter", # Random noise
#            pch = 1,          # Pch symbols
#            col = 'black',           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# e2 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.10, col = colors_mpha, ylim = c(0,4), outline = FALSE, xlab = '', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
# axis(1, labels = FALSE)
# mtext("10 meters", side=3, adj = 0, line = 0.7)
# 
# # Points
# stripchart(ic.dist.trimmed.10$Mean_IC ~ ic.dist.trimmed.10$Species_mpha, # Data
#            method = "jitter", # Random noise
#            pch = 1,          # Pch symbols
#            col = 'black',           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# e3 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.20, col = colors_mpha, ylim = c(0,4), outline = FALSE, xlab = 'Species', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
# axis(1, labels = FALSE)
# mtext("20 meters", side=3, adj = 0, line = 0.7)
# 
# # Points
# stripchart(ic.dist.trimmed.20$Mean_IC ~ ic.dist.trimmed.20$Species_mpha, # Data
#            method = "jitter", # Random noise
#            pch = 1,          # Pch symbols
#            col = 'black',           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# e4 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.40, col = colors_mpha, ylim = c(0,4), outline = FALSE, xlab = 'Species', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
# axis(1, labels = FALSE)
# mtext("40 meters", side=3, adj = 0, line = 0.7)
# 
# # Points
# stripchart(ic.dist.trimmed.40$Mean_IC ~ ic.dist.trimmed.40$Species_mpha, # Data
#            method = "jitter", # Random noise
#            pch = 1,          # Pch symbols
#            col = 'black',           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# par(mar = c(6, 4, 2, 1))
# e5 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.60, col = colors_mpha, ylim = c(0,4), outline = FALSE, xlab = '', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
# axis(1, labels = FALSE)
# mtext("60 meters", side=3, adj = 0, line = 0.7)
# 
# # Points
# stripchart(ic.dist.trimmed.60$Mean_IC ~ ic.dist.trimmed.60$Species_mpha, # Data
#            method = "jitter", # Random noise
#            pch = 1,          # Pch symbols
#            col = 'black',           # Color of the symbol
#            vertical = TRUE,   # Vertical mode
#            add = TRUE)        # Add it over
# 
# text(x = c(1:7),
#      ## Move labels to just below bottom of chart.
#      y = par("usr")[3] - 0.5,
#      ## Use names from the data list.
#      labels = c("MCAR", "MDEC","MPHA.ln1","MPHA.ln2", "MPHA.ln4","MSEN", "UNKN"),
#      ## Change the clipping region.
#      xpd = NA,
#      ## Rotate the labels by 35 degrees.
#      srt = 35,
#      ## Adjust the labels to almost 100% right-justified.
#      adj = 0.965)
# 
# mtext('Species', side = 1, line = 5)
# 
# dev.off()

# White quality only for IC dist
ic.dist.trimmed.white.quality$Species_mpha <- factor(ic.dist.trimmed.white.quality$Species_mpha, levels = c("MCAR", "MDEC", "MPHA.ln1", "MPHA.ln2", "MSEN.ln1", "MSEN.ln2"))

ic.dist.trimmed.white.quality.5 <- ic.dist.trimmed.white.quality[which(ic.dist.trimmed.white.quality$Depth == 5),]
ic.dist.trimmed.white.quality.10 <- ic.dist.trimmed.white.quality[which(ic.dist.trimmed.white.quality$Depth == 10),]
ic.dist.trimmed.white.quality.20 <- ic.dist.trimmed.white.quality[which(ic.dist.trimmed.white.quality$Depth == 20),]
ic.dist.trimmed.white.quality.40 <- ic.dist.trimmed.white.quality[which(ic.dist.trimmed.white.quality$Depth == 40),]
ic.dist.trimmed.white.quality.60 <- ic.dist.trimmed.white.quality[which(ic.dist.trimmed.white.quality$Depth == 60),]

png(file=here(here("morphology/figures/ic_dist_white_quality_by_depth_geno2.png")), width=14, height=4.25, res=300, units="in")

par(
  mar=c(2, 4, 2, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(1,5)
)

# Define a layout with equal panel heights but adjust space for the last plot
layout_matrix <- matrix(1:5, ncol = 5)
# layout(layout_matrix, heights = c(rep(1, 5 - 1), 1.2)) # Larger height for the last plot

f1 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality.5, col = colors_mpha, ylim = c(0.5,2.5), outline = FALSE, xlab = '', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("5 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(ic.dist.trimmed.white.quality.5$Mean_IC ~ ic.dist.trimmed.white.quality.5$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

f2 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality.10, col = colors_mpha, ylim = c(0.5,2.5), outline = FALSE, xlab = '', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("10 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(ic.dist.trimmed.white.quality.10$Mean_IC ~ ic.dist.trimmed.white.quality.10$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

f3 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality.20, col = colors_mpha, ylim = c(0.5,2.5), outline = FALSE, xlab = 'Species', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("20 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(ic.dist.trimmed.white.quality.20$Mean_IC ~ ic.dist.trimmed.white.quality.20$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

f4 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality.40, col = colors_mpha, ylim = c(0.5,2.5), outline = FALSE, xlab = 'Species', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("40 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(ic.dist.trimmed.white.quality.40$Mean_IC ~ ic.dist.trimmed.white.quality.40$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

par(mar = c(6, 4, 2, 1))
f5 <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.white.quality.60, col = colors_mpha, ylim = c(0.5,2.5), outline = FALSE, xlab = '', ylab = 'Average Intercorallite Distance', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("60 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(ic.dist.trimmed.white.quality.60$Mean_IC ~ ic.dist.trimmed.white.quality.60$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

text(x = c(1:6),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.3,
     ## Use names from the data list.
     labels = c("MCAR", "MDEC","MPHA.ln1","MPHA.ln2", "MPHA.ln4","MSEN"),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

mtext('Genetic lineage', side = 1, line = 5)

dev.off()

#### All septa data broken out by depth ####
septa.trimmed$Species_mpha <- factor(septa.trimmed$Species_mpha, levels = c("MCAR", "MDEC", "MPHA.ln1", "MPHA.ln2", "MSEN.ln1", "MSEN.ln2", "UNKN"))

septa.trimmed.5 <- septa.trimmed[which(septa.trimmed$Depth == 5),]
septa.trimmed.10 <- septa.trimmed[which(septa.trimmed$Depth == 10),]
septa.trimmed.20 <- septa.trimmed[which(septa.trimmed$Depth == 20),]
septa.trimmed.40 <- septa.trimmed[which(septa.trimmed$Depth == 40),]
septa.trimmed.60 <- septa.trimmed[which(septa.trimmed$Depth == 60),]

png(file=here(here("morphology/figures/septa_all_by_depth_geno.png")), width=4, height=11.25, res=300, units="in")

par(
  mar=c(2, 6, 2, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(5,1)
)

# Define a layout with equal panel heights but adjust space for the last plot
layout_matrix <- matrix(1:5, ncol = 1)
layout(layout_matrix, heights = c(rep(1, 5 - 1), 1.2)) # Larger height for the last plot

g1 <- boxplot(Mean_Septa ~ Species_mpha, data = septa.trimmed.5, col = colors_mpha, ylim = c(5,12), outline = FALSE, xlab = '', ylab = 'Number of septa protruding \n through tissue', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("5 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(septa.trimmed.5$Mean_Septa ~ septa.trimmed.5$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

g2 <- boxplot(Mean_Septa ~ Species_mpha, data = septa.trimmed.10, col = colors_mpha, ylim = c(5,12), outline = FALSE, xlab = '', ylab = 'Number of septa protruding \n through tissue', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("10 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(septa.trimmed.10$Mean_Septa ~ septa.trimmed.10$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

g3 <- boxplot(Mean_Septa ~ Species_mpha, data = septa.trimmed.20, col = colors_mpha, ylim = c(5,12), outline = FALSE, xlab = 'Species', ylab = 'Number of septa protruding \n through tissue', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("20 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(septa.trimmed.20$Mean_Septa ~ septa.trimmed.20$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

g4 <- boxplot(Mean_Septa ~ Species_mpha, data = septa.trimmed.40, col = colors_mpha, ylim = c(5,12), outline = FALSE, xlab = 'Species', ylab = 'Number of septa protruding \n through tissue', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("40 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(septa.trimmed.40$Mean_Septa ~ septa.trimmed.40$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

par(mar = c(6, 6, 2, 1))
g5 <- boxplot(Mean_Septa ~ Species_mpha, data = septa.trimmed.60, col = colors_mpha, ylim = c(5,12), outline = FALSE, xlab = '', ylab = 'Number of septa protruding \n through tissue', xaxt = 'n', las = 2)
axis(1, labels = FALSE)
mtext("60 meters", side=3, adj = 0, line = 0.7)

# Points
stripchart(septa.trimmed.60$Mean_Septa ~ septa.trimmed.60$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

text(x = c(1:7),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.7,
     ## Use names from the data list.
     labels = c("MCAR", "MDEC","MPHA.ln1","MPHA.ln2","MSEN.ln1", "MSEN.ln2", "UNKN"),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965)

mtext('Species', side = 1, line = 5)

dev.off()


# Plotting IC distance for each depth
plot_boxplot_depth_anova <- function(df, y_var = "Mean_IC") {
  
  # --- run Tukey or t-tests by depth ---
  tukey_results <- df %>%
    group_by(Depth) %>%
    group_split() %>%
    lapply(function(sub_df) {
      n_groups <- length(unique(sub_df$Species_mpha))
      if (n_groups < 2) return(NULL)  # skip single-group depths
      
      if (n_groups == 2) {
        # t-test case
        ttest <- t.test(sub_df[[y_var]] ~ sub_df$Species_mpha)
        sig <- ttest$p.value < 0.05
        labs <- if (sig) c("a", "b") else c("a", "a")
        names(labs) <- unique(sub_df$Species_mpha)
      } else {
        # ANOVA + Tukey case
        anova_model <- aov(sub_df[[y_var]] ~ sub_df$Species_mpha)
        tukey_test <- tryCatch(TukeyHSD(anova_model), error = function(e) NULL)
        if (is.null(tukey_test)) return(NULL)
        
        tukey_values <- tukey_test$`sub_df$Species_mpha`[, "p adj"]
        if (is.null(names(tukey_values)) || length(tukey_values) == 0) return(NULL)
        names(tukey_values) <- as.character(names(tukey_values))
        
        labs <- multcompLetters(tukey_values < 0.05)$Letters
      }
      
      data.frame(
        Species_mpha = names(labs),
        Letters = labs,
        Depth = unique(sub_df$Depth)
      )
    }) %>%
    bind_rows()
  
  # --- compute max y-values for label placement ---
  y_max <- df %>%
    group_by(Depth, Species_mpha) %>%
    summarize(max_y = max(.data[[y_var]], na.rm = TRUE), .groups = "drop")
  
  # --- only join if tukey_results has rows ---
  if (!is.null(tukey_results) && nrow(tukey_results) > 0) {
    tukey_results <- left_join(tukey_results, y_max,
                               by = c("Depth", "Species_mpha"))
  } else {
    tukey_results <- y_max %>% mutate(Letters = "")
  }
  
  # --- plot ---
  ggplot(df, aes(x = Species_mpha, y = .data[[y_var]], fill = Species_mpha)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    facet_grid(~ Depth, scales = "free", space = "free_x") +
    ylim(0.95,2.5) +
    geom_text(
      data = tukey_results, 
      aes(x = Species_mpha, y = max_y * 1.05, label = Letters),
      color = "black", size = 5, fontface = "bold", na.rm = TRUE
    ) +
    labs(x = "Genetic lineage", y = "Avg IC distance (mm)", fill = "lineage") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm")
    ) +
    scale_fill_manual(values = c("#ff6347", "#ff9a81", "#b6c590", "#6b8e23", "#b0e0e6", "cornflowerblue"))
}

#### Individual variation in septa count and IC dist ####
png(file=here("morphology/figures/sd_septa_and_ic_dist_quality_white_geno.png"), width=9, height=5, res=300, units="in")

par(
  mar=c(5, 4, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bty = 'n',
  mfrow = c(1,2)
)

hist(septa.trimmed$SD_Septa, main = 'Standard Deviation of \nSepta (n = 490)', xlab = 'Within-colony variation (count)')
hist(ic.dist.trimmed.white.quality$SD_IC, main = 'Standard Deviation of \nInter-columella Distance (n = 252)', xlab = 'Within-colony variation (mm)')

dev.off()


# Nice plot summarizing morphology across depths and genetic lineages
png(file=here("morphology/figures/morphology_anova.png"), width=10, height=8, res=300, units="in")

par(
  mar=c(4, 5, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

plot_boxplot_depth_anova(ic.dist.trimmed.white.quality)

dev.off()

















#### Sandbox ####
# ANOVA
ic.blurry.lm <- lm(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.blurry)

summary(ic.blurry.lm)
Anova(ic.blurry.lm, type = 3)

TUKEY100 <- TukeyHSD(aov(ic.blurry.lm))
plot(TUKEY100)

png(file=here("ic_dist_white_blurry.png"), width=8, height=6, res=300, units="in")

# Boxplot
g <- boxplot(Mean_IC ~ Species_mpha, data = ic.dist.trimmed.blurry, col = colors_mpha[-1], ylim = c(0.5,2.5), outline = FALSE, xlab = 'Species', ylab = 'Average Intercorallite Distance')

# Points
stripchart(ic.dist.trimmed.blurry$Mean_IC ~ ic.dist.trimmed.blurry$Species_mpha, # Data
           method = "jitter", # Random noise
           pch = 1,          # Pch symbols
           col = 'black',           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Generate Tukey groups
# Apply the function on my dataset
LABELS100 <- generate_label_df(TUKEY100 , "ic.dist.trimmed.blurry$Species_mpha")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.3*max(e$stats[nrow(e$stats),] )

#Add the labels
text( c(1:7) , e$stats[nrow(e$stats),]+over , LABELS100[,1]  , col= 'black' )

dev.off()

