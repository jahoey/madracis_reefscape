library(here)
library(stringr)
library(ade4)
library(marmap)

#### Load geographic distances ####
# geo <- as.matrix(read.csv(here("../../unicorn/jennifer/output/cur_sna_20m_20200303_updated_distmat.csv"), header = TRUE, row.names = 1)) # within a site

sites <- read.delim("~/Documents/CalAcademy_postdoc/Maps/site_coordinates.txt", header = FALSE)
colnames(sites) <- c('site', 'latitude', 'longitude')
sites$site2 <- c('kal', 'hul', 'est', 'sna', 'wat', 'sea', 'dir')


#### Load genetic distances ####
m <- as.matrix(read.table(here("angsd_results/coral_snps.ibsMat"), header = FALSE)) # angsd IBS matrix for 650 corals

# Read in IDs and subset IBS matrix to only individuals of interest
# Load all id names
ids <- read.csv(here("angsd_results/bam_lists/ids_recoded650.csv"), header = TRUE) # This is 650 Madracis colonies

# Load in only MDEC
ids_MDEC <- read.table(here("angsd_results/bam_lists/list_MDEC161.txt"), header = TRUE)

metadata <- as.data.frame(str_split_fixed(ids_MDEC$new_ID, '_', 6))
colnames(metadata) <- c('species', 'site', 'depth', 'collector', 'ID', 'extract')
metadata$depth <- as.numeric(metadata$depth)
metadata$id_original <- paste(metadata[,4], metadata[,5], sep = "_")

# Subset ids to only indices staying in analysis
ids2_sub_indices2keep <- as.numeric(rownames(ids[ids$new_ID %in% ids_MDEC$new_ID,]))

# Now subset genetic distances matrix to only MDEC
m_sub <- m[ids2_sub_indices2keep, ids2_sub_indices2keep]
colnames(m_sub) <- metadata$id_original
rownames(m_sub) <- metadata$id_original

#### Subset the MDEC genetic distance matrix to only colonies in sna_20 ####
# ids_MDEC_sna_20 <- metadata$id_original[which(metadata$site == 'sna' & metadata$depth == 20)]
# 
# m_sub_sna_20 <- m_sub[rownames(m_sub) %in% ids_MDEC_sna_20, colnames(m_sub) %in% ids_MDEC_sna_20]
# 
# #### Subset geographic distances to those in specific plot ####
# geo_sna_20 <- geo[rownames(geo) %in% ids_MDEC_sna_20, colnames(geo) %in% ids_MDEC_sna_20]
# 
# #### Order geographic distances matrix by genetic distances ####
# 
# geo_sna_20_ordered <- geo_sna_20[order(rownames(m_sub_sna_20)), order(colnames(m_sub_sna_20))]



#### Get lat/long for each individual ####
metadata_lonlat <- merge(metadata, sites, by.x = 'site', by.y = 'site2', all.x = TRUE)

MDEC_longlat <- metadata_lonlat[,c(7,10,9)]
MDEC_longlat_ordered <- MDEC_longlat[match(metadata$id_original, MDEC_longlat$id_original),]

MDEC_longlat_ordered$id_original == colnames(m_sub)

# Calculate pairwise distances using marmap
# Get bathymetry data
curacao <- getNOAA.bathy(lon1 = -69.4, lon2 = -68.5, lat1 = 11.95, lat2 = 12.6, resolution = 1)

plot(curacao, image = TRUE)
points(MDEC_longlat_ordered$longitude, MDEC_longlat_ordered$latitude, pch = 19)

# get.depth(curacao, locator = TRUE)
get.depth(curacao, MDEC_longlat_ordered[,2:3], locator = FALSE)

# Compute least cost paths
lonlat <- MDEC_longlat_ordered[,-1] # lon column then lat column for dist cost analysis
lonlat2 <- unique(lonlat) # only unique lon/lat for path cost analysis

tr <- trans.mat(curacao, min.depth = -0)
cost <- lc.dist(tr, MDEC_longlat_ordered[,2:3], res = 'path') # this is just for plotting on a map; probably want the 'dist' option for IBD; non-unique lat/lon seems to give errors, but maybe also a resolution thing, def a resolution thing

# Add least cost paths to map (if I can get it to work, seems to be computational)
paths <- lapply(cost, lines, col = 'black', lwd = 0.8, lty = 1)

cost2 <- lc.dist(tr, as.data.frame(MDEC_longlat_ordered[,2:3]), res = 'dist') # for IBD calculations in km
cost2 <- as.matrix(cost2)
colnames(cost2) <- MDEC_longlat_ordered$id_original
rownames(cost2) <- MDEC_longlat_ordered$id_original












#### And now plot! ####
# Within Snake Bay 20
# plot(m_sub_sna_20[lower.tri(m_sub_sna_20)] ~ geo_sna_20_ordered[lower.tri(geo_sna_20_ordered)], xlab = 'Geographic distance (m)', ylab = 'Genetic distance (IBS)', main = 'Snake Bay 20m MDEC')
# linreg <- lm(m_sub_sna_20[lower.tri(m_sub_sna_20)] ~ geo_sna_20_ordered[lower.tri(geo_sna_20_ordered)])
# summary(linreg)
# abline(linreg, col = "gray65")
# 
# # Mantel test for IBD using Rousset's individual distances and greater circle distances (works)
# class(geo_sna_20_ordered)
# class(m_sub_sna_20)
# geo <- as.dist(geo_sna_20_ordered) # b/c mantel.randtest requires objects of class dist
# gen <- as.dist(m_sub_sna_20)
# man <- mantel.randtest(geo, gen, nrepet=10000)
# man
# plot(man, nclass=30, main='Mantel Test', sub='p=0.463')

# For all MDEC
plot(m_sub[lower.tri(m_sub)] ~ jitter(cost2[lower.tri(cost2)]), xlab = 'Geographic distance (km)', ylab = 'Genetic distance (1-IBS)', main = 'MDEC')

linreg <- lm(m_sub[lower.tri(m_sub)] ~ cost2[lower.tri(cost2)])
summary(linreg)
abline(linreg, col = "gray65")

# There are probably some clones in here. This returns position in matrix: less_than_threshold <- which(m_sub < 0.07 & m_sub > 0, arr.ind = TRUE)
m_sub_noouts <- m_sub[-c(7,8,20,138,149,99), -c(7,8,20,138,149,99)] # leaving only 1 colony of potential (?) clone group
cost2_noouts <- cost2[-c(7,8,20,138,149,99), -c(7,8,20,138,149,99)]
colnames(m_sub_noouts) == colnames(cost2_noouts)
plot(m_sub_noouts[lower.tri(m_sub_noouts)] ~ jitter(cost2_noouts[lower.tri(cost2_noouts)]), xlab = 'Geographic distance (km)', ylab = 'Genetic distance (1-IBS)', main = 'MDEC no clones (?)')

linreg2 <- lm(m_sub_noouts[lower.tri(m_sub_noouts)] ~ cost2_noouts[lower.tri(cost2_noouts)])
summary(linreg2)
abline(linreg2, col = "gray65")

# Mantel test for IBD using Rousset's individual distances and greater circle distances (works)
# For all MDEC
class(cost2)
class(m_sub)
geo <- as.dist(cost2) # b/c mantel.randtest requires objects of class dist
gen <- as.dist(m_sub)
man <- mantel.randtest(geo, gen, nrepet=10000)
man
plot(man, nclass=30, main='Mantel Test', sub='p=0.358')

# For MDEC with potential clones removed
geo2 <- as.dist(cost2_noouts) # b/c mantel.randtest requires objects of class dist
gen2 <- as.dist(m_sub_noouts)
man2 <- mantel.randtest(geo2, gen2, nrepet=10000)
plot(man2, nclass=30, main='Mantel Test', sub='p=0.562')

