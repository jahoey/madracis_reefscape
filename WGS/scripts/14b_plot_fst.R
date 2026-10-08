## load packages
library(here)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(shadowtext)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#fst_colors <- brewer.pal(n = 8, name = "OrRd")
fst_colors <- rev(hcl.colors(10, palette = "OrRd")[-c(1:2)])
show_col(fst_colors)

#### Read in necessary files ####
# Load Fst matrix from unfolded SFS
# fst <- as.matrix(read.table(here("WGS/angsd_results/madracis606/fst/unfoldedsfs_pairwise_fst.txt"), header = TRUE, row.names = 1, sep = "\t", na.strings = c("", "NA"), check.names = FALSE))

# Load Fst matrix from folded SFS
fst2 <- as.matrix(read.table(here("WGS/angsd_results/madracis606/fst/foldedsfs_pairwise_fst.txt"), header = TRUE, row.names = 1, sep = "\t", na.strings = c("", "NA"), check.names = FALSE))

#### Plot ####
fst_long <- as.data.frame(as.table(fst2))
colnames(fst_long) <- c("Group1", "Group2", "FST")

# png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/fst/unfoldedsfs_fst.png", width=5.5, height=5, res=300, units="in")

png(file="~/Documents/CalAcademy_postdoc/Analyses/madracis_reefscape/WGS/angsd_results/madracis606/fst/foldedsfs_fst.png", width=5.5, height=5, res=300, units="in")

par(
  mar=c(4, 3, 3, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=16, # point size, which is the font size
  bty = 'n'
)

ggplot(fst_long, aes(x = Group1, y = Group2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = ifelse(is.na(FST), "", sprintf("%.4f", FST))),
    size = 3.5
  ) +
  scale_fill_gradientn(
    colors = fst_colors,
    na.value = "grey90"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    text = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    legend.margin = margin(l = -8)
  ) +
  labs(fill = expression(F[ST]))

dev.off()

# Plotting with shadows behind text
# ggplot(fst_long, aes(x = Group1, y = Group2, fill = FST)) +
#   geom_tile(color = "white") +
#   geom_shadowtext(
#     aes(label = ifelse(is.na(FST), "", sprintf("%.3f", FST))),
#     colour = "gray90",
#     bg.colour = "black",
#     bg.r = 0.08,
#     size = 3
#   ) +
#   scale_fill_gradientn(
#     colors = fst_colors,
#     na.value = "grey90"
#   ) +
#   coord_equal() +
#   theme_minimal() +
#   theme(
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   labs(fill = expression(F[ST]))
