library(maps)
library(mapdata)
library(prettymapr)

# Map show much of the Caribbean
png(file="~/Documents/CalAcademy_postdoc/Map/curacao_overview.png", width=10, height=8, res=300, units="in")
par(
  mar=c(5, 7, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=18, # point size, which is the font size
  bg=NA
)
map("worldHires", xlim=c(-100,-50), ylim=c(6,29), col="gray90", fill=TRUE) #plots the region of the USA & canada that I want

dev.off()

# Map just of Curacao
png(file="~/Documents/CalAcademy_postdoc/Maps/curacao.png", width=8, height=8, res=300, units="in")

par(
  mar=c(5, 7, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=18, # point size, which is the font size
  bg=NA
)

sites <- read.delim("~/Documents/CalAcademy_postdoc/Maps/site_coordinates.txt", header = FALSE)

map("worldHires", xlim=c(-69.25,-68.7), ylim=c(11.95,12.4), col="gray90", fill=TRUE)
points(sites$V3, sites$V2, col = "tomato", pch = 19, cex = 1.4)
addnortharrow(pos = 'topright', scale = 0.6, padin = c(0.5,0.75))
addscalebar(pos = 'bottomright')

dev.off()
