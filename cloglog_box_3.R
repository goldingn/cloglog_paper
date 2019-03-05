# plotting figure 2 (example) for the cloglog paper

# clear workspace
rm(list = ls())

# load packages and define a couple of small functions
library(RColorBrewer)
library(raster)
library(ggmap)

first <- function (x) x[1]
rescale <- function (x, ncell, range) {
  # map x from some range to [0, ncell]
  x <- x - min(range)
  x <- x / abs(diff(range))
  x <- x * ncell
  x
}

# set up colours
paired <- brewer.pal(12, 'Paired')
col_abs <- paired[1]
col_pres <- paired[2]

panel_col <- grey(0.3)
axis_col <- grey(0.4)
lab_col <- grey(0.5)
line_col <- grey(0.7)

# Start off using these data:
# http://datadryad.org/resource/doi:10.5061/dryad.v4p20
f <- 'http://datadryad.org/bitstream/handle/10255/dryad.63778/Banff_2012.csv'
df <- read.csv(f,
               stringsAsFactors = FALSE)

# define CRS's (UTM NAD zone 11 & wgs84)
nad11 <- CRS('+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
wgs84 <- CRS('+init=epsg:4326')

# trim down to lynx & drop segments, to shrink the map
df <- df[df$Species == 'lynx', ]
drop_segments <- c(61:130,
                   132:277,
                   286:305,
                   311:334,
                   335:391)
df <- df[!(df$segment %in% drop_segments), ]

# find the number of visits per segment & whether lynx ever seen
df_sub <- data.frame(visits = tapply(df$segment, df$segment, length),
                     ever_seen = tapply(df$y.seg, df$segment, first),
                     coord_x = tapply(df$x.centre, df$segment, first),
                     coord_y = tapply(df$y.centre, df$segment, first))

# transform coordinates to wgs84
coords <- SpatialPoints(df_sub[, c('coord_x', 'coord_y')],
                        proj4string = nad11)
coords <- spTransform(coords, wgs84)
df_sub[, c('coord_x', 'coord_y')] <- coords@coords

# get a base map
centroid <- colMeans(df_sub[, c('coord_x', 'coord_y')])
crop_extent <- as.vector(extent(coords))
diffs <- abs(diff(crop_extent))[c(1, 3)] * c(0.4, 0.05)
diffs <- rep(diffs, each = 2) * c(-1, 1, -1, 1)
crop_extent <- crop_extent + diffs
crop_extent <- round(crop_extent, 1)

tmp <- get_map(location = crop_extent[c(1, 3, 2, 4)],
               maptype = 'toner-lite',
               source = 'stamen')

# reorder to put smallest on top
o <- order(df_sub$visits, decreasing = TRUE)
df_sub <- df_sub[o, ]

png('~/temp/cloglog_box_3.png',
    width = 3000,
    height = 3000,
    pointsize = 80)

par(mar = c(4, 4, 2, 2) + 0.1,
    xpd = NA)

# set colour to reflect sampling effort
effort <- colorRampPalette(c(col_abs, col_pres))
colvec <- effort(14)[df_sub$visits]

raster::plot(tmp)

x_ax <- c(crop_extent[1],
          mean(crop_extent[1:2]),
          crop_extent[2])
y_ax <- c(crop_extent[3],
          mean(crop_extent[3:4]),
          crop_extent[4])

axis(side = 1,
     at = rescale(x_ax, ncol(tmp), crop_extent[1:2]),
     labels = round(x_ax, 1),
     col.axis = axis_col,
     cex.axis = 0.8,
     pos = 0)
axis(side = 2,
     at = rescale(y_ax, nrow(tmp), crop_extent[3:4]),
     labels = round(y_ax, 1),
     col.axis = axis_col,
     cex.axis = 0.8,
     las = 2,
     pos = 0)

mtext(text = 'longitude',
      side = 1,
      line = 2.5,
      col = lab_col,
      cex = 1)

mtext(text = 'latitude',
      side = 2,
      line = 1.5,
      col = lab_col,
      cex = 1)

rect(xleft = 0,
     ybottom = 0,
     xright = ncol(tmp),
     ytop = nrow(tmp),
     col = rgb(0.95, 0.95, 0.95, 0.7),
     border = line_col,
     xpd = NA)

points(x = rescale(df_sub$coord_x, ncol(tmp), crop_extent[1:2]),
       y = rescale(df_sub$coord_y, nrow(tmp), crop_extent[3:4]),
       pch = 21,
       bg = ifelse(df_sub$ever_seen, col_pres, col_abs),
       col = grey(0.4),
       lwd = 3,
       cex = log1p(df_sub$visits) * 1)


dev.off()

m <- glm(ever_seen ~ 1 + offset(log(df_sub$visits)),
         data = df_sub,
         family = binomial('cloglog'))

# get the expected number of lynx per 1km per visit
exp(m$coefficients)


