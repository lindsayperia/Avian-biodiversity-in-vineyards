library(tidyverse)
library(raster)
library(terra)
library(dplyr)
library(sf)
library(lidR)


##################################################################################################
# hand digitized land cover

# area of circle with 150m radius = 70685.83 m^2
digitized <- read_csv("Data/Preliminary Data/LandCoverbyHand.csv")

# create data frame with one column for each land cover type and one row for each point
unique(digitized$CoverType)

land <- data.frame(matrix(NA, nrow = 31, ncol = 9))
names(land) <- c("point", "developed", "woodland", "grassland", "shrubland", "vineyard", "orchard", "rowcrop", "water")
land$point <- 1:31

cover <- data.frame(t(land))
buff <- 70685.83

cover <- data.frame(matrix(NA, nrow = 8, ncol = 4))
names(cover) <- c("CoverType", "area", "pct", "point")
cover$CoverType <- c("developed", "woodland", "grassland", "shrubland", "vineyard", "orchard", "row crop", "water")

for(i in land$point) { 
  
  pt <- digitized %>%
    filter(point == i) # filtering to a specific point 
  sum <- pt %>%
    group_by(CoverType) %>% 
    summarize(area = sum(aream2)) %>%
    data.frame        # creating data frame of total area for each cover type
  sum$pct <- sum$area
  sum$pct <- sum$pct/buff*100 # area to percent cover per cover type
  sum$point <- i
  
# add variables to df with points as rows and cover types as columns
  try(land[i, 2] <- sum %>% 
        filter(CoverType == "developed") %>%
        dplyr::select(pct))
  try(land[i, 3] <- sum %>% 
        filter(CoverType == "woodland") %>%
        dplyr::select(pct))  
  try(land[i, 4] <- sum %>% 
        filter(CoverType == "grassland") %>%
        dplyr::select(pct))
  try(land[i, 5] <- sum %>% 
        filter(CoverType == "shrubland") %>%
        dplyr::select(pct))  
  try(land[i, 6] <- sum %>% 
        filter(CoverType == "vineyard") %>%
        dplyr::select(pct))
  try(land[i, 7] <- sum %>% 
        filter(CoverType == "orchard") %>%
        dplyr::select(pct))  
  try(land[i, 8] <- sum %>% 
        filter(CoverType == "row crop") %>%
        dplyr::select(pct))
  try(land[i, 9] <- sum %>% 
        filter(CoverType == "water") %>%
        dplyr::select(pct))  
  cover <- rbind(cover, sum)
  
}

landcover_final <- land %>% 
  replace(is.na(.), 0)

# read in distance to surface water rds
surface_water <- read_rds("Data/DistToWater.rds")

# LiDAR processing

# loading .LAZ's
las.cat <- readLAScatalog("4_LiDAR/files/")

# loading df of points
points <- st_read("Data/Final point count locations/PointCounts.shp")
points <- points  %>%
  st_transform("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

points <- points %>%
  arrange(Number) %>%
  dplyr::select("Number", "esrignss_l", "esrignss_1")
# empty data frame for lidar output
lidar <- data.frame(matrix(NA, ncol = 4, nrow = 31))
names(lidar) <- c("point", "canopysd", "heightsd","perc4m")
lidar$point <- 1:31

# calculating standard deviation of height and percent canopy cover <4m tall
for(i in points$Number) {
  roi <- clip_roi(las.cat, points[i,], radius = 150) # clipping to area of interest
  dtm <- rasterize_terrain(roi, res = 1, algorithm = tin()) # digital terrain model (terrain height in meters)
  dcm <- rasterize_canopy(roi, res = 1, algorithm = p2r(na.fill = tin())) # digital canopy model (canopy height in meters)
  height <- dcm - dtm # difference between canopy and terrain is canopy height
  lidar[i, 2] <- global(dcm, fun = "sd", na.rm = T) # calculating standard deviation of canopy height
  lidar[i, 3] <- global(height, fun = "sd", na.rm = T) # calculating standard deviation of height
  tall <- freq(height >= 4) # selecting canopy over 4 m tall (value=1)
  lidar[i,4] <- tall[2,3]/sum(tall$count)*100 # canopy over 4m / total raster 
}

names(lidar) <- c("point", "canopySD", "heightSD", "canopy")

# AudioMoth processing
AM <- read_csv("Data/splbyfile.csv")

AM.mean <- AM %>%
  group_by(Folder) %>%
  summarize(mean = mean(`A Weighted (Mean dB)`))

names(AM.mean) <- c("point", "all.dB")
AM.mean <- AM.mean[-32,]
AM.mean$point <- 1:31

# 
to_merge <- list(landcover_final, surface_water, lidar, AM.mean)
environmental.vars<- merge(to_merge, by = "point")

write_rds("Data/EnvironmentalVars.rds")



