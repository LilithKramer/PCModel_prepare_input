##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author:        Lilith Kramer
## Creation date: 25 May 2020
## Update:        17 Nov 2021 - minor bug fix, rasterize didn't work properly. Changed field value to 1, instead of character 'id'.
##
## Adapted from: Sven Teurlincx; 20151119_FetchCalc.r        
##
## Goal: Calculate fetch from a lake shapefile and weather data, which can be used as input for PCLake+
##
##
## Needed as input:
##   - Wind data from the Royal Netherlands Meteorological Institute (KNMI) 
##   - Shapefile with a polygon from a lake
##
## Process description:
##   - load wind data 
##   - calculate monthly avg wind direction from wind data
##   - load shapefile with polygon of lake
##   - rotate polygon (according to wind direction)
##   - rasterize polygon to 1m2 grid cells
##   - count cells (in wind direction) 
##   - calculate overall mean fetch
##
## By first calculating the monthly avg wind direction, then calculating the monthly avg fetch, and then calculating the mean of avg fetch 
## you will get a different mean than by first taking the avg wind dir and then calculating the fetch. 
## I think the current approach is more precise. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~
## Prepare environment ####
##~~~~~~~~~~~~~~~~~~~~~~~~~

## detach all packages
# lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE)

## load packages
library(data.table)
library(rgdal)      ## needs sp packages. NB: rgdal will be retired by the end of 2023. In the future this has to be changed to GDAL and PROJ. 
library(maptools)   ## elide function
library(raster)     ## crs function
library(sp)         ## part of raster and rgdal
library(openair)    ## to take a look at the wind data with a windrose

## set directories and filenames
dirWind  <- "..changetowinddatalocation/" ## directory to wind data
dirShape <- "..changetoshapefilelocation" ## directory to shapefile

fnWind  <- "etmgeg_215.txt" ## wind data from location Voorschoten
fnVM    <- "Valkenburgsemeer" ## shapefile from Lake Valkenburg



##~~~~~~~~~~~~~~~~~~
## Wind section ####
##~~~~~~~~~~~~~~~~~~

## direction from which the wind blows
## 360 degrees = north
##  90 degrees = east 
## 180 degrees = south
## 270 degrees = west

## wind speed and wind direction are measured 10m above ground level. Wind speed at ground level can be different, but is not available. 
## there are two options for wind speed available via the KNMI --> FG and FHVEC, daily mean and vector mean. I took the daily mean as this is more representative of the wind during the whole day. 

## load and format wind data
dtWind <- fread(file = paste0(dirWind, fnWindVo), skip = 47, fill = T)
dtWind <- dtWind[-1:-2,]
colnames(dtWind)[colnames(dtWind) %in% "# STN"] <- "STN"
dtWind[, Name := "Voorschoten"]
dtWind[, Year := as.numeric(substr(YYYYMMDD, 1, 4))] ## create year column (for averaging later on)
dtWind[, Month := as.numeric(substr(YYYYMMDD, 5, 6))] ## create month column (for averaging later on)
#dtWind <- dtWind[DDVEC > 0, ] ## 0 = no wind / variable wind, this does not occur in the example dataset
dtWind <- dtWind[, WS_ms := FG/10] ## convert units to m/s
dtWind <- dtWind[, list(Name, Year, Month, DDVEC, WS_ms)] ## keep only the columns you're interested in
rm(list = c("dirWind", "fnWind")) ## clean up environment

## (visual) inspection of wind data
windRose(dtWind[Year > 2013 & Year < 2016, ], type = "Name", ws = "WS_ms", wd = "DDVEC", angle = 20, width = 0.5, grid.line = 2.5)
dtWind[, median(DDVEC), by = list(Name)]
dtWind[, mean(DDVEC), by = list(Name, Year)]
dtWind[, mean(DDVEC), by = list(Name)]

## calculate the mean wind direction per year and month
dtWind_mean <- dtWind[Name %in% "Voorschoten", mean(DDVEC), by = list(Name, Year, Month)]


##~~~~~~~~~~~~~~~~~
## Map section ####
##~~~~~~~~~~~~~~~~~

## NB: Please make sure that you are using a cartesian coordinate projection in meters for your shapefile!
##     A cartesian projection assumes a flat surface, while for a geographic projection, the sphere of the earth is also part of the equation. 
##     e.g. Amersfoort / RD New projection is a cartesian projection in meters, while the WGS84 is a geographic projection with latitudes and longitudes. 
##     This is why Amersfoort / RD New works well with my script: you can put a shapefile with this projection type onto a grid with a resolution of 1 and then you'll immediately get the right distances in meters. 
##
##     The CRS defines the projection of the shapefile. (http://wiki.bk.tudelft.nl/toi-pedia/Coordinate_reference_systems)
##     There are different types of projection and you can transform one into another in QGIS if you wanted to (see: https://docs.qgis.org/2.8/en/docs/training_manual/vector_analysis/reproject_transform.html). 
##     But you can also use R (there's an example of a reprojection starting from line 106). 

## load shapefile with polygon into R
polyValkMeer <- rgdal::readOGR(dsn = dirShape, layer= fnVM)
polyValkMeer@data$AREA <- "Valkenburg" ## add area identifier to the water feature

## convert wind location data to a SpatialDataPointsFile and reproject it to the same coordinates as the shapefile
## 215	Voorschoten	lat 52.13333	lon 4.433333 (position of wind measurement Location)
# dfWindLoc <- data.frame(STN = "215", LON = 4.433333, LAT = 52.13333, ALT = 10)
# coordinates(dfWindLoc) <- c("LON", "LAT")
# proj4string(dfWindLoc) <- CRS("+proj=longlat +datum=WGS84")  
# spWindLoc = spTransform(dfWindLoc, proj4string(polyValkMeer))
# plot(polyValkMeer)
# plot(spWindLoc, add = T) ## is quite far away from lake, so needs a different zoom to show up, not included in this script



##~~~~~~~~~~~~~~~~~~~~~~~
## Fetch calculation ####
##~~~~~~~~~~~~~~~~~~~~~~~

## In my example the lake contains small island, I do not take their presence into account
## as I do a simple cumulative summation (every part of the grid that is water is included)
## a (probably slight) overestimation of the value is the result in this case

## this function rotates the polygon of the lake according to the wind direction. 
## if you plot the polygon including its rotation, it should seem as if the wind blows from left to right over the plot window
rotation_function <- function(rotation, polyChosen){
  if(!is.numeric(rotation)) {stop("Rotation value is not numeric.")}
  
  ## debug
  # rotation <- 184.0625
  # polyChosen <- polyValkMeer
  
  polyRotated <- elide(polyChosen, rotate=rotation+90, center = apply(bbox(polyChosen), 1, mean)) ## elide() from: maptools package
  crs(polyRotated) <- crs(polyChosen) ## comes from: raster package
  # plot(polyRotated)
  
  return(polyRotated)
}

## here we apply the rotation function of the mean wind directions
## we end up with a list of polygons that have been rotated according to the wind directions
lRotat <- lapply(dtWind_mean$V1, rotation_function, polyChosen = polyValkMeer)
plot(lRotat[[49]]) ## check if the polygons have been rotated


## this function converts polygons to rasters and then to matrices
## the resulting matrix contains the shape of the polygon as 0's and 1's (the 1's present the lake). 
## as we convert the polygon from an Amersfoort projection in meters with a resolution of one, each matrix cell
## represents 1 meter of the lake
polygon2raster2matrix <- function(polyChosen, resolution, selected_data_to_rasterize, background_value){
  
  # polyChosen <- lRotat[[1]]
  
  x_ext <- ceiling(abs(extent(polyChosen)[1]-extent(polyChosen)[2]))
  y_ext <- ceiling(abs(extent(polyChosen)[3]-extent(polyChosen)[4]))
  
  polyEmpty <- raster(ncol = x_ext, nrow = y_ext)
  extent(polyEmpty) <- extent(polyChosen)
  res(polyEmpty) <- resolution
  crs(polyEmpty) = crs(polyChosen)
  
  ## field = 1, stands for absence / presence with a value of 1 if present. 
  ## If you want to be more specific and chose a column name, you have to convert the data to a spatial*DataFrame, see lines 106+. 
  rasterChosen <- rasterize(polyChosen, polyEmpty, field = 1, background = background_value)
  #plot(rasterChosen)
  #polyChosen$area_ha
  
  matrixChosen <- as.matrix(rasterChosen)
  
  return(matrixChosen)
  
}

## here we apply the polygon2raster2matrix function to the list of rotated polygons
## we end up with a list of matrices with zero's and ones. 
lMatrices <- lapply(lRotat, polygon2raster2matrix, resolution = 1, background_value = 0)

## Fetch calculation
lCumsum_Mat <- lapply(lMatrices, FUN = function(x) apply(x, 2, cumsum)) ## adding up each row so we get the distance over the lake in meters 
lCumsm_Mat_lengths <- lapply(lCumsum_Mat, FUN = function(x) apply(x, 2, max)) ## getting the max number for each row (the actual distance the wind blows over the lake)
lCumsm_Mat_lengths_max <- lapply(lCumsm_Mat_lengths, FUN = max) ## defining the max distance the wind would blow over the lake 
max_each <- unlist(lCumsm_Mat_lengths_max)
lCumsm_Mat_lengths_mean <- lapply(lCumsm_Mat_lengths, FUN = mean) ## defining the mean distance the wind would blow over the lake 
mean_each <- unlist(lCumsm_Mat_lengths_mean)

fetch <- cbind(dtWind_mean, max = max_each, mean = mean_each) ## creating an overview of max and mean fetches
round(mean(fetch$mean), 0) ## final result --> fetch = 584 meter
