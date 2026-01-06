rm(list=ls())

## LOADING REQUIRED PACKAGES
library(adehabitatHR)
library(ks)
library(lme4)
library(rgdal)
library(raster)
library(dplyr)
library(rgeos)
library(amt)
library(tidyverse)
library(sf)
library(lubridate)
library(leaflet)
library(tiff)


# Guarei ##################
setwd("~/felipe.bufalo@unesp.br - Google Drive/My Drive/Unesp/Mestrado/Analises/Publication/Guarei/hr_guarei")
GUA01Data <- read.csv(file = "mov_guarei.csv", sep = ",", header = T); head(GUA01Data)

spdataGUA01 <- SpatialPointsDataFrame(coords = GUA01Data[,c("utm.x","utm.y")], data = GUA01Data, proj4string=CRS("+proj=utm +zone=22 +south +datum=WGS84"))

## creating grid for HR analysis

## define central reference for buffer based on the x/y min/max recorded locations for each fragment
xgua<-(xmax(spdataGUA01)+xmin(spdataGUA01))/2
ygua<-(ymax(spdataGUA01)+ymin(spdataGUA01))/2

# Define grid based on wanted extension and resolution
## GUA
cell_size <- 5*5 #define grid size
x <- seq(xgua-2000,xgua+2000,by=cell_size) #define extension - number is buffer in meters to be added to each side -> in this case we are creating a 4x4km grid
y <- seq(ygua-2000,ygua+2000,by=cell_size) ## the number takes the central point as reference
xyGUA <- expand.grid(x=x,y=y)
coordinates(xyGUA) <- ~x+y
gridded(xyGUA) <- TRUE
class(xyGUA)
projection(xyGUA) <- CRS("+proj=utm +zone=22 +south +datum=WGS84") ## CRS must be the same as the one used for the locations

# Check grid compared to data
plot(xyGUA)
plot(spdataGUA01, add=TRUE)

# now you will use the arguments 'grid = xyXXX, extent = 1' to run kernel function
# if you don't use the standardized grid, find the grid and extent values that best fit your data
# the grid is the number of pixels your extent will be divided so higher numbers = higher resolution
# the extent is how many times the extent of your data will be considered for the total extent of your analysis
# so, if the extent of your data changes for each individual, even if you use the same number here to run the
# analysis for different individual, the extent will be different and so the grid size.
# however, pay attention to the CRS of your grid and location points to avoid further errors
# ignore the warning messages!

## KernelUD Href 
kudGUA01 <- kernelUD(spdataGUA01, grid=xyGUA, extent=1, h="href")

hValsCorr <- seq(1.2, 0.3, by=-0.05)     
hValsA <- numeric(length(hValsCorr))
nPolygons <- numeric(length(hValsCorr))      
for (j in 1:length(hValsCorr)) {
  udKerHvals <- kernelUD( # calculate the kernel where different h values are given in each round of the loop
    spdataGUA01, grid=xyGUA, extent=1, h= kudGUA01@h$h*hValsCorr[j], kern = "bivnorm")
  #, same4all=F) # we feed to kernelUD() one individual at the time as a spatialpoint object
  hValsA[j] <- as.numeric(udKerHvals@h$h)
  hrKer95 <- getverticeshr(udKerHvals, percent = 95)
  nPolygons[j] <- length((hrKer95@polygons[[1]]@Polygons))
}

# inspect values
hValsA
nPolygons

#' We take the h value just one step above to estimate home ranges
minHval <- hValsA[max(which(nPolygons==1))]
# check it ...
#hValsA[11] ## change number to position of the last #1 from nPolygons
minHval

## minHval is, by definition, estimated using min h value, below which the HR polygon breaks up

#' And now we estimate and plot the 95% isopleth for each resident group using the h value that we just obtained. 
minGUA01 <- kernelUD(spdataGUA01, grid=xyGUA, extent=1, h = minHval, kern = c("bivnorm"))

# get HR area for each isopleth for HRef analysis
HRsizeHref <- data.frame(kernel.area(minGUA01, percent = c(1, seq(5,95,by=5), 99)))  
# to reshape into long format
HRsizeHrefLong <- reshape(HRsizeHref, varying = list(names(HRsizeHref)), times = names(HRsizeHref), ids = rownames(HRsizeHref), direction = "long")

names(HRsizeHrefLong)[3] <- "isopleth"
names(HRsizeHrefLong)[1] <- "group"
HRsizeHrefLong$group <- as.factor(rep("GUA01", nrow(HRsizeHrefLong)))
names(HRsizeHrefLong)[2] <- "area.href"
HRsizeHrefLong$isopleth <- as.numeric(HRsizeHrefLong$isopleth)
#HRsizeHrefLong$group <- as.factor(HRsizeHrefLong$group)
head(HRsizeHrefLong)


# add PA & IV
HRsizeHrefLong$PA <- HRsizeHrefLong$area.href/(max(HRsizeHrefLong$area.href))*100
HRsizeHrefLong$IV <- HRsizeHrefLong$isopleth

# identify core area as point of maximum deviation from area/%UD equality line
with(HRsizeHrefLong, which((IV-PA) == max(IV-PA)))
HRsizeHrefLong$IV[15] ## change with row number identified above
# 0.7, i.e. the 70% isopleth

# quantify intensity of use (Samuel, Pierce, Garton (1985) Identifying areas of concentrated use within the home range. Journal of Animal Ecology 54,711?719)

# plot area vs. UD% curves  
#tiff(file = "UD_curve.tiff", width = 1600, height = 1600, units = "px", res = 300)  ## if you want to save the plot, use this and the dev.off function in the end
plot(seq(0,100,length = 50), seq(0,100,length = 50), type = "l",lty = 3, lwd=2, xlab="Isopleth (%UD)", ylab="Proportion of use (%)")
lines(PA ~ IV, data = HRsizeHrefLong, col = 2, lwd=2)
points(x = 70, y = 43.5832274, type = "p", lwd=2) ## change X according to value identified above and Y to meet the line
#dev.off()

## getting the plot k95
k95 <- getverticeshr(minGUA01, percent = 95, unin = "m", unout = "ha")
as.data.frame(k95)

plot(k95, col = rainbow(7, alpha = 0.5))
plot(spdataGUA01, pch=16, cex=0.3, add=T)

## getting the plot k70
k70 <- getverticeshr(minGUA01, percent = 70, unin = "m", unout = "ha")
as.data.frame(k70)

plot(k70, col = rainbow(7, alpha = 0.5))
plot(spdataGUA01, pch=16, cex=0.3, add=T)

## getting the plot k50
k50 <- getverticeshr(minGUA01, percent = 50, unin = "m", unout = "ha")
as.data.frame(k50)

plot(k50, col = rainbow(7, alpha = 0.5))
plot(spdataGUA01, pch=16, cex=0.3, add=T)

datadir <- "k50.gpkg"
rgdal::writeOGR(k50,datadir, layer = "k50.gpkg", driver = "GPKG")


# Taquara ##################
setwd("~/felipe.bufalo@unesp.br - Google Drive/My Drive/Unesp/Mestrado/Analises/Publication/PEMD - Taquara/hr_taquara")

TAQData <- read.csv(file = "home_range_taquara.csv", sep = ",", header = T); head(TAQData)

spdataTAQ <- SpatialPointsDataFrame(coords = TAQData[,c("x","y")], data = TAQData, proj4string=CRS("+proj=utm +zone=22 +south +datum=WGS84"))

## creating grid for HR analysis

## define central reference for buffer based on the x/y min/max recorded locations for each fragment
xtaq<-(xmax(spdataTAQ)+xmin(spdataTAQ))/2
ytaq<-(ymax(spdataTAQ)+ymin(spdataTAQ))/2

# Define grid based on wanted extension and resolution
## Taquara
cell_size <- 5*5 #define grid size
x <- seq(xtaq-2000,xtaq+2000,by=cell_size) #define extension - number is buffer in meters to be added to each side -> in this case we are creating a 4x4km grid
y <- seq(ytaq-2000,ytaq+2000,by=cell_size) ## the number takes the central point as reference
xyTAQ <- expand.grid(x=x,y=y)
coordinates(xyTAQ) <- ~x+y
gridded(xyTAQ) <- TRUE
class(xyTAQ)
projection(xyTAQ) <- CRS("+proj=utm +zone=22 +south +datum=WGS84") ## CRS must be the same as the one used for the locations

# Check grid compared to data
plot(xyTAQ)
plot(spdataTAQ, add=TRUE)

# now you will use the arguments 'grid = xyXXX, extent = 1' to run kernel function
# if you don't use the standardized grid, find the grid and extent values that best fit your data
# the grid is the number of pixels your extent will be divided so higher numbers = higher resolution
# the extent is how many times the extent of your data will be considered for the total extent of your analysis
# so, if the extent of your data changes for each individual, even if you use the same number here to run the
# analysis for different individual, the extent will be different and so the grid size.
# however, pay attention to the CRS of your grid and location points to avoid further errors
# ignore the warning messages!

### KernelUD Href 
kudTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h="href")

hValsCorr <- seq(1.2, 0.3, by=-0.05)     
hValsA <- numeric(length(hValsCorr))
nPolygons <- numeric(length(hValsCorr))      
for (j in 1:length(hValsCorr)) {
  udKerHvals <- kernelUD( # calculate the kernel where different h values are given in each round of the loop
    spdataTAQ, grid=xyTAQ, extent=1, h= kudTAQ@h$h*hValsCorr[j], kern = "bivnorm")
  #, same4all=F) # we feed to kernelUD() one individual at the time as a spatialpoint object
  hValsA[j] <- as.numeric(udKerHvals@h$h)
  hrKer95 <- getverticeshr(udKerHvals, percent = 95)
  nPolygons[j] <- length((hrKer95@polygons[[1]]@Polygons))
}

# inspect values
hValsA
nPolygons

#' We take the h value just one step above to estimate home ranges
minHval <- hValsA[max(which(nPolygons==1))]
# check it ...
#hValsA[13] ## change number to position of the last #1 from nPolygons
minHval

## minHval is, by definition, estimated using min h value, below which the HR polygon breaks up

## Estimate Kernel UD again with minHval
kudTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h=minHval)

# get HR area for each isopleth for HRef analysis
HRsizeHref <- data.frame(kernel.area(kudTAQ, percent = c(1, seq(5,95,by=5), 99)))  
# to reshape into long format
HRsizeHrefLong <- reshape(HRsizeHref, varying = list(names(HRsizeHref)), times = names(HRsizeHref), ids = rownames(HRsizeHref), direction = "long")

names(HRsizeHrefLong)[3] <- "isopleth"
names(HRsizeHrefLong)[1] <- "group"
HRsizeHrefLong$group <- as.factor(rep("TAQ", nrow(HRsizeHrefLong)))
names(HRsizeHrefLong)[2] <- "area.href"
HRsizeHrefLong$isopleth <- as.numeric(HRsizeHrefLong$isopleth)
#HRsizeHrefLong$group <- as.factor(HRsizeHrefLong$group)
head(HRsizeHrefLong)

# add PA & IV
HRsizeHrefLong$PA <- HRsizeHrefLong$area.href/(max(HRsizeHrefLong$area.href))*100
HRsizeHrefLong$IV <- HRsizeHrefLong$isopleth

# identify core area as point of maximum deviation from area/%UD equality line
with(HRsizeHrefLong, which((IV-PA) == max(IV-PA)))
HRsizeHrefLong$IV[15] ## change with row number identified above
# 0.7, i.e. the 75% isopleth

# quantify intensity of use (Samuel, Pierce, Garton (1985) Identifying areas of concentrated use within the home range. Journal of Animal Ecology 54,711?719)

# plot area vs. UD% curves  
#tiff(file = "UD_curve.tiff", width = 1600, height = 1600, units = "px", res = 300)  ## if you want to save the plot, use this and the dev.off function in the end
plot(seq(0,100,length = 50), seq(0,100,length = 50), type = "l",lty = 3, lwd=2, xlab="Isopleth (%UD)", ylab="Proportion of use (%)")
lines(PA ~ IV, data = HRsizeHrefLong, col = 2, lwd=2)
points(x = 70, y = 36.9072788, type = "p", lwd=2) ## change X according to value identified above and Y to meet the line
#dev.off()

#' And now we estimate and plot the 95% isopleth for each resident group using the h value that we just obtained. 
minTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h = minHval, kern = c("bivnorm"))

## getting to the plot k95
k95 <- getverticeshr(minTAQ, percent = 95, unin = "m", unout = "ha")
as.data.frame(k95)

plot(k95, col = rainbow(7, alpha = 0.5))
plot(spdataTAQ, pch=16, cex=0.3, add=T)

## getting to the plot k70
k70 <- getverticeshr(minTAQ, percent = 70, unin = "m", unout = "ha")
as.data.frame(k70)

plot(k70, col = rainbow(7, alpha = 0.5))
plot(spdataTAQ, pch=16, cex=0.3, add=T)

datadir <- "k70.gpkg"
rgdal::writeOGR(k70,datadir, layer = "k70.gpkg", driver = "GPKG")

## getting to the plot k50
k50 <- getverticeshr(minTAQ, percent = 50, unin = "m", unout = "ha")
as.data.frame(k50)

# Suzano ##################
rm(list=ls())
setwd("~/felipe.bufalo@unesp.br - Google Drive/My Drive/Unesp/Mestrado/Analises/Publication/Suzano/hr_suzano")

TAQData <- read.csv(file = "home_range_suzano.csv", sep = ",", header = T); head(TAQData)

spdataTAQ <- SpatialPointsDataFrame(coords = TAQData[,c("x","y")], data = TAQData, proj4string=CRS("+proj=utm +zone=22 +south +datum=WGS84"))

## creating grid for HR analysis

## define central reference for buffer based on the x/y min/max recorded locations for each fragment
xtaq<-(xmax(spdataTAQ)+xmin(spdataTAQ))/2
ytaq<-(ymax(spdataTAQ)+ymin(spdataTAQ))/2

# Define grid based on wanted extension and resolution
## Taquara
cell_size <- 5*5 #define grid size
x <- seq(xtaq-2000,xtaq+2000,by=cell_size) #define extension - number is buffer in meters to be added to each side -> in this case we are creating a 4x4km grid
y <- seq(ytaq-2000,ytaq+2000,by=cell_size) ## the number takes the central point as reference
xyTAQ <- expand.grid(x=x,y=y)
coordinates(xyTAQ) <- ~x+y
gridded(xyTAQ) <- TRUE
class(xyTAQ)
projection(xyTAQ) <- CRS("+proj=utm +zone=22 +south +datum=WGS84") ## CRS must be the same as the one used for the locations

# Check grid compared to data
plot(xyTAQ)
plot(spdataTAQ, add=TRUE)

# now you will use the arguments 'grid = xyXXX, extent = 1' to run kernel function
# if you don't use the standardized grid, find the grid and extent values that best fit your data
# the grid is the number of pixels your extent will be divided so higher numbers = higher resolution
# the extent is how many times the extent of your data will be considered for the total extent of your analysis
# so, if the extent of your data changes for each individual, even if you use the same number here to run the
# analysis for different individual, the extent will be different and so the grid size.
# however, pay attention to the CRS of your grid and location points to avoid further errors
# ignore the warning messages!

### KernelUD Href 
kudTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h="href")

hValsCorr <- seq(1.2, 0.3, by=-0.05)     
hValsA <- numeric(length(hValsCorr))
nPolygons <- numeric(length(hValsCorr))      
for (j in 1:length(hValsCorr)) {
  udKerHvals <- kernelUD( # calculate the kernel where different h values are given in each round of the loop
    spdataTAQ, grid=xyTAQ, extent=1, h= kudTAQ@h$h*hValsCorr[j], kern = "bivnorm")
  #, same4all=F) # we feed to kernelUD() one individual at the time as a spatialpoint object
  hValsA[j] <- as.numeric(udKerHvals@h$h)
  hrKer95 <- getverticeshr(udKerHvals, percent = 95)
  nPolygons[j] <- length((hrKer95@polygons[[1]]@Polygons))
}

# inspect values
hValsA
nPolygons

#' We take the h value just one step above to estimate home ranges
minHval <- hValsA[max(which(nPolygons==1))]
# check it ...
#hValsA[13] ## change number to position of the last #1 from nPolygons
minHval

## minHval is, by definition, estimated using min h value, below which the HR polygon breaks up

## Estimate Kernel UD again with minHval
kudTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h=minHval)

# get HR area for each isopleth for HRef analysis
HRsizeHref <- data.frame(kernel.area(kudTAQ, percent = c(1, seq(5,95,by=5), 99)))  
# to reshape into long format
HRsizeHrefLong <- reshape(HRsizeHref, varying = list(names(HRsizeHref)), times = names(HRsizeHref), ids = rownames(HRsizeHref), direction = "long")

names(HRsizeHrefLong)[3] <- "isopleth"
names(HRsizeHrefLong)[1] <- "group"
HRsizeHrefLong$group <- as.factor(rep("TAQ", nrow(HRsizeHrefLong)))
names(HRsizeHrefLong)[2] <- "area.href"
HRsizeHrefLong$isopleth <- as.numeric(HRsizeHrefLong$isopleth)
#HRsizeHrefLong$group <- as.factor(HRsizeHrefLong$group)
head(HRsizeHrefLong)

# add PA & IV
HRsizeHrefLong$PA <- HRsizeHrefLong$area.href/(max(HRsizeHrefLong$area.href))*100
HRsizeHrefLong$IV <- HRsizeHrefLong$isopleth

# identify core area as point of maximum deviation from area/%UD equality line
with(HRsizeHrefLong, which((IV-PA) == max(IV-PA)))
HRsizeHrefLong$IV[15] ## change with row number identified above
# 0.7, i.e. the 70% isopleth

# quantify intensity of use (Samuel, Pierce, Garton (1985) Identifying areas of concentrated use within the home range. Journal of Animal Ecology 54,711?719)

# plot area vs. UD% curves  
#tiff(file = "UD_curve.tiff", width = 1600, height = 1600, units = "px", res = 300)  ## if you want to save the plot, use this and the dev.off function in the end
plot(seq(0,100,length = 50), seq(0,100,length = 50), type = "l",lty = 3, lwd=2, xlab="Isopleth (%UD)", ylab="Proportion of use (%)")
lines(PA ~ IV, data = HRsizeHrefLong, col = 2, lwd=2)
points(x = 70, y = 30.1346801, type = "p", lwd=2) ## change X according to value identified above and Y to meet the line
#dev.off()

#' And now we estimate and plot the 95% isopleth for each resident group using the h value that we just obtained. 
minTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h = minHval, kern = c("bivnorm"))

## getting to the plot k95
k95 <- getverticeshr(minTAQ, percent = 95, unin = "m", unout = "ha")
as.data.frame(k95)

plot(k95, col = rainbow(7, alpha = 0.5))
plot(spdataTAQ, pch=16, cex=0.3, add=T)

## getting to the plot k70
k70 <- getverticeshr(minTAQ, percent = 70, unin = "m", unout = "ha")
as.data.frame(k70)

plot(k70, col = rainbow(7, alpha = 0.5))
plot(spdataTAQ, pch=16, cex=0.3, add=T)

datadir <- "k95.gpkg"
rgdal::writeOGR(k95,datadir, layer = "k95.gpkg", driver = "GPKG")

## getting to the plot k50
k50 <- getverticeshr(minTAQ, percent = 50, unin = "m", unout = "ha")
as.data.frame(k50)

# Santa Maria ##################
rm(list=ls())
setwd("~/felipe.bufalo@unesp.br - Google Drive/My Drive/Unesp/Mestrado/Analises/Publication/Santa Maria Yness/hr_santamaria")

TAQData <- read.csv(file = "home_range_sma.csv", sep = ",", header = T); head(TAQData)

spdataTAQ <- SpatialPointsDataFrame(coords = TAQData[,c("x","y")], data = TAQData, proj4string=CRS("+proj=utm +zone=22 +south +datum=WGS84"))

## creating grid for HR analysis

## define central reference for buffer based on the x/y min/max recorded locations for each fragment
xtaq<-(xmax(spdataTAQ)+xmin(spdataTAQ))/2
ytaq<-(ymax(spdataTAQ)+ymin(spdataTAQ))/2

# Define grid based on wanted extension and resolution
## Taquara
cell_size <- 5*5 #define grid size
x <- seq(xtaq-2000,xtaq+2000,by=cell_size) #define extension - number is buffer in meters to be added to each side -> in this case we are creating a 4x4km grid
y <- seq(ytaq-2000,ytaq+2000,by=cell_size) ## the number takes the central point as reference
xyTAQ <- expand.grid(x=x,y=y)
coordinates(xyTAQ) <- ~x+y
gridded(xyTAQ) <- TRUE
class(xyTAQ)
projection(xyTAQ) <- CRS("+proj=utm +zone=22 +south +datum=WGS84") ## CRS must be the same as the one used for the locations

# Check grid compared to data
plot(xyTAQ)
plot(spdataTAQ, add=TRUE)

# now you will use the arguments 'grid = xyXXX, extent = 1' to run kernel function
# if you don't use the standardized grid, find the grid and extent values that best fit your data
# the grid is the number of pixels your extent will be divided so higher numbers = higher resolution
# the extent is how many times the extent of your data will be considered for the total extent of your analysis
# so, if the extent of your data changes for each individual, even if you use the same number here to run the
# analysis for different individual, the extent will be different and so the grid size.
# however, pay attention to the CRS of your grid and location points to avoid further errors
# ignore the warning messages!

### KernelUD Href 
kudTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h="href")

hValsCorr <- seq(1.2, 0.3, by=-0.05)     
hValsA <- numeric(length(hValsCorr))
nPolygons <- numeric(length(hValsCorr))      
for (j in 1:length(hValsCorr)) {
  udKerHvals <- kernelUD( # calculate the kernel where different h values are given in each round of the loop
    spdataTAQ, grid=xyTAQ, extent=1, h= kudTAQ@h$h*hValsCorr[j], kern = "bivnorm")
  #, same4all=F) # we feed to kernelUD() one individual at the time as a spatialpoint object
  hValsA[j] <- as.numeric(udKerHvals@h$h)
  hrKer95 <- getverticeshr(udKerHvals, percent = 95)
  nPolygons[j] <- length((hrKer95@polygons[[1]]@Polygons))
}

# inspect values
hValsA
nPolygons

#' We take the h value just one step above to estimate home ranges
minHval <- hValsA[max(which(nPolygons==1))]
# check it ...
#hValsA[13] ## change number to position of the last #1 from nPolygons
minHval

## minHval is, by definition, estimated using min h value, below which the HR polygon breaks up

## Estimate Kernel UD again with minHval
kudTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h=minHval)

# get HR area for each isopleth for HRef analysis
HRsizeHref <- data.frame(kernel.area(kudTAQ, percent = c(1, seq(5,95,by=5), 99)))  
# to reshape into long format
HRsizeHrefLong <- reshape(HRsizeHref, varying = list(names(HRsizeHref)), times = names(HRsizeHref), ids = rownames(HRsizeHref), direction = "long")

names(HRsizeHrefLong)[3] <- "isopleth"
names(HRsizeHrefLong)[1] <- "group"
HRsizeHrefLong$group <- as.factor(rep("TAQ", nrow(HRsizeHrefLong)))
names(HRsizeHrefLong)[2] <- "area.href"
HRsizeHrefLong$isopleth <- as.numeric(HRsizeHrefLong$isopleth)
#HRsizeHrefLong$group <- as.factor(HRsizeHrefLong$group)
head(HRsizeHrefLong)

# add PA & IV
HRsizeHrefLong$PA <- HRsizeHrefLong$area.href/(max(HRsizeHrefLong$area.href))*100
HRsizeHrefLong$IV <- HRsizeHrefLong$isopleth

# identify core area as point of maximum deviation from area/%UD equality line
with(HRsizeHrefLong, which((IV-PA) == max(IV-PA)))
HRsizeHrefLong$IV[15] ## change with row number identified above
# 0.7, i.e. the 70% isopleth

# quantify intensity of use (Samuel, Pierce, Garton (1985) Identifying areas of concentrated use within the home range. Journal of Animal Ecology 54,711?719)

# plot area vs. UD% curves  
#tiff(file = "UD_curve.tiff", width = 1600, height = 1600, units = "px", res = 300)  ## if you want to save the plot, use this and the dev.off function in the end
plot(seq(0,100,length = 50), seq(0,100,length = 50), type = "l",lty = 3, lwd=2, xlab="Isopleth (%UD)", ylab="Proportion of use (%)")
lines(PA ~ IV, data = HRsizeHrefLong, col = 2, lwd=2)
points(x = 70, y = 40.4102564, type = "p", lwd=2) ## change X according to value identified above and Y to meet the line
#dev.off()

#' And now we estimate and plot the 95% isopleth for each resident group using the h value that we just obtained. 
minTAQ <- kernelUD(spdataTAQ, grid=xyTAQ, extent=1, h = minHval, kern = c("bivnorm"))

## getting to the plot k95
k95 <- getverticeshr(minTAQ, percent = 95, unin = "m", unout = "ha")
as.data.frame(k95)

plot(k95, col = rainbow(7, alpha = 0.5))
plot(spdataTAQ, pch=16, cex=0.3, add=T)

## getting to the plot k70
k70 <- getverticeshr(minTAQ, percent = 70, unin = "m", unout = "ha")
as.data.frame(k70)

plot(k70, col = rainbow(7, alpha = 0.5))
plot(spdataTAQ, pch=16, cex=0.3, add=T)

datadir <- "k70.gpkg"
rgdal::writeOGR(k70,datadir, layer = "k70.gpkg", driver = "GPKG")

## getting to the plot k50
k50 <- getverticeshr(minTAQ, percent = 50, unin = "m", unout = "ha")
as.data.frame(k50)
