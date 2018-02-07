## University of Bayreuth - MSc Global Change Ecology 
## C7 Pattern of Land Use and Ecosystem Dynamics
## Project by Rebekka Riebl January/February 2018 
## based on previous work for the course M5 Remote Sensing
## Topic: Reforestation in Sicily, with the Commune Sperlinga in the Province Enna as case study.
## normal comments always refer to the following lines
## comments directly behind the command (in the same line) mean there's an issue!
## working directory: sicilyproject

### LIBRARIES needed - from old project, 1 or 2 might be redundant atm. (check again)
library(raster)
library(sp)
library(ggplot2)
library(RStoolbox)
library(rgdal) #maybe?
library(randomForest) #maybe??
library(lattice) #maybe??

### DATA preparation
## getting italy's administrative borders
italy <- getData("GADM", country="ITA", level=3)

## study commune - Province of Enna: lowest pop density in Sicily (wiki); 
## Sperlinga is commune with smallest population within Enna
sper <- italy[ which(italy$NAME_3=='Sperlinga'), ]
## match projections to Landsat data (UTM)
sper <-spTransform(sper,CRS="+proj=utm +zone=33 
                   +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### dataset from the 80s -  3rd May 1987
#windows users possible need to adapt the working directory paths using \ rather than /,
#as I have tried this solely under ubuntu MATE 17.04 and RStudio 1.0.136.
#usual working directory: sicilyproject
setwd("./raster_data/old")
allbandsoldsper <- stack(c("band1oldsper.tif", "band2oldsper.tif", 
                           "band3oldsper.tif", "band4oldsper.tif", 
                           "band5oldsper.tif", "band7oldsper.tif"))
## crop & mask dataset to study commune
cropoldsper <- crop(allbandsoldsper, sper)
oldsper <- mask(cropoldsper, sper)

### recent dataset from 2017 - 15 May 2017
setwd("../now")
allbandsnowsper <- stack(c("band1nowsper.tif", "band2nowsper.tif", 
                           "band3nowsper.tif", "band4nowsper.tif", 
                           "band5nowsper.tif", "band7nowsper.tif"))
cropnowsper <- crop(allbandsnowsper, sper)
nowsper <- mask(cropnowsper, sper)

### ANALYSIS - Land Cover Change
## Change Vector Analysis using RED and NIR bands
#Landsat 5: red=3, nir=4
#Landsat 8 has other band distributions: red=4, nir=5
cva_87_17 <- rasterCVA(oldsper[[3:4]], nowsper[[4:5]]) ##!!!make sure the right bands are used!!!!
#tasseledCap analysis producing greenness, wetness, brightness
#band 6 (thermal) needn't be excluded as it's not loaded in R in the first place
tc_87 <- tasseledCap(oldsper[[c(1:6)]], sat="Landsat5TM")
tc_17 <- tasseledCap(nowsper[[c(1:6)]], sat="Landsat8OLI")
#change vector analysis using greenness and wetness, threshold mimized in order to see change
#(default threshold (tmf=2) would give a yellow map)
cva_tc_87_17 <- rasterCVA(tc_87[[2:3]], tc_17[[2:3]], tmf=1.5)

#angle & magnitude RED/NIR
plot(cva_87_17) 
#angle & magnitude greenness/wetness
plot(cva_tc_87_17) ##worked so far##

### MULTIDATE ANALYSIS
# classification of change classes with QGIS
# trying to save layers in one tif (did not work for this data in M5), folder raster_data/now

x8717classes <- rgdal::readOGR("./vector_data", "changeclasses")
x8717classes <-spTransform(x8811classes,CRS="+proj=utm +zone=22 +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
library(lattice)
library(ggplot2)
library(randomForest)
x8717multidate <- superClass(x17, trainData = x8717classes, responseCol = "id")
arg <- list(at=seq(1,4,1), labels=c("for-for","for-nonfor", "nonfor-for", "nonfor-nonfor"))
#arg und col müssen noch angepasst werden
col <- c("green", "yellow", "darkgreen", "orange") #does darkgreen exist?? pick nice colours!
plot(x8717multidate$map, col=col, axis.arg=arg) #looks weird, redo classes in QGIS

#### SPATIAL DYNAMICS
## MOVING WINDOW ANALYSIS
nowsperndvi <- spectralIndices(nowsper, red="band4nowsper", 
                           nir="band5nowsper", indices ="NDVI")
plot(nowsperndvi)

focalnowsper_3 <- focal(nowsperndvi, w=matrix(1/9, ncol=3, nrow=3), fun=sd)
focalnowsper_7 <- focal(nowsperndvi, w=matrix(1/49, ncol=7, nrow=7), fun=sd)
fw <- focalWeight(nowsperndvi, 2, "Gauss")
focalnowsperGauss <- focal(nowsperndvi, w=fw)

plot(focalnowsper_3)
plot(focalnowsper_7)
plot(focalnowsperGauss)


#further packages allow different MWA - texture analysis - GLCM
library(glcm)
glcm.oldsper <- glcm(raster(oldsper, layer=1))
plot(glcm.oldsper) #mean, variance, homogeneity, contrast, dissimilarity, entropy, second moment, correlation
plot(glcm.oldsper_1$glcm_mean)
plot(glcm.oldsper_1$glcm_variance)
glcm.oldsper_3 <- glcm(raster(nowsper,layer=3))
plot(glcm.oldsper_3)
plot(glcm.oldsper_3$glcm_mean)
plot(glcm.oldsper_3$glcm_variance)

glcm.nowsper <- glcm(raster(nowsper, layer=1))
plot(glcm.nowsper) #mean, variance, homogeneity, contrast, dissimilarity, entropy, second moment, correlation
plot(glcm.nowsper_1$glcm_mean)
plot(glcm.nowsper_1$glcm_variance)
glcm.nowsper_3 <- glcm(raster(nowsper,layer=3))
plot(glcm.nowsper_3)
plot(glcm.nowsper_3$glcm_mean)
plot(glcm.nowsper_3$glcm_variance)

homogeneitydiff <- glcm.nowsper$glcm_homogeneity - glcm.oldsper$glcm_homogeneity
plot(homogeneitydiff)

### Local correlation

nowsper.ndvi <- raster("p224r63_2011")
nowsper.evi <- raster("p224r63_2011")

cor(nowsper.ndvi[], nowsper.evi[])

nowsper.correlation_ndvi_evi <- corLocal(nowsper.ndvi, nowsper.evi, ngb=11)
plot(nowsper.correlation_ndvi_evi)
