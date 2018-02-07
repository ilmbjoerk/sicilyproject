## started on 24 november 2017
## author: Rebekka Riebl
## remote sensing project for M5 in MSc Global Change Ecoloy, winter term 2017/18
## scientific topic: Reforestation in Sicily
## with the Commune Sperlinga in the Province Enna as case study.

#needed libraries
library(raster)
library(sp)
library(ggplot2)
library(RStoolbox)
library(rgdal)
library(randomForest)
library(lattice)

# italy administrative borders
italy <- getData("GADM", country="ITA", level=3)

# study commune
#Province of Enna: lowest pop density in Sicily (wiki); 
#Sperlinga is commune with smallest population
sper <- italy[ which(italy$NAME_3=='Sperlinga'), ]
#match projections to Landsat data
sper <-spTransform(sper,CRS="+proj=utm +zone=33 
                   +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### dataset from the 80s -  3rd May 1987
setwd("/media/rebekka/9016-4EF8/GCE/3rd/M5RemoteS/ownproject/raster_data/old")
allbandsoldsper <- stack(c("band1oldsper.tif", "band2oldsper.tif", 
                           "band3oldsper.tif", "band4oldsper.tif", 
                           "band5oldsper.tif", "band7oldsper.tif"))
#crop & mask dataset to study commune
cropoldsper <- crop(allbandsoldsper, sper)
oldsper <- mask(cropoldsper, sper)


## vegetation indices
#Landsat 5: red=3, nir=4
#NDVI
oldndvi <- spectralIndices(oldsper, red="band3oldsper", 
                           nir="band4oldsper", indices ="NDVI")
plot(oldndvi)
#MSAVI
oldmsavi <- spectralIndices(oldsper, red="band3oldsper",
                            nir="band4oldsper", indices ="MSAVI")
plot(oldmsavi)

### dataset from this year - 15 May 2017
#working directory
setwd("/media/rebekka/9016-4EF8/GCE/3rd/M5RemoteS/ownproject/raster_data/now")
allbandsnowsper <- stack(c("band1nowsper.tif", "band2nowsper.tif", "band3nowsper.tif", 
                    "band4nowsper.tif", "band5nowsper.tif", "band6nowsper.tif", 
                    "band7nowsper.tif"))
cropnowsper <- crop(allbandsnowsper, sper)
nowsper <- mask(cropnowsper, sper)

#landsat 8 has other band distributions: red=4, nir=5
#NDVI
nowndvi <- spectralIndices(nowsper, red="band4nowsper", nir="band5nowsper", 
                           indices ="NDVI")
plot(nowndvi)
#MSAVI
nowmsavi <- spectralIndices(nowsper, red="band4nowsper", nir="band5nowsper", 
                           indices ="MSAVI")
plot(nowmsavi)

## classifications
arg <- list(at=seq(1,3,1), labels=c("other","non-forest", "forest"))
col <- c("blue", "yellow", "green")

#unsupervised
#1987
ucold <- unsuperClass(oldsper, nClasses=3)
plot(ucold$map, col=col)

ucold2 <- unsuperClass(oldsper, nClasses=7)
plot(ucold$map)

#2017
ucnow <- unsuperClass(nowsper, nClasses=3)
plot(ucnow$map, col=col)

ucnow2 <- unsuperClass(nowsper, nClasses=7)
plot(ucnow2$map)


#supervised
#for both
setwd("/media/rebekka/9016-4EF8/GCE/3rd/M5RemoteS/ownproject/")

#1987
tdold <- rgdal::readOGR("./vector_data", "oldsperclasses")
scold <- superClass(oldsper, trainData = tdold, responseCol = "id") #installing randomForest needed
plot(scold$map, col=col, axis.arg=arg)

#2017
tdnow <- rgdal::readOGR("./vector_data", "nowsperclasses")
scnow <- superClass(nowsper, trainData = tdold, responseCol = "id")
plot(scnow$map, col=col, axis.arg=arg)

save.image("/media/rebekka/9016-4EF8/GCE/3rd/M5RemoteS/ownproject/sicilyworkspace.RData")
