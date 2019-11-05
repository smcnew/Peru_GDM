# Script to generate raster of bird diversity based on species distribution maps
#
library(sp)
library(rgdal)
library(geosphere)
library(maptools)
library(sf)
library(raster)
library(dplyr)
library(letsR)
library(viridis)

peru <- raster::getData("GADM", country="PER", level=0) #country shape
projection(peru) <- CRS("+init=epsg:4326")


setwd("/Users/sabrinamcnew/Documents")
species<-readOGR("./Bird_shapes")  ## shape files for all the species; get these from Birdlife
species <- spTransform (species, CRS("+init=epsg:4326"))

#using function letsR in order to make bird richness raster:

#find peru extent
minx <- NULL
maxx <- NULL
miny <- NULL
maxy <- NULL
for(i in 1:68){
  minx[i] <- peru@polygons[[1]]@Polygons[[i]]@coords[,1] %>%min()
  maxx[i] <- peru@polygons[[1]]@Polygons[[i]]@coords[,1] %>%max()
  miny[i] <- peru@polygons[[1]]@Polygons[[i]]@coords[,2] %>%min()
  maxy[i] <- peru@polygons[[1]]@Polygons[[i]]@coords[,2] %>%max()
}
minlong <- min(minx)
maxlong <- max(maxx)
minlat <- min(miny)
maxlat <- max(maxy)


#Function:
#At 0.01 degree this will take about 24h on desktop mac
PAM_birds <- lets.presab(species, res=0.01, count = T,
                         crs.grid = CRS("+init=epsg:4326"),
                         xmn = minlong, xmx = maxlong, ymn = minlat, ymx = maxlat)
summary(PAM_birds)

plot(PAM_birds$Richness_Raster,col=c("white",viridis(100)))

writeRaster(mask(PAM_birds$Richness_Raster, mask=peru), "birdrichness.grd",
            overwrite=T) #save
birdrast <- raster("birdrichness.grd")

write.csv(PAM_birds$Presence_and_Absence_Matrix, "bird_spp_richness_matrix.csv")

pdf("./output_plots/bird_richness.pdf")
plot(birdrast, col=viridis(100))
dev.off()
