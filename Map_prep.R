# This script pulls map and bioclim rasters for our climatic predictors
#
#
#
library(raster)
library(rgdal)
library(RStoolbox)
set.seed(1987)

# Load spatial data-------------------------------------------------------

# Load saved spatial data
peru <- raster::getData("GADM", country="PER", level=0) #download country shape
projection(peru) <- CRS("+init=epsg:4326")
perualt2 <- raster("./raw_data/peru_alt.grd") #read elevation raster
climate2 <- brick("./raw_data/peru_bioclim.grd") #read in bioclim brick

plot(climate2$bio1) #check and make sure they look as they should (bio1 is temp in C * 10)
plot(perualt2)
plot(peru,add=T)

# Metadata
metadata1 <- read.csv(file = "./formatted_data/GDM_metadata.csv")

# How to download data if we don't have it saved

#countries = c("BRA", "PER", "ECU", "BOL") #download more countries just for visualization sake
#brazil = do.call("bind", lapply(countries, function(x) raster::getData('GADM', country=x, level=0)))

# Elevation data
#peru_alt <- raster::getData('alt', country = "PER", mask=T) #elevation data
#proj4string(peru_alt) <- CRS("+init=epsg:4326")
#writeRaster(peru_alt, "peru_alt.grd") #save to dropbox peru elevation raster

# Climate data
#climate <- raster::getData("worldclim", var="bio", res = 2.5) # download bioclim data
#projection(climate) <- CRS("+init=epsg:4326")
#climate2 <- crop(climate, peru)
#climate <- mask(climate2, peru) #masked and cropped to shape of peru
#proj4string(climate) == proj4string(peru_alt) #check projections
#writeRaster(climate, "peru_bioclim.grd", overwrite=T) #write out bioclim brick



# Process spatial data ----------------------------------------------------

#Create a spatial data frame that has our community localities
CommunitySpatial <- SpatialPointsDataFrame(
  matrix(c(metadata1$community.Long,metadata1$community.Lat), ncol = 2),
  data.frame(ID=1:nrow(metadata1),community = metadata1$community,
             community.elev = metadata1$community.elev),
  proj4string=CRS("+init=epsg:4326"))

#Extract mean temp and mean precip for our communities
metadata1$temp <-raster::extract(climate2$bio1, CommunitySpatial)
metadata1$precip <-raster::extract(climate2$bio12, CommunitySpatial)
metadata1$temp_var <- raster::extract(climate2$bio4, CommunitySpatial)
metadata1$precip_var <- raster::extract(climate2$bio15, CommunitySpatial)

#PCA and scale raster data for temp variables and precip variables
#temp: using BioClim vars 1 - 11
tempPCA <- rasterPCA(subset(climate2, 1:11), spca = T) #creates a raster of PCAd temp vals
summary(tempPCA$model) #PC1 has 82% of variation
tempPCA$model$loadings #if you care about these things

#Pull PC1 and PC2 values for each community and add to metadata.
metadata1$tempPCA1 <- raster::extract(tempPCA$map$PC1, CommunitySpatial)
metadata1$tempPCA2 <- raster::extract(tempPCA$map$PC2, CommunitySpatial)

plot(temp ~ tempPCA1, metadata1) #PC1 is highly correlated with mean temp


#precip
precipPCA <- rasterPCA(subset(climate2, 12:19), spca = T)
summary(precipPCA$model) #PC1 has 84% of variation

metadata1$precipPCA1 <- raster::extract(precipPCA$map$PC1, CommunitySpatial)
metadata1$precipPCA2 <- raster::extract(precipPCA$map$PC2, CommunitySpatial)
with(metadata1, plot(precip, precipPCA1)) #highly correlated

#Write out our PCA rasters for later map generation:

writeRaster(tempPCA$map$PC1, "./formatted_data/tempPC1.grd", overwrite=T)
writeRaster(precipPCA$map$PC1, "./formatted_data/precipPC1.grd", overwrite=T)

