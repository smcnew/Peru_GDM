# This script formats spatial data, pulls climatic and npp values for each
# community, and creates maps of temperature, precip, elevation, and NPP.
# **Don't forget to download NPP file, it's not in the Github repo (too big).**

library(raster)
library(rgdal)
library(RStoolbox)
library(prettymapr)
library(viridis)
library(rasterVis)
library(dplyr)
set.seed(1987)

# Load data-------------------------------------------------------

# Spatial data
peru <- raster::getData("GADM", country="PER", level=0) #download country shape
projection(peru) <- CRS("+init=epsg:4326")
perualt2 <- raster("./raw_data/peru_alt.grd") #read elevation raster
climate2 <- brick("./raw_data/peru_bioclim.grd") #read in bioclim brick
npp <- raster("./formatted_data/npp.grd") # see below if you need to make this

tempPCA <- raster("./formatted_data/tempPC1.grd") #see below if you need to make this
precipPCA <- raster("./formatted_data/precipPC1.grd") #see below if you need to make this

#Bird richness from Birdlife Distributions
birdrast <- raster("formatted_data/birdrichness.grd")


plot(climate2$bio1) #check and make sure they look as they should (bio1 is temp in C * 10)
plot(perualt2)
plot(peru,add=T)

# Metadata
metadata1 <- read.csv(file = "./formatted_data/GDM_metadata.csv")

#Sampled hosts
sampledhosts <- read.csv("./raw_data/malaria-sabrinacopy.csv")

#Sampled parasites
inputhaplos <- read.csv("./raw_data/grouped-haplos.csv", row.names = 1)


# How to download spatial data if we don't have it saved

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


# Net Primary Productivity  -----------------------------------------------

### I tried many methods to create a raster of NPP that were all fairly
### complicated and a pain in the ass. I tried using the MODIS package but
### ran into GDAL programs (issues recognizing the .hdf format). I also tried to
### use the MODIS subset website but the extent was limited to downloading 100 km
### x 100 km blocks (too small!). I downloaded .hdf files by hand and processed
### them using HEG (also a pain!) but they looked super weird. The LAADS DAAC
### website is super unhelpful. Then I sort of realized that the NPP data
### maybe are not on the LAADS website because of "cloud contamination" issues.
### Eventually I found a good link to the U Montana http://files.ntsg.umt.edu/
### Website that has a different filtering process for NPP and managed to download
### average from 2000 through 2015 NPP.
### IMPORTANT: this file is too big to sync to github, so you will have to
### download it yourself:
### download.file("http://files.ntsg.umt.edu/data/NTSG_Products/MOD17/GeoTIFF/MOD17A3/GeoTIFF_30arcsec/MOD17A3_Science_NPP_mean_00_15.tif", destfile = "./raw_data/namethefile")

npp <- raster("./raw_data/MOD17A3_Science_NPP_mean_00_15.tif")
npp <- crop(npp, peru) # Mask and crop first because tif is of whole world
npp <- mask(npp, peru)
npp[is.na(npp[])] <- 0 # Turn Nas into 0s
npp <- crop(npp, peru) # Mask and crop again
npp <- mask(npp, peru)

plot(npp, col=viridis(100)) ### NPP are in kg carbon /m2 *10000


metadata1$npp.raster <- raster::extract(npp, CommunitySpatial) #pull vals for communities
writeRaster(npp, "./formatted_data/npp.grd", overwrite=T)
npp <- raster("./formatted_data/npp.grd")


# A bigger regional map ---------------------------------------------------
countries = c("PER", "BOL", "ECU", "COL", "VEN","GUY","SUR") #download more countries just for visualization sake
region = do.call("bind", lapply(countries,
                                function(x)
                                  raster::getData('GADM', country = x, level = 0)))

peru_alt <- raster::getData('alt', country = "PER", mask=F) #elevation data
bol_alt <- raster::getData('alt', country = "BOL", mask=F)
col_alt <- raster::getData('alt', country = "COL", mask=F)
vez_alt <- raster::getData('alt', country = "VEN", mask = T)
brazil_alt <- raster::getData('alt', country = "BRA", mask = T)
guy_alt <- raster::getData('alt', country = "GUY", mask = T)
sur_alt <- raster::getData('alt', country = "SUR", mask = T)
reg_alt <- merge(bol_alt, peru_alt, brazil_alt, col_alt, vez_alt, guy_alt, sur_alt)
proj4string(reg_alt) <- CRS("+init=epsg:4326")



# crop to create square based on rough google image
# xmin = -81.6
# xmax = -59.51
# y min = -18.9
# y max = 4.43

#crop raster and country lines to a smaller box
ex <- extent(c(-81.6, -49, -18.9, 4.43))
reg_alt <- crop(reg_alt, ex)
extent(reg_alt)
region <- crop(region, ex)

cuts <- seq(from = 0, to = 6000, by = 50) #for consistent elevation across maps
pdf("output_plots/region_map4.pdf", height = 8, width = 10.5)
plot(reg_alt, col = colorRampPalette(c("grey20","grey80"))(100),
     cex.axis=1.5, ext=ex)
plot(perualt2, breaks = cuts, col= viridis(length(cuts)), add=T)
plot(region, add=T, lwd = 2, border = "black")
plot(peru, add = T, border = "black", lwd = 2.8)
#addscalebar(plotepsg = 4326, label.cex = 2)
dev.off()


extent(reg_alt)


pdf("output_plots/peru_outline.pdf", height=5, width = 5, useDingbats = F)
plot(peru)
dev.off()

# Maps for figures --------------------------------------------------------

#Make a nice map of communities
pdf("./output_plots/samplingmap.pdf", useDingbats = F)
plot(perualt2, col= viridis(100))
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)
addscalebar(plotepsg = 4326)
#with(CommunitySpatial, text(CommunitySpatial,
#  labels=CommunitySpatial$ID, pos=4, cex=1.3))

dev.off()
## Need to make a consistent color scheme for elevation maps
cuts <- seq(from = 0, to = 6000, by = 50)


#plot our abiotic maps for export (including region map below)
pdf("./output_plots/rich_temp_precip_npp_raster_scaleb.pdf",
    useDingbats = F, height = 8, width = 5.5)
par(mfcol=c(3,2), mar = c(1,2,4,2))

plot(reg_alt, breaks = cuts, col= viridis(length(cuts)),
     cex.axis=1.5, ext=ex)
plot(region, add=T, lwd = 2, border = "grey41")
plot(peru, add = T, border = "grey10", lwd = 2.8)

plot(perualt2, breaks = cuts, col= viridis(length(cuts)), xaxt = "n", yaxt = "n", # cex.axis=1.5,
     box = FALSE, axes = F)
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)

plot(tempPCA$map$PC1, col=viridis(100),
     main= "Temperature", cex.main = 1.6, box = FALSE, axes = F)
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)


plot(precipPCA$map$PC1, col=viridis(100), xaxt = "n", yaxt = "n", #cex.axis=1.5,
     main="Precipitation", cex.main = 1.6, box = FALSE, axes = F)
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)


plot(npp/10000, col=viridis(100), xaxt = "n", yaxt = "n", #cex.axis=1.5,
     main="Net Primary Productivity", cex.main = 1.6, box = FALSE, axes = F)
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)


plot(birdrast_gam, col = viridis(100), xaxt = "n", yaxt = "n", #cex.axis=1.5,
     main = "Bird species richness", cex.main = 1.6, box = FALSE, axes = F)
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)
addscalebar(plotepsg = 4326)
dev.off()

# Shannon richness of host community (index of sampling effort)-----------------
# We want to use sampling effort as a predictor in parasite GDMs.
# Approach: take list of sampled hosts, assign them to community based on
# proximity to our community centers,

# Create unique lat long id, rounding both so we can match to community center
sampledhosts$latlong <- as.factor(paste(round(sampledhosts$neg.degrees.Lat,3),
                                        round(sampledhosts$neg.degrees.Long,3), sep="."))

community.lat.long.list <- unique(dplyr::select(inputhaplos, "latlong", "community.number"))

sampledhosts <- merge(sampledhosts, community.lat.long.list,
                      all.x=F, by="latlong") #add community numbers to host sampling

#community x host abundance  matrix
hostsurvey <- as.data.frame.matrix(table(sampledhosts$community.number, sampledhosts$Taxon))
metadata1$shannonD <- diversity(hostsurvey, "shannon")


# Write out metadata  -----------------------------------------------------

write.csv(metadata1, file = "./formatted_data/GDM_metadata.csv", row.names=F)


# Map of sampling localities ----------------------------------------------
sampling <- inputhaplos %>% select (latlong,
                     neg.degrees.Lat,
                     neg.degrees.Long,
                     Elev,
                     community.number,
                     Local,
                     Dept) %>%
  mutate(., latlongcomm = paste(round(neg.degrees.Lat, 2),
                             round(neg.degrees.Long, 2),
                             community.number, sep="."))


#-12.73.-73.96.2
Ns <- table(sampling$latlongcomm) %>% as.data.frame %>% rename(latlongcomm = Var1)
sampling <- sampling %>% distinct(latlongcomm, .keep_all=T)
filter(sampling, is.na(Elev))

#fix missing elevations by manually googling coordinates
sampling$Elev[sampling$latlongcomm == "-11.47.-74.79.2"] <- 1744
sampling$Elev[sampling$latlongcomm == "-12.73.-73.96.2"] <- 1837
sampling$Elev[sampling$latlongcomm == "-10.44.-75.35.2"] <- 2100
sampling$Elev[sampling$latlongcomm == "-11.51.-74.86.2"] <- 2640
sampling$Elev[sampling$latlongcomm == "-10.42.-75.32.2"] <- 1398
sampling$Elev[sampling$latlongcomm == "-12.79.-74.2"] <- 2750
sampling$Elev[sampling$latlongcomm == "-6.84.-77.49.4"] <- 1909
sampling$Elev[sampling$latlongcomm == "-6.72.-77.58.4"] <- 1730
sampling$Elev[sampling$latlongcomm == "-14.49.-69.28.9"] <- 3100
sampling$Elev[sampling$latlongcomm == "-13.1.-71.57.9"] <- 1963
sampling$Elev[sampling$latlongcomm == "-13.17.-71.6.9"] <- 2737
sampling$Elev[sampling$latlongcomm == "-13.16.-71.59.9"] <- 2600

#
filter(sampling, is.na(Elev))

sampling <- merge(sampling,Ns) %>%
  mutate( col = cut(
    Elev,
    breaks = seq(from = 0, to = 5100, by = 500),
    dig.lab = 5),
    size = cut(Freq, breaks = seq(from = 0, to = 120, by = 10)))
filter(sampling, community.number ==2)

max(sampling$Elev, na.rm=T)
local_spat <- SpatialPointsDataFrame(
  matrix(c(
    sampling$neg.degrees.Long, sampling$neg.degrees.Lat), ncol = 2),
  data.frame(
    Elev = sampling$Elev,
    N = sampling$Freq,
    community.number = sampling$community.number,
    col = sampling$col,
    size = sampling$size),
  proj4string = CRS("+init=epsg:4326")
)


length(levels(sampling$col))
cexs <- seq(from = 0.5, to = 1.5, length.out=length(levels(sampling$size)))
vircols <- viridis_pal(option = "B", end = 0.7)(length(levels(sampling$col)))

pdf("output_plots/locals_map2.pdf", useDingbats = F)
par(mfrow=c(5,4), mar=c(0,0,0,0))
for (i in 1:18) {
  plot(perualt2, legend = F, axes = F)
  comm <- local_spat[local_spat$community.number==i,]
  plot(comm, add=T, pch=1, cex = cexs[comm$size], lwd = 1,
       col = vircols[comm$col])
  mtext(i, side = 1, line = -2, adj = 0.1, cex = 1.5)
}
plot(NA,NA, axes=F)
legend("center", legend = levels(sampling$col), title = "Elevation (m)", pch = 19, pt.cex = 2, col = vircols, cex = 0.8)

plot(NA,NA, axes=F)
legend("center", legend = levels(sampling$size), title = "Samples (N)", pch = 21, pt.cex = cexs, cex = 0.7)
dev.off()

filter(sampling, community.number ==4 ) %>% select(Local, neg.degrees.Lat, neg.degrees.Long, latlongcomm) %>%
  select(latlongcomm) %>% distinct



