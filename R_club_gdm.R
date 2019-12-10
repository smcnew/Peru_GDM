

library(gdm)
library(dplyr)
library(stringr) #used to rename categories for plots
library(scales) #for alpha of colors
library(vegan) #distance matrices
library(gplots) #heatmaps
library(raster) # for map data


# Data ------------------------------------------------------------
community.number <- 1:18

# Parasite distance matrices
par_bc <- read.csv("./formatted_data/parasitebraycurtis.csv") %>% rename(., community.number=X) #Bray curtis
par_unifs <- read.csv("./formatted_data/malaria.generalized.unifracs.csv") %>%
  rename(., community.number=X) #Unifracs (generalized)

# Host distance matrices
host_ja <- read.csv("./formatted_data/hostjaccard.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Jaccard

host_unifs <- read.csv("./formatted_data/hostUnifrac.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Unifracs (unweighted)

# Metadata
metadata <- read.csv("./formatted_data/GDM_metadata.csv")

# Spatial data ------------------------------------------------------------

birdrast <- raster("./formatted_data/birdrichness.grd") #bird species richness raster
precipPC1 <- raster("./formatted_data/precipPC1.grd") #precip PCA raster
tempPC1 <- raster("./formatted_data/tempPC1.grd") #temp PCA raster
peru_alt <- raster("./raw_data/peru_alt.grd") #read elevation raster
npp <- raster("./formatted_data/npp.grd")

precip_re <- resample(precipPC1, peru_alt, "ngb") #resample to elev grid
temp_re <- resample(tempPC1, peru_alt, "ngb") #resample to elev grid
birdrast <- resample(birdrast, peru_alt, "ngb") #resample to elev grid
npp_re <- resample(npp, peru_alt, "ngb")

extent(peru_alt) == extent(precip_re) #should match
extent(peru_alt) == extent(temp_re) #should match

# Create spatial object of locations of communities
CommunitySpatial <- SpatialPointsDataFrame(
  matrix(c(
    metadata$community.Long, metadata$community.Lat
  ), ncol = 2),
  data.frame(
    ID = 1:nrow(metadata),
    community = metadata$community,
    community.elev = metadata$community.elev
  ),
  proj4string = CRS("+init=epsg:4326")
)

# Response variables  -----------------------------------------------------

par_unifs[,-1] %>% as.matrix %>% heatmap()
par_bc[,-1] %>% as.matrix %>% heatmap()


# Predictors --------------------------------------------------------------
par(mfrow=c(3,2), mar = c(2,4,5,1))
plot(birdrast, main = "Avian spp. richness")
plot(tempPC1, main = "Temperature PC1")
plot(precipPC1, main = "Precip PC1")
plot(peru_alt, main = "Elevation")
plot(npp, main = "Net primary productivity")

head(metadata)


# GDM FUN -----------------------------------------------------------------

# Base function to run GDM
gdm.fun <- function (r, e, d) {
  dplyr::select(metadata, e) %>% as.data.frame(.) %>%
    mutate(., community.number = 1:18) -> e.formatted  # select predictor vars
  e.r.sitepair <-
    formatsitepair(
      bioData = r,
      bioFormat = 3,
      #response = matrix
      XColumn = "community.Long",
      YColumn = "community.Lat",
      #required: lat and long columns (in one or both response/predictors)
      siteColumn = "community.number",
      #required: community (site) ID, should be same in both matrices
      distPreds = d,
      predData = e.formatted,
      weightType = "equal"
    ) #predictor data = scaled metadata
  dplyr::select(e.r.sitepair,-s1.matrix_2,-s2.matrix_2) -> e.r.sitepair #remove second matrix
  gdm(e.r.sitepair, geo = T)
}

# Function to permute gdm many times in order to estimate variable importance
p.gdm.fun <- function (r, e, d) {
  dplyr::select(metadata, e) %>% as.data.frame(.) %>%
    mutate(., community.number = 1:18) -> e.formatted  # select predictor vars
  e.r.sitepair <-
    formatsitepair(
      bioData = r,
      bioFormat = 3,
      #response = matrix
      XColumn = "community.Long",
      YColumn = "community.Lat",
      #required: lat and long columns (in one or both response/predictors)
      siteColumn = "community.number",
      #required: community (site) ID, should be same in both matrices
      distPreds = d,
      predData = e.formatted,
      weightType = "equal"
    ) #predictor data = scaled metadata
  dplyr::select(e.r.sitepair,-s1.matrix_2,-s2.matrix_2) -> e.r.sitepair #remove second matrix  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm = 100)
  gdm.varImp(e.r.sitepair,
             geo = T,
             parallel = T,
             nPerm = 100) #select number of permutations
}


# Create GDM models  ------------------------------------------------------

pe <- c("community.number", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","community.elev","total.host",
        "npp.raster")
pd <- list(host_ja, host_unifs)

mod_spp <- gdm.fun(r = par_bc, e = pe, d = pd)
mod_phy <- gdm.fun(r = par_unifs, e = pe, d = pd)

dev.off()
plot(mod_spp, plot.layout = c(2,4)) #visualize full species model
plot(mod_phy, plot.layout = c(2,4)) #visualize full phylo model

summary(mod_phy)

#Permute models to calculate variable importance
gdm_par_bc <- p.gdm.fun(r = par_bc, e = pe, d = pd) #full model Bray-Curtis
gdm_par_uf <- p.gdm.fun(r = par_unifs, e = pe, d = pd) #full model generalized unifracs

#Output
gdm_par_uf

dev.off()
gdm_par_uf[[2]][,1] %>% sort(decreasing = TRUE) %>% barplot()


# Create some maps  -------------------------------------------------------

rasterdata_p <- brick(addLayer(precip_re, peru_alt, birdrast)) #rasters we want to predict onto

enviro_meta_p_p <- # assemble variables we want to use in the model
  dplyr::select(
    metadata,
    community.number,
    community.Lat,
    community.Long,
    precipPCA1,
    community.elev,
    total.host
  ) %>%
  as.data.frame(.) %>%
  mutate(., community.number = 1:18) %>% as.matrix(.)

enviro_table_p <- #format site by site table
  formatsitepair(
    bioData = par_unifs,
    bioFormat = 3,
    XColumn = "community.Long",
    YColumn = "community.Lat",
    siteColumn = "community.number",
    predData = enviro_meta_p_p,
    weightType = "equal"
  )

 map_gdm_par_phy <- gdm(enviro_table_p, geo=T) #create the model

  rastTrans <- gdm.transform( map_gdm_par_phy, rasterdata_p) #transform rasters from model
  plot(rastTrans)
  rastDat <- na.omit(getValues(rastTrans)) #get rid of NAs
  pcaSamp <- prcomp(rastDat) # Put the rasters into a PCA so we can reduce down to three dimensions
  pcaRast <- predict(rastTrans, pcaSamp, index = 1:3) # Predict model vals onto pca rasters
  pcaRast[[1]] <- (pcaRast[[1]] - pcaRast[[1]]@data@min) /
    (pcaRast[[1]]@data@max - pcaRast[[1]]@data@min) * 255
  pcaRast[[2]] <- (pcaRast[[2]] - pcaRast[[2]]@data@min) /
    (pcaRast[[2]]@data@max - pcaRast[[2]]@data@min) * 255
  pcaRast[[3]] <- (pcaRast[[3]] - pcaRast[[3]]@data@min) /
    (pcaRast[[3]]@data@max - pcaRast[[3]]@data@min) * 255
  par(mfrow = c(1, 1))
  plotRGB(pcaRast, r = 1, g = 2, b = 3)


