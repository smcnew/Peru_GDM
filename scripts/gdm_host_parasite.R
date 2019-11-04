#### R script to run GDM analysis for host and parasite turnover. Includes:
#### 1) loading of distance matrices and raster data
#### 2) using GDM functions to analyze
####

#Packages
#**Code to install GDM from archive:**
#Download package gdm from CRAN archive
#url <-   "https://cran.r-project.org/src/contrib/gdm_1.3.11.tar.gz"
#url <-"https://cran.r-project.org/src/contrib/Archive/gdm/gdm_1.3.10.tar.gz"
#gdm <- "gdm_1.3.11.tar.gz"
#download.file(url = url, destfile = gdm)
#install.packages(pkgs=gdm, type="source", repos=NULL)

library(gdm)
library(dplyr) #we're a little tidy
library(stringr) #used to rename categories for plots
library(scales) #for alpha of colors
library(vegan) #distance matrices
set.seed(1987)

# Read in data ------------------------------------------------------------
community.number <- 1:18

#Metadata
metadata <- read.csv("./formatted_data/GDM_metadata.csv")

#Parasite distance matrices
parBC <- read.csv("./formatted_data/parasitebraycurtis.csv") %>% rename(., community.number=X)
parGunifracs <- read.csv("./formatted_data/malaria.generalized.unifracs.csv") %>%
  rename(., community.number=X)


#HPL distance matrices
hplWunifs <- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/wunifs", x, "csv",sep="."))))

hplUunifs <- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/uunifs", x, "csv",sep="."))))

hplGunifs <- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/gunifs", x, "csv",sep="."))))

hplJaccard<- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/jaccards", x, "csv",sep="."))))

hplBC<- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/braycurtis", x, "csv",sep="."))))


#Host distance matrices
hostJ <- read.csv("./formatted_data/hostjaccard.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Jaccard

hostunifs <- read.csv("./formatted_data/hostUnifrac.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Unifracs

#Spatial data
#birdrast <- raster("./formatted_data/birdrichness.grd") #bird species richness raster
precipPC1 <- raster("./formatted_data/precipPC1.grd") #precip PCA raster
tempPC1 <- raster("./formatted_data/tempPC1.grd") #temp PCA raster
peru_alt <- raster("./raw_data/peru_alt.grd") #read elevation raster
npp <- raster("./formatted_data/npp.grd")


# Functions ---------------------------------------------------------------
# Create functions to run GDMs because we'll need to run many different models.
# Basic form of GDM:
# turnover ~ rain + temp + elev + npp + parasite richness + symbiont turnover
# where "symbiont" means predicting host turnover as function of parasite turnover
# and vice versa.


# Function Arguments:
# r, response dissimilarity matrix, either host or parasite,
# e, vector of abiotic variables we care about
# d, predictor distance matrices, i.e. host turnover and parasite turnover

# Base function to run GDM
gdm.fun <- function (r, e, d) {
  dplyr::select(metadata, e) %>% as.data.frame(.) %>%
    mutate(., community.number=1:18) -> e.formatted  # select predictor vars
  e.r.sitepair <- formatsitepair(bioData = r, bioFormat = 3, #response = matrix
                                 XColumn= "community.Long", YColumn = "community.Lat", #required: lat and long columns (in one or both response/predictors)
                                 siteColumn = "community.number", #required: community (site) ID, should be same in both matrices
                                 distPreds = d,
                                 predData = e.formatted, weightType="equal") #predictor data = scaled metadata
  dplyr::select(e.r.sitepair, -s1.matrix_2, -s2.matrix_2) -> e.r.sitepair #remove second matrix
  gdm(e.r.sitepair, geo=T)
  #plot(gdm(e.r.sitepair, geo=T))
}

# Function to permute gdm many times in order to estimate variable importance
p.gdm.fun <- function (r, e, d) {
  dplyr::select(metadata, e) %>% as.data.frame(.) %>%
    mutate(., community.number=1:18) -> e.formatted  # select predictor vars
  e.r.sitepair <- formatsitepair(bioData = r, bioFormat = 3, #response = matrix
                                 XColumn= "community.Long", YColumn = "community.Lat", #required: lat and long columns (in one or both response/predictors)
                                 siteColumn = "community.number", #required: community (site) ID, should be same in both matrices
                                 distPreds = d,
                                 predData = e.formatted, weightType="equal") #predictor data = scaled metadata
  dplyr::select(e.r.sitepair, -s1.matrix_2, -s2.matrix_2) -> e.r.sitepair #remove second matrix  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm = 100)
  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm=100)
}

# Environmental function, only takes climate/npp data, (no distance matrices)
gdm.fun.e <- function (r, e) {
  dplyr::select(metadata, e) %>% as.data.frame(.) %>%
    mutate(., community.number=1:18) -> e.formatted  # select predictor vars
  e.r.sitepair <- formatsitepair(bioData = r, bioFormat = 3, #response = matrix
                                 XColumn= "community.Long", YColumn = "community.Lat", #required: lat and long columns (in one or both response/predictors)
                                 siteColumn = "community.number", #required: community (site) ID, should be same in both matrices
                                 predData = e.formatted, weightType="equal") #predictor data = scaled metadata
  gdm(e.r.sitepair, geo=T)
}

# A function to select the best model from a permutated bunch,
# based on most variance explained with fewest predictors
modsel <- function(model){
  best <- (model[[1]][2,] - sort(1:ncol(model[[1]]), decreasing=T)) %>% which.max #best model
  na.omit(model[[2]][,best])
}
modsel(pb) %>% names()
modsel(gdmgU) %>% round(digits=2) %>% scale() %>% rownames()




# Parasite GDMs -----------------------------------------------------------

# Distance matrices: parWunifracs (phylo turnover) parBC (species turnover)
# Variables of interest: elevation, distance (lat/long), temp, precip, NPP,
# host richness and host turnover
# format input vars
pe <- c("community.number", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","community.elev","total.host",
        "npp.raster")
pd <- list(hostJ, hostunifs)

plot(gdm.fun(parBC, e=pe, d=pd), plot.layout = c(2,4)) #visualize full species model
plot(gdm.fun(parGunifracs, e=pe, d=pd), plot.layout = c(2,4)) #visualize full phylo model

# Permute model to test importance of predictor variables

gdm_p_bc <- p.gdm.fun(parBC, pe, pd) #full model Bray-Curtis
gdm_p_uf <- p.gdm.fun(parGunifracs, pe, pd) #full model generalized unifracs

# Save models as R objects because they took a long ass time to run
sapply(c("gdm_p_bc", "gdm_p_uf"), function (x)
  saveRDS(get(x), file=paste("./GDM_results/", x, ".Rdata", sep="")))

# Select best model
modsel(gdm_p_bc) #Distance, Precip, Host turnover
modsel(gdm_p_uf) #Distance, Precip, Host richness, Elevation


# Best models; model select and then just run the GDM with the predictors we care about.
gdm_b_par_spp <- gdm.fun(parBC, e= c("community.number", "community.Lat", "community.Long", "precipPCA1"), d=pd)

# Host turnover not important for phylogenetic turnover so we use gdm.fun.e (no distance matrix predictor)
gdm_b_par_phy <- gdm.fun.e(parGunifracs, e= c("community.number", "community.Lat", "community.Long", names(modsel(gdm_p_uf))[-1]))

summary(gdm_b_par_phy)
