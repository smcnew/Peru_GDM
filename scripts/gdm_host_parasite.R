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

# Data ------------------------------------------------------------
community.number <- 1:18

#Metadata
metadata <- read.csv("./formatted_data/GDM_metadata.csv")

#Host distance matrices
host_ja <- read.csv("./formatted_data/hostjaccard.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Jaccard

host_unifs <- read.csv("./formatted_data/hostUnifrac.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Unifracs (unweighted)

#Parasite distance matrices
par_bc <- read.csv("./formatted_data/parasitebraycurtis.csv") %>% rename(., community.number=X) #Bray curtis
par_unifs <- read.csv("./formatted_data/malaria.generalized.unifracs.csv") %>%
  rename(., community.number=X) #Unifracs (generalized)

#HPL distance matrices
hpl_unifs <- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/gunifs", x, "csv",sep="."))))

hpl_bc<- as.list(sapply(1:3, function(x)
  read.csv(file=paste("./formatted_data/braycurtis", x, "csv",sep="."))))


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
  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm=1000)
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
  best <- (model[[1]][2,] - sort(1:ncol(model[[1]]), decreasing = T)) %>% which.max #best model
  na.omit(model[[2]][,best])
}


#
# Parasite GDMs -----------------------------------------------------------

# Distance matrices: par_unifs (phylo turnover) par_bc (species turnover)
# Variables of interest: elevation, distance (lat/long), temp, precip, NPP,
# host richness and host turnover
# Function Arguments:
# r, response dissimilarity matrix, either host or parasite,
# e, vector of abiotic variables we care about
# d, predictor distance matrices, i.e. host turnover and parasite turnover

pe <- c("community.number", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","community.elev","total.host",
        "npp.raster")
pd <- list(host_ja, host_unifs)

plot(gdm.fun(r = par_bc, e = pe, d = pd), plot.layout = c(2,4)) #visualize full species model
plot(gdm.fun(r = par_unifs, e = pe, d = pd), plot.layout = c(2,4)) #visualize full phylo model

# Permute model to test importance of predictor variables

gdm_par_bc <- p.gdm.fun(r = par_bc, e = pe, d = pd) #full model Bray-Curtis
gdm_par_uf <- p.gdm.fun(r = par_unifs, e = pe, d = pd) #full model generalized unifracs

# Save models as R objects because they took a long ass time to run
sapply(c("gdm_par_bc", "gdm_par_uf"), function (x)
  saveRDS(get(x), file=paste("./GDM_results/", x, ".rds", sep="")))

# What factors are in the best models?
modsel(gdm_par_bc) #Distance, Precip, Host turnover
modsel(gdm_par_uf) #Distance, Precip, Host richness, Elevation


# Best models; model select and then just run the GDM with the predictors we care about.
gdm_b_par_spp <- gdm.fun(r = par_bc,
                         e = c("community.number", "community.Lat", "community.Long", "precipPCA1"), d=pd)

# Host turnover not important for phylogenetic turnover so we use gdm.fun.e (no distance matrix predictor)
gdm_b_par_phy <- gdm.fun.e(r = par_unifs,
                           e = c("community.number", "community.Lat", "community.Long",
                                 names(modsel(gdm_p_uf))[-1]))

saveRDS(gdm_b_par_spp, file = "GDM_results/gdm_b_par_spp.rds")
saveRDS(gdm_b_par_phy, file = "GDM_results/gdm_b_par_phy.rds")

#
# GDM host turnover -------------------------------------------------------------
# Distance matrices: host_unifs (phylo turnover) and host_ja (species turnover)
# Function Arguments:
# r, response dissimilarity matrix, either host or parasite,
# e, vector of abiotic variables we care about
# d, predictor distance matrices, i.e. host turnover and parasite turnover

hd <- list(par_unifs, par_bc)
he <- c("community.number", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","community.elev","parasite.richness",
        "npp.raster")

plot(gdm.fun(host_ja, e=he, d=hd), plot.layout = c(3,4)) #visualize full models
plot(gdm.fun(host_unifs, e=he, d=hd), plot.layout = c(3,4))

#Permute models to id most important variables
gdm_host_ja  <- p.gdm.fun(host_ja, he, hd)
gdm_host_uf <- p.gdm.fun(host_unifs, he, hd)

sapply(c("gdm_host_uf", "gdm_host_ja"), function (x) #save models bc permutations take so long
  saveRDS(get(x), file=paste("./GDM_results/", x, ".rds", sep="")))


#Best models (after selection)
gdm_b_host_spp <-  gdm.fun.e(host_ja, e = c("community.number",
                                           "community.Lat",
                                           "community.Long",
                                           names(modsel(gdm_host_ja))[-1]))
gdm_b_host_phy <-  gdm.fun(host_unifs, e = c("community.number",
                                            "community.Lat",
                                            "community.Long",
                                            names(modsel(gdm_host_uf))[-c(1,6)]), d=hd)

#Write out best models for spline plotting later
saveRDS(gdm_b_host_spp, file = "GDM_results/gdm_b_host_spp.rds")
saveRDS(gdm_b_host_phy, file = "GDM_results/gdm_b_host_phy.rds")

#

# GDM parasites by genus  -------------------------------------------------
#
# Response vars: unifracs (phylo turnover) and bray curtis (species turnover)
# Minor challenge: because not all communities are represented in each genus
# we have to subset metadata for just the communities we're modeling.

# Create new functions that include subsetting metadata
# Args:
# r = response matrix, e = enviro variables
# d = distance (host) matrix. g = genus: 1 = haemo, 2 = leuco, 3 = plasmo

# Basic gdm function
gdm.fun.hlp <- function (r, e, d, g) {
  dplyr::select(metadata, e) %>%
    filter(community.number %in% community.numbers[[g]]) %>%
    as.data.frame(.) %>%
    mutate(., community.number=community.numbers[[g]]) -> e.formatted  # select predictor vars
  d.formatted <- lapply(d, function(x) filter(x, community.number %in% community.numbers[[g]]))
  e.r.sitepair <- formatsitepair(bioData = r, bioFormat = 3, #response = matrix
                                 XColumn= "community.Long", YColumn = "community.Lat", #required: lat and long columns (in one or both response/predictors)
                                 siteColumn = "community.number", #required: community (site) ID, should be same in both matrices
                                 distPreds = d.formatted,
                                 predData = e.formatted, weightType="equal") #predictor data = scaled metadata
  dplyr::select(e.r.sitepair, -s1.matrix_2, -s2.matrix_2) -> e.r.sitepair #remove second matrix
  gdm(e.r.sitepair, geo=T)
}

# Permutation function
p.gdm.fun.hlp <- function (r, e, d, g) { #function takes r, response variables i.e. dissim matrix, e, vector of ecol variables we care about.
  dplyr::select(metadata, e) %>% filter(community.number %in%
                                          community.numbers[[g]]) %>% as.data.frame(.) %>%
    mutate(., community.number=community.numbers[[g]]) -> e.formatted # select predictor vars
  d.formatted <- lapply(d, function(x) filter(x, community.number %in% community.numbers[[g]]))
  e.r.sitepair <- formatsitepair(bioData = r, bioFormat = 3, #response = matrix
                                 XColumn= "community.Long", YColumn = "community.Lat", #required: lat and long columns (in one or both response/predictors)
                                 siteColumn = "community.number", #required: community (site) ID, should be same in both matrices
                                 distPreds = d.formatted,
                                 predData = e.formatted, weightType="equal") #predictor data = scaled metadata
  dplyr::select(e.r.sitepair, -s1.matrix_2, -s2.matrix_2) -> e.r.sitepair #remove second matrix  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm = 100)
  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm=1000)
}

#format data to have community number in front of all columns
h <- 1:18
l <- c(4,5,6,7,9,10,15,16)
p <- c(1,3,4,8,9,11,13)
community.numbers <- list(h,l,p)
community.numbers
names(community.numbers) <-c(rep("community.number",3))

# Response variables
gunifs_genera <- mapply(cbind, community.numbers, hpl_unifs)
gunifs_genera <- lapply(gunifs_genera, setNames, "community.number")

braycurtis_genera <- mapply(cbind, community.numbers, hpl_bc)
braycurtis_genera <- lapply(braycurtis_genera, setNames,"community.number")


# HAEMOPROTEUS
# Predictor variables
pe <- c("community.number", "community.elev", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","total.host", "npp.raster", "shannonD")
pd <- list(host_ja, host_unifs)

# Full models
gdm.fun.hlp(braycurtis_genera[[1]], e=pe, d=pd, g=1) %>% plot(., plot.layout = c(3,4))
gdm.fun.hlp(gunifs_genera[[1]], e=pe, d=pd, g=1) %>% plot(., plot.layout = c(3,4))

# Permuted models

gdm_h_bc <- p.gdm.fun.hlp(braycurtis_genera[[1]], pe, pd, 1) #full model bc
gdm_h_uf <- p.gdm.fun.hlp(gunifs_genera[[1]], pe, pd, 1) #full model gUnifs

haemo_results <- list(gdm_h_bc, gdm_h_uf)
names(haemo_results) <- c("h_gunifs", "h_braycurtis")

#save results
sapply(c("gdm_h_bc", "gdm_h_uf"), function (x)
  saveRDS(get(x), file=paste("./GDM_results/",x,".rds", sep="")))

#LEUCO (permuted)
gdm_l_bc <- p.gdm.fun.hlp(braycurtis_genera[[2]], pe, pd, 2) #full model bc
gdm_l_uf <- p.gdm.fun.hlp(gunifs_genera[[2]], pe, pd, 2) #full model gUnifs

leuco_results <- list(gdm_l_bc, gdm_l_uf)
names(leuco_results) <- c("l_gunifs", "l_braycurtis")

#save results
sapply(c("gdm_l_bc", "gdm_l_uf"), function (x)
  saveRDS(get(x), file=paste("./GDM_results/",x,".rds", sep="")))


#PLASMO
gdm_p_bc <- p.gdm.fun.hlp(braycurtis_genera[[3]], pe, pd, 3) #full model bc
gdm_p_uf <- p.gdm.fun.hlp(gunifs_genera[[3]], pe, pd, 3) #full model gUnifs

plasmo_results <- list(gdm_p_bc, gdm_p_uf)
names(plasmo_results) <- c("p_gunifs", "p_braycurtis")

#save results
sapply(c("gdm_p_bc", "gdm_p_uf"), function (x)
  saveRDS(get(x), file=paste("./GDM_results/",x,".rds", sep="")))

gdm_p_bc[[1]]


