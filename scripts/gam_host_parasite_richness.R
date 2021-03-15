# This script models host and parasite richness using Generalized
# Additive Models (GAMS).
#
#
library(dplyr)
library(mgcv) # the GAM package
library(raster) # to deal with spatial data
library(viridis) # map cols
library(sjPlot)
set.seed(1987)

# Data --------------------------------------------------------------------

#Metadata
metadata <- read.csv("./formatted_data/GDM_metadata.csv")

#Spatial data
precipPC1 <- raster("./formatted_data/precipPC1.grd") #precip PCA raster
tempPC1 <- raster("./formatted_data/tempPC1.grd") #temp PCA raster
peru_alt <- raster("./raw_data/peru_alt.grd") #read elevation raster
npp <- raster("./formatted_data/npp.grd")
birdrast <- raster("./formatted_data/birdrichness.grd") #bird species richness raster from birdlife dist

precipRe <- resample(precipPC1, peru_alt, "ngb") #resample to elev grid
tempRe <- resample(tempPC1, peru_alt, "ngb") #resample to elev grid
birdrast <- resample(birdrast, peru_alt, "ngb") #resample to elev grid
nppRe <- resample(npp, peru_alt, "ngb")

extent(peru_alt)==extent(precipRe) #should match
extent(peru_alt)==extent(tempRe) #should match

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



# Models ------------------------------------------------------------------

# Parasite richness
par(mfrow=c(3,3))
gam_pNULL <- gam(parasite.richness ~  s(shannonD, k=3) + s(npp.raster, k=3) +
                   s(precipPCA1, k=3) + s(tempPCA1, k=3) + s(total.host),
                 data = metadata,  method = "REML", family="nb")
summary(gam_pNULL)
tab_model(gam_pNULL)


gam_p <- gam(parasite.richness ~  s(shannonD, k=3) ,
             data = metadata,  method = "REML", family="nb")
tab_model(gam_p)
summary(gam_p)
pdf("./output_plots/parasite_gam.pdf", useDingbats = F)
par(mfrow=c(3,2))
gam.check(gam_p)
plot(gam_p, xlab="Sampled host community diversity (Shannon's D)")
plot(metadata$shannonD, fitted(gam_p), xlab="Sampled host community diversity (Shannon's D)" )
dev.off()



# Host richness
par(mfrow=c(2,2))

gam_h <- gam(total.host ~ s(community.elev, k=3) + #s(community.Lat) +
               #               s(tempPCA1, k=3) + s(precipPCA1, k=3) +
               s(npp.raster, k=3),
             data = metadata,  method = "REML", family="nb") # use a negative binomial because var > mean
tab_model(gam_h)
summary(gam_h)
pdf("./output_plots/host_gam.pdf", useDingbats = F, height=9, width=7)
par(mfrow=c(4,2), mar=c(4,5,3,5))
gam.check(gam_h)
plot(gam_h)
plot(metadata$community.elev, fitted(gam_h), xlab="Elevation (m)")
plot(metadata$npp.raster/10000, fitted(gam_h), xlab="Net primary productivity Kg C/m2")
dev.off()
?plot.gam

#logit
# inverse logit formula: exp(x)/(1+exp(x))
x <- -0.5
exp(x)/(1+exp(x))

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
logit2prob(0)
logit2odds <- function(logit){
  odds <- exp(logit)
  return(odds)
}
logit2odds(0)
logit2odds(-.5)
logit2odds(-1)

# Projected maps  ---------------------------------------------------------
#can we predict gam onto raster?
gamhrast <- brick(addLayer(peru_alt, nppRe))
names(gamhrast) <- c("community.elev", "npp.raster")
rastTrans <-
  predict(gamhrast, gam_h) #predictions come as log(richness) because of neg binom model


pdf("./output_plots/host_map_GAM_communities.pdf", useDingbats = F)
par(mfrow = c(1, 1))
plot(exp(rastTrans), col = viridis(100)) #exp. to back transform from log-link (negative binomial)
plot(
  CommunitySpatial,
  add = T,
  pch = 21,
  bg = "white",
  cex = 1.3
)
dev.off()

plot(birdrast, col = viridis(100))
head(gam_h)

cor(exp(raster::extract(rastTrans, CommunitySpatial)), metadata$total.host) # highly corrleated with our estimation of richness

# Save raster
writeRaster(exp(rastTrans), "./formatted_data/gam_bird_rich.grd", overwrite=T)






# scratch -----------------------------------------------------------------

plot(parasite.richness ~ parasite.observed, metadata, ylab = "estimated parasite richness", xlab = "observed parasite richness")
metadata
infections <- peruhaploabun %>% rowSums()

plot(x = infections, y = metadata$parasite.richness, ylab = "estimated parasite richness")

test <- matrix(nrow  = 100, ncol = 18)
for ( i in 1:100){
peruspec.rar <- rrarefy(peruhaploabun, 44)
test[i,] <- peruspec.rar %>% apply(., 1, function(x) sum(x > 0))
}
rare_rich <- apply(test, 2, mean)

plot(x = rare_rich, y = metadata$parasite.richness, ylab = "estimated parasite richness",
     xlab = "Rarefied richness")

gam(rare_rich ~  s(shannonD, k=3),
    data = metadata,  method = "REML", family="nb") %>% summary()

gam(rare_rich ~  s(total.host),
    data = metadata,  method = "REML", family="nb") %>% summary()
gam(rare_rich ~  s(shannonD, k=3) + s(npp.raster, k=3) +
                   s(precipPCA1, k=3) + s(tempPCA1, k=3) + s(total.host),
                 data = metadata,  method = "REML", family="nb") %>% summary()
