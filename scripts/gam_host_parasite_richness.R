# This script models host and parasite richness using Generalized
# Additive Models (GAMS).
#
#
library(dplyr)
library(mgcv) # the GAM package
library(raster)

set.seed(1987)

# Data --------------------------------------------------------------------

#Metadata
metadata <- read.csv("./formatted_data/GDM_metadata.csv")

#Spatial data
precipPC1 <- raster("./formatted_data/precipPC1.grd") #precip PCA raster
tempPC1 <- raster("./formatted_data/tempPC1.grd") #temp PCA raster
peru_alt <- raster("./raw_data/peru_alt.grd") #read elevation raster
npp <- raster("./formatted_data/npp.grd")

precipRe <- resample(precipPC1, peru_alt, "ngb") #resample to elev grid
tempRe <- resample(tempPC1, peru_alt, "ngb") #resample to elev grid
birdrast <- resample(birdrast, peru_alt, "ngb") #resample to elev grid
nppRe <- resample(npp, peru_alt, "ngb")

extent(peru_alt)==extent(precipRe) #should match
extent(peru_alt)==extent(tempRe) #should match



# Models ------------------------------------------------------------------

# Parasite richness
par(mfrow=c(3,3))
gam_pNULL <- gam(parasite.richness ~  s(shannonD, k=3) + s(npp.raster, k=3) +
                   s(precipPCA1, k=3) + s(tempPCA1, k=3),
                 data = metadata,  method = "REML", family="nb")
summary(gam_pNULL)

gam_p <- gam(parasite.richness ~  s(shannonD, k=3) ,
             data = metadata,  method = "REML", family="nb")
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
             data = metadata,  method = "REML", family="nb")
summary(gam_h)

pdf("./output_plots/host_gam.pdf", useDingbats = F, height=9, width=7)
par(mfrow=c(4,2), mar=c(4,5,3,5))
gam.check(gam_h)
plot(gam_h)
plot(metadata$community.elev, fitted(gam_h), xlab="Elevation")
plot(metadata$npp.raster/1000, fitted(gam_h), xlab="Net primary productivity")
dev.off()


# Projected maps  ---------------------------------------------------------
#can we predict gam onto raster?
gamhrast <- brick(addLayer(peru_alt, nppRe))
names(gamhrast) <- c("community.elev", "npp.raster")
rastTrans <- predict(gamhrast, gam_h)


pdf("./output_plots/host_map_GAM_communities.pdf", useDingbats = F)
par(mfrow=c(1,1))
plot(exp(rastTrans), col=viridis(100))
plot(CommunitySpatial, add=T, pch=21, bg="white", cex=1.3)
dev.off()

plot(birdrast, col=viridis(100))
head(gam_h)
cor(exp(raster::extract(rastTrans, CommunitySpatial)), metadata$total.host)



