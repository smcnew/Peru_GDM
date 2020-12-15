library(gdm)
library(geosphere)
library(dplyr) #we're a little tidy
library(stringr) #used to rename categories for plots
library(scales) #for alpha of colors
library(vegan) #distance matrices
library(gplots) #heatmaps
library(ape) #tree
library(GUniFrac) #for calculating unifrac distance
library(emmeans) #host traits
library(lwgeom)
library(FSA) #for dunn/kruskal wallace test
set.seed(1987)


# Data --------------------------------------------------------------------


metadata <- read.csv("./formatted_data/GDM_metadata.csv")
inputhaplos <- read.csv("./raw_data/grouped-haplos.csv", row.names = 1)
malariatree <- read.nexus("./raw_data/perumalariatree.nex")
sampledhosts <- read.csv("./raw_data/malaria-sabrinacopy.csv")

CommunitySpatial <- SpatialPointsDataFrame(
  matrix(c(metadata$community.Long,metadata$community.Lat), ncol = 2),
  data.frame(ID=1:nrow(metadata),community.number = metadata$community.number,
             community.elev = metadata$community.elev),
  proj4string=CRS("+init=epsg:4326"))


#Distance matrices
par_bc <- read.csv("./formatted_data/parasitebraycurtis.csv") %>% rename(., community.number=X) #Bray curtis
par_unifs <- read.csv("./formatted_data/malaria.generalized.unifracs.csv") %>%
  rename(., community.number=X) #Unifracs (generalized)
host_ja <- read.csv("./formatted_data/hostjaccard.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Jaccard

host_unifs <- read.csv("./formatted_data/hostUnifrac.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Unifracs (unweighted)

# Functions ---------------------------------------------------------------
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
} #base GDM

p.gdm.fun <- function (r, e, d) { #permuted GDM
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

# Make a plot
varcols = c(
  "#332288",
  "#BBBBBB",
  "#66CCEE",
  "#AA3377",
  "#CCBB44",
  "#228833",
  "#4477AA",
  "#EE6677",
  "#882255"
)

# Assign colors to predictor variables. Last color gets repeated because it should
# represent both host richness and parasite richness (i.e. effect of "symbiont richness"
# on turnover)
colmatch <-
  data.frame(cols = c(varcols, varcols[9]), factors = c(
    c(
      "precipPCA1",
      "Geographic",
      "community.elev",
      "tempPCA1",
      "matrix_1",
      "npp.raster",
      "shannonD",
      "total.host",
      "parasite.richness",
      "host_mass"
    )
  ))



gbarplot2 <- function(model, title) {
  par(mar = c(5, 5, 4, 2))
  # Choose best model based on most variance explained using fewest predictors
  modsel <- function(model) {
    best <-
      (model[[1]][2, ] - sort(1:ncol(model[[1]]), decreasing = T)) %>% which.max #best model
    na.omit(model[[2]][, best])
  }

  vals <-
    modsel(model) # Extract variable importance from best model
  factor_cols <-
    colmatch$cols[match(names(sort(vals, decreasing = T)),
                        colmatch$factors)] %>% as.character # Assign col
  # Need to rename factor names to make them consistent and clean
  namestring <- str_replace_all(
    names(sort(vals, decreasing = T)),
    c(
      "Geographic" = "Dist",
      "tempPCA1" = "Temp",
      "precipPCA1" = "Precip",
      "community.elev" = "Elev",
      "total.host" = "Symb rich",
      "matrix_1" = "Symb turn",
      "shannonD" = "Samp Effort",
      "parasite.richness" = "Symb rich",
      "npp.raster" = "NPP"
    )
  )
  # Create a barplot of variable importance from most to least important.
  mp <-
    barplot(
      as.numeric(scale(sort(vals, decreasing = T), center = F)), #scale variables
      #as.numeric(sort(vals, decreasing = T)),
      col = factor_cols,
      main = title,
      axis = F,
      axisnames = F,
      cex.axis = 1.1,
      ylab = "Scaled importance",
      cex.lab = 1.3
    )
  mtext(
    paste("Variance explained = ", round(model[[1]][2], 2), "%", sep = ""),
    side = 3,
    line = 0.4,
    cex = 0.7
  )
  text(
    mp,
    par("usr")[3],
    labels = namestring,
    srt = 45,
    #par("usr")[3] = ymin of plot
    adj = c(1.1, 1.1),
    xpd = T,
    cex = 1.3
  )
}

# Map illustrating distance -----------------------------------------------

library(sf)
library(geosphere)

km_dist <- distm(CommunitySpatial, fun = distHaversine)/1000
metadata1

host_ja

pdf("output_plots/example_distance.pdf", useDingbats = F)
plot(perualt2)
points(CommunitySpatial[c(1,7,10),], pch = 19)
text(CommunitySpatial[c(1,7,10),], labels = c("    1", "    7", "      10"))
dev.off()

inputhaplos %>% filter(community.number ==10)



# Would excluding rare parasites include model fit ------------------------

peruhaploabun <- as.data.frame.matrix(
  table(inputhaplos$community.number, inputhaplos$Lineage.name))
t_peruhaploabun <- t(peruhaploabun) %>% as.data.frame %>% mutate(n_inf = rowSums(.[,1:18]))
filter(t_peruhaploabun, n_inf < 5) %>% nrow()
nrow(t_peruhaploabun)

t_peruhaploabun <- t(peruhaploabun) %>% as.data.frame %>% mutate(n_inf = rowSums(.[,1:18])) %>%
  filter(n_inf > 4) %>% select (-n_inf) %>% t()

rowSums(t_peruhaploabun) #min infection = 25
infection.min <- min(rowSums(t_peruhaploabun))

common_bc <- avgdist(t_peruhaploabun, sample=infection.min, dmethod="bray",
                  iterations = 1000) %>% as.matrix %>% as.data.frame()
pdf("output_plots/common_parasite_bc.pdf")
heatmap.2(as.matrix(common_bc[,-1]),  main="Parasite Bray Curtis", trace="none")
dev.off()
heatmap.2(as.matrix(par_bc[,-1]),  main="Parasite Bray Curtis", trace="none")

#do a gdm of the new matrix

pe <- c("community.number", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","community.elev","total.host",
        "npp.raster")
pd <- list(host_ja, host_unifs)

common_bc <- mutate(common_bc, community.number = 1:18) %>% select(community.number, everything())

gdm.fun(r = common_bc, e = pe, d = pd)#visualize full species model

# Permute model to test importance of predictor variables

gdm_par_bc_common <- p.gdm.fun(r = common_bc, e = pe, d = pd) #full model Bray-Curtis

pdf("output_plots/common_parasites_barplot.pdf", width = 5, height = 11)
par(mfrow = c(5, 2))
gbarplot2(gdm_par_bc_common, "Common parasites")
dev.off()
#

# Rarefy less -------------------------------------------------------------
peruhaploabun <- as.data.frame.matrix(
  table(inputhaplos$community.number, inputhaplos$Lineage.name))
rowSums(peruhaploabun) #min infection = 25
rare_less_bc <- avgdist(peruhaploabun, sample=70, dmethod="bray",
                     iterations = 1000) %>% as.matrix %>% as.data.frame()


str(rare_less_bc)
# Communities dropped because of low n: 1, 12, 14, 16, 17, 18

coms <- 1:18
coms <- coms[!coms %in% c(1, 12, 14, 16, 17, 18)]

rare_less_bc <- mutate(rare_less_bc, community.number = coms) %>% select(community.number, everything())


# Run model rarefying to 70. Long format because we need to subset predictors to remove dropped communities
dplyr::select(metadata, pe) %>% filter(community.number %in%
                                         rare_less_bc$community.number) %>% as.data.frame(.) %>%
  mutate(., community.number = rare_less_bc$community.number) -> e.formatted # select predictor vars
d.formatted <-
  lapply(pd, function(x)
    filter(x, community.number %in% rare_less_bc$community.number))
e.r.sitepair <-
  formatsitepair(
    bioData = rare_less_bc,
    bioFormat = 3,
    #response = matrix
    XColumn = "community.Long",
    YColumn = "community.Lat",
    #required: lat and long columns (in one or both response/predictors)
    siteColumn = "community.number",
    #required: community (site) ID, should be same in both matrices
    distPreds = d.formatted,
    predData = e.formatted,
    weightType = "equal"
  ) #predictor data = scaled metadata
dplyr::select(e.r.sitepair,-s1.matrix_2,-s2.matrix_2) -> e.r.sitepair #remove second matrix  gdm.varImp(e.r.sitepair, geo=T, parallel=T, nPerm = 100)
rare_70_parasite <- gdm.varImp(e.r.sitepair,
                               geo = T,
                               parallel = T,
                               nPerm = 100)



# Rarefy unifs
peruspec.rar <- rrarefy(peruhaploabun, infection.min) #this rarefies dataset once.
infection.min = 70
gUnifs <- function () {
  rrarefy(peruhaploabun, infection.min) %>%
    GUniFrac(., malariatree, alpha=c(0, 0.5, 1)) %>%
    .$unifracs %>% .[, , "d_0.5"] %>% return(.)
}

# Make a holding dataframe, so we can rarefy and calculate unifracs 1000 times
meangUnifs <- data.frame(matrix(0,nrow=18,ncol=18))
sapply(1:1000, function(x) {x; meangUnifs <<- meangUnifs + gUnifs()})
meangUnifs <- meangUnifs/1000

# Run phylo model

meangUnifs <- as.data.frame(meangUnifs) %>% mutate(community.number = 1:18) %>% select(community.number, everything())
gdm_par_unifs_common <- p.gdm.fun(r = meangUnifs, e = pe, d = pd)



# Rarefy without dropping communities
rrarefy(peruhaploabun, 70) %>% vegdist(., method = "bray", binary = T, diag = T) %>% as.matrix
b <- as.matrix(vegdist(peruspec.rar, method="bray", binary=T, diag=T))
gUnifs <- function () {
  rrarefy(peruhaploabun, infection.min) %>%
    GUniFrac(., malariatree, alpha=c(0, 0.5, 1)) %>%
    .$unifracs %>% .[, , "d_0.5"] %>% return(.)
}

rarefy_bc_70 <- function() {
  rrarefy(peruhaploabun, 70) %>% vegdist(., method = "bray", binary = T, diag = T) %>% as.matrix %>%
    return(.)
}
rarefy_bc_70()
meanbc <- data.frame(matrix(0,nrow=18,ncol=18))
sapply(1:1000, function(x) {x; meanbc <<- meanbc + rarefy_bc_70() })
meanbc <- meanbc/1000

meanbc <- as.data.frame(meanbc) %>% mutate(community.number = 1:18) %>% select (community.number, everything())
gdm_par_bc_common2 <- p.gdm.fun(r = meanbc, e = pe, d = pd)

pdf("output_plots/rarefy_70_barplot.pdf", width = 5, height = 11)
par(mfrow = c(5, 2))
gbarplot2(rare_70_parasite, "Parasite species turnover rarefying to 70 infections")
gbarplot2(gdm_par_unifs_common, "Parasite phylogenetic turnover rarefy = 70")
gbarplot2(gdm_par_bc_common2, "Parasite species rarefy = 70, all communities")
dev.off()


#
# Correlations among predictors  ------------------------------------------
pdf("output_plots/predictor_correlations.pdf", useDingbats = F)
par(mfrow=c(3,2), mar = c(5,5,3,3))
plot(tempPCA1 ~ precipPCA1, metadata, xlab = "Precipitation (PC1)", ylab = "Temperature (PC1)")
text(0, -4, "R = 0.35, VIF = 0.149")
with(metadata, cor(tempPCA1, precipPCA1)) #0.35
lm(tempPCA1 ~ precipPCA1, metadata) %>% summary()
1/(1-.13)

plot(x = metadata$npp.raster/10000, y = metadata$tempPCA1, xlab = "NPP (Kg C/m2)", ylab = "Temperature (PC1)")
text(1.5, -4, "R = 0.51, VIF = 1.27")
with(metadata, cor(npp.raster, tempPCA1)) #0.51
lm(tempPCA1 ~ npp.raster, metadata) %>% summary()
1/(1-0.21) #VIF

plot(x = metadata$npp.raster/10000, y = metadata$precipPCA1, xlab = "NPP (Kg C/m2)", ylab = "Precipitation (PC1)")
text(1.5, -3.5, "R = 0.62, VIF = 1.64")
lm(precipPCA1 ~ npp.raster, metadata) %>% summary()
with(metadata, cor(npp.raster, precipPCA1)) #0.62
1/(1-0.39)#VIF

plot(community.elev ~ tempPCA1, metadata, xlab = "Temperature (PC1)", ylab = "Elevation (m)")
text(-3.7, 500, "R = 0.86, VIF = 4")
lm(community.elev ~ tempPCA1, metadata) %>% summary()
1/(1-0.75) #VIF
with(metadata, cor(community.elev, tempPCA1))

plot(x = metadata$npp.raster/10000, y = metadata$community.elev, xlab = "NPP (Kg C/m2)", ylab = "Elev (m)")
text(.75, 1000, "R = -0.2, VIF = 1.05")
lm(community.elev ~ npp.raster, metadata) %>% summary()
with(metadata, cor(npp.raster, community.elev)) #-0.2
1/(1-0.05)#VIF

dev.off()


# What happens if we use either temp or elevation in models:
# Parasites
# Only elev
pe_e <- c("community.number", "community.Lat", "community.Long",
        "precipPCA1","community.elev","total.host",
        "npp.raster")
pe_t <- c("community.number", "community.Lat", "community.Long",
          "precipPCA1","tempPCA1","total.host",
          "npp.raster")
pd <- list(host_ja, host_unifs)

spp_par_elev_gdm <- p.gdm.fun(par_bc, e = pe_e, d = pd)
spp_par_temp_gdm <- p.gdm.fun(par_bc, e = pe_t, d = pd)

phy_par_elev_gdm <- p.gdm.fun(par_unifs, e = pe_e, d = pd)
phy_par_temp_gdm <- p.gdm.fun(par_unifs, e = pe_t, d = pd)

# Hosts

he_e <- c("community.number", "community.Lat", "community.Long",
          "precipPCA1","community.elev","parasite.richness",
          "npp.raster")
he_t <- c("community.number", "community.Lat", "community.Long",
          "precipPCA1","tempPCA1", "parasite.richness",
          "npp.raster")

hd <- list(par_unifs, par_bc)

spp_host_elev_gdm <- p.gdm.fun(host_ja, e = he_e, d = hd)
spp_host_temp_gdm <- p.gdm.fun(host_ja, e = he_t, d = hd)

phy_host_elev_gdm <- p.gdm.fun(host_unifs, e = he_e, d = hd)
phy_host_temp_gdm <- p.gdm.fun(host_unifs, e = he_t, d = hd)

gdm_host_ja <- readRDS("GDM_results/gdm_host_ja.rds")
gdm_host_uf <- readRDS("GDM_results/gdm_host_uf.rds")


# parasites
gdm_par_bc <- readRDS("GDM_results/gdm_par_bc.rds")
gdm_par_uf <- readRDS("GDM_results/gdm_par_uf.rds")

pdf("output_plots/dropping_temp_elev.pdf", height = 7, width = 8)
par(mfrow=c(4,4))
plot.new(); mtext("Host spp turnover"); mtext(line = -2, "Best model includes both,\nAlone, temp and elev gain importance", cex = 0.6)
gbarplot2(gdm_host_ja, "Both temp and elev")
gbarplot2(spp_host_temp_gdm, "Temp only")
gbarplot2(spp_host_elev_gdm, "Elev only")

plot.new(); mtext("Host phylo turnover"); mtext(line = -2, "Best model includes both,\nWithout temp model explains less var", cex = 0.6)
gbarplot2(gdm_host_uf, "")
gbarplot2(phy_host_temp_gdm, "")
gbarplot2(phy_host_elev_gdm, "")

plot.new(); mtext("Parasite spp turnover"); mtext(line = -2, "Neither temp nor elev very important", cex = 0.6)
gbarplot2(gdm_par_bc, "")
gbarplot2(spp_par_temp_gdm, "")
gbarplot2(spp_par_elev_gdm, "")

plot.new(); mtext("Parasite phylo turnover"); mtext(line = -2, "Effect of elev is distinct from effect of temp;\nModel without elev explains less var", cex = 0.6)
gbarplot2(gdm_par_uf, "")
gbarplot2(phy_par_temp_gdm, "")
gbarplot2(phy_par_elev_gdm, "")
dev.off()


# Analysis using kilometers -----------------------------------------------

km_dist <- distm(CommunitySpatial, fun = distHaversine) / 1000
km_dist <- as.data.frame(km_dist) %>% mutate(community.number = 1:18) %>% select(community.number, everything()) #distance matrix in km

gdm.fun.km <- function (r, e, d) {
  dplyr::select(metadata, e) %>% as.data.frame(.) %>%
    mutate(., community.number=1:18) -> e.formatted  # select predictor vars
  e.r.sitepair <- formatsitepair(
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
  gdm(e.r.sitepair, geo=F)
}

CommunitySpatial

coordinates(CommunitySpatial)
select(metadata, community.Lat, community.Long)

load(system.file("./data/gdm.RData", package="gdm"))
gdmExpData[1:3,] %>% str()
sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
str(envTab)
summary(envTab$Lat)



pointDistance(CommunitySpatial, allpairs = TRUE, lonlat= T)
euc_dist <- pointDistance(select(metadata, community.Lat, community.Long), allpairs = T, lonlat= F)
haversine_dist <- distm(CommunitySpatial, fun = distHaversine)
greatc_dist <- st_as_sf(metadata, coords = c("community.Long", "community.Lat"), crs = "+init=epsg:4326") %>% st_distance()


cor(as.vector(haversine_dist)[as.vector(haversine_dist)>0], as.vector(greatc_dist)[as.vector(greatc_dist) > 0 ])
lm (as.vector(euc_dist)[as.vector(euc_dist)>0] ~ as.vector(haversine_dist)[as.vector(haversine_dist)>0]) %>% summary()


plot (unique(as.vector(euc_dist))*111, unique(as.vector(greatc_dist))/1000,
      xlab = "Euclidean distance (km; 111 km = 1 degree)", ylab = "Haversine distance (km)")
abline(a = 0, b = 1)
dev.off()

# What's the degree to km conversion at the extremes of our latitudinal range?
min(metadata$community.Lat) #-14.05
max(metadata$community.Lat) #-5.1
median(metadata$community.Long) #-77

pointDistance(p1 = c(-77, -14), p2 = c(-77, -13), lonlat = TRUE) # 110.63 km
pointDistance(p1 = c(-77, -4), p2 = c(-77, -5), lonlat = TRUE) # 110.58 km

# Plot distance vs. difference in precipitation
pdf("output_plots/precipxdistance.pdf")
precip_differences <- dist(metadata$precipPCA1, method = "manhattan", diag = T, upper = T) %>% as.matrix() %>% c()
plot(x = c(km_dist), y = precip_differences, xlab = "Distance between communities (km)", ylab = "Difference in precipitation (PC1)")
dev.off()

points_sf <- st_as_sf(metadata, coords = c("community.Long", "community.Lat"), crs = "+init=epsg:4326")

st_distance(points_sf)/1000
plot(as.vector(distm(CommunitySpatial, fun = distHaversine) / 1000), as.vector(st_distance(points_sf)/1000))
abline(a = 0, b = 1)
dev.off()
metadata

test <- formatsitepair(
  bioData = host_ja,
  bioFormat = 3,
  XColumn = "community.Long",
  YColumn = "community.Lat",
  siteColumn = "community.number",
  #distPreds = d,
  predData = metadata[,c("community.number", "community.Lat", "community.Long",
                      "tempPCA1", "precipPCA1","community.elev", "npp.raster")],
  weightType = "equal",
  distPreds = hd
) #predictor data = scaled metadata
# Distance between point 1 and 2 = 0.979

sqrt((-69.2128 - -76.10601)^2 + (-11.711 - -9.727500)^2)

summary(test$distance)
head(test)
CommunitySpatial@coords %>% str()



# Hosts
he <- c("community.number", "community.Lat", "community.Long",
        "tempPCA1", "precipPCA1","community.elev", "npp.raster")
hd <- list(km_dist, par_unifs)

# Parasites
pe <- c("community.number", "community.Lat", "community.Long",
"precipPCA1","community.elev","total.host",
"npp.raster")
pd <- list(host_ja, km_dist)

gdm_km_host_spp <-  gdm.fun.km(host_ja, e = he, d = hd)
gdm_km_host_phy <-  gdm.fun.km(host_unifs, e = he, d = hd)
gdm.fun(host_ja, e = he, d = hd) %>% plot(plot.layout = c(3,4))
gdm.fun.km(host_ja, e = he, d = hd) %>% plot(plot.layout = c(3,4))

gdm_km_par_spp <-  gdm.fun.km(par_bc, e = c("community.number", "community.Lat", "community.Long",
                                            "precipPCA1"), d = pd)
gdm_km_par_phy <-  gdm.fun.km(par_unifs, e = c("community.number", "community.Lat", "community.Long",
                                                       "precipPCA1","community.elev","total.host"), d = pd)


plot(gdm_km_host_phy, plot.layout = c(4,2))
plot(gdm_km_host_spp, plot.layout = c(4,2))
plot(gdm_km_par_phy, plot.layout = c(4,3))
plot(gdm_km_par_spp, plot.layout = c(4,3))


# Best models from original analysis
gdm_b_host_phy <- readRDS("GDM_results/gdm_b_host_phy.rds")
gdm_b_host_spp <- readRDS("GDM_results/gdm_b_host_spp.rds")
gdm_b_par_phy <- readRDS("GDM_results/gdm_b_par_phy.rds")
gdm_b_par_spp <- readRDS("GDM_results/gdm_b_par_spp.rds")


plot(gdm_b_host_phy, plot.layout = c(4,2))
plot(gdm_b_host_spp , plot.layout = c(4,2))
plot(gdm_b_par_phy , plot.layout = c(4,3))
plot(gdm_b_par_spp, plot.layout = c(4,3))

# Make some comparison plots: host spp, host phylo, par spp, par phylo

isplineExtract(gdm_km_host_spp) %>% as.data.frame %>% select(x.matrix_2, y.matrix_2) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (km)", ylab = "Ecol. distance", main = "Host species turnover")

pdf("output_plots/distance_inkm.pdf")
par(mfrow = c(4,2))

isplineExtract(gdm_km_host_spp) %>% as.data.frame %>% select(x.matrix_2, y.matrix_2) %>%
plot(type = "l", lwd = 3, xlab = "Geographic distance (km)", ylab = "Ecol. distance", main = "Host species turnover")
isplineExtract(gdm_b_host_spp) %>% as.data.frame %>% select(x.Geographic, y.Geographic) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (degrees)", ylab = "Ecol. distance", main = "Host species turnover")

isplineExtract(gdm_km_host_phy) %>% as.data.frame %>% select(x.matrix_2, y.matrix_2) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (km)", ylab = "Ecol. distance", main = "Host phylo turnover")
isplineExtract(gdm_b_host_phy) %>% as.data.frame %>% select(x.Geographic, y.Geographic) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (degrees)", ylab = "Ecol. distance", main = "Host phylo turnover")


isplineExtract(gdm_km_par_spp) %>% as.data.frame %>% select(x.matrix_2, y.matrix_2) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (km)", ylab = "Ecol. distance", main = "Parasite species turnover")
isplineExtract(gdm_b_par_spp) %>% as.data.frame %>% select(x.Geographic, y.Geographic) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (degrees)", ylab = "Ecol. distance", main = "Parasite species turnover")

isplineExtract(gdm_km_par_phy) %>% as.data.frame %>% select(x.matrix_2, y.matrix_2) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (km)", ylab = "Ecol. distance", main = "Parasite phylo turnover")
isplineExtract(gdm_b_par_phy) %>% as.data.frame %>% select(x.Geographic, y.Geographic) %>%
  plot(type = "l", lwd = 3, xlab = "Geographic distance (degrees)", ylab = "Ecol. distance", main = "Parasite phylo turnover")
dev.off()


# Host characteristics ----------------------------------------------------

sampledhosts %>% head()
# Create unique lat long id, rounding both so we can match to community center
sampledhosts$latlong <- as.factor(paste(round(sampledhosts$neg.degrees.Lat,3),
                                        round(sampledhosts$neg.degrees.Long,3), sep="."))
community.lat.long.list <- unique(dplyr::select(inputhaplos, "latlong", "community.number"))
sampledhosts <- merge(sampledhosts, community.lat.long.list, all.x=F, by="latlong") #add community numbers to host sampling
sampledhosts$community.number <- as.factor(sampledhosts$community.number)
head(sampledhosts)

pdf("output_plots/screened_bird_mass.pdf", width = 9, height = 5)
boxplot(log10(body.mass..g.) ~ community.number, sampled_hosts_one, ylab = "Log body mass (g)", xlab = "Community")
dev.off()
dunnTest(body.mass..g. ~community.number, sampledhosts, method = "bonferroni") %>% .$res %>% filter(P.adj < 0.05)

filter(sampledhosts, community.number == 2)

metadata$host_mass <- aggregate (log(body.mass..g.) ~ community.number, sampledhosts[sampledhosts$family!="Trochilidae",], median)[,2]
metadata$host_mass <- aggregate ((body.mass..g.) ~ community.number, sampledhosts, median)[,2]

pe <- c("community.number", "community.Lat", "community.Long",
        "precipPCA1","community.elev","total.host",
        "npp.raster", "host_mass", "shannonD")
pd <- list(host_ja, host_unifs)
gdm.fun(par_bc, e = pe, d = pd) %>% plot(plot.layout  = c(3,4))
gdm_par_spp_mass <-  gdm.fun(par_bc, e = c("community.number", "community.Lat", "community.Long", "precipPCA1", "host.mass"), d = pd)

spp_mass_gdm <- p.gdm.fun(par_bc, e = pe, d = pd)
phy_mass_gdm <- p.gdm.fun(par_unifs, e = pe, d = pd)

gbarplot2(spp_mass_gdm, "parasite spp")
gbarplot2(phy_mass_gdm, "parasite phylo")


head(sampledhosts)
# Barplot of full model

gbarplot_full(phy_mass_gdm, "Parasite phylo turnover, full model")
gbarplot_full(spp_mass_gdm, "Parasite spp turnover, full model")
dev.off()
modsel(spp_mass_gdm)

gbarplot_full <- function(model, title) {
  par(mar = c(5, 5, 4, 2))
  vals <-model[[2]][,1]
  factor_cols <-
    colmatch$cols[match(names(sort(vals, decreasing = T)),
                        colmatch$factors)] %>% as.character # Assign col
  # Need to rename factor names to make them consistent and clean
  namestring <- str_replace_all(
    names(sort(vals, decreasing = T)),
    c(
      "Geographic" = "Dist",
      "tempPCA1" = "Temp",
      "precipPCA1" = "Precip",
      "community.elev" = "Elev",
      "total.host" = "Host rich",
      "matrix_1" = "Host turn",
      "shannonD" = "Samp Effort",
      "parasite.richness" = "Symb rich",
      "npp.raster" = "NPP",
      "host_mass_log" = "Host mass"
    )
  )
  # Create a barplot of variable importance from most to least important.
  mp <-
    barplot(
      as.numeric(scale(sort(vals, decreasing = T), center = F)), #scale variables
      #as.numeric(sort(vals, decreasing = T)),
      col = factor_cols,
      main = title,
      axis = F,
      axisnames = F,
      cex.axis = 1.1,
      ylab = "Scaled importance",
      cex.lab = 1.3
    )
  mtext(
    paste("Variance explained = ", round(model[[1]][2], 2), "%", sep = ""),
    side = 3,
    line = 0.4,
    cex = 0.7
  )
  text(
    mp,
    par("usr")[3],
    labels = namestring,
    srt = 45,
    #par("usr")[3] = ymin of plot
    adj = c(1.1, 1.1),
    xpd = T,
    cex = 1.3
  )
}



#
#









