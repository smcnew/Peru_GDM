# Script to create result figures based on GDM models.
# Loads models that were saved as .rds from gdm_host_parasite.R script

library(raster) # for map data
library(tools) # for quick reading in of models
library(stringr) # for replacing string names
library(gdm) #to extract splines
library(scales) #for alpha of colors
library(dplyr)

set.seed(1987)

#
# Load models  ------------------------------------------------------------

# Permuted models (to generate barplots)

# hosts
gdm_host_ja <- readRDS("GDM_results/gdm_host_ja.rds")
gdm_host_uf <- readRDS("GDM_results/gdm_host_uf.rds")


# parasites
gdm_par_bc <- readRDS("GDM_results/gdm_par_bc.rds")
gdm_par_uf <- readRDS("GDM_results/gdm_par_uf.rds")

#parasite x genus
gdm_h_bc <- readRDS("GDM_results/gdm_h_bc.rds")
gdm_h_uf <- readRDS("GDM_results/gdm_h_uf.rds")
gdm_l_bc <- readRDS("GDM_results/gdm_l_bc.rds")
gdm_l_uf <- readRDS("GDM_results/gdm_l_uf.rds")
gdm_p_bc <- readRDS("GDM_results/gdm_p_bc.rds")
gdm_p_uf <- readRDS("GDM_results/gdm_p_uf.rds")

# Best models of hosts and parasites (for splines)
gdm_b_host_phy <- readRDS("GDM_results/gdm_b_host_phy.rds")
gdm_b_host_spp <- readRDS("GDM_results/gdm_b_host_spp.rds")
gdm_b_par_phy <- readRDS("GDM_results/gdm_b_par_phy.rds")
gdm_b_par_spp <- readRDS("GDM_results/gdm_b_par_spp.rds")

# Shortcut code to load everything in GDM_results file (if one's familiar with
# what's in there, and ignores advice to avoid assigning objects through loops)
models <- list.files("GDM_results/")
for (i in 1:length(models)) {
  assign(file_path_sans_ext(models)[i], readRDS(paste("GDM_results/",
  models[i], sep = "")))
}

# Load spatial data -------------------------------------------------------

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

# Community metadata
metadata <- read.csv("./formatted_data/GDM_metadata.csv")

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

# Load distance matrices (for maps) -------------------------------------------------------
#Host distance matrices
host_ja <- read.csv("./formatted_data/hostjaccard.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Jaccard

host_unifs <- read.csv("./formatted_data/hostUnifrac.csv") %>%
  rename(community.number = X) %>% mutate(community.number = 1:18) #Unifracs (unweighted)

#Parasite distance matrices
par_bc <- read.csv("./formatted_data/parasitebraycurtis.csv") %>% rename(., community.number=X) #Bray curtis
par_unifs <- read.csv("./formatted_data/malaria.generalized.unifracs.csv") %>%
  rename(., community.number=X) #Unifracs (generalized)

#
# Barplots ----------------------------------------------------------------
# Plots of variable importance based on best model of permuted models

# Fucntion to create a barplot from the best model based on backwards selection
# from full model. Arguments: model = full permuted model, title = title

# Choose some good categorical colors
varcols = c(
  "#332288",
  "#BBBBBB",
  "#66CCEE",
  "#AA3377",
  "#CCBB44",
  "#228833",
  "#4477AA",
  "#EE6677"
)

# Assign colors to predictor variables. Last color gets repeated because it should
# represent both host richness and parasite richness (i.e. effect of "symbiont richness"
# on turnover)
colmatch <-
  data.frame(cols = c(varcols, varcols[8]), factors = c(
    c(
      "precipPCA1",
      "Geographic",
      "community.elev",
      "tempPCA1",
      "matrix_1",
      "npp.raster",
      "shannonD",
      "total.host",
      "parasite.richness"
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
      as.numeric(scale(sort(vals, decreasing = T), center = F)),
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

# Probably could clean code by putting models in list

pdf("./output_plots/colored_bar_scaled.pdf",
    width = 5,
    height = 11)
par(mfrow = c(5, 2))
gbarplot2(gdm_host_ja, "Host species")
gbarplot2(gdm_host_uf, "Host phylogenetic")
gbarplot2(gdm_par_bc, "All parasite species")
gbarplot2(gdm_par_uf, "All parasite phylogenetic")
gbarplot2(gdm_h_bc, "Haemoproteus species")
gbarplot2(gdm_h_uf, "Haemoproteus phylogenetic")
gbarplot2(gdm_l_bc, "Leucocytozoon species")
gbarplot2(gdm_l_uf, "Leucocytozoon phylogenetic")
gbarplot2(gdm_p_bc, "Plasmodium species")
gbarplot2(gdm_p_uf, "Plasmodium phylogenetic")
dev.off()


# Spline plot  ------------------------------------------------------------
# Code to create splines from best models
# Challenge: Best models for different measures of turnover include different
# combos of variables. So we need to plot splines sometimes, but not always.

haplo.splines1 <- isplineExtract(gdm_b_par_spp) # extract x and y val splines
haplo.splines2 <- isplineExtract(gdm_b_par_phy)
haplo.splines3 <- isplineExtract(gdm_b_host_spp)
haplo.splines4 <- isplineExtract(gdm_b_host_phy)

modlist <-
  list(haplo.splines1,
       haplo.splines2,
       haplo.splines3,
       haplo.splines4)

# Function that rename columns to match display names (i.e. "tempPCA" -> "Temp" )
# and turns x and y into dataframes. Takes and returns a gdm model.
rename_fun <- function (mod) {
  colnames(mod$x) <- str_replace_all(
    colnames(mod$x),
    c(
      "Geographic" = "Dist",
      "tempPCA1" = "Temp",
      "precipPCA1" = "Precip",
      "community.elev" = "Elev",
      "total.host" = "Symb_rich",
      "matrix_1" =
        "Symb_turn",
      "shannonD" =
        "Samp_Effort",
      "parasite.richness" =
        "Symb_rich",
      "npp.raster" =
        "NPP"
    )
  )
  mod$x <- as.data.frame(mod$x)
  colnames(mod$y) <- str_replace_all(
    colnames(mod$y),
    c(
      "Geographic" = "Dist",
      "tempPCA1" = "Temp",
      "precipPCA1" = "Precip",
      "community.elev" = "Elev",
      "total.host" = "Symb_rich",
      "matrix_1" =
        "Symb_turn",
      "shannonD" =
        "Samp_Effort",
      "parasite.richness" =
        "Symb_rich",
      "npp.raster" =
        "NPP"
    )
  )
  mod$y <- as.data.frame(mod$y)
  mod
}

# Apply function to model list
modlist <- lapply(modlist, rename_fun)

# This is a nightmare. Basically, we need to go through each data frame and add
# dummy columns of 0s for all the variables that weren't in the best model.
# That way we can plot all variables in the same order, with the same color.
# If they have values in the best model, great. If they weren't, they'll just be
# plotted as 0s and will be invisible.
# ETA this is actually pretty clever, I just need to comment better.

fncols <- function(data, cname) {
  add <-
    cname[!cname %in% names(data$x)] # identify variables not in best model
  if (length(add) != 0)
    data$x[add] <- 0 # add 0s as vals for dummy vars
  data$x <-
    select(data$x, variables) # Reorder columns so they're the same for each model
  add <- cname[!cname %in% names(data$y)] # do same for y vals.
  if (length(add) != 0)
    data$y[add] <- 0
  data$y <- select(data$y, variables)
  data
}
variables <-
  c("Dist",
    "Temp",
    "Precip",
    "Elev",
    "Symb_rich",
    "Symb_turn",
    "NPP")
# Apply data-modifying function to our model list.

modlist <-lapply (modlist, function (x) fncols(data = x, cname = variables))

lapply(modlist, function (i) colnames(i$x)) # should all be the same now
lapply(modlist, function (i) (i$x)[1, ]) # factors not in best model should have 0s

# Plot it all up
pdf("./output_plots/spline_bestmodels.pdf",
    height = 11,
    width = 7)
namestring <- colnames(modlist[[1]]$x)
par(mfrow = c(4, 2), mar = c(4.5, 4.5, 2, 2))
for (i in 1:ncol(modlist[[1]]$x)) {
  ymaxs <-
    sapply(modlist, function (i)
      apply(i$y, 2, max)) %>% t() # y max for each predictor
  xmaxs <-
    sapply(modlist, function (i)
      apply(i$x, 2, max)) %>% t() # x max for each predictor
  xmins <-
    sapply(modlist, function (i)
      apply(i$x, 2, min)) %>% t() # x min for each predictor
  plot(
    NA,
    NA, #Create blank plots for each factor
    xlab = paste("Dissimilarity in", namestring[i]),
    ylab = paste("f(", namestring[i], ")", sep = ""),
    ylim = c(0, 1.1 * max(ymaxs[, i])),
    xlim = c(min(xmins[, i]), max(xmaxs[, i])),
    cex.axis = 1.2,
    cex.lab = 1.4
  )
  # now plot lines from each model
  # 10,000 points for using mapply
    mapply (
    function(j, k, l)
      lines(
        j$x[, i],
        j$y[, i],
        col = k,
        lty = l,
        lwd = 3
      ),
    j = modlist,  # Plot one line per model
    k = c(varcols[c(7, 7, 4, 4)]),  # color = host or parasite
    l = c(1, 2, 1, 2) # line type = species or phylo turnover
  )
}

# Legend
plot(
  NA,
  NA,
  xlim = c(1, 1),
  ylim = c(1, 1),
  axes = F,
  ylab = "",
  xlab = ""
)
legend(
  "center",
  legend = c(
    "Parasite species",
    "Parasite phylogenetic",
    "Host species",
    "Host phylogenetic"
  ),
  col = varcols[c(7, 7, 4, 4)],
  lwd = 3,
  cex = 1.5,
  title = "Turnover model",
  lty = c(1, 2, 1, 2)
)
dev.off()


#

# Host vs. parasite plots  ------------------------------------------------
# Plot total turnover as a function of "ecological distance" and compare
# hosts vs. parasites and spp. vs. phylo turnover

#plots of best models (post model selection)
pdf("./output_plots/parasites_vs_hosts_bestmodel.pdf", width=10, height=5, useDingbats = F)
par(mfrow=c(1,2), mar=c(5,5,4,4))
alph=0.7
plot(gdm_b_par_spp$ecological, gdm_b_par_spp$observed, pch=24, cex=1.2, bg=alpha(varcols[7], alph),
     xlab="Predicted ecological distance", cex.axis=1.3, cex.lab=1.5,
     ylab="Observed compositional dissimilarity", ylim=c(0,1),
     main="Species dissimilarity", cex.main=1.8, xlim=c(0.5, 4.5))
points(gdm_b_host_spp$ecological, gdm_b_host_spp$observed, pch=21, cex=1.2, bg=alpha(varcols[4],alph))
lines(smooth.spline(gdm_b_par_spp$ecological, gdm_b_par_spp$observed, spar=1.2), lwd=3, col=varcols[7])
lines(smooth.spline(gdm_b_host_spp$ecological, gdm_b_host_spp$observed, spar=1.1), lwd=3, col=varcols[4], lty=2)
legend("bottomright", pch=c(21,24), pt.bg=varcols[c(4,7)], cex=1,
       legend=c("Hosts", "Parasites"))


plot(gdm_b_par_phy$ecological, gdm_b_par_phy$observed, pch=21, cex=1.2, bg=alpha(varcols[4], alph),
     xlab="Predicted ecological distance", cex.axis=1.3, cex.lab=1.5,
     ylab="Observed compositional dissimilarity", ylim=c(0,1),
     main="Phylogenetic dissimilarity", cex.main=1.8, xlim=c(0.2,1.8))
points(gdm_b_host_phy$ecological, gdm_b_host_phy$observed, pch=24, cex=1.2, bg=alpha(varcols[7], alph))
lines(smooth.spline(gdm_b_host_phy$ecological, gdm_b_host_phy$observed), col=varcols[4], lwd=3, lty=2)
lines(smooth.spline(gdm_b_par_phy$ecological, gdm_b_par_phy$observed, spar=1.2), lwd=3, col=varcols[7])
legend("bottomright", pch=c(21,24), pt.bg=varcols[c(4,7)], cex=1,
       legend=c("Hosts", "Parasites"))

dev.off()
#
# Maps --------------------------------------------------------------------

# Create visualization function, takes a GDM model and a stack of raster layers
# and projects the model data back onto the layers to create a map. Similar colors
# on the map indicate similar communities. Instead of using saved models
# we'll make them again because it's *very important* that the predictors
# and the raster layers are in the same order.


mapfun <- function(model, rasterdata) {
  rastTrans <- gdm.transform(model, rasterdata)
  plot(rastTrans)
  rastDat <- na.omit(getValues(rastTrans))
  pcaSamp <- prcomp(rastDat)
  pcaRast <- predict(rastTrans, pcaSamp, index = 1:3)
  pcaRast[[1]] <- (pcaRast[[1]] - pcaRast[[1]]@data@min) /
    (pcaRast[[1]]@data@max - pcaRast[[1]]@data@min) * 255
  pcaRast[[2]] <- (pcaRast[[2]] - pcaRast[[2]]@data@min) /
    (pcaRast[[2]]@data@max - pcaRast[[2]]@data@min) * 255
  pcaRast[[3]] <- (pcaRast[[3]] - pcaRast[[3]]@data@min) /
    (pcaRast[[3]]@data@max - pcaRast[[3]]@data@min) * 255
  par(mfrow = c(1, 1))
  plotRGB(pcaRast, r = 1, g = 2, b = 3)
}

#HOSTS PHYLOGENETIC TURNOVER: GDMs suggest best predictors of
#host turnover are precip, temp, elev and distance and npp.

rasterdata_h <- brick(addLayer(precip_re, temp_re, peru_alt, npp_re))

enviro_meta_h <-
  dplyr::select(
    metadata,
    community.number,
    community.Lat,
    community.Long,
    precipPCA1,
    tempPCA1,
    community.elev,
    npp.raster
  ) %>%
  as.data.frame(.) %>%
  mutate(., community.number = 1:18) %>% as.matrix(.)

enviro_table_host_phy <-
  formatsitepair(
    bioData = host_unifs,
    bioFormat = 3, #bioFormat = 3 means distance matrix
    XColumn = "community.Long",
    YColumn = "community.Lat",
    siteColumn = "community.number",
    predData = enviro_meta_h,
    weightType = "equal"
  )

map_gdm_host_phy <- gdm(enviro_table_host_phy, geo=T)


#results:
pdf("./output_plots/host_phylo_turnover_map.pdf")
mapfun(map_gdm_host_phy, rasterdata_h)
dev.off()

#HOSTS SPECIES TURNOVER: best predictors of
#host turnover are precip, temp, elev and distance and npp.
enviro_table_host_spp <-
  formatsitepair(
    bioData = host_ja,
    bioFormat = 3,
    XColumn = "community.Long",
    YColumn = "community.Lat",
    siteColumn = "community.number",
    predData = enviro_meta_h,
    weightType = "equal"
  )

map_gdm_host_spp <- gdm(enviro_table_host_spp, geo = T)

pdf("./output_plots/host_spp_turnover_map.pdf")
mapfun(map_gdm_host_spp, rasterdata_h)
dev.off()

# Host turnover was a significant predictor of parasite species turnover
# To create a map of parasite species turnover we will need a raster of
# host turnover. Approach: transform environmental spatial predictors using the
# host turnover model. Extract values from those rasters and PCA them to produce
# principle components


rast_trans_host <- gdm.transform(map_gdm_host_spp, rasterdata_h)
rast_dat_h <- na.omit(getValues(rast_trans_host)) #pulls values for each square of each raster
pca_samp_h <- prcomp(rast_dat_h) # reduce dimensions to distill 6 predictors to PCs
pca_rast_h <- predict(rast_trans_host, pca_samp_h, index = 1:4) # predict vals based on PC xs
metadata <- cbind(metadata, extract(pca_rast_h, CommunitySpatial)) # extract values to use as predictors
pca_samp_h %>% summary()

#PARASITES PHYLO: based on GDMs of parasite unifs. Most important variables include
#elevation (by a lot), precipitation, host species richness, and geographic distance.

#raster stack
rasterdata_p <- brick(addLayer(precip_re, peru_alt, birdrast))

#select metadata
enviro_meta_p_p <-
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

#create site pair table
enviro_table_p <-
  formatsitepair(
    bioData = par_unifs,
    bioFormat = 3,
    XColumn = "community.Long",
    YColumn = "community.Lat",
    siteColumn = "community.number",
    predData = enviro_meta_p_p,
    weightType = "equal"
  )

map_gdm_par_phy <- gdm(enviro_table_p, geo=T)

pdf("./output_plots/parasite_phylo_turnover_map.pdf")
mapfun(map_gdm_par_phy, rasterdata_p) #plot results
dev.off()

#PARASITE SPP: Precip, elevation, host spp turnover (PC1 and PC2), not distance

rasterdata_p_s <- brick(addLayer(precip_re, peru_alt, pca_rast_h[[1:2]]))

#select metadata
enviro_meta_p_s <-
  dplyr::select(
    metadata,
    community.number,
    community.Lat,
    community.Long,
    precipPCA1,
    community.elev,
    layer.1,
    layer.2
    #layer.3,
    #layer.4
  ) %>%
  as.data.frame(.) %>%
  mutate(., community.number = 1:18) %>% as.matrix(.)

#create site pair table
enviro_table_p_s <- formatsitepair(
  bioData = par_bc,
  bioFormat = 3,
  XColumn = "community.Long",
  YColumn = "community.Lat",
  siteColumn = "community.number",
  predData = enviro_meta_p_s,
  weightType = "equal"
)

map_gdm_par_spp <- gdm(enviro_table_p_s, geo=F) #distance not an important predictor

pdf("./output_plots/parasite_species_turnover_map.pdf")
mapfun(map_gdm_par_spp, rasterdata_p_s) #plot results
dev.off()
