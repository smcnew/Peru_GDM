# Script to create result figures based on GDM models.
# Loads models that were saved as .rds from gdm_host_parasite.R script

library(raster) # for map data
library(tools) # for quick reading in of models
library(stringr) # for replacing string names
library(gdm) #to extract splines
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
#models <- list.files("GDM_results/")
#for (i in 1:length(models)) {
#  assign(file_path_sans_ext(models)[i], readRDS(paste("GDM_results/",
#  models[i], sep = "")))
#}

# Load spatial data -------------------------------------------------------

#birdrast <- raster("./formatted_data/birdrichness.grd") #bird species richness raster
precipPC1 <- raster("./formatted_data/precipPC1.grd") #precip PCA raster
tempPC1 <- raster("./formatted_data/tempPC1.grd") #temp PCA raster
peru_alt <- raster("./raw_data/peru_alt.grd") #read elevation raster
npp <- raster("./formatted_data/npp.grd")

precipRe <- resample(precipPC1, peru_alt, "ngb") #resample to elev grid
tempRe <- resample(tempPC1, peru_alt, "ngb") #resample to elev grid
birdrast <- resample(birdrast, peru_alt, "ngb") #resample to elev grid
nppRe <- resample(npp, peru_alt, "ngb")

extent(peru_alt) == extent(precipRe) #should match
extent(peru_alt) == extent(tempRe) #should match

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
# Stupid long and not very elegant code to create splines from best models
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
# Maps --------------------------------------------------------------------
