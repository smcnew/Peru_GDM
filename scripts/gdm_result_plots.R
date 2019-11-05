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

extent(peru_alt)==extent(precipRe) #should match
extent(peru_alt)==extent(tempRe) #should match

#
# Barplots ----------------------------------------------------------------
# Plots of variable importance based on best model of permuted models

# Fucntion to create a barplot from the best model based on backwards selection
# from full model. Arguments: model = full permuted model, title = title

# Choose some good categorical colors
varcols = c("#332288", "#BBBBBB", "#66CCEE",  "#AA3377",
            "#CCBB44", "#228833", "#4477AA", "#EE6677")

# Assign colors to predictor variables. Last color gets repeated because it should
# represent both host richness and parasite richness (i.e. effect of "symbiont richness"
# on turnover)
colmatch <- data.frame( cols = c(varcols, varcols[8]), factors = c(c("precipPCA1",
                                                                     "Geographic",
                                                                     "community.elev",
                                                                     "tempPCA1",
                                                                     "matrix_1",
                                                                     "npp.raster",
                                                                     "shannonD",
                                                                     "total.host",
                                                                     "parasite.richness"
)))



gbarplot2<- function(model, title) {
  par(mar=c(5,5,4,2))
  # Choose best model based on most variance explained using fewest predictors
  modsel <- function(model){
    best <- (model[[1]][2,] - sort(1:ncol(model[[1]]), decreasing=T)) %>% which.max #best model
    na.omit(model[[2]][,best])
  }

  vals <- modsel(model) # Extract variable importance from best model
  factor_cols <- colmatch$cols[match(names(sort(vals, decreasing=T)),
                                     colmatch$factors)] %>% as.character # Assign col
  # Need to rename factor names to make them consistent and clean
  namestring <- str_replace_all(
    names(sort(vals, decreasing=T)),
    c("Geographic" = "Dist",
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
  mp <- barplot(as.numeric(scale(sort(vals, decreasing = T), center=F)), col=factor_cols,
                main=title, axis=F, axisnames=F, cex.axis=1.1,
                ylab="Scaled importance", cex.lab=1.3)
  mtext(paste("Variance explained = ", round(model[[1]][2], 2), "%", sep=""),
        side = 3, line = 0.4, cex=0.7)
  text(mp, par("usr")[3], labels=namestring, srt=45, #par("usr")[3] = ymin of plot
       adj=c(1.1,1.1), xpd=T, cex=1.3)
}

# Probably could clean code by putting models in list

pdf("./output_plots/colored_bar_scaled.pdf", width=5, height=11)
par(mfrow=c(5,2))
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


# Rename columns to be consistent and clean
rename_fun <- function (frame) {
  str_replace_all(
    colnames(frame),
    c("Geographic" = "Dist",
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
}

# Rename columns and then turn into data frames
colnames(haplo.splines1$x) <- rename_fun(haplo.splines1$x)
colnames(haplo.splines1$y) <- rename_fun(haplo.splines1$y)
colnames(haplo.splines2$x) <- rename_fun(haplo.splines2$x)
colnames(haplo.splines2$y) <- rename_fun(haplo.splines2$y)
colnames(haplo.splines3$x) <- rename_fun(haplo.splines3$x)
colnames(haplo.splines3$y) <- rename_fun(haplo.splines3$y)
colnames(haplo.splines4$x) <- rename_fun(haplo.splines4$x)
colnames(haplo.splines4$y) <- rename_fun(haplo.splines4$y)
haplo.splines1$x <-as.data.frame(haplo.splines1$x)
haplo.splines1$y <-as.data.frame(haplo.splines1$y)
haplo.splines2$x <-as.data.frame(haplo.splines2$x)
haplo.splines2$y <-as.data.frame(haplo.splines2$y)
haplo.splines3$x <-as.data.frame(haplo.splines3$x)
haplo.splines3$y <-as.data.frame(haplo.splines3$y)
haplo.splines4$x <-as.data.frame(haplo.splines4$x)
haplo.splines4$y <-as.data.frame(haplo.splines4$y)


variables <- c("Dist", "Temp", "Precip", "Elev", "Symb_rich", "Symb_turn", "NPP")

# This is a nightmare. Basically, we need to go through each data frame and add
# dummy columns of 0s for all the variables that weren't in the best model.
# That way we can plot all variables in the same order, with the same color.
# If they have values in the best model, great. If they weren't, they'll just be
# plotted as 0s and will be invisible.
# ETA this is actually pretty clever, I just need to comment better.

fncols <- function(data, cname) {
  add <-cname[!cname %in% names(data$x)] # variables not in best model
  if(length(add)!= 0) data$x[add] <- 0 # add 0s as vals for dummy vars
  data$x <- select(data$x, variables) # Reorder everything so it's consistent
  add <-cname[!cname %in% names(data$y)] # do same for y vals.
  if(length(add)!= 0) data$y[add] <- 0
  data$y <- select(data$y, variables)
  data
}

haplo.splines1 <- fncols(haplo.splines1, variables)
haplo.splines2 <- fncols(haplo.splines2, variables)
haplo.splines3 <- fncols(haplo.splines3, variables)
haplo.splines4 <- fncols(haplo.splines4, variables)

colnames(haplo.splines3$x) #colnames should all be the same now for all dfs
head(haplo.splines4$x) # factors not in best model will have 0s


# oo boy.


pdf("./output_plots/spline_bestmodels.pdf", height=11, width=7)
namestring <- colnames(haplo.splines1$x)
par(mfrow=c(4,2), mar=c(4.5,4.5,2,2))
for (i in 1:ncol(haplo.splines1$x)) {
  ymaxs <- lapply(list(haplo.splines1,haplo.splines2,haplo.splines3,haplo.splines4),
                  function (x) apply(x[[2]], 2, max)) %>% as.data.frame() %>% t()
  xmaxs <- lapply(list(haplo.splines1,haplo.splines2,haplo.splines3,haplo.splines4),
                  function (x) apply(x[[1]], 2, max)) %>% as.data.frame() %>% t()
  xmins <- lapply(list(haplo.splines1,haplo.splines2,haplo.splines3,haplo.splines4),
                  function (x) apply(x[[1]], 2, min)) %>% as.data.frame() %>% t()
  plot(haplo.splines1$x[,i], haplo.splines1$y[,i], type="l", lwd=3, col=varcols[7],
       xlab=paste("Dissimilarity in", namestring[i]), ylab=paste("f(",namestring[i],")",sep=""),
       ylim=c(0, 1.1*max(ymaxs[,i])), xlim=c(min(xmins[,i]), max(xmaxs[,i])),
       cex.axis=1.2, cex.lab=1.4)
  lines(haplo.splines2$x[,i], haplo.splines2$y[,i], lwd=3, col=varcols[7], lty=2)
  lines(haplo.splines3$x[,i], haplo.splines3$y[,i], lwd=3, col=varcols[4])
  lines(haplo.splines4$x[,i], haplo.splines4$y[,i], lwd=3, col=varcols[4], lty=2)
  #for (j in 1:4){
  # modlist <- list(haplo.splines1,haplo.splines2,haplo.splines3, haplo.splines4)
  #text(x = 1.1*max(modlist[[j]]$x[,i]) , y = max(modlist[[j]]$y[,i]),
  #    labels = names[j], cex=1)
}
plot(NA, NA, xlim=c(1,1), ylim=c(1,1), axes=F, ylab="", xlab="")
legend("center", legend = c("Parasite species", "Parasite phylogenetic",
                            "Host species","Host phylogenetic"),
       col = varcols[c(7,7,4,4)], lwd=3, cex=1.5,
       title="Turnover model", lty=c(1,2,1,2))
dev.off()


#
# Maps --------------------------------------------------------------------


