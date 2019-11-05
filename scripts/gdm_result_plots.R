

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


# Barplots ----------------------------------------------------------------
# Plots of variable importance based on best model of permuted models
# Model names: gdm_h_ja, gdm_h_
##Model selected barplot
gbarplot2<- function(model, title) {
  par(mar=c(5,5,4,2))
  modsel <- function(model){
    best <- (model[[1]][2,] - sort(1:ncol(model[[1]]), decreasing=T)) %>% which.max #best model
    na.omit(model[[2]][,best])
  }
  vals <- modsel(model)
  piecols <- colmatch$cols[match(names(sort(vals, decreasing=T)), colmatch$factors)] %>% as.character
  namestring <- str_replace_all(
    names(sort(vals, decreasing=T)),
    c("Geographic" = "Dist",
      "tempPCA1" = "Temp",
      "precipPCA1" = "Precip",
      "community.elev" = "Elev",
      "total.host" = "Symb rich",
      "matrix_1" =
        "Symb turn",
      "shannonD" =
        "Samp Effort",
      "parasite.richness" =
        "Symb rich",
      "npp.raster" =
        "NPP"
    )
  )
  mp <- barplot(as.numeric(scale(sort(vals, decreasing = T), center=F)), col=piecols,
                main=title, axis=F, axisnames=F, cex.axis=1.4,
                ylab="Scaled importance", cex.lab=1.5)
  mtext(paste("Variance explained = ", round(model[[1]][2], 2), "%", sep=""),
        side = 3, line = 0.4, cex=0.7)
  text(mp, par("usr")[3], labels=namestring, srt=45, #par("usr")[3] = ymin of plot
       adj=c(1.1,1.1), xpd=T, cex=1.3)
}

gbarplot2(gdmHUW, "Host species")

pdf("./output_plots/colored_bar_scaled.pdf", width=5, height=11)
par(mfrow=c(5,2))
gbarplot2(gdm_h_ja, "Host species")
gbarplot2(gdmHUW, "Host phylogenetic")
gbarplot2(gdmBC, "All parasite species")
gbarplot2(gdmgU, "All parasite phylogenetic")
gbarplot2(hb, "Haemoproteus species")
gbarplot2(hg, "Haemoproteus phylogenetic")
gbarplot2(lb, "Leucocytozoon species")
gbarplot2(lg, "Leucocytozoon phylogenetic")
gbarplot2(pb, "Plasmodium species")
gbarplot2(pg, "Plasmodium phylogenetic")
dev.off()

par(mfrow=c(6,2))

gbarplot2(hj, "h")
gbarplot2(hb, "hb")
gbarplot2(hu, "hu")
gbarplot2(hg, "hg")
gbarplot2(lj, "lj")
gbarplot2(lb, "lb")
gbarplot2(lu, "lu")
gbarplot2(lg, "lg")
gbarplot2(pj, "pj")
gbarplot2(pb, "pb")
gbarplot2(pu, "pu")
gbarplot2(pg, "pg")



# Spline plot  ------------------------------------------------------------
# Stupid long and not very elegant code to create splines from best models

haplo.splines1 <- isplineExtract(modbp1) #extract x and y val splines
haplo.splines2 <- isplineExtract(modbp2)
haplo.splines3 <- isplineExtract(modbh1)
haplo.splines4 <- isplineExtract(modbh2)

head(haplo.splines1$x)

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


variables <- c("Dist", "Temp", "Precip", "Elev", "Symb_rich", "Symb_turn",
               "NPP")

fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data$x)]
  if(length(add)!=0) data$x[add] <- 0 #data$x["Dist"]
  data$x <- select(data$x, variables)
  add <-cname[!cname%in%names(data$y)]
  if(length(add)!=0) data$y[add] <- 0
  data$y <- select(data$y, variables)
  data
}

haplo.splines1 <- fncols(haplo.splines1, variables)
haplo.splines2 <- fncols(haplo.splines2, variables)
haplo.splines3 <- fncols(haplo.splines3, variables)
haplo.splines4 <- fncols(haplo.splines4, variables)

colnames(haplo.splines4$x)


pdf("./output_plots/spline_bestmodels", height=11, width=7)
par(mfrow=c(4,2), mar=c(4.5,4.5,2,2))
namestring <- colnames(haplo.splines1$x)
for (i in 1:ncol(haplo.splines1$x)) {
  ymaxs <- lapply(list(haplo.splines1,haplo.splines2,haplo.splines3,haplo.splines4),
                  function (x) apply(x[[2]], 2, max)) %>% as.data.frame() %>% t()
  xmaxs <- lapply(list(haplo.splines1,haplo.splines2,haplo.splines3,haplo.splines4),
                  function (x) apply(x[[1]], 2, max)) %>% as.data.frame() %>% t()
  xmins <- lapply(list(haplo.splines1,haplo.splines2,haplo.splines3,haplo.splines4),
                  function (x) apply(x[[1]], 2, min)) %>% as.data.frame() %>% t()
  plot(haplo.splines1$x[,i], haplo.splines1$y[,i], type="l", lwd=3, col=varcols[1],
       xlab=paste(namestring[i], "dissimilarity"), ylab=paste("f(",namestring[i],")",sep=""),
       ylim=c(0, 1.1*max(ymaxs[,i])), xlim=c(min(xmins[,i]), max(xmaxs[,i])),
       cex.axis=1.2, cex.lab=1.4)
  lines(haplo.splines2$x[,i], haplo.splines2$y[,i], lwd=3, col=varcols[1], lty=2)
  lines(haplo.splines3$x[,i], haplo.splines3$y[,i], lwd=3, col=varcols[4])
  lines(haplo.splines4$x[,i], haplo.splines4$y[,i], lwd=3, col=varcols[4], lty=2)
  #for (j in 1:4){
  # modlist <- list(haplo.splines1,haplo.splines2,haplo.splines3, haplo.splines4)
  #text(x = 1.1*max(modlist[[j]]$x[,i]) , y = max(modlist[[j]]$y[,i]),
  #    labels = names[j], cex=1)
}
plot(NA, NA, xlim=c(1,1), ylim=c(1,1), axes=F, ylab="", xlab="")
legend("center", legend = c("Parasite species", "Parasite phylogenetic", "Host species","Host phylogenetic"), col = varcols[c(1,1,4,4)], lwd=3, cex=1.5,
       title="Turnover model", lty=c(1,2,1,2))
dev.off()



#
# Maps --------------------------------------------------------------------


