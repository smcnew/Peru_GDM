# Analysis McNew et al. 2019 in prep
## Author: Sabrina McNew
***
***
#### Code and data for Generalized Dissimilarity Modeling 

### File structure 
The code and files are set up to be run as a [project](https://r4ds.had.co.nz/workflow-projects.html)
in RStudio. This means file paths are relative within the script (for ex. "raw_data" will be a subdirectory within 
whatever working directory you place the project in).  

Two directories were not synced to GitHub: "GDM_results" and "output_plots", because their contents are 
too large. You should create them before starting. 

### A fairly complete list of packages needed (check beginning of each script
### for complete list): 
```

# Data wrangling 
library(dplyr) 
library(stringr) 
library(scales) 
library(tools) #base 

# Spatial data tools 
library(raster) 
library(rgdal)
library(RStoolbox) 

# Community ecology tools  
library(iNEXT) 
library(gplots)
library(vegan)
library(gdm)
library(GUniFrac)
library(ape)

# Generalized dissimilartiy modeling. Check documentation to make sure best 
# version of package is still on CRAN.

library(gdm) 


```
### Guide to scripts (more or less in the order they should be run)

1. compiling_hosts.R
+ code to assemble lists of hosts based on BirdLife shapefiles and museum records. 
2. prep_distance_matrices.R
+ generates distance matrices (i.e. quantifies turnover) based on presence/absence or abundance data.
3. prep_spatial_dat.R
  + formats spatial data for our climatic predictors (i.e. temp/precip)
4. bird_species_richness_map.R  
  + creates a raster of bird abundance over study area. You will need to procure
shapefiles from BirdLife; they're too big to upload and store easily.   
5. gdm_host_parasite.R
  + Run and save GDM models.  
6. gdm_result_plots.R
  + Plot and format GDM results.  
7. gam_host_parasite_richness.R 
+ Model spatial variation in host and parasite richness using generalized additive models (GAMs)

### Output information
The folder turnover_rasters contains rasters of species and phylogenetic turnover
for both hosts and parasites. The rasters can be read into R using the stack()
function from the Raster package. Each stack contains 3 layers, corresponding
to the first 3 PC axes of the dimension-reduced best model (see gdm_result_plots
or GDM vignette for more information). Maps are created using the plotRGB()
function where the red, green, and blue channels correspond to the first 3  
pc axes. 