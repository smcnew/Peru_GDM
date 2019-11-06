# Analysis McNew et al. 2019 in prep
## Author: Sabrina McNew
***
***
### Code and data for Generalized Dissimilarity Modeling 

#### File structure 
The code and files are set up to be run as a [project](https://r4ds.had.co.nz/workflow-projects.html)
in RStudio. This means file paths are relative within the script (for ex. "raw_data" will be a subdirectory within 
whatever working directory you place the project in).  

Two directories were not synced to GitHub: "formatted_data" and "output_plots", because their contents are 
too large. You should create them before starting. 

#### A fairly complete list of packages needed: 
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
#### Guide to scripts (more or less in the order they should be run)

1. prep_distance_matrices.R
+ generates distance matrices (i.e. quantifies turnover) based on presence/absence or abundance data.
2. map_prep.R
+ formats spatial data for our climatic predictors (i.e. temp/precip)
3. bird_species_richness_map.R 
+ creates a raster of bird abundance over study area. You will need to procure
shapefiles from BirdLife; they're too big to upload and store easily. 
4. gdm_host_parasite.R
+ Run and save GDM models 
5. gdm_result_plots.R
+ Plot and format GDM results 
6. gam_host_parasite_richness.R 
+ Model spatial variation in host and parasite richness using generalized additive models (GAMs)