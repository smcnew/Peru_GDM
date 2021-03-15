# Maps of turnover and richness
Packages needed: raster, rgdal
## Turnover maps
Maps are uploaded as raster files (.grd and .gri) which encode contain spatial data and can be read into R (or arcGIS or other program of your choice). Each file contains three layers, corresponding to the first three PC axes of the dimension-reduced best model. The rasters can be read into R using the raster::stack() function (e.g. raster::stack("Downloads/gdm_host_phy_raster")). Maps as in figure 4 can be reproduced using the raster::plotRGB() function. See code and comments in gdm_result_plots.R and/or GDM vignette for more information. 

Included are: 
1. gdm_host_phy_raster (bird phylogenetic turnover) 
2. gdm_host_spp_raster (bird species turnover) 
3. gdm_par_phy_raster (parasite phylogenetic turnover)
4. gdm_par_spp_raster (parasite species turnover)

Pdfs of the maps (as seen in figure 4) are uploaded as 
1. host_phylo_turnover_map.pdf
2. host_species_turnover_map.pdf
3. parasite_phylo_turnover_map.pdf
4. parasite_species_turnover_map.pdf

## Species richness map
Raster of bird species richness is included as a single layer (birdrichness.grd/gri)
It can be read into R using raster::raster("path/birdrichness")
It can be plotted using raster::plot(birdrichness)
It can be visualized ouside of R by opening the pdf (bird_richness.pdf). 

Occasionally, pdfs look blurry in mac preview. If this is true for you, try opening in another viewer (e.g. Adobe Illustrator/Reader). 
