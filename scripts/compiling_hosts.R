## Script to compile the list of hosts at each local
##

library(sp)
library(rgdal)
library(geosphere)
library(maptools)
library(sf)
library(raster)
library(letsR)
library(plyr)
library(dplyr)
library(tidyr)

# Data --------------------------------------------------------------------


# Spatial data
#peru <- raster::getData("GADM", country="PER", level=0) #download country shape
#projection(peru) <- CRS("+init=epsg:4326")
peru_alt <- raster("raw_data/peru_alt.grd") #elev raster
metadata1 <- read.csv(file = "raw_data/bare_metadata.csv") #locals of communities

CommunitySpatial <- SpatialPointsDataFrame(
  matrix(c(metadata1$community.Long,metadata1$community.Lat), ncol = 2),
  data.frame(ID=1:nrow(metadata1),community = metadata1$community,
             community.elev = metadata1$community.elev),
  proj4string=CRS("+init=epsg:4326"))


haplos <- read.csv("raw_data/grouped-haplos.csv") # list of infections
observedhosts <- read.csv("raw_data/observedhosts.csv") # list of all specimens collected by MSB
taxonomy <- read.csv("raw_data/taxonomy-synonyms.csv") #list of taxon name synonymized
birdelevscomplete <- read.csv("raw_data/birdelevscomplete.csv") #csv of elevation

# Host shape files
species <- readOGR("Bird_shapes") #2087 species
species <- spTransform (species, CRS("+init=epsg:4326"))
species2 <- lapply(species@polygons,
                   function(x) SpatialPolygons(list(x),
                  proj4string = CRS("+init=epsg:4326"))) #create a list of each species indpendently
setNames(species2, species@data$SCINAME)

head(species2)


# Visual Tests -------------------------------------------------------------------
dev.off()
plot(peru_alt)
plot(CommunitySpatial[4,], add=T, pch= 19) #visualise a community
plot(species[species@data$SCINAME=="Campylorhynchus fasciatus",], add=T) #plot a distribution
which("Campylorhynchus fasciatus" == species@data$SCINAME) #find index for species
plot(species2[[71]], add=T, lty = 2, col = "red") #make sure that the list of spatial polygons lines up as well
taxonomy[taxonomy$MSB.taxon=="Campylorhynchus fasciatus",] #check names; synonyms will be listed in "taxonomy"
birdelevscomplete %>% filter(taxon == "Campylorhynchus fasciatus") #check elevation
metadata1 %>% filter(community.number == 4) #check whether elevation ranges match

dev.off()

# Overlap shape files -----------------------------------------------------
# overlap shape files with community centers to create initial list of BirdLife hosts
communitydf <- lapply(species2, function(x) over(CommunitySpatial, x)) %>% #intersect communities over each species
  plyr::ldply() %>% #turn into dataframe
  replace(is.na(.), 0) %>% # make binary 0/1
  mutate(., taxon = species@data$SCINAME) %>% #add species scientific name
  ddply(., "taxon", numcolwise(sum)) #Each row is a polygon- so there are duplicates for spp w disjointed distributions. Sum so each row is 1 spp

head(communitydf)


communitydf$totals <- apply(communitydf[,2:19], 1, sum) #sum up number of spp. that intersected
communitydf <- communitydf[communitydf$totals > 0,] #remove species that don't intersect
communitydf <- merge(communitydf, birdelevscomplete[,1:3], all.x=T) #merge dataframe with taxonomic name list
sum(is.na(communitydf$Minimum.elevation)) #check that each spp has elevation
head(communitydf)

communitydf$taxon[duplicated(communitydf$taxon)] %>% length #check to make sure no spp are duplicated


# Elevation filter --------------------------------------------------------
## Now filter by elevation

sppelev <- mapply (c, communitydf$Minimum.elevation, communitydf$Maximum.elevation, SIMPLIFY = F) # create a list of min and max elev for each spp.
commelev <- metadata1$community.elev # a vector of the community elevations

#### intersect ranges of community elevation; keeps all species that could occur in elevation range of community.
commelevrange <- mapply(c, metadata1$minelev, metadata1$maxelev, SIMPLIFY=F) ## create a list of the elevation range of each community
longcommelevrange <-rep(commelevrange, each=length(sppelev)) #create a really long list, each comm repeated nspecies times,

dim(communitydf)

#elev.binary2, a dataframe with the presence or absence of each host in community elevation range.
elev.binary2 <- matrix(mapply(function (x, y)
  ifelse(sum(x[1]:x[2] %in% y[1]:y[2]) > 0, 1, 0),
  longcommelevrange, sppelev ), nrow=1362, ncol=18) ## intersect comm elev ranges and spp elev ranges. turn into matrix

elev.binary2 <- as.data.frame(elev.binary2) %>%  #clean up
  mutate(taxon = communitydf$taxon) %>%  #add taxon
  select(taxon, everything()) #make taxon #1 column


# Now combine shapefile dataframe with elevation dataframe
communitydf[,2:19] <- elev.binary2[,2:19]*communitydf[,2:19] # multiply distribution intersects x elevation intersects
communitydf$finaltotal <- apply(communitydf[2:19], 1, sum) #final total of presences in any community
communitydf <- communitydf[communitydf$finaltotal>0,] # filter dataframe to only include species present
dim(communitydf)  ##1274 host species


# Check community list against taxa we actually sampled ------------------------------------------------

#First, determine how many of our real sampled hosts (from collection) are on
#the list generated from Birdlife shapefiles

communityhostlist <-
  sapply(2:19, function(x)
    as.character(subset(communitydf$taxon, subset = communitydf[, x] == 1))) ## create a character list of taxa for ea. community
samplinghostlist <-
  aggregate(Taxon ~ community.number, data = haplos, function(x)
    as.character(unique(x)))$Taxon

hostsnotincommunity <-
  mapply(function(x, y)
    x[!x %in% y], samplinghostlist, communityhostlist) #list of sampled hosts not in BirdLife community list; includes quite a few
hostsincommunity <-
  mapply(function(x, y)
    x[x %in% y], samplinghostlist, communityhostlist) #list of sampled hosts in BirdLife community

sapply(1:18, function (x)
  length(hostsnotincommunity[[x]]) - sum(hostsnotincommunity[[x]] %in% as.character(species$SCINAME))) ##number of hosts not in distribution list
sapply(1:18, function (x)
  length(hostsincommunity[[x]]) / (length(hostsincommunity[[x]]) + length(hostsnotincommunity[[x]]))) ##Proportion of hosts in distribution list
sapply(1:18, function (x)
  paste(
    length(hostsincommunity[[x]]),
    length(hostsnotincommunity[[x]]),
    sep = ":"
  )) #ratio of hosts in to hosts outside

# Many hosts sampled by MSB not in BirdLife community lists, due to coarseness
# of shapefiles and taxonomy issues.

# Add MSB list  -----------------------------------------------------------
observedhosts$latlong <-
  with(observedhosts, paste(round(neg.degrees.Lat, 3),
                            round(neg.degrees.Long, 3),
                            sep =".")) #create a lat-long-elev ID within 3 decimal points, to identify sampling locals

haplos$latlong <-
  paste(round(haplos$neg.degrees.Lat, 3),
        round(haplos$neg.degrees.Long, 3),
        sep = ".") #create a lat-long ID within 3 decimal places

samplingsites <- select(haplos, community.number, latlong, Local, Elev, Dept, community) #subset haplos to just have id info
observedhosts2 <- merge(observedhosts, samplingsites[,c(1:2)], all.x=F, all.y=F) #Add community ID to MSB observed host list
length(unique(observedhosts2$latlong)) ### 252 sites in MSB list id'd out of 269 in our haplotype sampled list

sapply(1:18, function(x) unique(observedhosts2$Department[observedhosts2$community.number==x])) #Depts included in ea. community
sapply(1:18, function(x) length(observedhosts2$Department[observedhosts2$community.number==x])) #Samples included in ea. community
observedhosts2$latlong[observedhosts2$community.number==1] %>% unique()


# apply an elevation filter
observedhosts2<- select(metadata1, community.number, minelev, maxelev) %>% merge(., observedhosts2, all.y=T)

#add a column to tell whether the observation is within community elevation bounds
observedhosts2$inrange <- mapply(function (a,b,c) ifelse(c %in% a:b, 1,0),
                                 a = observedhosts2$minelev,
                                 b = observedhosts2$maxelev,
                                 c = observedhosts2$Elevation )
observedhosts2 <- observedhosts2[observedhosts2$inrange==1,] #filter out about 750 observations that aren't in range out of 100k

observedhostexport <- select(observedhosts2, Scientific.name, community.number) %>% #create smaller dataframe for export, with just taxa and community #
                      distinct() %>% #remove dups
                      droplevels()
table(observedhostexport$community) #number of taxa per community
write.csv(unique(observedhostexport$Scientific.name), "observedhostsinstudy.csv")

# Combine Birdlife and MSB lists ------------------------------------------
#turn MSB hosts into a species x community table

observedhoststable <- table(observedhostexport$Scientific.name,observedhostexport$community.number) %>%
  as.data.frame.matrix()

observedhoststable$MSB.taxon <- rownames(observedhoststable)
rownames(observedhoststable) <- c()

observedhoststable <- merge(observedhoststable, taxonomy, all.x=T, all.y=F) #synonimize taxonomic names
observedhoststable<- observedhoststable[!is.na(observedhoststable$BirdTree.taxon),] #remove taxa that aren't in the master list (i.e. bats n lizards)

#Merge taxonomy of BirdLife dataframe

communitydf <- dplyr::rename(communitydf, Birdlife.taxon = taxon) %>% merge(., taxonomy, all.x=T, all.y=F)

##combine observed and theoretical lists.

df1 <- select(communitydf, hostphylogeny.taxon, "1","2", "3", "4", "5", "6", "7", "8", "9",
              "10", "11", "12", "13", "14", "15", "16", "17", "18")
df2 <- select(observedhoststable, hostphylogeny.taxon, "1","2", "3", "4", "5", "6", "7", "8", "9",
              "10", "11", "12", "13", "14", "15", "16", "17", "18")


allhostlist <- rbind(df1, df2)
allhostlist <- ddply(allhostlist, "hostphylogeny.taxon", numcolwise(sum))


allhostlist[allhostlist > 1] <- 1 #Some cells have 2 where MSB and bird life agree. Change to 1
head(allhostlist)

apply(allhostlist[,2:19], 2, sum) #host richness
metadata1$total.host <- apply(allhostlist[,2:19], 2, sum) #add host richness to metadata

write.csv(allhostlist, "formatted_data/allhostlist.csv", row.names=FALSE)
write.csv(metadata1, file = "raw_data/bare_metadata.csv", row.names = F) #updated metadata

# Sanity checks -----------------------------------------------------------

duplicated(allhostlist$hostphylogeny.taxon) %>% sum() #see if any taxa are duplicated
check_all_hosts <- merge(taxonomy, allhostlist)
colnames(check_all_hosts) [5:22] <- paste("com_", 1:18, sep="")

# Check 1: Make sure that taxa with range maps overlap the point
head(check_all_hosts)

for (j in 1:18){
pdf(paste("community", j, ".pdf", sep=""), useDingbats = F)
  community <- paste("com_", j, sep="") %>% print()
taxa <- filter(check_all_hosts, get(community) ==1) %>% select("Birdlife.taxon") %>%
  pull() %>% as.character %>% stri_remove_empty() # All hosts in this community
random_taxa <- taxa[sample(length(taxa), 49)] #pull 49 random species
par(mfrow=c(7,7), mar = c(0,0,0,0))
for (i in 1:49) {
plot(peru_alt, legend = F, axes =F )
plot(CommunitySpatial[j,], add=T, pch= 19) #visualise a community
plot(species[species@data$SCINAME==random_taxa[i],], add=T) #plot a distribution
mtext(line = -1, adj=-.1, random_taxa[i], cex = 0.6)
}
dev.off()
}

## Double check species whose distributions weren't quite right
observedhosts2[observedhosts2$Scientific.name=="Phaeomyias murina"&observedhosts2$community.number==18,]

filter(taxonomy, Birdlife.taxon =="Thamnophilus shumbae" )
plot(peru_alt)
plot(species[species@data$SCINAME=="Grallaria nuchalis",], add=T)

# Check 2 add elevations and filter species out whose elevations don't
# intersect with community elevations
dim(check_all_hosts)
check_all_hosts$taxon <- check_all_hosts$hostphylogeny.taxon
check_all_hosts <- merge(check_all_hosts, birdelevscomplete, all.x=TRUE)

# for each community, match the species' elevational range against the
# community elevation range. Returns a list of the species that do not match
# elevation range for that community
for (j in 1:18) {
community <- paste("com_", j, sep="")
echeck <-NULL
subset <- filter(check_all_hosts, get(community) ==1)
subset <- na.omit(subset)
for (i in 1:nrow(subset)){
  echeck[i] <- ifelse(sum(metadata1$minelev[j]:metadata1$maxelev[j] %in%
                            subset$Minimum.elevation[i]:subset$Maximum.elevation[i])>0, TRUE, FALSE)
}
print(community)
which(echeck ==F) %>% print()
}

# Community 4
#Com elevation range: 1800 - 2500
#Elaenia chiriquensis collected < 350 m lower
#Myornis senilis collected > 350 m higher
#Asthenes fuliginosa 3150

# Community 13 problem: many species collected by MSB outside of "official spp range"
# but within elevation of community. Conclusion: Published elevational ranges not totally accurate.
echeck <-NULL
subset <- filter(check_all_hosts, com_13 ==1)
subset <- na.omit(subset)
for (i in 1:nrow(subset)){
  echeck[i] <- ifelse(sum(metadata1$minelev[13]:metadata1$maxelev[13] %in%
                            subset$Minimum.elevation[i]:subset$Maximum.elevation[i])>0, TRUE, FALSE)
}
which(echeck==F)
subset[185,]
observedhosts2[observedhosts2$Scientific.name=="Xenops minutus"&observedhosts2$community.number==13,]

haplos %>% filter(community.number ==13) %>% select(Local, Elev, Dept, neg.degrees.Lat, neg.degrees.Long)

#


