# Script to create distance matrices (response variables) for GDM.
# Distance matrices quantify turnover between communities
# We will use Jaccard / Bray Curtis to quantify species turnover
# We will use Unifracs to quantify phylogenetic turnover.
# Created using R version 3.5.1 in R Studio November 2019

set.seed(1987)
library(vegan) #for rarefying
library(iNEXT) #for rarefying and estimating parasite community richness
library(gplots) #heatmaps
library(GUniFrac) #for calculating unifrac distance
library(ape) #for reading phylogeny
library(dplyr) #for data organization
library(stringr) #for renaming strings

# Input data ----------------------------------------------------------------

# Parasite data: csv list of parasites, each row is an infection
# includes lat long for each record, the community number it was collected from
# and "Lineage.name", the parasite species.

inputhaplos <- read.csv("./raw_data/grouped-haplos.csv", row.names = 1)
filter(inputhaplos, community.number ==7)
#tree of all malaria lineages
malariatree <- read.nexus("./raw_data/perumalariatree.nex")

# Host data: Each row is a bird species, each column is a community. Values are
# presence or absence of that species in that community.
# Clean data to make binary and then transform and put into df.Data are in
# formatted data folder because host lists were compiled in "compiling_hosts.R"
inputhosts <- read.csv("./formatted_data/allhostlist.csv", row.names = 1,
                       col.names=paste("X", 1:18)) %>%
                        t(.) %>% as.data.frame(.)
colnames(inputhosts) <- colnames(inputhosts) %>% str_replace_all(., " ", "_") #replace space with underscore to match phylo
inputhosts[c(1,7),] %>% apply(., 2, sum) %>% as.data.frame() -> tempdat
names(tempdat) <- "col1"
filter(tempdat, col1 == 2) %>% nrow()



# Tree of hosts from Bird Trees
hosttree <- read.nexus("./raw_data/host-tree-oct12c.nex")

# Community metadata: each row is a community.
metadata1 <- read.csv("./raw_data/bare_metadata.csv", row.names=1)


# Host species and phylogenetic turnover ---------------------------------------
# We only have abundance for our host data,
# so we will use Jaccard and unweighted Unifracs to characterize differences
# between communities. We don't have to rarefy because our communites are sampled equally
# Key dataframes: inputhosts and hosttree

#Are there any tree tips not in host list or vice versa?
hosttree$tip.label[!(hosttree$tip.label %in% colnames(inputhosts))]
colnames(inputhosts)[!(colnames(inputhosts) %in% hosttree$tip.label)]
#"Colaptes_cinereicappillus" "Piezorina_cinerea" spelling need to be changed
colnames(inputhosts)[colnames(inputhosts) =="Colaptes_cinereicappillus"] <- "Colaptes_cinereicappilus"
colnames(inputhosts)[colnames(inputhosts) =="Piezorina_cinerea"] <- "Piezorhina_cinerea"

#Unweighted Unifrac:
hostUnifrac <- GUniFrac(inputhosts, hosttree, alpha=c(0, 0.5, 1)) %>%
  .$unifracs %>% .[, , "d_UW"]

#Jaccard
hostJ <- as.matrix(vegdist(inputhosts, method="jaccard", binary=T, diag=T)) %>%
  as.data.frame()

#write out files
write.csv(hostJ, "formatted_data/hostjaccard.csv")
write.csv(hostUnifrac, "formatted_data/hostUnifrac.csv")

#
# Rarefy parasites --------------------------------------------------------
# Sampling effort and therefore number of parasites we detected varied among
# communities. To control for this variation, we'll rarefy (subsample) infections
# to the minimum number detected in any one community.

# 1. Create site x haplotype abundance matrix; row = communities, cols = lineages
# Values are the number of infections detected per community.
peruhaploabun <- as.data.frame.matrix(
  table(inputhaplos$community.number, inputhaplos$Lineage.name))

#peruspec <- apply(peruhaploabun,1, function(x) sum(x!=0)) # observed number of haplotypes per community

# 2. Figure out what minimum number of infections is
rowSums(peruhaploabun) # observed number of infections per community
infection.min <- min(rowSums(peruhaploabun))
#peruspec.rar <- rrarefy(peruhaploabun, infection.min) #this rarefies dataset once.

# Parasite species turnover  --------------------------
# Species turnover can be calculated using Jaccard (presence absence) or
# Bray-Curtis (includes abundance info) distance. Downstream we only used BC
# because it includes more info and Jaccard vals were all very close to 1

# Basic function. vegdist returns a vector, needs to be put into matrix.
#b <- as.matrix(vegdist(peruspec.rar, method="bray", binary=T, diag=T))

# avdist: calculate distance many times from different rarefied samples (fast)

meanBC <- avgdist(peruhaploabun, sample=infection.min, dmethod="bray",
                  iterations = 1000) %>% as.matrix %>% as.data.frame()
write.csv(meanBC, "formatted_data/parasitebraycurtis.csv")

#meanJ <- avgdist(peruhaploabun, sample = infection.min, dmethod="jaccard",
#                 interations = 1000, binary=T) %>% as.matrix %>% as.data.frame()
#write.csv(meanJ, "parasitejaccard.csv")



# Parasite phylogenetic turnover  -----------------------------------------
# Unifracs can be unweighted (presence/absence), weighted (takes abundance into account)
# or "generalized" some combo of the both that appears to be the best
# Unweighted is "d_UW", weighted is "d_1", generalized is "d_0.5"
# Create unifrac distances using our rarefied dataset (peruspec.rar) and the lineage tree (perutree)
# Check and make sure our lineags from the tree and from the matrix match,
# should be 0 for both
malariatree$tip.label[!(malariatree$tip.label %in% colnames(peruspec.rar))]
colnames(peruspec.rar)[!(colnames(peruspec.rar) %in% malariatree$tip.label)]

# Make a function to extract the unifrac matrix
gUnifs <- function () {
  rrarefy(peruhaploabun, infection.min) %>%
    GUniFrac(., malariatree, alpha=c(0, 0.5, 1)) %>%
    .$unifracs %>% .[, , "d_0.5"] %>% return(.)
}

# Make a holding dataframe, so we can rarefy and calculate unifracs 1000 times
meangUnifs <- data.frame(matrix(0,nrow=18,ncol=18))
sapply(1:1000, function(x) {x; meangUnifs <<- meangUnifs + gUnifs()})
meangUnifs <- meangUnifs/1000

heatmap.2(as.matrix(meangUnifs)) #heat map to plot distances
write.csv(meangUnifs, file = "formatted_data/malaria.generalized.unifracs.csv")

#
# Turnover each parasite genus -------------------------------
# The code calculates all distance matrices just in case
#create vectors of lineages in each genus
filter(inputhaplos, Genus=="H"|Genus=="C") %>%.$Lineage.name %>%
  unique() %>% droplevels()->h.lineages
filter(inputhaplos, Genus=="L") %>% .$Lineage.name %>% unique() %>% droplevels() ->l.lineages
filter(inputhaplos, Genus=="P") %>% .$Lineage.name %>% unique() %>% droplevels ()->p.lineages
length(h.lineages) + length(l.lineages) + length(p.lineages) #check to see equals 400

#subset trees. setdiff() drops all tips where the tip labels aren't x.lineages
htree <- drop.tip(malariatree, setdiff(malariatree$tip.label, h.lineages))
ptree <- drop.tip(malariatree, setdiff(malariatree$tip.label, p.lineages))
ltree <- drop.tip(malariatree, setdiff(malariatree$tip.label, l.lineages))

#subset and create abundance dataframes, communities = rows, cols = lineages
filter(inputhaplos, Genus == "H"|Genus =="C") %>% droplevels() %>%
  with(., table(community.number, Lineage.name)) %>%
  as.data.frame.matrix()-> h.abund

filter(inputhaplos, Genus == "P") %>% droplevels() %>%
  with(., table(community.number, Lineage.name)) %>%
  as.data.frame.matrix()-> p.abund

filter(inputhaplos, Genus == "L") %>% droplevels() %>%
  with(., table(community.number, Lineage.name)) %>%
  as.data.frame.matrix()-> l.abund


#put everything in a list because lists are great
parasitetrees <- list(htree, ltree, ptree)
parasiteabund <- list(h.abund, l.abund, p.abund)


#Problem: by subsetting data for some communities we have very few infections
#Solution (?) drop communities with < N infections
minsampling <- 19

#figure out which commnities are left
#Creates a list for H, P, L, of which communities satisfy minimum sampling
community.numbers <-
  sapply(parasiteabund, function(x)
    as.numeric(rownames(x)[rowSums(x) > minsampling]))
names(community.numbers) <-c(rep("community.number",3))
lapply(community.numbers, length) #number of communities for each genera



#Create presence absence matrices for each genus of parasites
parasiteabund <- lapply(parasiteabund, function(x) x %>%
                          mutate(total = rowSums(.)) %>%  #add row sums
                          filter(total > minsampling) %>% #select communities w > min inf
                          dplyr::select(-total)) #remove rowsums

infectionmin <- lapply(parasiteabund, function(x) rowSums(x)) #infs per comm.


#Unifracs: create some functions to extract our weighted and unweighted unifs

wUnifs <- function (abund, tree) {
  rrarefy(abund, min(rowSums(abund))) %>%       #rarefy dataframe
    GUniFrac(., tree, alpha=c(0, 0.5, 1)) %>% #gunifrac fun, default wght
    .$unifracs %>% .[, , "d_1"] %>% return(.) #extract just weighted unif df
}

uUnifs <- function (abund, tree) {
  rrarefy(abund, min(rowSums(abund))) %>%       #rarefy dataframe
    GUniFrac(., tree, alpha=c(0, 0.5, 1)) %>% #gunifrac fun, default wght
    .$unifracs %>% .[, , "d_UW"] %>% return(.) #extract just weighted unif df
}

gUnifs <- function (abund, tree) {
  rrarefy(abund, min(rowSums(abund))) %>%       #rarefy dataframe
    GUniFrac(., tree, alpha=c(0, 0.5, 1)) %>% #gunifrac fun, default wght
    .$unifracs %>% .[, , "d_0.5"] %>% return(.) #extract just generalized unif df
}

#This is confusing but good code: For 1000 iterations calculate the wUnifs fun
#on our three abundance matrices using 3 different trees. For each wUnif frame
#add to existing matrix. After we'll divide by n interations to get the mean.

hplstore <- lapply(parasiteabund, function (x)
  data.frame(matrix(0,nrow=nrow(x), ncol=nrow(x)))) #1) create a set of dataframes

sapply(1:1000, function(x) {x; hplstore <<- mapply( "+",
                                                    mapply(wUnifs, parasiteabund, parasitetrees) , hplstore)}) # 2) run "loop"
hplstore <- lapply(hplstore, function(x) x/1000) # 3) calculate av.

sapply(1:3, function(x) write.csv(hplstore[[x]],
                                  file=paste("./formatted_data/wunifs", x, "csv",
                                             sep="."), row.names = F)) #write out results

#do the same for unweighted unifs
hplstoreU <- lapply(parasiteabund, function (x)
  data.frame(matrix(0,nrow=nrow(x), ncol=nrow(x))))
sapply(1:1000, function(x) {x; hplstoreU <<- mapply( "+",
                                                     mapply(uUnifs, parasiteabund, parasitetrees) , hplstoreU)})
hplstoreU <- lapply(hplstoreU, function(x) x/1000) # 3) calculate av.
sapply(1:3, function(x) write.csv(hplstoreU[[x]],
                                  file=paste("./formatted_data/uunifs", x,
                                             "csv", sep="."), row.names = F))

#do the same for generalized unifs
hplstoreG <- lapply(parasiteabund, function (x)
  data.frame(matrix(0,nrow=nrow(x), ncol=nrow(x))))
sapply(1:1000, function(x) {x; hplstoreG <<- mapply( "+",
                                                     mapply(gUnifs, parasiteabund, parasitetrees) , hplstoreG)})
hplstoreG <- lapply(hplstoreG, function(x) x/1000) # 3) calculate av.
sapply(1:3, function(x) write.csv(hplstoreG[[x]],
                                  file=paste("./formatted_data/gunifs", x,
                                             "csv", sep="."), row.names = F))



#Weighted and unweighted unifs are written out as "wunifs-1" etc. 1=H, 2=L, 3=P


#Non-phylogenetic methods
jaccards <- lapply(parasiteabund, function (x) avgdist(x, sample = min(rowSums(x)),
                                                       dmethod="jaccard", interations = 1000, binary=T))
jaccards <- lapply(jaccards, function (x) as.data.frame(as.matrix(x)))
sapply(1:3, function (x) write.csv(jaccards[[x]],
                                   file=paste("./formatted_data/jaccards", x,
                                              "csv", sep="."), row.names = F))

braycurtis <- lapply(parasiteabund, function (x) avgdist(x, sample = min(rowSums(x)),
                                                         dmethod="bray", interations = 1000, binary=T))
braycurtis <- lapply(braycurtis, function (x) as.data.frame(as.matrix(x)))
sapply(1:3, function (x) write.csv(braycurtis[[x]],
                                   file=paste("./formatted_data/braycurtis", x,
                                              "csv", sep="."), row.names = F))


# Estimate parasite community richness & EXPORT METADATA ------------
# Use the parasite abundance matrix to estimate parasite species richness
# for our communities (one of our GDM predictors) using iNEXT package.

outA <- iNEXT(t(peruhaploabun), q=0, datatype="abundance") #transform matrix
metadata1$parasite.richness <- outA$AsyEst$Estimator[outA$AsyEst$Diversity=="Species richness"]
metadata1$parasite.observed <- outA$AsyEst$Observed[outA$AsyEst$Diversity=="Species richness"]

#cleanup and writeout metdata; remove bad bioclim columns.
write.csv(metadata1, file = "./formatted_data/GDM_metadata.csv",
          row.names=F)

# Heatmaps of distance matrices -----------------------------------------------------------

pdf("./output_plots/heat_map_HostU.pdf")
heatmap.2(as.matrix(hostUnifrac[,-1]),  main="Host UniFracs", trace="none")
dev.off()
pdf("./output_plots/heat_map_HostJ.pdf")
heatmap.2(as.matrix(hostJ[,-1]),  main="Host Jaccard", trace="none")
dev.off()
pdf("./output_plots/heat_map_ParGd.pdf")
heatmap.2(as.matrix(parGunifracs[,-1]),  main="Parasite Unifracs", trace="none")
dev.off()
pdf("./output_plots/heat_map_ParBC.pdf")
heatmap.2(as.matrix(parBC[,-1]),  main="Parasite Bray Curtis", trace="none")
dev.off()



