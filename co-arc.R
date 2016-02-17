setwd("~/OneDrive/R Projects/Banksia ARC/PolyEnvSES/")
require(ape)
source("age.range.correlation.SES.R")

tree <- read.nexus(file="Concatenated alignments 20.trees.124.mcc.txt")
overlap <- read.csv("banksia_cooc_cij.csv", header=TRUE, row.names=1)
speciesnames <- row.names(overlap)
keptspecies <- vector()
keptspecies <- speciesnames[speciesnames %in% tree$tip.label]
lostspecies <- speciesnames[!speciesnames %in% tree$tip.label]
#ecodata <- read.csv("data_banksia_allspp_flowermonths.csv", header=TRUE)
#rownames(ecodata) <- ecodata[,1]
#keptspecies <- keptspecies[keptspecies %in% rownames(ecodata)]
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keptspecies])
geo.for.arc <- overlap[keptspecies, keptspecies]
colnames(geo.for.arc) <- keptspecies
rownames(geo.for.arc) <- keptspecies
cooc_geoarc <- age.range.correlation.SES(tree, geo.for.arc, n=500)

save(cooc_geoarc, file="cooc_geoarc")
