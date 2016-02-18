setwd("~/Dropbox/R Projects/Banksia ARC/MultiArc")
library(ape)
library(phyloclim)
source("age.range.correlation.2.R")

trees <- read.nexus(file="Concatenated alignments 20.trees.124.txt")
overlap <- read.csv("geo_overlap.csv", header=TRUE, row.names=1)
rownames(overlap) <- colnames(overlap)
speciesnames <- row.names(overlap)
keptspecies <- vector()
keptspecies <- speciesnames[speciesnames %in% trees[[1]]$tip.label]
ecodata <- read.csv("data_banksia_allspp_flowermonths.csv", header=TRUE)
rownames(ecodata) <- ecodata[,1]
keptspecies <- keptspecies[keptspecies %in% rownames(ecodata)]
trees <- lapply(trees, drop.tip, tip=as.vector(unlist(trees$tip.label))[!as.vector(unlist(trees$tip.label[[1]])) %in% keptspecies])
geo.for.arc <- overlap[keptspecies, keptspecies]
colnames(geo.for.arc) <- keptspecies
rownames(geo.for.arc) <- keptspecies
geo.for.arc[lower.tri(geo.for.arc)] <- t(geo.for.arc)[lower.tri(geo.for.arc)]

multiarc <- function(overlap, trees, burnin=0, n=length(trees) - burnin, tri = "upper", reps.per.tree = 1){
    trees <- trees[(burnin + 1):(length(trees))]
    trees <- trees[seq(from = 1, to = length(trees), length.out = n)]
    #print(trees)
    empirical.x <-array(dim=c(0,2))
    random.x <-array(dim=c(0,2))

    for(i in 1:length(trees)){
        print(paste("Replicate", i))
        x <- age.range.correlation(trees[[i]], overlap, tri, reps.per.tree)
        #print(random.x)
        print(x)
        random.x <- rbind(random.x, x$MonteCarlo.replicates)
        empirical.x <- rbind(empirical.x, x$linear.regression$coefficients)
        print(dim(random.x))
        print(dim(empirical.x))
    }
    list(empirical.replicates = empirical.x, MonteCarlo.replicates = random.x)
}

