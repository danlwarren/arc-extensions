setwd("~/Dropbox/R Projects/Banksia ARC/Co-occurrence/")
require(ape)
source("age.range.correlation.2.R")
require(phyloclim)
require(picante)

coocurrence.arc <- function (phy, overlap, tri = "upper", n = 1000, null.model = "independentswap") 
{
  age <- branching.times(phy)
  ovlap <- overlap
  if (tri == "upper") 
    ovlap[lower.tri(ovlap)] <- t(ovlap)[lower.tri(ovlap)]
  if (tri == "lower") 
    ovlap[upper.tri(ovlap)] <- t(ovlap)[upper.tri(ovlap)]
  id <- match(phy$tip.label, rownames(ovlap))
  #print(phy$tip.label)
  #print(rownames(ovlap))
  #print(id)
  ovlap <- ovlap[id, id]
  overlap <- sapply(names(age), nested.mean.overlap, phy = phy, 
                    olap = ovlap)
  x <- cbind(age, overlap)
  x.lm <- lm(overlap ~ age)
  randomization <- function(phy, o, n, age) {
    o <- randomizeMatrix(o , null.model=null.model)
    o <- sapply(names(age), nested.mean.overlap, phy = phy, 
                olap = o)
    o <- cbind(age, o)
    o <- lm(o[, 2] ~ age)
    o$coefficients
  }
  ### MODIFIED BIT ###
  random.x <-array(dim=c(n,2))
  
  for (i in 1:n){
    print(paste(i,"of",n,"iterations",sep=" "))
    print(Sys.time())
    random.x[i,] <- randomization(o = ovlap, phy = phy,age = age)
  }
  ### END MODIFIED BIT ###
  
  f.intercept <- length(which(random.x[,1] > x.lm$coefficients[1]))/n
  f.slope <- length(which(random.x[,2] > x.lm$coefficients[2]))/n
  f <- c(f.intercept, f.slope)
  p <- sapply(f, function(x) 2 * min(x, 1 - x))
  sig <- cbind(f, p)
  rownames(sig) <- c("intercept", "slope")
  list(age.range.correlation = x, linear.regression = x.lm, 
       sig = sig, MonteCarlo.replicates = t(random.x))
}


fixp <- function(x){
  f.intercept <- length(which(x$MonteCarlo.replicates[1, ] > x$linear.regression$coefficients[1]))/length(x$MonteCarlo.replicates[1,])
  f.slope <- length(which(x$MonteCarlo.replicates[2, ] > x$linear.regression$coefficients[2]))/length(x$MonteCarlo.replicates[2,])
  f <- c(f.intercept, f.slope)
  p <- sapply(f, function(x) 2 * min(x, 1 - x))
  sig <- cbind(f, p)
  rownames(sig) <- c("intercept", "slope")
  list(age.range.correlation = x$age.range.correlation, linear.regression = x$linear.regression, 
       sig = sig, MonteCarlo.replicates = x$MonteCarlo.replicates)
}

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
geoarc <- fixp(coocurrence.arc(tree, geo.for.arc, n=5, null.model="frequency"))

save(geoarc, file="geoarc")

plot.clustarc <- function(x, ...){
  p <- plot(x$age.range.correlation, type="n", ...)
  for(i in 1:length(x$MonteCarlo.replicates[1,])){
    #print((geoarc$MonteCarlo.replicates[i,2], geoarc$MonteCarlo.replicates[i,1]))
    abline(a = x$MonteCarlo.replicates[1,i], b = x$MonteCarlo.replicates[2,i], col="gray90")
  }
  abline(a = x$linear.regression$coefficients[1], b = x$linear.regression$coefficients[2], col = "black")
  points(x$age.range.correlation)
}

png(file="Co-occurrence_ARC.png",  res=200, width=1440, height=1440)
plot.clustarc(geoarc, xlab="Age", ylab="Cij", main="Co-Occurrence",  las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
dev.off()