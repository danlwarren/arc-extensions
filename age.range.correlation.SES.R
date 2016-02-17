# modified from the function in phyloclim
# uses a loop instead of sapply() for the randomizations
# this is just to allow a counter to be printed to the screen
require("phyloclim")

age.range.correlation.SES <- function (phy, overlap, tri = "upper", n = 1000)
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
        id <- sample(seq(along = o[, 1]))
        rownames(o) <- colnames(o) <- colnames(o)[id]
        rand.o <- sapply(names(age), nested.mean.overlap, phy = phy,
            olap = o)
        o <- cbind(age, rand.o)
        o <- lm(o[, 2] ~ age)
        list(coef = o$coefficients, fitted.values = o$fitted.values, random.overlap = rand.o)
    }
### MODIFIED BIT ###
    random.x <-array(dim=c(n,2))
    random.overlap <-array(dim=c(n,length(age)))

    for (i in 1:n){
      print(paste(i,"of",n,"iterations",sep=" "))
      print(Sys.time())
      thisrand  <- randomization(o = ovlap, phy = phy,age = age)
      random.x[i,] <- thisrand$coef
      random.overlap[i,] <- thisrand$random.overlap
    }
### END MODIFIED BIT ###

    f.intercept <- length(which(random.x[,1] > x.lm$coefficients[1]))/n
    f.slope <- length(which(random.x[,2] > x.lm$coefficients[2]))/n
    f <- c(f.intercept, f.slope)
    p <- sapply(f, function(x) 2 * min(x, 1 - x))
    sig <- cbind(f, p)
    rownames(sig) <- c("intercept", "slope")
    sds <- apply(random.overlap, 2, sd)
    means <- apply(random.overlap, 2, mean)
    print(means)
    print(sds)
    print(x[,"overlap"])
    ses <- (x[,"overlap"] - means)/sds
    random.overlap <- cbind(age, t(random.overlap))
    colnames(random.overlap) <- c("age", paste("rep", seq(1,n)))
    list(age.range.correlation = x, linear.regression = x.lm,
        sig = sig, sim.overlaps = random.overlap,
        MonteCarlo.replicates = t(random.x), standard.effects = ses, means = means, sds = sds)
}
