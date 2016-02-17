setwd("~/OneDrive/R Projects/Banksia ARC/PolyEnvSES")
library(ggplot2)
library(RColorBrewer)
load("cluster_geoarc")
load("cooc_geoarc")
load("poly_geoarc")

fixp <- function(x){
  f.intercept <- length(which(x$MonteCarlo.replicates[1, ] > x$linear.regression$coefficients[1]))/length(x$MonteCarlo.replicates[1,])
  f.slope <- length(which(x$MonteCarlo.replicates[2, ] > x$linear.regression$coefficients[2]))/length(x$MonteCarlo.replicates[2,])
  f <- c(f.intercept, f.slope)
  p <- sapply(f, function(x) 2 * min(x, 1 - x))
  sig <- cbind(f, p)
  rownames(sig) <- c("intercept", "slope")

  list(age.range.correlation = x$age.range.correlation, linear.regression = x$linear.regression,
       sig = sig, MonteCarlo.replicates = x$MonteCarlo.replicates, sim.residuals = x$sim.residuals,
       empirical.residuals = x$empirical.residualls, standard.effects = x$standard.effects)
}

fixed_cluster <- fixp(cluster_geoarc)
fixed_cooc <- fixp(cooc_geoarc)
fixed_poly <- fixp(poly_geoarc)

fixed_cluster$sig
fixed_cooc$sig
fixed_poly$sig



plot.clustarc.ses <- function(x, ...){
  ses.df <- data.frame(cbind(x$age.range.correlation, abs(x$standard.effects)))
  names(ses.df) <- c("age", "overlap", "ses")
  my.colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  ses.plot <- ggplot(ses.df, aes(x=age, y=overlap, colour=ses))

  # Add lines from reps.  Intercepts are in row 1, slopes in row 2
  for(i in 1:length(x$MonteCarlo.replicates[1,])){
    ses.plot <- ses.plot + geom_abline(slope=x$MonteCarlo.replicates[2,i],
                                       intercept=x$MonteCarlo.replicates[1,i],
                                       color="grey86")
  }

  ses.plot <- ses.plot + geom_point(aes(size=5)) +
    scale_color_gradientn(colours = my.colors(100)) +
    geom_abline(slope=x$linear.regression$coefficients[2], intercept=x$linear.regression$coefficients[1]) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())


  return(ses.plot)
}
plot.clustarc.ses(fixed_cooc)
plot.clustarc.ses(fixed_cluster)
plot.clustarc.ses(fixed_poly)

#
#
# png(file="Geo_ARC.png",  res=200, width=1440, height=1440)
# geoplot <- plot.clustarc(geoarc, xlab="Age", ylab="Overlap", main="Geography", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Temp_ARC.png",  res=200, width=1440, height=1440)
# geoplot <- plot.clustarc(temparc, xlab="Age", ylab="Overlap", main="Temperature", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Growth_ARC.png",  res=200, width=1440, height=1440)
# growthplot <- plot.clustarc(growtharc, xlab="Age", ylab="Overlap", main="Growth", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Moisture_ARC.png",  res=200, width=1440, height=1440)
# moistplot <- plot.clustarc(moistarc, xlab="Age", ylab="Overlap", main="Moisture", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Rad_ARC.png",  res=200, width=1440, height=1440)
# radplot <- plot.clustarc(radarc, xlab="Age", ylab="Overlap", main="Radiation", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Soil_ARC.png",  res=200, width=1440, height=1440)
# soilplot <- plot.clustarc(soilarc, xlab="Age", ylab="Overlap", main="Soil Type", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Wind_ARC.png",  res=200, width=1440, height=1440)
# windplot <- plot.clustarc(windarc, xlab="Age", ylab="Overlap", main="Wind", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Precip_ARC.png",  res=200, width=1440, height=1440)
# precplot <- plot.clustarc(precarc, xlab="Age", ylab="Overlap", main="Precipitation", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
# png(file="Prod_ARC.png",  res=200, width=1440, height=1440)
# prodplot <- plot.clustarc(prodarc, xlab="Age", ylab="Overlap", main="Productivity", las=1, cex.lab=1.5, cex.axis=1.5, cex.main=2)
# dev.off()
#
#
# fixed_geoarc$sig
# fixed_geoarc$linear.regression
# fixed_growtharc$sig
# fixed_growtharc$linear.regression
# fixed_heightarc$sig
# fixed_heightarc$linear.regression
# fixed_moistarc$sig
# fixed_moistarc$linear.regression
# fixed_seedarc$sig
# fixed_seedarc$linear.regression
# fixed_windarc$sig
# fixed_windarc$linear.regression
# fixed_precarc$sig
# fixed_precarc$linear.regression
# fixed_timearc$sig
# fixed_timearc$linear.regression
# fixed_prodarc$sig
# fixed_prodarc$linear.regression
# fixed_sproutarc$sig
# fixed_sproutarc$linear.regression
# fixed_radarc$sig
# fixed_radarc$linear.regression
# fixed_soilarc$sig
# fixed_soilarc$linear.regression

