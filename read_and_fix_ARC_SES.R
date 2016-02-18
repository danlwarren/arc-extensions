library(ggplot2)
library(RColorBrewer)

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
