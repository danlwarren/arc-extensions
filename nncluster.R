require(fields)  # Used for Euclidean distance matrices

# xy is a two-columnmatrix of x and y values
# species is a vector of species names of the same length as ncol(xy)
nncluster <- function(xy, species){

  # Figure out how many species we've got, initialize a data frame
  species <- as.factor(species)
  levs <- levels(species)
  df <- as.data.frame(matrix (nrow=length(levels(species)), ncol=length(levels(species))))
  rownames(df) <- levels(species)
  colnames(df) <- levels(species)

  # Do each pair of species
  for(i in 1:length(levs)){
    for(j in i:length(levs)){
      if(i != j){

        # Need at least two points to calculate
        if(nrow(xy[species == levs[i],]) > 1 && nrow(xy[species == levs[j],]) > 1){

          print(paste("Calculating", levs[i], "vs", levs[j]))
          within1 <- rdist(xy[species == levs[i],],xy[species == levs[i],])
          within2 <- rdist(xy[species == levs[j],],xy[species == levs[j],])
          between <- rdist(xy[species == levs[i],],xy[species == levs[j],])
          score1 <- rep(NA,length(within1[1,]))
          score2 <- rep(NA,length(within1[2,]))

          for(k in 1:length(within1[1,])){
            thisscore <- min(within1[k,-(k)])/min(between[k,])
            score1[k] <- thisscore
          }

          for(k in 1:length(within2[1,])){
            thisscore <- min(within2[k,-(k)])/min(between[,k])
            score2[k] <- thisscore
          }

          score1 <- length(which(score1>1))/length(score1)
          score2 <- length(which(score2>1))/length(score2)
          w <- mean(c(score1, score2))

        }
        else{w <- NA}
        df[i,j] <- w
      }
    }

  }
  return( df)
}
