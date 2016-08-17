require(fields)  # Used for Euclidean distance matrices

# xy is a two-columnmatrix of x and y values
# species is a vector of species names of the same length as ncol(xy)
nncluster <- function(xy, species){

  # Figure out how many species we've got, initialize a data frame
  species.names <- unique(species)
  df <- as.data.frame(matrix (nrow=length(species.names), ncol=length(species.names)))
  rownames(df) <- species.names
  colnames(df) <- species.names

  # Do each pair of species
  for(i in 1:length(species.names)){
    for(j in i:length(species.names)){
      if(i != j){

        # Need at least two points to calculate
        if(nrow(xy[species == species.names[i],]) > 1 && nrow(xy[species == species.names[j],]) > 1){

          print(paste("Calculating", species.names[i], "vs", species.names[j]))
          within1 <- rdist(xy[species == species.names[i],],xy[species == species.names[i],])
          within2 <- rdist(xy[species == species.names[j],],xy[species == species.names[j],])
          between <- rdist(xy[species == species.names[i],],xy[species == species.names[j],])
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
