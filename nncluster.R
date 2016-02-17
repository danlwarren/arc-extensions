require(fields)  # Used for Euclidean distance matrices
nncluster <- function(x, y){
  y <- as.factor(y)
  levs <- levels(y)
  df <- as.data.frame(matrix (nrow=length(levels(y)), ncol=length(levels(y))))
  rownames(df) <- levels(y)
  colnames(df) <- levels(y)
  for(i in 1:length(levs)){
   for(j in i:length(levs)){
     if(i != j){
        if(nrow(x[y == levs[i],]) > 1 && nrow(x[y == levs[j],]) > 1){
            print(paste("Calculating", levs[i], "vs", levs[j]))
            within1 <- rdist(x[y == levs[i],],x[y == levs[i],])
            within2 <- rdist(x[y == levs[j],],x[y == levs[j],])
            between <- rdist(x[y == levs[i],],x[y == levs[j],])
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