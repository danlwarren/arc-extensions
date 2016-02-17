setwd("~/Dropbox/R Projects/Banksia ARC/PolyEnvArc")
require(adehabitat)
require(maptools)
require(spatstat)
require(sp) # the classes and methods that make up spatial ops in R
require(maptools) # tools for reading and manipulating spatial objects
require(rgeos)
require(rgdal)
require(gpclib)
require(ggplot2)
require(scales)
source("NNCH.R")

pointdata <- read.csv("Banksia PC.csv", header=TRUE)
speciesnames <- levels(pointdata$Species)
keptspecies <- vector()
for(i in 1:length(speciesnames)){ #Looping over all species
  print(paste("Generating ranges for", speciesnames[i]))
  speciespoints <- subset(pointdata, Species == speciesnames[i])
  print(paste("Generating range for", speciesnames[i],"using", length(speciespoints[,1]), "points"))
  if(length(speciespoints[,1]) > 4){ #If there are enough points for a range map
    thisrange<-NNCH(speciespoints[3:4], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("geo_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[5:6], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("prec_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[7:8], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("growth_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[9:10], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("rad_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[11:12], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("temp_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[13:14], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("wind_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[15:16], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("soil_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[17:18], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("moist_",speciesnames[i], sep=""))
    
    thisrange<-NNCH(speciespoints[19:20], k=length(speciespoints[,1]))
    NNCH.export.shapefile(thisrange, paste("prod_",speciesnames[i], sep=""))
    
    keptspecies <- c(keptspecies, speciesnames[i])
  }
  else{ #Not enough points, remove species from matrix so it is not used for further calculations
    print (paste("Removing",speciesnames[i]))
  }
}

geo <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(geo) <- keptspecies
rownames(geo) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("geo_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("geo_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    geo[i,j] <- FT
  }
}
write.csv(geo, file="geo_overlap.csv")

prec <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(prec) <- keptspecies
rownames(prec) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("prec_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("prec_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    prec[i,j] <- FT
  }
}
write.csv(prec, file="prec_overlap.csv")

growth <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(growth) <- keptspecies
rownames(growth) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("growth_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("growth_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    growth[i,j] <- FT
  }
}
write.csv(growth, file="growth_overlap.csv")

rad <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(rad) <- keptspecies
rownames(rad) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("rad_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("rad_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    rad[i,j] <- FT
  }
}
write.csv(rad, file="rad_overlap.csv")

temp <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(temp) <- keptspecies
rownames(temp) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("temp_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("temp_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    temp[i,j] <- FT
  }
}
write.csv(temp, file="temp_overlap.csv")

wind <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(wind) <- keptspecies
rownames(wind) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("wind_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("wind_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    wind[i,j] <- FT
  }
}
write.csv(wind, file="wind_overlap.csv")

soil <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(soil) <- keptspecies
rownames(soil) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("soil_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("soil_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    soil[i,j] <- FT
  }
}
write.csv(soil, file="soil_overlap.csv")

moist <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(moist) <- keptspecies
rownames(moist) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("moist_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("moist_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    moist[i,j] <- FT
  }
}
write.csv(moist, file="moist_overlap.csv")


prod <- matrix(nrow = length(keptspecies), ncol = length(keptspecies))
colnames(prod) <- keptspecies
rownames(prod) <- keptspecies
for(i in 1:length(keptspecies)){
  shapename1 <- paste("prod_", keptspecies[i], ".shp", sep = "")
  thisshape1 <- readShapeSpatial(shapename1)
  #print (paste("Reading" ,shapename1))
  range1 <- as(thisshape1[1,], "gpc.poly")
  for(j in 1:length(thisshape1)){
    range1 <- union(range1, as(thisshape1[j,], "gpc.poly"))
  }
  #plot(range1, xlab = shapename1)
  for(j in i:length(keptspecies)){
    shapename2 <- paste("prod_", keptspecies[j], ".shp", sep = "")
    thisshape2 <- readShapeSpatial(shapename2)
    #print (paste("Reading" ,shapename2))
    range2 <- as(thisshape2[1,], "gpc.poly")
    for(k in 1:length(thisshape2)){
      range2 <- union(range2, as(thisshape2[k,], "gpc.poly"))
    }
    #range1 now contains range of species i, range2 contains species j
    FT <- area.poly(intersect(range1, range2))/min(area.poly(range1), area.poly(range2)) 
    pair <- paste(keptspecies[i], "vs", keptspecies[j])
    print(paste("Species", i, "and", j, ":", pair))
    #range3 <- intersect(range1,range2)
    plot(append.poly(range1, range2), xlab = pair)
    prod[i,j] <- FT
  }
}
write.csv(prod, file="prod_overlap.csv")