for (k in 1:length(set.choice)){
  
  set.x <- set.choice[[k]][1,]; set.y <- set.choice[[k]][2,]
  
  index.y <- rep(F,length(train.label)); index.x <- rep(F,length(train.label))
  for(i in 1:length(train.label)){
    if (train.label[i] %in% set.y){
      index.y[i] <- T
    }
    if (train.label[i] %in% set.x){
      index.x[i] <- T
    }
  }
  
  data.X <- train.x[index.x,]
  data.Y <- train.x[index.y,]
  
  m.x <- dim(data.X)[1]
  n.y <- dim(data.Y)[1]
  
  
  pdf(paste0("Set",k,".pdf"),         # File name
      width = 6, height = 6, # Width and height in inches
      bg = "black")          # Background color
  
  par(mfrow = c(8,8),mai = c(0,0,0,0), oma = c(0,0,0,0))
  for(i in 1:32){
    plot(as.raster(matrix(data.X[sample(1:m.x,1),],28,28, byrow = T)))
  }
  
  for(i in 1:32){
    plot(as.raster(matrix(data.Y[sample(1:n.y,1),],28,28, byrow = T)))
  }
  
  # Closing the graphical device
  dev.off() 
  
  
}
