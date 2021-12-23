# Prediction (binary exposure)
Mean.Parameter_pred <- function(dt, ind=NULL){
  
  if(ind==1){
    xpred.temp <- xpred.Mult <- xpred.mult
    xpred.temp[[1]] <- xpred.Mult[[1]] <- rep(1, 2*n)
    xcut <- Xcut1
  } else {
    xpred.temp <- xpred.Mult <- xpred.mult
    xpred.temp[[1]] <- xpred.Mult[[1]] <- rep(0, 2*n)
    xcut <- Xcut1
  }
  
  Terminal <- which(dt$Terminal==1)
  Terminal.len <- length(Terminal)
  
  if(Terminal.len>1){
    Split.hist <- list()
    Value.hist <- list()
    Side.hist <- list()
    for(i in 1:Terminal.len){
      parent.node <- dt$parent[Terminal[i]]
      current.node <- dt$position[Terminal[i]]
      Split.hist[[i]] <- dt$Split[dt$position==parent.node]
      Value.hist[[i]] <- dt$Value[dt$position==parent.node]
      Side.hist[[i]] <- (current.node)%%2
      while(parent.node != 1){
        current.node <- dt$position[dt$position==parent.node]
        parent.node <- dt$parent[dt$position==parent.node] 
        Split.hist[[i]] <- c(Split.hist[[i]], dt$Split[dt$position==parent.node])
        Value.hist[[i]] <- c(Value.hist[[i]], dt$Value[dt$position==parent.node]) 
        Side.hist[[i]] <- c(Side.hist[[i]], (current.node)%%2)
      }
      Split.hist[[i]] <- rev(Split.hist[[i]])
      Value.hist[[i]] <- rev(Value.hist[[i]])
      Side.hist[[i]] <- rev(Side.hist[[i]])
    }
    
    ind.list <- NULL
    T <- rep(NA, n*2)
    for(i in 1:Terminal.len){
      count = 0
      while(count < length(Split.hist[[i]])){
        count <- count + 1
        if(Side.hist[[i]][count] == 0){
          sub.ind <- which(xpred.temp[[Split.hist[[i]][count]]] < xcut[[Split.hist[[i]][count]]][Value.hist[[i]][count]])
          xpred.temp <- lapply(xpred.temp, function(x) x[sub.ind])
        }else{
          sub.ind <- which(xpred.temp[[Split.hist[[i]][count]]] >= xcut[[Split.hist[[i]][count]]][Value.hist[[i]][count]])
          xpred.temp <- lapply(xpred.temp, function(x) x[sub.ind])
        }
      }
      if(length(xpred.temp[[P+2]])==0){
       xpred.temp <-xpred.Mult
      }else{
      ind.list <- c(ind.list, xpred.temp[[(P+2)]])
      T[xpred.temp[[(P+2)]]] <- dt$MU[Terminal[i]]
      xpred.temp <- lapply(xpred.Mult, function(x) x[-ind.list])
      }
    }
  }else{
    T <- rep(dt$MU[1], 2*n)
  }
  return(list(T=T))
}



# Prediction (continuous exposure)
Mean.Parameter_pred_cont <- function(dt, ind=NULL){
  
  xpred.temp <- xpred.Mult <- xpred.mult
  xpred.temp[[1]] <- xpred.Mult[[1]] <- rep(ind, 2*n)
  xcut <- Xcut1

  Terminal <- which(dt$Terminal==1)
  Terminal.len <- length(Terminal)
  
  if(Terminal.len>1){
    Split.hist <- list()
    Value.hist <- list()
    Side.hist <- list()
    for(i in 1:Terminal.len){
      parent.node <- dt$parent[Terminal[i]]
      current.node <- dt$position[Terminal[i]]
      Split.hist[[i]] <- dt$Split[dt$position==parent.node]
      Value.hist[[i]] <- dt$Value[dt$position==parent.node]
      Side.hist[[i]] <- (current.node)%%2
      while(parent.node != 1){
        current.node <- dt$position[dt$position==parent.node]
        parent.node <- dt$parent[dt$position==parent.node] 
        Split.hist[[i]] <- c(Split.hist[[i]], dt$Split[dt$position==parent.node])
        Value.hist[[i]] <- c(Value.hist[[i]], dt$Value[dt$position==parent.node]) 
        Side.hist[[i]] <- c(Side.hist[[i]], (current.node)%%2)
      }
      Split.hist[[i]] <- rev(Split.hist[[i]])
      Value.hist[[i]] <- rev(Value.hist[[i]])
      Side.hist[[i]] <- rev(Side.hist[[i]])
    }
    
    ind.list <- NULL
    T <- rep(NA, n*2)
    for(i in 1:Terminal.len){
      count = 0
      while(count < length(Split.hist[[i]])){
        count <- count + 1
        if(Side.hist[[i]][count] == 0){
          sub.ind <- which(xpred.temp[[Split.hist[[i]][count]]] < xcut[[Split.hist[[i]][count]]][Value.hist[[i]][count]])
          xpred.temp <- lapply(xpred.temp, function(x) x[sub.ind])
        }else{
          sub.ind <- which(xpred.temp[[Split.hist[[i]][count]]] >= xcut[[Split.hist[[i]][count]]][Value.hist[[i]][count]])
          xpred.temp <- lapply(xpred.temp, function(x) x[sub.ind])
        }
      }
      if(length(xpred.temp[[P+2]])==0){
        xpred.temp <-xpred.Mult
      }else{
        ind.list <- c(ind.list, xpred.temp[[(P+2)]])
        T[xpred.temp[[(P+2)]]] <- dt$MU[Terminal[i]]
        xpred.temp <- lapply(xpred.Mult, function(x) x[-ind.list])
      }
    }
  }else{
    T <- rep(dt$MU[1], 2*n)
  }
  return(list(T=T))
}


