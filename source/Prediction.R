# Prediction (binary exposure)
Mean.Parameter_pred <- function(dt, ind=NULL){
  
  if(ind==1){
   xcut <- Xcut1;
  }else{
   xcut <- Xcut0;
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
  
  ind.list <- list()
  T <- rep(0, n*2)
  xpred.temp <- xpred.mult
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
      ind.list[[i]] <- xpred.temp[[(P+1)]]
      xpred.temp <- lapply(xpred.mult, function(x) x[-ind.list[[i]]])
      T[ind.list[[i]]] <- dt$MU[Terminal[i]]
  }
  }else{
    T <- rep(dt$MU, 2*n)
  }
  return(list(T=T))
}



