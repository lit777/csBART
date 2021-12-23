# Fun. of PRUNE alteration
PRUNE <-  function(sigma2, sigma_mu, dt,  R, prop.prob,  Obs, ind=NULL){
    
  if(ind==1){
    xpred <- Xpred1.list; xcut <- Xcut1; n <- n1
  }else{
    if(ind==0){
      xpred <- Xpred0.list; xcut <- Xcut0; n <- n0
    }else{
      xpred <- Xpred.list; xcut <- Xcut; n <- n
    }
  }

    singly.position <- as.numeric(names(which(table(dt$parent[which(dt$Terminal==1)])==2)))
    singly.inode <- ifelse(length(singly.position)==1, singly.position, sample(singly.position, 1)) # pick a singly internal parent node

    subset.ind <- which(dt$position==singly.inode)
    prop.pred <- dt$Split[subset.ind]
    prop.rule <- dt$Value[subset.ind]
    begin <- dt$begin[subset.ind]
    end <- dt$end[subset.ind]
    subset.ind <- which(dt$position==2*singly.inode)
    begin.L <- dt$begin[subset.ind]
    end.L <- dt$end[subset.ind]
    temp.L <- Obs[begin.L:end.L]
    subset.ind <- which(dt$position==(2*singly.inode+1))
    begin.R <- dt$begin[subset.ind]
    end.R <- dt$end[subset.ind]
    temp.R <- Obs[begin.R:end.R]

    enough.unique <- which(mapply(function(x) length(unique(xpred[[x]][Obs[begin:end]])), 1:P)>=2)
    prop.prob <- prop.prob[prop.pred]/sum(prop.prob[enough.unique]) # P(selecting the jth attribute to split on)
    unique.len <- length(unique(xpred[[prop.pred]][Obs[begin:end]]))
    
    # transition
    transition_forward <- 0.28 * 1/ length(singly.position)
    transition_back <- 0.28 * 1 / (length(which(dt$Terminal==1)) - 1) * (prop.prob) * (1 / unique.len)
  
    # Transition ratio
    TRANS <- log(0.28) - (log(length(which(dt$Terminal==1)) - 1)) + log(max(prop.prob, 0)) - log(unique.len) - log(0.28) + log(length(singly.position))
  
    # Likelihood ratio
    nlL <- length(temp.L)
    nlR <- length(temp.R)
  
    LH <- log(sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))/sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))) + (sigma_mu / (2*sigma2) * ( -(sum(R[temp.L]))^2/(sigma2+nlL*sigma_mu) -(sum(R[temp.R]))^2/(sigma2+nlR*sigma_mu) + (sum(R[union(temp.R, temp.L)]))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
  
    # Structure ratio
    d <- 1
    while(singly.inode >= 2^d){
      d <- d+1
    }
    d <- d-1
#    d <- dt$level[which(dt$position==singly.inode)]-1
    STR <- -log(alpha) - 2*log((1- alpha / (2 + d)^beta )) + log((1 + d)^beta - alpha) - log(max(prop.prob, 0)) + log(unique.len)
  
    r <- TRANS+LH+STR
    if(r > log(runif(1))){
      dt.new <- dt
      subset.ind <- which(dt.new$parent==singly.inode)
      dt.new <- lapply(dt.new, function(x) return(x[-subset.ind]))
      subset.ind <- which(dt.new$position==singly.inode)
      dt.new$Split[subset.ind] <- NA
      dt.new$Value[subset.ind] <- NA
      dt.new$Terminal[subset.ind] <- 1
      Obs[begin:end] <- sort(Obs[begin:end]) 
      
      return(list(dt=dt.new, Obs=Obs))
    }
  return(list(dt=dt, Obs=Obs))
}
