# Fun. of grow (first) alteration
GROW.first <- function(sigma2, sigma_mu, dt,  R, prop.prob, Obs, ind=NULL){
  
  if(ind==1){
    xpred <- Xpred1.list; xcut <- Xcut1; P <- P+1
  } else {
    if(ind==0){
      xpred <- Xpred0.list; xcut <- Xcut0
    } else {
      xpred <- Xpred.list; xcut <- Xcut; prop.prob <- prop.prob[-1]/sum(prop.prob[-1])
    }
  }
  
  prop.pred <- sample(1:P, 1, replace=FALSE,  prob = prop.prob) # pick a predictor
  if(ind == 1 & prop.pred == 1){
      prop.rule <- 2
      value <- 1
  }else{
      prop.rule <- sample(2:length(xcut[[prop.pred]]), 1)
      value <- xcut[[prop.pred]][prop.rule]
  }
  
  R.L <- which(xpred[[prop.pred]] < value) # Observations < prop.rule
  R.R <- setdiff(1:n, R.L)   # Observations >= prop.rule

  # Transition ratio (log scale)
  TRANS <- log(p.prune) - log(max(prop.prob[prop.pred],0)) + log(length(xcut[[prop.pred]])-1) - log(p.grow)
  
  # Likelihood ratio (log scale)
  nlL <- length(R.L)
  nlR <- length(R.R)
  
  LH <- log(sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))/sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))) + (sigma_mu / (2*sigma2) * ( (sum(R[R.L]))^2/(sigma2+nlL*sigma_mu) +(sum(R[R.R]))^2/(sigma2+nlR*sigma_mu) - (sum(R))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
  
  # Structure ratio (log scale)
  d <- 0
  STR <- log(alpha) + 2*log((1- alpha / (2 + d)^beta )) - log((1 + d)^beta - alpha) + log(max(prop.prob[prop.pred],0)) - log(length(xcut[[prop.pred]])-1)
  
  r <- TRANS+LH+STR
  
  if(r > log(runif(1))){
    # New tree structure
    dt.new <- dt
    dt.new$Split <- prop.pred
    dt.new$Terminal <- 0
    dt.new$Value <- prop.rule
    dt.add <- list(position=c(2,3), parent=c(1,1),  Terminal=rep(1,2), Split=rep(NA,2), Value=rep(NA,2), MU=rep(NA,2), begin=c(1,length(R.L)+1), end=c(length(R.L), n))
    dt.new <- mapply(c, dt.new, dt.add, SIMPLIFY = F)
    
    Obs.new <- Obs
    Obs.new[1:length(R.L)] <- sort(R.L)
    Obs.new[(length(R.L)+1):n] <- sort(R.R)
    return(list(dt=dt.new, Obs=Obs.new))
  }
  return(list(dt=dt, Obs=Obs))
}


# Fun. of grow alteration
GROW <- function(sigma2, sigma_mu, dt,  R, prop.prob, Obs, ind=NULL){
  
  if(ind==1){
    xpred <- Xpred1.list; xcut <- Xcut1; P <- P+1
  } else {
    if(ind==0){
      xpred <- Xpred0.list; xcut <- Xcut0
    } else {
      xpred <- Xpred.list; xcut <- Xcut; prop.prob <- prop.prob[-1]/sum(prop.prob[-1])
    }
  }

  prop.tnode <- sample(which(dt$Terminal==1), 1, replace=FALSE) # pick a terminal node
  
  begin <- dt$begin[prop.tnode]
  end <- dt$end[prop.tnode]
  Obs.ind <- Obs[begin:end]
  if(length(Obs.ind) < 2){
    return(list(dt=dt, Obs=Obs)) # return the current tree if there is no covariate with enough unique values
  }
  enough.unique <- which(mapply(function(x) length(unique(xpred[[x]][Obs[begin:end]])), 1:P)>=2)

  prop.pred <- sample(enough.unique, 1, replace=FALSE,  prob = prop.prob[enough.unique]) # pick a predictor
  
  xpred.prop.pred <- xpred[[prop.pred]][Obs.ind]
  xcut.prop.pred <- xcut[[prop.pred]]
  unique.len <- length(unique(xpred.prop.pred)) # Num. of unique values
  if(unique.len == 1){
    return(list(dt=dt, Obs=Obs)) # return the current tree if there is no covariate with enough unique values
  }
  prop.rule <- sort(unique(xpred.prop.pred))[sample(1:(unique.len-1), 1)+1]
  prop.rule <- which(xcut.prop.pred==prop.rule)
  
  R.L <- Obs.ind[which(xpred.prop.pred < xcut.prop.pred[prop.rule])] # Observations < prop.rule
  R.R <- setdiff(Obs.ind, R.L)   # Observations >= prop.rule
  
  # New tree structure
  dt.new <- dt
  dt.new.position <- dt.new$position[prop.tnode]
  dt.new$Terminal[prop.tnode] <- 0

  dt.add <- list( position=c(2*dt.new.position,2*dt.new.position+1), parent=rep(dt.new.position,2),  Terminal=rep(1,2), Split=rep(NA,2), Value=rep(NA,2), MU=rep(NA,2), begin=c(begin,(begin+length(R.L))), end=c((begin+length(R.L)-1),(end)))
  dt.new <- mapply(c, dt.new, dt.add, SIMPLIFY = F)
  
  ### Prune step
  # find internal nodes with two terminal child nodes
  count <- length(which(table(dt.new$parent[which(dt.new$Terminal==1)])==2))

  # Transition ratio
  TRANS <- log(0.28) + log(length(which(dt$Terminal==1))) - log(max(prop.prob[prop.pred]/sum(prop.prob[enough.unique]),0)) + log((unique.len-1)) - log(0.28) - log(count)
  
  # Likelihood ratio
  nlL <- length(R.L)
  nlR <- length(R.R)
  
  LH <- log(sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))/sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))) + (sigma_mu / (2*sigma2) * ( (sum(R[R.L]))^2/(sigma2+nlL*sigma_mu) +(sum(R[R.R]))^2/(sigma2+nlR*sigma_mu) - (sum(R[union(R.R,R.L)]))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
  
  # Structure ratio
  d <- 1
  while(dt$position[prop.tnode] >= 2^d){
    d <- d+1
  }
  d < d-1
  STR <- log(alpha) + 2*log((1- alpha / (2 + d)^beta )) - log((1 + d)^beta - alpha) + log(max(prop.prob[prop.pred]/sum(prop.prob[enough.unique]),0)) - log((unique.len-1))
  
  r <- TRANS+LH+STR
  if(r > log(runif(1))){

    dt.new$Split[prop.tnode] <- prop.pred
    dt.new$Value[prop.tnode] <- prop.rule
    
    Obs.new <- Obs
    Obs.new[begin:(begin+length(R.L)-1)] <- sort(R.L)
    Obs.new[(begin+length(R.L)):(end)] <- sort(R.R)
    return(list(dt = dt.new, Obs = Obs.new))
  }
  return(list(dt = dt, Obs = Obs))
}

# Fun. of grow (continuous) alteration
GROW.first.cont <- function(sigma2, sigma_mu, dt,  R, prop.prob, Obs, ind=NULL){
  
  xpred <- Xpred1.list; xcut <- Xcut1; P <- P+1

  prop.pred <- sample(1:P, 1, replace=FALSE,  prob = prop.prob) # pick a predictor
    prop.rule <- sample(2:length(xcut[[prop.pred]]), 1)
    value <- xcut[[prop.pred]][prop.rule]
  R.L <- which(xpred[[prop.pred]] < value) # Observations < prop.rule
  R.R <- setdiff(1:n, R.L)   # Observations >= prop.rule
    
  # Transition ratio (log scale)
  TRANS <- log(p.prune) - log(max(prop.prob[prop.pred],0)) + log(length(xcut[[prop.pred]])-1) - log(p.grow)
  
  # Likelihood ratio (log scale)
  nlL <- length(R.L)
  nlR <- length(R.R)
  
  LH <- log(sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))/sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))) + (sigma_mu / (2*sigma2) * ( (sum(R[R.L]))^2/(sigma2+nlL*sigma_mu) +(sum(R[R.R]))^2/(sigma2+nlR*sigma_mu) - (sum(R))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
  
  # Structure ratio (log scale)
  d <- 0
  STR <- log(alpha) + 2*log((1- alpha / (2 + d)^beta )) - log((1 + d)^beta - alpha) + log(max(prop.prob[prop.pred],0)) - log(length(xcut[[prop.pred]])-1)
  
  r <- TRANS+LH+STR
  
  if(r > log(runif(1))){
    # New tree structure
    dt.new <- dt
    dt.new$Split <- prop.pred
    dt.new$Terminal <- 0
    dt.new$Value <- prop.rule
    dt.add <- list(position=c(2,3), parent=c(1,1),  Terminal=rep(1,2), Split=rep(NA,2), Value=rep(NA,2), MU=rep(NA,2), begin=c(1,length(R.L)+1), end=c(length(R.L), n))
    dt.new <- mapply(c, dt.new, dt.add, SIMPLIFY = F)
    
    Obs.new <- Obs
    Obs.new[1:length(R.L)] <- sort(R.L)
    Obs.new[(length(R.L)+1):n] <- sort(R.R)
    return(list(dt=dt.new, Obs=Obs.new))
  }
  return(list(dt=dt, Obs=Obs))
}

