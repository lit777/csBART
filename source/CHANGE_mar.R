# Fun. of CHANGE alteration  
CHANGE <- function(sigma2, sigma_mu, dt,  R, prop.prob, Obs, ind=NULL){
  
  if(ind==1){
    xpred <- Xpred1.list; xcut <- Xcut1; P <- P+1
  } else {
    if(ind==0){
      xpred <- Xpred0.list; xcut <- Xcut0;
    } else {
      xpred <- Xpred.list; xcut <- Xcut; prop.prob <- prop.prob[-1] / sum(prop.prob[-1])
    }
  }

  # find nodes with two terminal child nodes (singly internal parent nodes)
  singly.position <- as.numeric(names(which(table(dt$parent[which(dt$Terminal==1)])==2)))
  singly.inode <- ifelse(length(singly.position)==1, singly.position, sample(singly.position, 1)) # pick a singly internal parent node
  
  subset.ind <- which(dt$position==singly.inode)
  begin <- dt$begin[subset.ind]
  end <- dt$end[subset.ind]
  Obs.ind <- Obs[begin:end]
  enough.unique <- which(mapply(function(x) length(unique(xpred[[x]][Obs[begin:end]])), 1:P)>=2)
  prop.pred <- sample(enough.unique, 1, replace=FALSE,  prob = prop.prob[enough.unique])

  unique.len <- length(unique(xpred[[prop.pred]][Obs.ind])) # Num. of unique values
  prop.rule <- value <- sort(unique(xpred[[prop.pred]][Obs.ind]))[sample(1:(unique.len-1), 1)+1]

  RL.star <- Obs.ind[which(xpred[[prop.pred]][Obs.ind] < value)]
  RR.star <- setdiff(Obs.ind, RL.star)
  
  RL <- Obs[dt$begin[which(dt$position==2*singly.inode)]:dt$end[which(dt$position==2*singly.inode)]]
  RR <- Obs[dt$begin[which(dt$position==(2*singly.inode+1))]:dt$end[which(dt$position==(2*singly.inode+1))]]

  # Likelihood ratio
  nlL <- length(RL)
  nlR <- length(RR)
  nlL_star <- length(RL.star)
  nlR_star <- length(RR.star)
  
  LH <- log(sqrt((sigma2/sigma_mu+nlL)*(sigma2/sigma_mu+nlR))) - log(sqrt((sigma2/sigma_mu+nlL_star)*(sigma2/sigma_mu+nlR_star))) + (0.5 / sigma2 * ( (sum(R[RL.star]))^2/(nlL_star + sigma2/sigma_mu) + (sum(R[RR.star]))^2/(nlR_star + sigma2/sigma_mu) -  (sum(R[RL]))^2/(nlL + sigma2/sigma_mu) - (sum(R[RR]))^2/(nlR + sigma2/sigma_mu)))
  
  if(LH > log(runif(1))){
    # New tree structure
    dt.new <- dt
    subset.ind1 <- which(dt.new$position==singly.inode)
    dt.new$Terminal[subset.ind1] <- 0
    dt.new$Split[subset.ind1] <- prop.pred
    dt.new$Value[subset.ind1] <- which(xcut[[prop.pred]]==prop.rule)
    
    subset.ind2 <- which(dt.new$position==2*singly.inode)
    dt.new$begin[subset.ind2] <- begin 
    dt.new$end[subset.ind2] <- begin+length(RL.star)-1
    subset.ind3 <- which(dt.new$position==(2*singly.inode+1))
    dt.new$begin[subset.ind3] <- begin+length(RL.star) 
    dt.new$end[subset.ind3] <- end
    
    Obs.new <- Obs
    Obs.new[begin:(begin+length(RL.star)-1)] <- sort(RL.star)
    Obs.new[(begin+length(RL.star)):(end)] <- sort(RR.star)
    return(list(dt=dt.new, Obs=Obs.new))
  }  
  return(list(dt=dt, Obs=Obs))
}
