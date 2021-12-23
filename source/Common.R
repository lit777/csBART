# Sampling mean parameters
Mean.Parameter <- function(sigma2, sigma_mu, dt,  Obs=Obs, R, ind=NULL){
  
  if(ind==1){
    xpred <- Xpred1; xcut <- Xcut1; nn <- n1
  }else{
    if(ind==0){
      xpred <- Xpred0; xcut <- Xcut0; nn <- n0
    }else{
      xpred <- Xpred; xcut <- Xcut; nn <- n
    }
  }
  
  T <- rep(0, nn)
  terminal <- which(dt$Terminal==1)
   for(i in 1:length(terminal)){
     Obs.ind <- Obs[(dt$begin[terminal[i]]):(dt$end[terminal[i]])]
     Var <- 1/(1/sigma_mu+length(Obs.ind)/sigma2)
     Mean <- Var * (sum(R[Obs.ind])/sigma2)
     T[Obs.ind] <- dt$MU[terminal[i]] <- rnorm(1,Mean, sqrt(Var))
  }
  return(list(T=T, dt=dt))
}
