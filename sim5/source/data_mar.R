#n=300
#P=100
cov <- list()
for(i in 1:P){
  cov[[i]] <- rnorm(n,0,1)
}
Xpred <- do.call(cbind, cov)
Xcut <- lapply(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t]))) # e.g. unique values of predictors

reg1 <- ifelse(Xpred[,1]< 0, 1, -1 )
reg2 <- ifelse(Xpred[,2]< 0, -1, 1 )

Y_trt <- rbinom(n, 1, pnorm(0.5+reg1+reg2-0.5*abs(Xpred[,3]-1)+1.5*Xpred[,4]*Xpred[,5])) # exposure model
Y_1 <- rnorm(n, 1*reg1+1.5*reg2-1+2*abs(Xpred[,3]+1) + 2*Xpred[,4]+ exp(0.5*Xpred[,5])-0.5*1*abs(Xpred[,6]) - 1*1*abs(Xpred[,7]+1) + 1*Xpred[,8]+ 1*Xpred[,9] + 1*Xpred[,10] + 0.5*Xpred[,11]+ 0.5*Xpred[,12] + 0.5*Xpred[,13] - 0.5*Xpred[,14] - 0.5*Xpred[,15] - 0.5*Xpred[,16] - exp(0.2*Xpred[,17]), 0.3)

Y_0 <- rnorm(n, 1*reg1+1.5*reg2-0+2*abs(Xpred[,3]+1) + 2*Xpred[,4]+ exp(0.5*Xpred[,5])-0.5*0*abs(Xpred[,6]) - 1*0*abs(Xpred[,7]+1) + 1*Xpred[,8]+ 1*Xpred[,9] + 1*Xpred[,10] + 0.5*Xpred[,11]+ 0.5*Xpred[,12] + 0.5*Xpred[,13] - 0.5*Xpred[,14] - 0.5*Xpred[,15] - 0.5*Xpred[,16] - exp(0.2*Xpred[,17]), 0.3)

Y_out <- Y_trt*Y_1 + (1-Y_trt)*Y_0

shift <- mean(Y_out)
Y_out <- Y_out - shift


Xpred1 <- cbind(Y_trt, Xpred)
Xcut1 <- lapply(1:dim(Xpred1)[2], function(t) sort(unique(Xpred1[,t]))) # e.g. unique values of predictors


Y1 <- Y_out
n1 <- length(Y1)

Z <- rnorm(n, qnorm(mean(Y_trt)), 1)    # latent variable N(0,1)


#------ Initial Setup (priors, initial values and hyper-parameters)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE
sigma2_1 <- var(Y_out)    # Initial value of SD^2
nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_1)
lambda_1 <- uniroot.all(f, c(0.1^5,10))
sigma_mu_1 <- max((min(Y_out)/(-2*sqrt(m)))^2, (max(Y_out)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

sigma2 <- 1               # sigma for the exposure model
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- uniroot.all(f, c(0.1^5,10))
sigma_mu <- max((min(Z)/(-2*sqrt(m)))^2, (max(Z)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

dir.alpha  <- dir.alpha1 <- 0.3         # Hyper-parameter on selection probabilities
alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2                 # default setting (alpha, beta) = (0.95, 2)


