##### Simulation: Three part incidental truncation (dependence between MVP) model ########
############## Estimation assuming endogenous access #############
# Andrés Ramírez Hassan
# Febraury 1st, 2023
rm(list = ls())
set.seed(010101)
library(doParallel)
library(snow)
library(MCMCpack)
N <- 2500 # 20000
J <- 1
h1 <- 3 # Dim(Tj), see below
H <- h1 # Dim(T), see below
k1 <- 3 # Dim(Bj), See below
K <- k1 # Dim(B), See below
l1 <- 3 #  Dim(Dj), See below
L <- l1 # Dim(D), See below

##### Multivariate probit: Access #####
a1 <- c(1, -1, 1)

rho <- 1
SIGMA <- rho*matrix(c(1,0.7,0.6,0.7,1, 0.8, 0.6, 0.8, 1), 3*J, 3*J)
SIGMA <- rho*matrix(c(1,0.0,0.0,0.0,1, 0.8, 0.0, 0.8, 1), 3*J, 3*J)
# isSymmetric.matrix(SIGMA)
# matrixcalc::is.positive.definite(SIGMA)
##### Multivariate probit: Selection #####
# b1 <- c(1, -0.5, 0.5) # This set of parameters does not achieve good identification in the selection equation without exclusion restrictions. 
# This is because a12 and a13 have a linear dependence to b12 and b13
b1 <- c(1, -0.5, -1.2) # This set of parameters does not achieve good identification in the selection equation without exclusion restrictions. 

# Groups: 3^J
Comb <- matrix(c(0, 0, 1, 0, 1, 1), byrow = TRUE, 3, 2) # Combination access/use: First column is access and second is use 
CombNew <- Comb

##### SUR: Outcome #####
d1 <- c(1, -0.5, 1)

THETA <- c(a1, b1, d1)
#### Data preparation ####
IJ <- diag(J)
CJ <- matrix(0, J, J)
ZJ <- matrix(0, J, L)
XJ <- matrix(0, J, K)

#### Hyperparameters ####
b0 <- rep(0, K)
d0 <- rep(0, L)
t0 <- c(b0, d0)
T0 <- 1000*diag(K+L)
T0i <- solve(T0)
r0 <- 2*J + 2
R0 <- diag(K+L)

#### Gibbs sampler: Functions ####
PostTheta <- function(Sigma, Cl, Y, ZX){
  #Z : regressors use equation
  #X : regressors quantity equation
  #Y : quantity
  #A : Binary access
  #C : Binary use
  #Sigma: Covariance matrix (parameter)
  #Cl: Latent use (parameter)
  # 
  # Sigma = SIGMA[2:3,2:3]
  # Cl = Clp
  # A = A
  # C = C
  # Y = Y
  # ZX = ZX
  WY <- cbind(t(Cl), Y)
  XtX <- matrix(0, K+L, K+L)
  Xtwy <- matrix(0, K+L, 1)
  ids <- cumsum(c(1,c(K,L)))
  for(m in 2:3^J){
    idGood <- which(CombNew[m,]==1)
    columns <- sum(c(H,K,L)*c(rep(1,length(idGood)),rep(0,(3-length(idGood)))))
    JJ <- matrix(0, K+L, columns)
    idcov <- NULL
    for(l in 1:length(idGood)){
      JJ[ids[idGood[l]]:(ids[idGood[l]+1]-1),ids[l]:(ids[l+1]-1)] <- diag(3)
      idcovl <- ids[idGood[l]]:(ids[idGood[l]+1]-1)
      idcov <- c(idcov, idcovl)
    }
    G <- Groups[[m]]
    for(i in G){
      ZXm <- ZX[[i]][idGood,]
      if(length(idGood) == 1){
        ZXmi <- ZXm[idcov]
        XtXi <- JJ%*%ZXmi%*%solve(Sigma[idGood,idGood])%*%ZXmi%*%t(JJ)
        XtX <- XtX + XtXi
        Xtwyi <- JJ%*%ZXmi%*%solve(Sigma[idGood,idGood])%*%WY[i,idGood]
        Xtwy <- Xtwy + Xtwyi
      }else{
        ZXmi <- ZXm[,idcov]
        XtXi <- JJ%*%t(ZXmi)%*%solve(Sigma[idGood,idGood])%*%ZXmi%*%t(JJ)
        XtX <- XtX + XtXi
        Xtwyi <- JJ%*%t(ZXmi)%*%solve(Sigma[idGood,idGood])%*%WY[i,idGood]
        Xtwy <- Xtwy + Xtwyi
      }
    }
  }
  Tn <- solve(XtX + T0i)
  tn <- Tn%*%(Xtwy + T0i%*%t0)
  Tpost <- MASS::mvrnorm(1, tn, Tn)
  return(Tpost)
}

# THETAPost <- PostTheta(Sigma = SIGMA, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
# cbind(THETAPost,THETA)

PostSig <- function(theta, Cl, Y, C, A, ZZ, ZX){
  # theta = THETAp; Cl = Clp; Y = Y; C = C; A = A; ZZ = ZZ; ZX = ZX
  CY <- cbind(t(Cl), Y)
  eqs <- J
  WtW <- matrix(0, eqs, eqs)
  G0i <- which(A[,1] == 1)
  for(i in G0i){
    WtWi <- (CY[i, 1:eqs] - sum(ZX[[i]][1:eqs,]*theta))^2
    WtW <- WtW + WtWi
  }
  sigma11r <- LaplacesDemon::rinvwishart(r0-2*J+N, WtW + diag(eqs))
  
  eqs <- 2
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(C[,1] == 1)
  for(i in G1i){
    WtWgi <- (CY[i, 1:eqs] - ZX[[i]][1:eqs,]%*%theta)%*%t(CY[i, 1:eqs] - ZX[[i]][1:eqs,]%*%theta)
    WtWg1 <- WtWg1 + WtWgi
  }
  r11 <- WtWg1 + diag(eqs)
  r11i <- solve(r11[1:(eqs-1),1:(eqs-1)])
  r21 <- r11[eqs,1:(eqs-1)]
  r22 <- r11[eqs,eqs]
  r22.1 <- r22-t(r21)%*%r11i%*%r21
  sigma22.1r <- LaplacesDemon::rinvwishart(r0+length(G1i), r22.1)
  mr <- t(r21)%*%r11i
  sigma21.1r <- LaplacesDemon::rmatrixnorm(t(mr), as.matrix(Matrix::forceSymmetric(r11i)), sigma22.1r)
  sigma21r <- t(sigma21.1r)%*%sigma11r
  sigma22r <- sigma22.1r + sigma21r%*%sigma21.1r
  sigma11_1rns <- cbind(rbind(sigma11r, sigma21r), c(sigma21r, sigma22r))
  sigma11_1r <- sigma11_1rns
  
  return(sigma11_1r)
  }

# Sig11_1 <- PostSig(Theta = THETA, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
# cbind(matrixcalc::vech(Sig11_1), matrixcalc::vech(SIGMA))

PostCl <- function(theta, Sigma, Ai, Ci, Yi, Zi, Xi, Cli){
#  i <- 1
#  Ai <- A[i]; Ci <- C[i]; Yi <- Y[i]; Zi <- Z[i,]; Xi <- X[i,]; Cli <- Cl[i]
  if(Ai == 0){
    Cli <- Cli
  }else{
    if(Ci == 0){
      lb <- -Inf
      ub <- 0
      mij <- Zi%*%theta[1:3]
      sdij <- Sigma[1,1]
    }else{
      lb <- 0 
      ub <- Inf
      mij <- Zi%*%theta[1:3] + Sigma[2,1]/Sigma[2,2]*(Yi-Xi%*%theta[4:6])
      sdij <- Sigma[1,1]-Sigma[2,1]^2/Sigma[2,2]
    }

    Cli <- EnvStats::rnormTrunc(1, mean = mij, sd = sdij^0.5, min = lb, max = ub)
    if(Cli <= -20){ # Computational issues
      Cli <- -20
    }
    if(Cli >= 20){ # Computational issues
      Cli <- 20
    }
  }
  
  return(Cli)
}

# i <- 13
# m = Gs[i]; theta = THETA; Sigma = SIGMA; ACli = c(Al[i,],Cl[i,]);
# Yi = Y[i,]; Ai = A[i,]; Ci = C[i,]; WZXi = WZX[[i]]
# ACliPost <- PostACl(m = 1, theta = THETA, Sigma = SIGMA, ACli = c(Al[i,],Cl[i,]), Yi = Y[i,], Ai = A[i,], Ci = C[i,], WZXi = WZX[[i]])
# cbind(ACli, ACliPost)

#### Gibbs sampler: Implementation ####
S <- 5000 # 1100
thin <- 10
burnin <- 1000 + thin
ThetaPost <- matrix(NA, S, K + L)
SigmaPost <- array(NA, c(2*J, 2*J, S))
ThetaPostNOst <- matrix(NA, S, K + L)
SigmaPostNOst <- array(NA, c(2*J, 2*J, S))
AClPost <- array(NA, c(N, J, S))
findraws <- seq(burnin, S, thin)

Rep <- 100 

ALPHApost <- array(0, c(H, 3, Rep))

THETApost <- array(0, c(K+L, 3, Rep))
SIGMApost <- array(0, c(3, 3, Rep))
THETApostNOst <- array(0, c(K+L, 3, Rep))
SIGMApostNOst <- array(0, c(3, 3, Rep))
MultiVarMod <- vector(mode='list', length=Rep)

rep <- 1

cn <- detectCores() # 6
cl <- makeCluster(cn, type = "SOCK")
registerDoParallel(cl)

wzx <- rnorm(N)
w1 <- rnorm(N); z1 <- rnorm(N); x1 <- rnorm(N)
W <- cbind(1, w1, wzx)
Z <- cbind(1, z1, wzx)
X <- cbind(1, x1, wzx)

#### Data preparation ####
WW <- lapply(1:N, function(i){kronecker(IJ, t(W[i, ]))}) 
ZZ <- lapply(1:N, function(i){cbind(kronecker(IJ, t(Z[i, ])), ZJ)})
XX <- lapply(1:N, function(i){cbind(XJ, kronecker(IJ, t(X[i, ])))})
ZX <- lapply(1:N, function(i){rbind(ZZ[[i]], XX[[i]])})

while(rep <= Rep){
  tick <- Sys.time()
  
  U <- MASS::mvrnorm(N, mu = rep(0, 3*J), Sigma = SIGMA)
  Al1 <- W%*%a1 + U[, 1]
  Al <- Al1
  A1 <- Al1 > 0
  A <- A1
  
  Cl1 <- Z%*%b1 + U[, 2]
  Cl <- Cl1
  C1o <- Cl1 > 0
  C1 <- C1o
  C1[which(A1==0)] <- 0
  C <- C1
  AC <- cbind(A1, C1)

  Groups <- list()
  Gs <- rep(0, N)
  for(g in 1:dim(Comb)[1]){
    gg <- NULL
    for(i in 1:N){
      ggi <- sum(Comb[g,]==AC[i,]) == 2*J
      gg <- c(gg, ggi) 
      if(ggi == TRUE){
        Gs[i] <- g
      }
    }
    Groups[[g]] <- which(gg == TRUE)
  }
  
  Yl <- X%*%d1 + U[, 3]
  Y1 <- Yl
  Y1[which(C1==0)] <- NA
  Y <- Y1
  CY <- cbind(C, Y)
  ClY <- cbind(Cl, Y)
  # plot(wzx, Yl)
  # abline(lm(Yl~wzx), col = "blue")
  # points(wzx, Y1, col = "red")
  # abline(lm(Y1~wzx), col = "green")
  # 
  # la <- dnorm(W%*%a1)/pnorm(W%*%a1)
  # s2la <- 1 - SIGMA[2,1]^2*la*(W%*%a1+la)
  # cov(Z[,3]/s2la^0.5,la/s2la^0.5)*SIGMA[2,1]
  # 
  # lc <- dnorm(Z%*%b1)/pnorm(Z%*%b1)
  # cov(X[,3],lc)*SIGMA[3,2]
  # 
  RegA1 <- MCMCprobit(A1 ~ W - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  A1Hat <- summary(RegA1)
  ResA1 <- cbind(A1Hat$statistics[,1], A1Hat$quantiles[,1], A1Hat$quantiles[,5])
  # 
  # RegC1 <- MCMCprobit(C1[which(A1==1)] ~ Z[which(A1==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
  #                     b0 = 0, B0 = 1000^-1*diag(3))
  # C1Hat <- summary(RegC1)
  # ResC1 <- cbind(C1Hat$statistics[,1], C1Hat$quantiles[,1], C1Hat$quantiles[,5])
  # 
  # RegY1 <- MCMCregress(Y1[which(C1==1)] ~ X[which(C1==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
  #                      b0 = 0, B0 = 1000^-1*diag(3))
  # Y1Hat <- summary(RegY1)
  # ResY1 <- cbind(Y1Hat$statistics[-4,1], Y1Hat$quantiles[-4,1], Y1Hat$quantiles[-4,5])
  # 
  # ResUniVar <- rbind(ResA1, ResC1, ResY1)
  ALPHApost[,,rep] <- ResA1
  
  SIGMAp <- SIGMA[2:3,2:3]
  THETAp <- THETA[4:9]
  Clp <- Cl
  
  
  #### Parallel code ####
  clusterExport(cl, list("J", "K", "L", "k1", "PostCl", "ZZ", "C", "THETAp",
                         "SIGMAp", "H", "h1", "ZX", "A", "CombNew", "Gs", "N",
                         "Y", "Clp", "W", "Z", "X", "XX", "WW", "Cl"))
  
  for(s in 1:S){
    ticks <- Sys.time()
    clusterExport(cl, list("THETAp", "SIGMAp", "Clp"))
    idA <- which(A==1)
    Clp <- t(parSapply(cl, 1:N, function(i){PostCl(theta = THETAp, Sigma = SIGMAp, Ai = A[i], Ci = C[i], Yi = Y[i], Zi = Z[i], Xi = X[i], Cli = Cl[i])}))
    # Clp <- t(Cl)
    SIGMAp <- PostSig(theta = THETAp, Cl = Clp, Y = Y, C = C, A = A, ZZ = ZZ, ZX = ZX)
    SIGMApNew <- SIGMAp
    SIGMApNew[1,1] <- 1
    THETAp <- PostTheta(Sigma = SIGMAp, Cl = Clp, Y = Y, ZX = ZX)
    THETASp12 <- THETAp[1:H]/SIGMAp[1,1]^0.5

    ThetaPost[s,] <- c(THETASp12, THETAp[4:6])
    SigmaPost[,,s] <- SIGMApNew
    ThetaPostNOst[s,] <- THETAp
    SigmaPostNOst[,,s] <- SIGMAp 
    AClPost[,,s] <- Clp
    tocks <- Sys.time()
    print(tocks-ticks)
    print(s)
  }
  
  thetaHat <- coda::mcmc(ThetaPost[findraws,])
  RestTheta <- summary(thetaHat)
  THETApost[,,rep] <- cbind(RestTheta$statistics[,1], RestTheta$quantiles[,c(1,5)]) 
  
  thetaHatNOst <- coda::mcmc(ThetaPostNOst[findraws,])
  RestThetaNOst <- summary(thetaHatNOst)
  THETApostNOst[,,rep] <- cbind(RestThetaNOst$statistics[,1], RestThetaNOst$quantiles[,c(1,5)]) 
  
  SigmaHatVech <- coda::mcmc(t(sapply(1:S, function(s){matrixcalc::vech(SigmaPost[,,s])})))
  RestSigma <- summary(coda::mcmc(SigmaHatVech[findraws,]))
  SIGMApost[,,rep] <- cbind(RestSigma$statistics[,1], RestSigma$quantiles[,c(1,5)]) 
  
  SigmaHatVechNOst <- coda::mcmc(t(sapply(1:S, function(s){matrixcalc::vech(SigmaPostNOst[,,s])})))
  RestSigmaNOst <- summary(coda::mcmc(SigmaHatVechNOst[findraws,]))
  SIGMApostNOst[,,rep] <- cbind(RestSigmaNOst$statistics[,1], RestSigmaNOst$quantiles[,c(1,5)]) 
  
  PostResults <- list(ALPHApost = ALPHApost, THETApost = THETApost, SIGMApost = SIGMApost, THETApostNOst = THETApostNOst, SIGMApostNOst = SIGMApostNOst)
  
  save(PostResults, file = "PostResultsOneProduct3StagesExAccess.RData")
  
  MultiVarMod[[rep]] <- list(ThetaPost = ThetaPost[findraws,],
                             ThetaPostNOst = ThetaPostNOst[findraws,],
                             SigmaPost = t(sapply(findraws, function(s){matrixcalc::vech(SigmaPost[,,s])})),
                             SigmaPostNOst = t(sapply(findraws, function(s){matrixcalc::vech(SigmaPostNOst[,,s])})))
  save(MultiVarMod, file = "PostDrawsOneProduct3StagesExAccess.RData")
  
  print(rep)
  tock <- Sys.time()
  print(tock-tick)
  rep <- rep + 1
}

stopCluster(cl)

# load("PostResults.RData")
