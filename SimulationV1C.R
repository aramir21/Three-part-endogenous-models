##### Simulation: Three part incidental truncation (dependence between MVP) model ########
############## Estimation assuming endogenous access #############
# Andrés Ramírez Hassan
# Febraury 1st, 2023

##### Multivariate probit: Selection #####
rm(list = ls())
set.seed(010101)
library(doParallel)
library(snow)
N <- 2500 # 20000
J <- 1
Jdgp <- 3
h1 <- 2 # Dim(Tj), see below
H <- h1  # Dim(T), see below
k1 <- 3 # Dim(Bj), See below
K <- k1  # Dim(B), See below
l1 <- 3 # Dim(Dj), See below
L <- l1  # Dim(D), See below

##### Multivariate probit: Access #####
a1 <- c(1, -1)
a2 <- c(0.8, -1.2)
a3 <- c(1.1, -0.7)
rho <- 1
SIGMAdgp <- rho*matrix(c(1,0.6,0.4,0.5,0.4,0.2,0,0,0,0.6,1,0.5,0.4,0.3,0.4,0,0,0,0.4,0.5,1,
                      0.5,0.3,0.5,0,0,0,0.5,0.4,0.5,1,0.4,0.5,0.3,0.4,0.2,0.4,0.3,0.3,0.4,1,
                      0.4,0.3,0.3,0.1,0.2,0.4,0.5,0.5,0.4,1,0.5,0.3,0.1,0,0,0,0.3,0.3,0.5,1,0.6,
                      0.4,0,0,0,0.4,0.3,0.3,0.6,1,0.3,0,0,0,0.2,0.1,0.1,0.4,0.3,1), 3*Jdgp, 3*Jdgp)


SIGMA <- rho*matrix(c(1,0.5,0.0,0.5, 1, 0.3, 0.0, 0.3, 1), 3*J, 3*J)
# isSymmetric.matrix(SIGMA)
# matrixcalc::is.positive.definite(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)
b2 <- c(0.5, 1.5, -1)
b3 <- c(1, 1, -1)

# Groups: 3^J
Comb <- matrix(c(0, 0, 1, 0, 1, 1), byrow = TRUE, 3, 2) # Combination access/use: First column is access and second is use 
CombNew <- Comb

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)
d2 <- c(1, 2, 2)
d3 <- c(1, 1.5, 1.8)

THETA <- c(a1, b1, d1)
#### Data preparation ####
IJ <- diag(J)
CJ <- matrix(0, J, J)
WJ <- matrix(0, J, K+L)
ZJ1 <- matrix(0, J, H)
ZJ2 <- matrix(0, J, L)
XJ <- matrix(0, J, H+K)

#### Hyperparameters ####
a0 <- rep(0, H)
b0 <- rep(0, K)
d0 <- rep(0, L)
t0 <- c(a0, b0, d0)
T0 <- 1000*diag(H+K+L)
T0i <- solve(T0)
r0 <- J + 2
R0 <- diag(H+K+L)
R0A <- diag(J)

#### Gibbs sampler: Functions ####
PostTheta <- function(Sigma, Al, Cl, A, C, Y, WZX){
  # Sigma = SIGMA
  # Al = Al
  # Cl = Cl
  # A = A
  # C = C
  # Y = Y
  # WZX = WZX
  
  WY <- cbind(Al, Cl, Y)
  XtX <- matrix(0, H+K+L, H+K+L)
  Xtwy <- matrix(0, H+K+L, 1)
  ids <- cumsum(c(1,c(h1,k1,l1)))
  for(m in 1:3^J){
    idGood <- c(1,which(CombNew[m,]==1)+1)
    JJ <- matrix(0, H+K+L, H+(length(idGood)-1)*K) # This assumes K = L. Be carefull in application
    idcov <- NULL
    for(l in 1:length(idGood)){
      idscol <- ids[l]:(ids[l+1]-1)
      JJ[ids[idGood[l]]:(ids[idGood[l]+1]-1),idscol] <- diag(length(idscol))
      idcovl <- ids[idGood[l]]:(ids[idGood[l]+1]-1)
      idcov <- c(idcov, idcovl)
    }
    G <- Groups[[m]]
    for(i in G){
      ZXm <- WZX[[i]][idGood,]
      if(m==1){
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

PostSig <- function(Theta, Al, Cl, A, C, Y, WZX){
  # Theta = THETA
  # Al = Al
  # Cl = Cl
  # A = A
  # C = C
  # Y = Y
  # WZX = WZX
  
  ACY <- cbind(Al, Cl, Y)
  eqs <- J
  WtW <- matrix(0, eqs, eqs)
  for(i in 1:N){
    WtWi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
    WtW <- WtW + WtWi
  }
  sigma11r <- LaplacesDemon::rinvwishart(r0-2*J+N, WtW + diag(eqs))
  
  eqs <- 2
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(A[,1] == 1)
  for(i in G1i){
    WtWgi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
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
  
  eqs <- 3
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(A[,1] == 1 & C[,1] == 1)
  for(i in G1i){
    WtWgi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
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
  sigma21r <- t(sigma21.1r)%*%sigma11_1r
  sigma22r <- sigma22.1r + sigma21r%*%sigma21.1r
  sigma11_1rns <- cbind(rbind(sigma11_1r, sigma21r), c(sigma21r, sigma22r))
  sigma11_1r <- sigma11_1rns 
  return(sigma11_1r)
}

# Sig11_1 <- PostSig(Theta = THETA, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
# cbind(matrixcalc::vech(Sig11_1), matrixcalc::vech(SIGMA))

PostACl <- function(m, theta, Sigma, ACli, Yi, Ai, Ci, WZXi){
  # i <- 9
  # m = Gs[i]; theta = THETA; Sigma = SIGMA; ACli = c(Al[i,],Cl[i,]);
  # Yi = Y[i,]; Ai = A[i,]; Ci = C[i,]; WZXi = WZX[[i]]
  
  if(m == 1){
    idGood1 <- 1
    idGoodn <- c(2, 3) 
  }else{
    if(m == 2){
      idGood1 <- c(1, 2)
      idGoodn <- 3
    }else{
      idGood1 <- c(1, 2)
      idGoodn <- NULL 
    }
  }
  
  ACYli <- c(ACli,Yi)
  ACi <- c(Ai, Ci)
  if(length(idGood1) == 1){
    for(j in idGood1){
      if(ACi[j] == 0){
        lb <- -Inf
        ub <- 0
      }else{
        lb <- 0 
        ub <- Inf
      }
      mij <- WZXi[j,]%*%theta
      sdij <- Sigma[j,j]
      ACli[j] <- EnvStats::rnormTrunc(1, mean = mij, sd = sdij^0.5, min = lb, max = ub)
      if(ACli[j] <= -20){ # Computational issues
        ACli[j] <- -20
      }
      if(ACli[j] >= 20){ # Computational issues
        ACli[j] <- 20
      }
    }
  }else{
    for(j in idGood1){
      if(ACi[j] == 0){
        lb <- -Inf
        ub <- 0
      }else{
        lb <- 0 
        ub <- Inf
      }
      mij <- WZXi[j,]%*%theta + Sigma[j,-c(j,idGoodn)]%*%solve(Sigma[-c(j,idGoodn),-c(j,idGoodn)])%*%(ACYli[-c(j,idGoodn)] - WZXi[-c(j,idGoodn),]%*%theta)
      sdij <- Sigma[j,j]-Sigma[j,-c(j,idGoodn)]%*%solve(Sigma[-c(j,idGoodn),-c(j,idGoodn)])%*%(Sigma[-c(j,idGoodn),j])
      ACli[j] <- EnvStats::rnormTrunc(1, mean = mij, sd = sdij^0.5, min = lb, max = ub)
      if(ACli[j] <= -20){ # Computational issues
        ACli[j] <- -20
      }
      if(ACli[j] >= 20){ # Computational issues
        ACli[j] <- 20
      }
    }
  }
  return(ACli)
}
# i <- 1
# m = Gs[i]; theta = THETA; Sigma = SIGMA; ACli = c(Al[i,],Cl[i,]);
# Yi = Y[i,]; Ai = A[i,]; Ci = C[i,]; WZXi = WZX[[i]]
# ACliPost <- PostACl(m = 1, theta = THETA, Sigma = SIGMA, ACli = c(Al[i,],Cl[i,]), Yi = Y[i,], Ai = A[i,], Ci = C[i,], WZXi = WZX[[i]])
# cbind(ACli, ACliPost)

#### Gibbs sampler: Implementation ####
S <- 5000
thin <- 10
burnin <- 1000
ThetaPost <- matrix(NA, S, H + K + L)
SigmaPost <- array(NA, c(3*J, 3*J, S))
ThetaPostNOst <- matrix(NA, S, H + K + L)
SigmaPostNOst <- array(NA, c(3*J, 3*J, S))
AClPost <- array(NA, c(N, 2*J, S))
findraws <- seq(burnin, S, thin)

Rep <- 100 

THETApost <- array(0, c(H+K+L, 3, Rep))
SIGMApost <- array(0, c(3*J*(3*J+1)/2, 3, Rep))
THETApostNOst <- array(0, c(H+K+L, 3, Rep))
SIGMApostNOst <- array(0, c(3*J*(3*J+1)/2, 3, Rep))

rep <- 1

cn <- detectCores() 
cl <- makeCluster(cn, type = "SOCK")
registerDoParallel(cl)

W <- cbind(1, rnorm(N))
# Z <- cbind(W, rnorm(N))
Z <- cbind(1, rnorm(N), rnorm(N))
# X <- cbind(W, rnorm(N, 0, 1))
X <- cbind(1, rnorm(N), rnorm(N, 0, 1))

WW <- lapply(1:N, function(i){cbind(kronecker(IJ, t(W[i, ])), WJ)})
ZZ <- lapply(1:N, function(i){cbind(ZJ1, kronecker(IJ, t(Z[i, ])), ZJ2)})
XX <- lapply(1:N, function(i){cbind(XJ, kronecker(IJ, t(X[i, ])))})
WZ <- lapply(1:N, function(i){rbind(WW[[i]], ZZ[[i]])})
WZX <- lapply(1:N, function(i){rbind(WW[[i]], ZZ[[i]], XX[[i]])})

while(rep <= Rep){
  tick <- Sys.time()
  U <- MASS::mvrnorm(N, mu = rep(0, 3*Jdgp), Sigma = SIGMAdgp)
  Al1 <- W%*%a1 + U[, 1]
  Al2 <- W%*%a2 + U[, 2]
  Al3 <- W%*%a3 + U[, 3]
  Al <- Al1
  A1 <- Al1 > 0
  A2 <- Al2 > 0
  A3 <- Al3 > 0
  A <- A1
  
  Cl1 <- Z%*%b1 + U[, 4]
  Cl2 <- Z%*%b2 + U[, 5]
  Cl3 <- Z%*%b3 + U[, 6]
  Cl <- Cl1
  C1 <- Cl1 > 0
  C2 <- Cl2 > 0
  C3 <- Cl3 > 0
  C1[which(A1==0)] <- 0
  C2[which(A2==0)] <- 0
  C3[which(A3==0)] <- 0
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
  
  Y1 <- X%*%d1 + U[, 7]
  Y2 <- X%*%d2 + U[, 8]
  Y3 <- X%*%d3 + U[, 9]
  Y1[which(C1==0)] <- 0
  Y2[which(C2==0)] <- 0
  Y3[which(C3==0)] <- 0
  Y <- Y1
  CY <- cbind(C, Y)
  ClY <- cbind(Cl, Y)
  
  SIGMAp <- SIGMA
  THETAp <- THETA
  Clp <- Cl
  Alp <- Al
  AClp <- cbind(Alp, Clp)
  
  #### Parallel code ####
  clusterExport(cl, list("J", "K", "L", "k1", "PostACl", "ZZ", "C", "THETAp",
                         "SIGMAp", "H", "h1", "WZX", "A", "CombNew", "Gs", "N",
                         "Y", "AClp"))
  
  for(s in 1:S){
    # ticks <- Sys.time()
    clusterExport(cl, list("THETAp", "SIGMAp", "AClp"))
    AClp <- t(parSapply(cl, 1:N, function(i){PostACl(m = Gs[i], theta = THETAp, Sigma = SIGMAp, ACli = AClp[i,], Yi = Y[i,], WZXi = WZX[[i]], Ci = C[i,], Ai = A[i,])}))
    SIGMAp <- PostSig(Theta = THETAp, Al = AClp[,1], Cl = AClp[,2], A = A, C = C, Y = Y, WZX = WZX)
    SIGMASpProb <- diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))%*%SIGMAp[1:(2*J),1:(2*J)]%*%diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))
    SIGMApNew <- SIGMAp
    SIGMApNew[1:(2*J),1:(2*J)] <- SIGMASpProb 
    THETAp <- PostTheta(Sigma = SIGMAp, Al = AClp[,1], Cl = AClp[,2], A = A, C = C, Y = Y, WZX = WZX)
    THETASp12 <- THETAp[1:2]/SIGMAp[1,1]^0.5
    THETASp34 <- THETAp[3:5]/SIGMAp[2,2]^0.5
    
    
    ThetaPost[s,] <- c(THETASp12, THETASp34,THETAp[-c(1:5)])
    SigmaPost[,,s] <- SIGMApNew
    ThetaPostNOst[s,] <- THETAp
    SigmaPostNOst[,,s] <- SIGMAp
    AClPost[,,s] <- AClp
    # tocks <- Sys.time()
    # print(tocks-ticks)
    # print(s)
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
  
  PostResults <- list(THETApost = THETApost, SIGMApost = SIGMApost,
                      THETApostNOst = THETApostNOst, SIGMApostNOst = SIGMApostNOst)
  
  save(PostResults, file = "PostResultsV1CnewV1.RData")
  print(rep)
  tock <- Sys.time()
  print(tock-tick)
  rep <- rep + 1
}

stopCluster(cl)

