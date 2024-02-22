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
J <- 3
h1 <- 2; h2 <- 2; h3 <- 2 # Dim(Tj), see below
H <- h1 + h2 + h3 # Dim(T), see below
k1 <- 3; k2 <- 3; k3 <- 3 # Dim(Bj), See below
K <- k1 + k2 + k3 # Dim(B), See below
l1 <- 3; l2 <- 3; l3 <- 3 # Dim(Dj), See below
L <- l1 + l2 + l3 # Dim(D), See below

##### Multivariate probit: Access #####
a1 <- c(1, -1)
a2 <- c(0.8, -1.2)
a3 <- c(1.1, -0.7)
rho <- 1
SIGMA <- rho*matrix(c(1,0.6,0.4,0.5,0.4,0.2,0,0,0,0.6,1,0.5,0.4,0.3,0.4,0,0,0,0.4,0.5,1,
                      0.5,0.3,0.5,0,0,0,0.5,0.4,0.5,1,0.4,0.5,0.3,0.4,0.2,0.4,0.3,0.3,0.4,1,
                      0.4,0.3,0.3,0.1,0.2,0.4,0.5,0.5,0.4,1,0.5,0.3,0.1,0,0,0,0.3,0.3,0.5,1,0.6,
                      0.4,0,0,0,0.4,0.3,0.3,0.6,1,0.3,0,0,0,0.2,0.1,0.1,0.4,0.3,1), 3*J, 3*J)
SIGMA13 <- matrix(c(0.3, 0.2, 0.3, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2), 3, 3) 
SIGMA[1:3,7:9] <- SIGMA13
SIGMA[7:9,1:3] <- t(SIGMA13)
# isSymmetric.matrix(SIGMA)
# matrixcalc::is.positive.definite(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)
b2 <- c(0.5, 1.5, -1)
b3 <- c(1, 1, -1)
# Groups: 3^J
CombId <- rbind(c(0, 0), c(1, 0), c(1, 1))
id1Comb <- rep(1:J, 1, each = J^(J-1))
id2Comb <- rep(1:J, J, each = J)
id3Comb <- rep(1:J, J^(J-1))
idsComb <- cbind(id1Comb, id2Comb, id3Comb)
Comb1 <- matrix(0, 3^J, 2)
Comb2 <- matrix(0, 3^J, 2)
Comb3 <- matrix(0, 3^J, 2)
for(id in 1:3^J){
  Comb1[id, ] <- CombId[id1Comb[id],]
  Comb2[id, ] <- CombId[id2Comb[id],]
  Comb3[id, ] <- CombId[id3Comb[id],]
}
Comb <- cbind(Comb1, Comb2, Comb3) # Combination access/use per drug: A1 C1 A2 C2 A3 C3
CombNew <- cbind(Comb[,c(1,3,5)],Comb[,c(2,4,6)]) # First three columns is access and the three last ones are use: A1 A2 A3 C1 C2 C3
# n <- NULL
# for(g in 1:27){
#   ni <- length(Groups[[g]])
#   n <- c(n, ni)
# }

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)
d2 <- c(1, 2, 2)
d3 <- c(1, 1.5, 1.8)

THETA <- c(a1, a2, a3, b1, b2, b3, d1, d2, d3)
#### Data preparation ####
IJ <- diag(J)
CJ1 <- matrix(0, J, 2*J)
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
r0 <- 3*J + 2
R0 <- diag(H+K+L)
R0A <- diag(J)

#### Gibbs sampler: Implementation ####
S <- 5000 # 1100
thin <- 10
burnin <- 1000 + thin
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

cn <- detectCores() # 6
cl <- makeCluster(cn, type = "SOCK")
registerDoParallel(cl)

W <- cbind(1, rnorm(N))
z1 <- rnorm(N) 
# Z <- cbind(W, z1)
Z <- cbind(1, rnorm(N), rnorm(N))
# X <- Z
X <- cbind(1, rnorm(N), rnorm(N, 0, 1))

WW <- lapply(1:N, function(i){cbind(kronecker(IJ, t(W[i, ])), WJ)})
ZZ <- lapply(1:N, function(i){cbind(ZJ1, kronecker(IJ, t(Z[i, ])), ZJ2)})
XX <- lapply(1:N, function(i){cbind(XJ, kronecker(IJ, t(X[i, ])))})
WZ <- lapply(1:N, function(i){rbind(WW[[i]], ZZ[[i]])})
WZX <- lapply(1:N, function(i){rbind(WW[[i]], ZZ[[i]], XX[[i]])})


while(rep <= Rep){
  tick <- Sys.time()
  
  U <- MASS::mvrnorm(N, mu = rep(0, 3*J), Sigma = SIGMA)
  Al1 <- W%*%a1 + U[, 1]
  Al2 <- W%*%a2 + U[, 2]
  Al3 <- W%*%a3 + U[, 3]
  IRMA1 <- dnorm(W%*%a1)/pnorm(W%*%a1)
  plot(IRMA1)
  cov(IRMA1,cbind(W[,2],Z[,3]))*SIGMA[4,1]*SIGMA[1,1]/0.84^0.5
  IRMC1 <- dnorm(Z%*%b1)/pnorm(Z%*%b1)
  plot(IRMC1)
  cov(IRMC1,cbind(Z[,2],X[,3]))*SIGMA[4,7]*SIGMA[4,4]
  Al <- cbind(Al1, Al2, Al3)
  A1 <- Al1 > 0
  A2 <- Al2 > 0
  A3 <- Al3 > 0
  A <- cbind(A1, A2, A3)
  
  Cl1 <- Z%*%b1 + U[, 4]
  Cl2 <- Z%*%b2 + U[, 5]
  Cl3 <- Z%*%b3 + U[, 6]
  
  Cl <- cbind(Cl1, Cl2, Cl3)
  C1 <- Cl1 > 0
  C2 <- Cl2 > 0
  C3 <- Cl3 > 0
  C1[which(A1==0)] <- 0
  C2[which(A2==0)] <- 0
  C3[which(A3==0)] <- 0
  C <- cbind(C1, C2, C3)
  AC <- cbind(A1, C1, A2, C2, A3, C3)
  
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
  Y <- cbind(Y1, Y2, Y3)
  CY <- cbind(C, Y)
  ClY <- cbind(Cl, Y)
  
  RegA1 <- MCMCprobit(A1 ~ W - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(2))
  RegA2 <- MCMCprobit(A2 ~ W - 1, mcmc = S - burnin, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(2))
  RegA3 <- MCMCprobit(A3 ~ W - 1, mcmc = S - burnin, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(2))
  A1Hat <- summary(RegA1)
  A2Hat <- summary(RegA2)
  A3Hat <- summary(RegA3)
  Ahat <- cbind(c(A1Hat$statistics[,1], A2Hat$statistics[,1], A3Hat$statistics[,1]),
                c(A1Hat$quantiles[,1], A2Hat$quantiles[,1], A3Hat$quantiles[,1]),
                c(A1Hat$quantiles[,5], A2Hat$quantiles[,5], A3Hat$quantiles[,5]))
  
  RegC1 <- MCMCprobit(C1[which(A1==1)] ~ Z[which(A1==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  RegC2 <- MCMCprobit(C2[which(A2==1)] ~ Z[which(A2==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  RegC3 <- MCMCprobit(C3[which(A3==1)] ~ Z[which(A3==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  C1Hat <- summary(RegC1)
  C2Hat <- summary(RegC2)
  C3Hat <- summary(RegC3)
  Chat <- cbind(c(C1Hat$statistics[,1], C2Hat$statistics[,1], C3Hat$statistics[,1]),
                c(C1Hat$quantiles[,1], C2Hat$quantiles[,1], C3Hat$quantiles[,1]),
                c(C1Hat$quantiles[,5], C2Hat$quantiles[,5], C3Hat$quantiles[,5]))
  
  RegY1 <- MCMCregress(Y1[which(C1==1)] ~ X[which(C1==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  RegY2 <- MCMCregress(Y2[which(C2==1)] ~ X[which(C2==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  RegY3 <- MCMCregress(Y3[which(C3==1)] ~ X[which(C3==1),] - 1, mcmc = S - burnin + 10, burnin = burnin, thin = thin,
                      b0 = 0, B0 = 1000^-1*diag(3))
  Y1Hat <- summary(RegY1)
  Y2Hat <- summary(RegY2)
  Y3Hat <- summary(RegY3)
  Yhat <- cbind(c(Y1Hat$statistics[-4,1], Y2Hat$statistics[-4,1], Y3Hat$statistics[-4,1]),
                c(Y1Hat$quantiles[-4,1], Y2Hat$quantiles[-4,1], Y3Hat$quantiles[-4,1]),
                c(Y1Hat$quantiles[-4,5], Y2Hat$quantiles[-4,5], Y3Hat$quantiles[-4,5]))
  
  THETApost[,,rep] <- 
  
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
    # SIGMAp <- PostSig(Theta = THETAp, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
    SIGMAp <- PostSig(Theta = THETAp, Al = AClp[,1:3], Cl = AClp[,4:6], A = A, C = C, Y = Y, WZX = WZX)
    SIGMASpProbV1 <- diag(1/(diag(SIGMAp)^0.5))%*%SIGMAp%*%diag(1/(diag(SIGMAp)^0.5))
    SIGMASpProb <- diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))%*%SIGMAp[1:(2*J),1:(2*J)]%*%diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))
    SIGMApNew <- SIGMAp
    SIGMApNew[1:(2*J),1:(2*J)] <- SIGMASpProb 
    # THETAp <- PostTheta(Sigma = SIGMAp, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
    THETAp <- PostTheta(Sigma = SIGMAp, Al = AClp[,1:3], Cl = AClp[,4:6], A = A, C = C, Y = Y, WZX = WZX)
    THETASp12 <- THETAp[1:2]/SIGMAp[1,1]^0.5
    THETASp34 <- THETAp[3:4]/SIGMAp[2,2]^0.5
    THETASp56 <- THETAp[5:6]/SIGMAp[3,3]^0.5
    THETASp79 <- THETAp[7:9]/SIGMAp[4,4]^0.5
    THETASp02 <- THETAp[10:12]/SIGMAp[5,5]^0.5
    THETASp35 <- THETAp[13:15]/SIGMAp[6,6]^0.5
    
    ThetaPost[s,] <- c(THETASp12, THETASp34, THETASp56, THETASp79, THETASp02, THETASp35,THETAp[-c(1:15)])
    SigmaPost[,,s] <- SIGMApNew
    ThetaPostNOst[s,] <- THETAp
    SigmaPostNOst[,,s] <- SIGMAp #SIGMASpProbV1
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
  
  PostResults <- list(THETApost = THETApost, SIGMApost = SIGMApost, THETApostNOst = THETApostNOst, SIGMApostNOst = SIGMApostNOst)
  
  save(PostResults, file = "PostResultsV3BnewV2.RData")
  print(rep)
  tock <- Sys.time()
  print(tock-tick)
  rep <- rep + 1
}

stopCluster(cl)

# load("PostResults.RData")
