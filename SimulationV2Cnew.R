##### Simulation: Three part incidental truncation model ########
############## Estimation assuming exogenous access #############
# Andrés Ramírez Hassan
# Febraury 1st, 2023
rm(list = ls())
set.seed(010101)
library(doParallel)
library(snow)
N <- 2500 # 20000
J <- 3
h1 <- 3; h2 <- 3; h3 <- 3 # Dim(Tj), see below
H <- h1 + h2 + h3 # Dim(T), see below
k1 <- 3; k2 <- 3; k3 <- 3 # Dim(Bj), See below
K <- k1 + k2 + k3 # Dim(B), See below
l1 <- 3; l2 <- 3; l3 <- 3 # Dim(Dj), See below
L <- l1 + l2 + l3 # Dim(D), See below

##### Multivariate probit: Access #####
a1 <- c(1, -1, 1)
a2 <- c(0.8, -1.2, 0.5)
a3 <- c(1.1, -0.7, 0.8)
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
b1 <- c(1, -0.5, 0.5)
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
Comb <- cbind(Comb1, Comb2, Comb3) # Combination access/consumption
CombNew <- cbind(Comb[,c(1,3,5)],Comb[,c(2,4,6)])
# n <- NULL
# for(g in 1:27){
#   ni <- length(Groups[[g]])
#   n <- c(n, ni)
# }

##### SUR: Outcome #####
d1 <- c(1, -0.5, 1)
d2 <- c(1, 2, 2)
d3 <- c(1, 1.5, 1.8)

THETA <- c(a1, a2, a3, b1, b2, b3, d1, d2, d3)
#### Data preparation ####
IJ <- diag(J)
CJ <- matrix(0, J, J)
ZJ <- matrix(0, J, L)
XJ <- matrix(0, J, K)

#### Hyperparameters ####
a0 <- rep(0, H)
A0 <- 1000*diag(H)
A0i <- solve(A0)
om0 <- J + 2
Om0 <- diag(J)

b0 <- rep(0, K)
d0 <- rep(0, L)
t0 <- c(b0, d0)
T0 <- 1000*diag(K+L)
T0i <- solve(T0)
r0 <- 2*J + 2
R0 <- diag(K+L)
#### Posterior parameters #####
omn <- om0 + N

#### Gibbs sampler: Functions ####
PostAlpha <- function(Omega, Al, W){
  WtW <- matrix(0, H, H)
  Wta <- matrix(0, H, 1)
  for(i in 1:N){
    WtWi <- t(W[[i]])%*%solve(Omega)%*%W[[i]]
    WtW <- WtW + WtWi
    Wtai <- t(W[[i]])%*%solve(Omega)%*%Al[i,]
    Wta <- Wta + Wtai
  }
  Tn <- solve(WtW + A0i)
  tn <- Tn%*%(Wta + A0i%*%a0)
  tpost <- MASS::mvrnorm(1, mu = tn, Sigma = Tn)
  return(tpost)
}

# AlphaPost <- PostAlpha(Omega=OMEGA, Al = Al, W = WW)
# cbind(AlphaPost, c(a1, a2, a3))

PostOmega <- function(alpha, Al, W){
  WtW <- matrix(0, J, J)
  for(i in 1:N){
    WtWi <- (Al[i, ] - W[[i]]%*%alpha)%*%t(Al[i, ] - W[[i]]%*%alpha)
    WtW <- WtW + WtWi
  }
  Ompost <- LaplacesDemon::rinvwishart(omn, Matrix::nearPD((WtW + Om0))$mat)
  return(Ompost)
}
# OmegaPost <- PostOmega(alpha = c(a1, a2, a3), Al = Al, W = WW)

PostAl <- function(alpha, Omega, Ali, Ai, Wi){
  for(j in 1:J){
    if(Ai[j] == 0){
      lb <- -Inf
      ub <- 0
      mij <- Wi[j,]%*%alpha + Omega[j,-j]%*%solve(Omega[-j,-j])%*%(Ali[-j] - Wi[-j,]%*%alpha)
      sdij <- Omega[j,j]-Omega[j,-j]%*%solve(Omega[-j,-j])%*%(Omega[j,-j])
      Ali[j] <- EnvStats::rnormTrunc(1, mean = mij, sd = sdij^0.5, min = lb, max = ub)
      if(Ali[j] <= -20){ # Computational issues
        Ali[j] <- -20
      }
    }else{
      lb <- 0 
      ub <- Inf
      mij <- Wi[j,]%*%alpha + Omega[j,-j]%*%solve(Omega[-j,-j])%*%(Ali[-j]- Wi[-j,]%*%alpha)
      sdij <- Omega[j,j]-Omega[j,-j]%*%solve(Omega[-j,-j])%*%(Omega[j,-j])
      Ali[j] <- EnvStats::rnormTrunc(1, mean = mij, sd = sdij^0.5, min = lb, max = ub)
      if(Ali[j] >= 20){ # Computational issues
        Ali[j] <- 20
      }
    }
  }
  return(Ali)
}
# i <- 15
# PostAl(alpha = c(a1, a2, a3), Omega = OMEGA, Ali = Al[i,], Ai = A[i,], Wi = WW[[i]])
# Al[i,]

PostTheta <- function(Sigma, Cl, Y, ZX){
  WY <- cbind(Cl, Y)
  XtX <- matrix(0, K+L, K+L)
  Xtwy <- matrix(0, K+L, 1)
  ids <- cumsum(c(1,c(k1,k2,k3,l1,l2,l3)))
  for(m in 2:3^J){
    idGood <- which(CombNew[m,]==1)
    JJ <- matrix(0, K+L, length(idGood)*3)
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

# ThetaPost <- PostTheta(Sigma = SIGMA, Cl = Cl, Y = Y, ZX = ZX)
# theta <- c(b1, b2, b3, d1, d2, d3)
# cbind(ThetaPost, theta)


PostSig <- function(theta, Cl, Y, C, A, ZZ, ZX){
  WtWg0 <- 0
  G0i <- which(A[,1] == 1)
  for(i in G0i){
    WtWgi <- (Cl[i, 1] - ZZ[[i]][1,]%*%theta)^2
    WtWg0 <- WtWg0 + WtWgi
  }
  rr <- WtWg0 + 1
  sigma11r <- LaplacesDemon::rinvwishart(r0-4+length(G0i), rr)
  
  WtWg1 <- matrix(0, 2, 2)
  G1i <- which(A[,1] == 1 & A[,2] == 1)
  for(i in G1i){
    WtWgi <- (Cl[i, 1:2] - ZZ[[i]][1:2,]%*%theta)%*%t(Cl[i, 1:2] - ZZ[[i]][1:2,]%*%theta)
    WtWg1 <- WtWg1 + WtWgi
  }
  r11 <- WtWg1 + diag(2)
  r11i <- 1/(r11[1,1])
  r21 <- r11[2,1]
  r22 <- r11[2,2]
  r22.1 <- r22-r21^2*r11i
  sigma22.1r <- LaplacesDemon::rinvwishart(r0+length(G1i), r22.1)
  mr <- r21*r11i
  sigma21.1r <- LaplacesDemon::rmatrixnorm(mr, r11i, sigma22.1r)
  sigma21r <- sigma21.1r*sigma11r
  sigma22r <- sigma22.1r + sigma21.1r*sigma21r
  sigma11_1rns <- rbind(c(sigma11r, sigma21r), c(sigma21r, sigma22r))
  sigma11_1r <- sigma11_1rns 
  
  WtW11f <- matrix(0, J, J)
  G20 <- which(A[,1]==1 & A[,2]==1 & A[,3]==1)
  for(i in G20){
    WtW11fi <- (Cl[i,] - ZZ[[i]]%*%theta)%*%t(Cl[i,] - ZZ[[i]]%*%theta)
    WtW11f <- WtW11f + WtW11fi
  }
  r11Gf <- WtW11f  + diag(J)
  r11fi <- solve(r11Gf[1:2,1:2])
  r21f <- r11Gf[1:2,3]
  r22f <- r11Gf[3,3]
  r22.1f <- r22f-t(r21f)%*%r11fi%*%r21f
  sigma22.1rf <- LaplacesDemon::rinvwishart(r0+length(G20), r22.1f)
  mf <- t(r21f%*%r11fi)
  sigma21.1rf <- LaplacesDemon::rmatrixnorm(mf, as.matrix(Matrix::forceSymmetric(r11fi)), as.matrix(Matrix::forceSymmetric(sigma22.1rf)))
  sigma21rf <- t(sigma21.1rf)%*%sigma11_1r
  sigma22rf <- sigma22.1rf + t(sigma21.1rf)%*%t(sigma21rf)
  sigma11_1rfns <- cbind(rbind(sigma11_1r, sigma21rf), c(sigma21rf, sigma22rf))
  Sigma11 <- sigma11_1rfns 
  
  ClY1 <- cbind(Cl, Y[,1])
  RtR11 <- matrix(0, J+1, J+1)
  G <- which(C[,1]==1)
  for(i in G){
    RtR11i <- (ClY1[i,] - ZX[[i]][1:(J+1),]%*%theta)%*%t(ClY1[i,] - ZX[[i]][1:(J+1),]%*%theta)
    RtR11 <- RtR11 + RtR11i
  }
  R11G <- RtR11 + diag(J+1)
  R11i <- solve(R11G[1:J,1:J])
  R21 <- R11G[1:J,(J+1)]
  R22 <- R11G[J+1,J+1] 
  R22.1 <- R22-t(R21)%*%R11i%*%R21
  Sigma22.1 <- LaplacesDemon::rinvwishart(r0+length(G), R22.1)
  M <- t(R21%*%R11i) 
  Sigma21.1 <- LaplacesDemon::rmatrixnorm(M, as.matrix(Matrix::forceSymmetric(R11i)), as.matrix(Matrix::forceSymmetric(Sigma22.1)))
  Sigma21 <- t(Sigma21.1)%*%Sigma11
  Sigma22 <- Sigma22.1 + t(Sigma21.1)%*%t(Sigma21)
  Sigma11_1 <- cbind(rbind(Sigma11, Sigma21), c(Sigma21, Sigma22))
  
  ClY12 <- cbind(Cl, Y[,1:2])
  RtR11n <- matrix(0, J+2, J+2)
  Gn <- which(C[,1]==1 & C[,2]==1)
  for(i in Gn){
    RtR11ni <- (ClY12[i,] - ZX[[i]][1:(J+2),]%*%theta)%*%t(ClY12[i,] - ZX[[i]][1:(J+2),]%*%theta)
    RtR11n <- RtR11n + RtR11ni
  }
  R11Gn <- RtR11n + diag(J+2)
  R11ni <- solve(R11Gn[1:(J+1),1:(J+1)])
  R21n <- R11Gn[1:(J+1),J+2]
  R22n <- R11Gn[J+2,J+2] 
  R22.1n <- R22n-t(R21n)%*%R11ni%*%R21n
  Sigma22.1n <- LaplacesDemon::rinvwishart(r0+length(Gn), R22.1n)
  Mn <- t(R21n%*%R11ni) 
  Sigma21.1n <- LaplacesDemon::rmatrixnorm(Mn, as.matrix(Matrix::forceSymmetric(R11ni)), as.matrix(Matrix::forceSymmetric(Sigma22.1n)))
  Sigma21n <- t(Sigma21.1n)%*%Sigma11_1
  Sigma22n <- Sigma22.1n + t(Sigma21.1n)%*%t(Sigma21n)
  Sigma11_1n <- cbind(rbind(Sigma11_1, Sigma21n), c(Sigma21n, Sigma22n))
  
  ClY123 <- cbind(Cl, Y)
  RtR11f <- matrix(0, J+3, J+3)
  Gf <- which(C[,1]==1 & C[,2]==1 & C[,3]==1)
  for(i in Gf){
    RtR11fi <- (ClY123[i,] - ZX[[i]]%*%theta)%*%t(ClY123[i,] - ZX[[i]]%*%theta)
    RtR11f <- RtR11f + RtR11fi
  }
  R11Gf <- RtR11f + diag(J+3)
  R11fi <- solve(R11Gf[1:(J+2),1:(J+2)])
  R21f <- R11Gf[1:(J+2),J+3]
  R22f <- R11Gf[J+3,J+3] 
  R22.1f <- R22f-t(R21f)%*%R11fi%*%R21f
  Sigma22.1f <- LaplacesDemon::rinvwishart(r0+length(Gf), R22.1f)
  Mf <- t(R21f%*%R11fi) 
  Sigma21.1f <- LaplacesDemon::rmatrixnorm(Mf, as.matrix(Matrix::forceSymmetric(R11fi)), as.matrix(Matrix::forceSymmetric(Sigma22.1f)))
  Sigma21f <- t(Sigma21.1f)%*%Sigma11_1n
  Sigma22f <- Sigma22.1f + t(Sigma21.1f)%*%t(Sigma21f)
  Sigma11_1f <- cbind(rbind(Sigma11_1n, Sigma21f), c(Sigma21f, Sigma22f))
  
  return(Sigma11_1f)
}
# Sig11_1 <- PostSig(theta = theta, Cl = Cl, A = A, Y = Y, C = C, ZZ = ZZ, ZX = ZX)
# Sig11_1
# Sigma11_1 <- SIGMA[1:(J+3),1:(J+3)]
# Sigma11_1
# cbind(matrixcalc::vech(Sig11_1), matrixcalc::vech(Sigma11_1))

PostCl <- function(m, theta, Sigma, Cli, Yi, Xi, Zi, ZXi, Ci){
  idGood1 <- which(CombNew[m,1:3]==1)
  idGoodn <- which(CombNew[m,1:6]==0)
  CYli <- c(Cli,Yi)
  for(j in idGood1){
    if(Ci[j] == 0){
      lb <- -Inf
      ub <- 0
    }else{
      lb <- 0 
      ub <- Inf
    }
    if(length(idGood1)==1){
      mij <- Zi[j,]%*%theta
      sdij <- Sigma[j,j]
    }else{
      mij <- Zi[j,]%*%theta + Sigma[j,-c(j,idGoodn)]%*%solve(Sigma[-c(j,idGoodn),-c(j,idGoodn)])%*%(CYli[-c(j,idGoodn)] - ZXi[-c(j,idGoodn),]%*%theta)
      sdij <- Sigma[j,j]-Sigma[j,-c(j,idGoodn)]%*%solve(Sigma[-c(j,idGoodn),-c(j,idGoodn)])%*%(Sigma[-c(j,idGoodn),j])
    }
    Cli[j] <- EnvStats::rnormTrunc(1, mean = mij, sd = sdij^0.5, min = lb, max = ub)
    if(Cli[j] <= -20){ # Computational issues
      Cli[j] <- -20
    }
    if(Cli[j] >= 20){ # Computational issues
      Cli[j] <- 20
    }
  }
  return(Cli)
}
# i <- 36
# m <- Gs[i]
# theta = theta; Sigma = SIGMA; Cli = Cl[i,]; Yi = Y[i,]; Xi = XX[[i]]
# Zi = ZZ[[i]]; ZXi = ZX[[i]]; Ci = C[i,]
# EnsCl <- replicate(500, PostCl(m = m, theta = theta, Sigma = SIGMA, Cli = Cl[i,], Yi = Y[i,], Xi = XX[[i]], Zi = ZZ[[i]], ZXi = ZX[[i]], Ci = C[i,]))
# l <- 3
# plot(EnsCl[l,], type = "l")
# abline(h = Cli[l], col = "red")
# Cl[i,]

S <- 5000 # 1100
thin <- 10
burnin <- 1000 + thin
AlphaPost <- matrix(NA, S, H)
OmegaPost <- array(NA, c(J, J, S))
AlPost <- array(NA, c(N, J, S))
ThetaPost <- matrix(NA, S, K + L)
SigmaPost <- array(NA, c(2*J, 2*J, S))
ThetaPostNOst <- matrix(NA, S, K + L)
SigmaPostNOst <- array(NA, c(2*J, 2*J, S))
OmegaPostNOst <- array(NA, c(J, J, S))
ClPost <- array(NA, c(N, J, S))
findraws <- seq(burnin, S, thin)

Rep <- 100 # 2

ALPHApost <- array(0, c(H, 3, Rep))
THETApost <- array(0, c(K+L, 3, Rep))
OMEGApost <- array(0, c(J*(J+1)/2, 3, Rep))
SIGMApost <- array(0, c(2*J*(2*J+1)/2, 3, Rep))

THETApostNOst <- array(0, c(K+L, 3, Rep))
SIGMApostNOst <- array(0, c(2*J*(2*J+1)/2, 3, Rep))
OMEGApostNOst <- array(0, c(J*(J+1)/2, 3, Rep))


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
  Al2 <- W%*%a2 + U[, 2]
  Al3 <- W%*%a3 + U[, 3]
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
  
  #### Gibbs sampler: Implementation ####
  ALPHAp <- c(a1, a2, a3)
  OMEGA <- SIGMA[1:J,1:J]
  OMEGAp <- OMEGA
  OMEGASp <- OMEGA
  SIGMAp <- SIGMA[(J+1):(3*J),(J+1):(3*J)]
  THETAp <- c(b1, b2, b3, d1, d2, d3)
  Clp <- Cl
  Alp <- Al
  

  clusterExport(cl, list("J", "K", "L", "k1", "PostCl", "ZZ", "C", "THETAp",
                         "SIGMAp", "XX", "Y", "ZX", "X", "Z", "H", "h1", "ALPHAp",
                         "OMEGAp", "PostAl", "WW", "W", "A", "CombNew", "Gs", "N"))
  
  Wp <- W
  for(s in 1:S){
    # ticks <- Sys.time()
    clusterExport(cl, list("THETAp", "SIGMAp", "Clp", "Alp", "OMEGAp", "ALPHAp", "OMEGASp"))
    Alp <- t(parSapply(cl, 1:N, function(i){PostAl(alpha = ALPHAp, Omega = OMEGAp, Ali = Alp[i,], Ai = A[i,], Wi = WW[[i]])}))
    OMEGAp <- PostOmega(alpha = ALPHAp, Al = Alp, W = WW)
    OMEGASp <- diag(1/(diag(OMEGAp)^0.5))%*%OMEGAp%*%diag(1/(diag(OMEGAp)^0.5))
    ALPHAp <- PostAlpha(Omega=OMEGAp, Al = Alp, W = WW)
    ALPHASp12 <- ALPHAp[1:3]/OMEGAp[1,1]^0.5
    ALPHASp34 <- ALPHAp[4:6]/OMEGAp[2,2]^0.5
    ALPHASp56 <- ALPHAp[7:9]/OMEGAp[3,3]^0.5
    
    Clp <- t(parSapply(cl, 1:N, function(i){PostCl(m = Gs[i], theta = THETAp, Sigma = SIGMAp, Cli = Clp[i,], Yi = Y[i,], Xi = XX[[i]], Zi = ZZ[[i]], ZXi = ZX[[i]], Ci = C[i,])}))
    SIGMAp <- PostSig(theta = THETAp, Cl = Clp, A = A, C = C, Y = Y, ZZ = ZZ, ZX = ZX)
    SIGMASpProb <- diag(1/(diag(SIGMAp[1:J,1:J])^0.5))%*%SIGMAp[1:J,1:J]%*%diag(1/(diag(SIGMAp[1:J,1:J])^0.5))
    SIGMApNew <- SIGMAp
    SIGMApNew[1:J,1:J] <- SIGMASpProb 
    THETAp <- PostTheta(Sigma = SIGMAp, Cl = Clp, Y = Y, ZX = ZX)
    THETASp13 <- THETAp[1:3]/SIGMAp[1,1]^0.5
    THETASp46 <- THETAp[4:6]/SIGMAp[2,2]^0.5
    THETASp79 <- THETAp[7:9]/SIGMAp[3,3]^0.5
    
    AlphaPost[s,] <- c(ALPHASp12,ALPHASp34,ALPHASp56)
    OmegaPost[,,s] <- OMEGASp
    OmegaPostNOst[,,s] <- OMEGAp
    AlPost[,,s] <- Alp
    ThetaPost[s,] <- c(THETASp13, THETASp46, THETASp79, THETAp[-c(1:9)])
    SigmaPost[,,s] <- SIGMApNew
    ThetaPostNOst[s,] <- THETAp
    SigmaPostNOst[,,s] <- SIGMAp
    ClPost[,,s] <- Clp
    # tocks <- Sys.time()
    # print(tocks-ticks)
    # print(s)
  }
  
  alphaHat <- coda::mcmc(AlphaPost[findraws,])
  RestAlpha <- summary(alphaHat)
  ALPHApost[,,rep] <- cbind(RestAlpha$statistics[,1], RestAlpha$quantiles[,c(1,5)])  
  
  OmegaHatVech <- coda::mcmc(t(sapply(1:S, function(s){matrixcalc::vech(OmegaPost[,,s])})))
  RestOmega <- summary(coda::mcmc(OmegaHatVech[findraws,]))
  OMEGApost[,,rep] <- cbind(RestOmega$statistics[,1], RestOmega$quantiles[,c(1,5)]) 
  
  OmegaHatVechNOst <- coda::mcmc(t(sapply(1:S, function(s){matrixcalc::vech(OmegaPostNOst[,,s])})))
  RestOmegaNOst <- summary(coda::mcmc(OmegaHatVechNOst[findraws,]))
  OMEGApostNOst[,,rep] <- cbind(RestOmegaNOst$statistics[,1], RestOmegaNOst$quantiles[,c(1,5)]) 
  
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
  
  PostResults <- list(ALPHApost = ALPHApost, THETApost = THETApost,
                      OMEGApost = OMEGApost, SIGMApost = SIGMApost,
                      OMEGApostNOst = OMEGApostNOst,
                      THETApostNOst = THETApostNOst, SIGMApostNOst = SIGMApostNOst)
  
  save(PostResults, file = "PostResultsJ3ExAccess.RData")
  print(rep)
  tock <- Sys.time()
  print(tock-tick)
  rep <- rep + 1
}
stopCluster(cl)



