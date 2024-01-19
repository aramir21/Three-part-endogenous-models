##### Simulation: Three part incidental truncation (dependence between MVP) model ########
############## Estimation assuming endogenous access #############
# Andrés Ramírez Hassan
# Febraury 1st, 2023
rm(list = ls())
set.seed(010101)
library(doParallel)
library(snow)
Data <- read.csv("SSPP456.csv", header = TRUE, sep = ",")
attach(Data)
str(Data)
N <- dim(Data)[1]
J <- 3

## Create this objects as matrices ## 
# CatMunp <- as.matrix(fastDummies::dummy_cols(cat_mpio)[,-1])
# write.csv(CatMunp, file = "CatMunp.csv")
CatMunp <- read.csv("CatMunp.csv", header = TRUE)
CatMunpA <- CatMunp[,2] + CatMunp[,7] # High category
CatMunpB <- CatMunp[,3] + CatMunp[,4] # Medium category
CatMunpC <- CatMunp[,5] + CatMunp[,6] # Low category
# Use as reference CatMunpA
W <- cbind(1, CatMunpB, CatMunpC, est4, est5, urban) #regressors access including intercept

A1 <- sitieneacueducto
A2 <- sitieneenergia
A3 <- sitienegas
A <- cbind(A1, A2, A3) #binary access
summary(A)

Z <- cbind(1, est4, est5, lnincome, lnpfullacu, lnpfullener,	lnpfullgas) #regressors use including intercept 
C1 <- siconsumeacu
C2 <- siconsumeener
C3 <- siconsumegas
C <- cbind(C1, C2, C3)
AC <- cbind(A1, C1, A2, C2, A3, C3)

male <- gene_jefe
X <- cbind(1, est4, est5, altitude1, edu1, edu2, edu3, edu4, male, edadjefe, ncuartos, npersonashogar, tienehijos, lnincome, 
           lnpfullacu, lnpfullener,	lnpfullgas) #regressors quantity including intercept
Y1 <- log(xobsacu)
Y2 <- log(xobsener)
Y3 <- log(xobsgas)
Y <- cbind(Y1, Y2, Y3)

CY <- cbind(C, Y)

h1 <- dim(W)[2]; h2 <- dim(W)[2]; h3 <- dim(W)[2] # Dim(Tj), see below
H <- h1 + h2 + h3 # Dim(T), see below
k1 <- dim(Z)[2]; k2 <- dim(Z)[2]; k3 <- dim(Z)[2] # Dim(Bj), See below
K <- k1 + k2 + k3 # Dim(B), See below
l1 <- dim(X)[2]; l2 <- dim(X)[2]; l3 <- dim(X)[2] # Dim(Dj), See below
L <- l1 + l2 + l3 # Dim(D), See below

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

#### Data preparation ####
IJ <- diag(J)
CJ1 <- matrix(0, J, 2*J)
CJ <- matrix(0, J, J)
WJ <- matrix(0, J, K+L)
ZJ1 <- matrix(0, J, H)
ZJ2 <- matrix(0, J, L)
XJ <- matrix(0, J, H+K)

WW <- lapply(1:N, function(i){cbind(kronecker(IJ, t(W[i, ])), WJ)})
ZZ <- lapply(1:N, function(i){cbind(ZJ1, kronecker(IJ, t(Z[i, ])), ZJ2)})
XX <- lapply(1:N, function(i){cbind(XJ, kronecker(IJ, t(X[i, ])))})
WZ <- lapply(1:N, function(i){rbind(WW[[i]], ZZ[[i]])})
WZX <- lapply(1:N, function(i){rbind(WW[[i]], ZZ[[i]], XX[[i]])})


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

#### Gibbs sampler: Functions ####
PostTheta <- function(Sigma, Al, Cl, A, C, Y, WZX){
  WY <- cbind(Al, Cl, Y)
  XtX <- matrix(0, H+K+L, H+K+L)
  Xtwy <- matrix(0, H+K+L, 1)
  ids <- cumsum(c(1,c(h1,h2,h3,k1,k2,k3,l1,l2,l3)))
  for(m in 1:3^J){
    idGood <- c(1:3,which(CombNew[m,]==1)+3)
    idGood1 <- which(CombNew[m,c(1:3)]==1)
    idGood2 <- which(CombNew[m,c(4:6)]==1)
    JJ <- matrix(0, H+K+L, H+length(idGood1)*k1+length(idGood2)*l1)
    idsCol <- cumsum(c(1,c(h1,h2,h3, rep(k1,length(idGood1)), rep(l1, length(idGood2)))))
    idcov <- NULL
    for(l in 1:length(idGood)){
      idscol <- idsCol[l]:(idsCol[l+1]-1)
      JJ[ids[idGood[l]]:(ids[idGood[l]+1]-1),idscol] <- diag(length(idscol))
      idcovl <- ids[idGood[l]]:(ids[idGood[l]+1]-1)
      idcov <- c(idcov, idcovl)
    }
    G <- Groups[[m]]
    for(i in G){
      ZXm <- WZX[[i]][idGood,]
      ZXmi <- ZXm[,idcov]
      XtXi <- JJ%*%t(ZXmi)%*%solve(Sigma[idGood,idGood])%*%ZXmi%*%t(JJ)
      XtX <- XtX + XtXi
      Xtwyi <- JJ%*%t(ZXmi)%*%solve(Sigma[idGood,idGood])%*%WY[i,idGood]
      Xtwy <- Xtwy + Xtwyi
    }
  }
  Tn <- solve(XtX + T0i)
  tn <- Tn%*%(Xtwy + T0i%*%t0)
  Tpost <- MASS::mvrnorm(1, tn, Tn)
  return(Tpost)
}
# THETAPost <- PostTheta(Sigma = SIGMAp, Al = Alp, Cl = Clp, A = A, C = C, Y = Y, WZX = WZX)
# cbind(THETAPost,THETAp)

PostSig <- function(Theta, Al, Cl, A, C, Y, WZX){
  ACY <- cbind(Al, Cl, Y)
  eqs <- J
  WtW <- matrix(0, eqs, eqs)
  for(i in 1:N){
    WtWi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
    WtW <- WtW + WtWi
  }
  sigma11r <- LaplacesDemon::rinvwishart(r0-2*J+N, WtW + diag(eqs))
  
  eqs <- 4
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
  
  eqs <- 5
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(A[,1] == 1 & A[,2] == 1)
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
  
  eqs <- 6
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(A[,1] == 1 & A[,2] == 1 & A[,3] == 1)
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
  
  eqs <- 7
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(C[,1] == 1)
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
  
  eqs <- 8
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(C[,1] == 1 & C[,2] == 1)
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
  
  eqs <- 9
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(C[,1] == 1 & C[,2] == 1 & C[,3] == 1)
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
  idGood1 <- c(1:3, which(CombNew[m,1:3]==1) + 3)
  idGoodn <- which(CombNew[m,1:6]==0) + 3
  ACYli <- c(ACli,Yi)
  ACi <- c(Ai, Ci)
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
  return(ACli)
}

# i <- 400
# m = Gs[i]; theta = THETAp; Sigma = SIGMAp; ACli = c(Alp[i,],Clp[i,]);
# Yi = Y[i,]; Ai = A[i,]; Ci = C[i,]; WZXi = WZX[[i]]
# ACliPost <- PostACl(m = 1, theta = THETAp, Sigma = SIGMAp, ACli = c(Alp[i,],Clp[i,]), Yi = Y[i,], Ai = A[i,], Ci = C[i,], WZXi = WZX[[i]])
# cbind(ACli, ACliPost)

#### Gibbs sampler: Implementation ####
S <- 20000 # 6000 # 1100
thin <- 10 # 5
burnin <- 10000 + thin
ThetaPost <- matrix(NA, S, H + K + L)
SigmaPost <- array(NA, c(3*J, 3*J, S))
ThetaPostNOst <- matrix(NA, S, H + K + L)
SigmaPostNOst <- array(NA, c(3*J, 3*J, S))
AClPost <- array(NA, c(N, 2*J, S))
findraws <- seq(burnin, S, thin)

SIGMAp <- diag(3*J) # Set an initial value for covariance matrix
THETAp <- rep(1, H+K+L) # Set an initial value for location parameters
Al1 <- matrix(c(ifelse(A1==1, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Al2 <- matrix(c(ifelse(A2==1, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Al3 <- matrix(c(ifelse(A3==1, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Alp <- cbind(Al1, Al2, Al3)
Cl1 <- matrix(c(ifelse(C1==0, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Cl2 <- matrix(c(ifelse(C2==0, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Cl3 <- matrix(c(ifelse(C3==0, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Clp <- cbind(Cl1, Cl2, Cl3)
AClp <- cbind(Alp, Clp)

# THETApost <- array(0, c(H+K+L, 3, Rep))
# SIGMApost <- array(0, c(3*J*(3*J+1)/2, 3, Rep))

cn <- 6 # detectCores() # 6
cl <- makeCluster(cn, type = "SOCK")
registerDoParallel(cl)


#### Parallel code ####
clusterExport(cl, list("J", "K", "L", "k1", "PostACl", "ZZ", "C", "THETAp",
                       "SIGMAp", "H", "h1", "WZX", "A", "CombNew", "Gs", "N",
                       "Y", "AClp"))

for(s in 1:S){
  ticks <- Sys.time()
  clusterExport(cl, list("THETAp", "SIGMAp", "AClp"))
  AClp <- t(parSapply(cl, 1:N, function(i){PostACl(m = Gs[i], theta = THETAp, Sigma = SIGMAp, ACli = AClp[i,], Yi = Y[i,], WZXi = WZX[[i]], Ci = C[i,], Ai = A[i,])}))
  # SIGMAp <- PostSig(Theta = THETAp, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
  SIGMAp <- PostSig(Theta = THETAp, Al = AClp[,1:3], Cl = AClp[,4:6], A = A, C = C, Y = Y, WZX = WZX)
  SIGMASpProb <- diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))%*%SIGMAp[1:(2*J),1:(2*J)]%*%diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))
  SIGMApNew <- SIGMAp
  SIGMApNew[1:(2*J),1:(2*J)] <- SIGMASpProb 
  # THETAp <- PostTheta(Sigma = SIGMAp, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
  THETAp <- PostTheta(Sigma = SIGMAp, Al = AClp[,1:3], Cl = AClp[,4:6], A = A, C = C, Y = Y, WZX = WZX)
  THETASp12 <- THETAp[1:h1]/SIGMAp[1,1]^0.5
  THETASp34 <- THETAp[(h1+1):(h1+h2)]/SIGMAp[2,2]^0.5
  THETASp56 <- THETAp[(h1+h2+1):(h1+h2+h3)]/SIGMAp[3,3]^0.5
  THETASp79 <- THETAp[(h1+h2+h3+1):(h1+h2+h3+k1)]/SIGMAp[4,4]^0.5
  THETASp02 <- THETAp[(h1+h2+h3+k1+1):(h1+h2+h3+k1+k2)]/SIGMAp[5,5]^0.5
  THETASp35 <- THETAp[(h1+h2+h3+k1+k2+1):(h1+h2+h3+k1+k2+k3)]/SIGMAp[6,6]^0.5
  
  ThetaPost[s,] <- c(THETASp12, THETASp34, THETASp56, THETASp79, THETASp02, THETASp35,THETAp[-c(1:(H+K))])
  SigmaPost[,,s] <- SIGMApNew
  ThetaPostNOst[s,] <- THETAp
  SigmaPostNOst[,,s] <- SIGMAp
  AClPost[,,s] <- AClp
  PostResults <- list(ThetaPost = ThetaPost, SigmaPost = SigmaPost, ThetaPostNOst = ThetaPostNOst, SigmaPostNOst = SigmaPostNOst)
  save(PostResults, file = "PostResultsAppSSPPV1.RData")
  tocks <- Sys.time()
  print(tocks-ticks)
  print(s)
}

load("PostResultsAppSSPPV1.RData")
S <- 20000 # 6000 # 1100
thin <- 2 # 5
burnin <- 18000 + thin
findraws <- seq(burnin, S, thin)
ThetaPost <- PostResults[["ThetaPost"]]
thetaHat <- coda::mcmc(ThetaPost[findraws,])
plot(thetaHat)
RestTheta <- summary(thetaHat)
RestThetaQuant <- RestTheta$quantiles
RestThetaStats <- RestTheta$statistics[,1:2]
rownames(RestThetaQuant) <- c(rep(colnames(W), J), rep(colnames(Z), J), rep(colnames(X), J))
write.csv(RestThetaQuant, file = "LocationQuantSSPPV2.csv")
rownames(RestThetaStats) <- c(rep(colnames(W), J), rep(colnames(Z), J), rep(colnames(X), J))
write.csv(RestThetaStats, file = "LocationStatsSSPPV2.csv")

COVAR <- PostResults[["SigmaPostNOst"]]
COVAR1 <- t(sapply(findraws, function(s){diag(1/(diag(COVAR[,,s])^0.5))%*%COVAR[,,s]%*%diag(1/(diag(COVAR[,,s])^0.5))})) 
SigmaHat <- coda::mcmc(COVAR1[,c(2:9,12:18,22:27,32:36,42:45,52:54,62:63,72)])
plot(SigmaHat)
RestSigma <- summary(SigmaHat)
RestSigmaQuant <- RestSigma$quantiles
write.csv(RestSigmaQuant, file = "ScaleQuantSSPPV2.csv")
RestSigmaStats <- RestSigma$statistics[,1:2]
write.csv(RestSigmaStats, file = "ScaleStatsSSPPV2.csv")

library(coda)
autocorr.plot(thetaHat)
raftery.diag(thetaHat, r=0.01)
geweke.diag(thetaHat)
heidel.diag(thetaHat)

autocorr.plot(SigmaHat)
raftery.diag(SigmaHat, r=0.01)
geweke.diag(SigmaHat)
heidel.diag(SigmaHat)

stopCluster(cl)

