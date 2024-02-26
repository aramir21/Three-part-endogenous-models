########### R Script - Three part incidental truncation model ########## 
# 
# Creation date:      ?
# Modification date:  08/15/2023
# 
# Description:
#   
# Notes:
#     -
#     -
#     -
#   
# Author: Andres Ramirez & Santiago Velasquez

rm(list = ls())
set.seed(0101)

# Libraries
library(doParallel)
library(snow)
library(dplyr)
library(matrixcalc)

# Folder directory
# setwd("C:/Users/tatia/OneDrive - Universidad EAFIT/IEFIC_2018")
#setwd("C:/Users/svelasquez1/Dropbox (BFI)/SVB_4_Consumo/Data/process/MC/Final")
# Import database
Data <- read.csv("Base_finalCREDITO.csv", header = TRUE, sep = ",")
# Data <- read.csv("Data/Base_final.csv", header = TRUE, sep = ",") %>%
#   mutate(usaLI=ifelse(accesoLI==0,0,usaLI),
#          usaTC=ifelse(accesoTC==0,0,usaTC),
#          cuanto_usa_TC=ifelse(accesoTC==0,0,cuanto_usa_TC),
#          cuanto_usa_LI=ifelse(accesoLI==0,0,cuanto_usa_LI),
#          cuanto_usa_TC=ifelse(accesoTC==1 & usaTC==0,0,cuanto_usa_TC),
#          cuanto_usa_LI=ifelse(accesoLI==1 & usaLI==0,0,cuanto_usa_LI))
# Data<-Data[!(Data$accesoTC==1 & is.na(Data$usaTC)==T),]
# Data<-Data[!(Data$accesoLI==1 & is.na(Data$usaLI)==T),]
# Data <- Data[!(is.na(Data$cuenta_ahorro)==T),] #38.856
# Data <- Data[!(is.na(Data$edu_2)==T),] #37375
# Data <- Data[!(is.na(Data$INGTOTOB)==T),] #37287
# Data <- Data[!(is.na(Data$str_1)==T),] #37.017
# readr::write_csv(Data, "Data/Base_final.csv")
attach(Data)
str(Data)
N <- dim(Data)[1]
J <- 2

## Create this objects as matrices ## 

## Access
# Reference categories: 
W <- cbind(1, oci_pensionado, cuenta_ahorro, req_ingreso, edad,
           cot_salud_pension) #regressors access including intercept

A1 <- accesoTC  
A2 <- accesoLI
A <- cbind(A1, A2) #binary access
summary(A)

## Extensive margin
# quite capacidad de deuda,
# Reference categories: str_1, edu_1
Z <- cbind(1, genero, jefe_hogar,edad, edad^2, edu_2, edu_3, edu_4, edu_5, 
             haymenores, vivienda_propia,tiene_negocio,ind_riqueza_bienes,
             tiene_carro,tiene_moto,oci_pensionado, log(INGTOTOB),
             str_2, str_3,str_4, str_5, str_6,
             otroscreditos , LIhogar , TChogar) #regressors use including intercept
#cor(Z)
C1 <- usaTC
C2 <- usaLI 
C <- cbind(C1, C2)
AC <- cbind(A1, C1, A2, C2)

## Intensive margin
# Reference categories: edu_1, vivienda_otros, INI, str_1
X <- cbind(1, genero, jefe_hogar, edad, edad^2, edu_2, edu_3, edu_4, edu_5 ,
             haymenores, vivienda_propia, tiene_negocio, ind_riqueza_bienes,
             tiene_carro, tiene_moto, oci_pensionado,log(INGTOTOB),
             str_2, str_3, str_4, str_5, str_6) #regressors quantity including intercept
Y1 <- log(cuanto_usa_TC)
Y2 <- log(cuanto_usa_LI)
Y <- cbind(Y1, Y2)

CY <- cbind(C, Y)


h1 <- dim(W)[2]; h2 <- dim(W)[2] # Dim(Tj), see below
H <- h1 + h2 # Dim(T), see below
k1 <- dim(Z)[2]; k2 <- dim(Z)[2] # Dim(Bj), See below
K <- k1 + k2 # Dim(B), See below
l1 <- dim(X)[2]; l2 <- dim(X)[2] # Dim(Dj), See below
L <- l1 + l2 # Dim(D), See below

# Groups: 3^J
CombId <- rbind(c(0, 0), c(1, 0), c(1, 1))
id1Comb <- rep(1:3, 1, each = 3^(J-1))
id2Comb <- rep(1:3, 3^(J-1))
idsComb <- cbind(id1Comb, id2Comb)
Comb1 <- matrix(0, 3^J, 2)
Comb2 <- matrix(0, 3^J, 2)
for(id in 1:3^J){
  Comb1[id, ] <- CombId[id1Comb[id],]
  Comb2[id, ] <- CombId[id2Comb[id],]
}
Comb <- cbind(Comb1, Comb2) # Combination access/use per drug: A1 C1 A2 C2 A3 C3
CombNew <- cbind(Comb[,c(1,3)],Comb[,c(2,4)]) # First two columns is access and the two last ones are use: A1 A2  C1 C2 


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
  ids <- cumsum(c(1,c(h1,h2,k1,k2,l1,l2)))
  for(m in 1:3^J){
    idGood <- c(1:2,which(CombNew[m,]==1)+2) # Que ecuacion estoy estimando
    idGood1 <- which(CombNew[m,c(1:2)]==1)
    idGood2 <- which(CombNew[m,c(3:4)]==1)
    JJ <- matrix(0, H+K+L, H+length(idGood1)*k1+length(idGood2)*l1)
    idsCol <- cumsum(c(1,c(h1,h2, rep(k1,length(idGood1)), rep(l1, length(idGood2)))))
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


### Standard errors
PostSig <- function(Theta, Al, Cl, A, C, Y, WZX){
  ACY <- cbind(Al, Cl, Y)
  eqs <- J
  WtW <- matrix(0, eqs, eqs)
  for(i in 1:N){
    WtWi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
    WtW <- WtW + WtWi
  }
  sigma11r <- LaplacesDemon::rinvwishart(r0-2*J+N, WtW + diag(eqs))
  
  eqs <- 3
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
  
  eqs <- 4
  WtWg1 <- matrix(0, eqs, eqs)
  G1i <- which(A[,1] == 1 & A[,2] == 1)
  for(i in G1i){
    WtWgi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
    WtWg1 <- WtWg1 + WtWgi
  }
  r11 <- WtWg1 + diag(eqs)
  r11i <- solve(r11[1:(eqs-1),1:(eqs-1)], tol =-0.0000000000117)
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
  
  # eqs <- 6
  # WtWg1 <- matrix(0, eqs, eqs)
  # G1i <- which(A[,1] == 1 & A[,2] == 1 & A[,3] == 1)
  # for(i in G1i){
  #   WtWgi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
  #   WtWg1 <- WtWg1 + WtWgi
  # }
  # r11 <- WtWg1 + diag(eqs)
  # r11i <- solve(r11[1:(eqs-1),1:(eqs-1)])
  # r21 <- r11[eqs,1:(eqs-1)]
  # r22 <- r11[eqs,eqs]
  # r22.1 <- r22-t(r21)%*%r11i%*%r21
  # sigma22.1r <- LaplacesDemon::rinvwishart(r0+length(G1i), r22.1)
  # mr <- t(r21)%*%r11i
  # sigma21.1r <- LaplacesDemon::rmatrixnorm(t(mr), as.matrix(Matrix::forceSymmetric(r11i)), sigma22.1r)
  # sigma21r <- t(sigma21.1r)%*%sigma11_1r
  # sigma22r <- sigma22.1r + sigma21r%*%sigma21.1r
  # sigma11_1rns <- cbind(rbind(sigma11_1r, sigma21r), c(sigma21r, sigma22r))
  # sigma11_1r <- sigma11_1rns 
  # 
  eqs <- 5
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
  
  eqs <- 6
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
  
  # eqs <- 9
  # WtWg1 <- matrix(0, eqs, eqs)
  # G1i <- which(C[,1] == 1 & C[,2] == 1 & C[,3] == 1)
  # for(i in G1i){
  #   WtWgi <- (ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)%*%t(ACY[i, 1:eqs] - WZX[[i]][1:eqs,]%*%Theta)
  #   WtWg1 <- WtWg1 + WtWgi
  # }
  # r11 <- WtWg1 + diag(eqs)
  # r11i <- solve(r11[1:(eqs-1),1:(eqs-1)])
  # r21 <- r11[eqs,1:(eqs-1)]
  # r22 <- r11[eqs,eqs]
  # r22.1 <- r22-t(r21)%*%r11i%*%r21
  # sigma22.1r <- LaplacesDemon::rinvwishart(r0+length(G1i), r22.1)
  # mr <- t(r21)%*%r11i
  # sigma21.1r <- LaplacesDemon::rmatrixnorm(t(mr), as.matrix(Matrix::forceSymmetric(r11i)), sigma22.1r)
  # sigma21r <- t(sigma21.1r)%*%sigma11_1r
  # sigma22r <- sigma22.1r + sigma21r%*%sigma21.1r
  # sigma11_1rns <- cbind(rbind(sigma11_1r, sigma21r), c(sigma21r, sigma22r))
  # sigma11_1r <- sigma11_1rns
  
  return(sigma11_1r)
}

# Sig11_1 <- PostSig(Theta = THETAp, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
# cbind(matrixcalc::vech(Sig11_1), matrixcalc::vech(SIGMAp))


PostACl <- function(m, theta, Sigma, ACli, Yi, Ai, Ci, WZXi){
  idGood1 <- c(1:2, which(CombNew[m,1:2]==1) + 2)
  idGoodn <- which(CombNew[m,1:4]==0) + 2
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
S <- 1000 # 6000 # 1100
thin <- 5
burnin <- S*0.2
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
Alp <- cbind(Al1, Al2)
Cl1 <- matrix(c(ifelse(C1==0, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Cl2 <- matrix(c(ifelse(C2==0, EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0), EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf))),N,1)
Clp <- cbind(Cl1, Cl2)
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
  SIGMAp <- PostSig(Theta = THETAp, Al = AClp[,1:2], Cl = AClp[,3:4], A = A, C = C, Y = Y, WZX = WZX)
  SIGMASpProb <- diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))%*%SIGMAp[1:(2*J),1:(2*J)]%*%diag(1/(diag(SIGMAp[1:(2*J),1:(2*J)])^0.5))
  SIGMApNew <- SIGMAp
  SIGMApNew[1:(2*J),1:(2*J)] <- SIGMASpProb 
  # THETAp <- PostTheta(Sigma = SIGMAp, Al = Al, Cl = Cl, A = A, C = C, Y = Y, WZX = WZX)
  THETAp <- PostTheta(Sigma = SIGMAp, Al = AClp[,1:2], Cl = AClp[,3:4], A = A, C = C, Y = Y, WZX = WZX)
  THETASp12 <- THETAp[1:h1]/SIGMAp[1,1]^0.5 # Estandarizar por la desviacion estandar de acceso a marijuana
  THETASp34 <- THETAp[(h1+1):(h1+h2)]/SIGMAp[2,2]^0.5
  THETASp79 <- THETAp[(h1+h2+1):(h1+h2+k1)]/SIGMAp[3,3]^0.5
  THETASp02 <- THETAp[(h1+h2+k1+1):(h1+h2+k1+k2)]/SIGMAp[4,4]^0.5
  
  ThetaPost[s,] <- c(THETASp12, THETASp34, THETASp79, THETASp02,THETAp[-c(1:(H+K))])
  SigmaPost[,,s] <- SIGMApNew
  ThetaPostNOst[s,] <- THETAp
  SigmaPostNOst[,,s] <- SIGMAp
  AClPost[,,s] <- AClp
  PostResults <- list(ThetaPost = ThetaPost, SigmaPost = SigmaPost, ThetaPostNOst = ThetaPostNOst, SigmaPostNOst = SigmaPostNOst)
  save(PostResults, file = "Results/PostResultsAppCREDV2.RData")
  tocks <- Sys.time()
  print(tocks-ticks)
  print(s)
}



load("Results/PostResultsAppCREDV2.RData")
S <- 100 # 6000 # 1100
thin <- 5
burnin <- 10 + thin
findraws <- seq(burnin, S, thin)
ThetaPost <- PostResults[["ThetaPost"]]
thetaHat <- coda::mcmc(ThetaPost[findraws,])
# plot(thetaHat)
RestTheta <- summary(thetaHat)
RestThetaQuant <- RestTheta$quantiles
RestThetaStats <- RestTheta$statistics[,1:2]
rownames(RestThetaQuant) <- c(rep(colnames(W), J), rep(colnames(Z), J), rep(colnames(X), J))
write.csv(RestThetaQuant, file = "Results/LocationQuantCRED2.csv")
rownames(RestThetaStats) <- c(rep(colnames(W), J), rep(colnames(Z), J), rep(colnames(X), J))
write.csv(RestThetaStats, file = "Results/LocationStatsCRED2.csv")
#....................
# guardar en excel para overleaf (OPCIONAL)
base1 <- readr::read_csv("Results/LocationQuantCRED2.csv")
base2 <- readr::read_csv("Results/LocationStatsCRED2.csv")
writexl::write_xlsx(base1,"Results/LocationQuantCRED2.xlsx")
writexl::write_xlsx(base2,"Results/LocationStatsCRED2.xlsx")
#....................

COVAR <- PostResults[["SigmaPostNOst"]]
COVAR1 <- t(sapply(findraws, function(s){diag(1/(diag(COVAR[,,s])^0.5))%*%COVAR[,,s]%*%diag(1/(diag(COVAR[,,s])^0.5))})) 
SigmaHat <- coda::mcmc(COVAR1[,c(2:6,9:12,16:18,23:24,30)]) # Matriz de varianzas y covarianzas
par(mar = c(1,1,1,1))
plot(SigmaHat)
RestSigma <- summary(SigmaHat)
RestSigmaQuant <- RestSigma$quantiles
write.csv(RestSigmaQuant, file = "Results/ScaleQuantCREDV2.csv")
RestSigmaStats <- RestSigma$statistics[,1:2]
write.csv(RestSigmaStats, file = "Results/ScaleStatsCREDV2.csv")

#....................
# guardar en excel para overleaf (OPCIONAL)
base1 <- readr::read_csv("Results/ScaleQuantCREDV2.csv")
base2 <- readr::read_csv("Results/ScaleStatsCREDV2.csv")
writexl::write_xlsx(base1,"Results/ScaleQuantCREDV2.xlsx")
writexl::write_xlsx(base2,"Results/ScaleStatsCREDV2.xlsx")
#....................

library(coda)
autocorr.plot(thetaHat)
raftery.diag(thetaHat, r=0.01) #este es el que no hace por tamanio de la muestra
geweke.diag(thetaHat)
heidel.diag(thetaHat)

autocorr.plot(SigmaHat)
raftery.diag(SigmaHat, r=0.01)
geweke.diag(SigmaHat)
heidel.diag(SigmaHat)

stopCluster(cl)










