### Analysis simulation ###
rm(list = ls())
library(MCMCpack)
library(matrixcalc)
load("PostResultsV3BnewV2.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
a2 <- c(0.8, -1.2)
a3 <- c(1.1, -0.7)
J <- 3
rho <- 1
SIGMA <- rho*matrix(c(1,0.6,0.4,0.5,0.4,0.2,0,0,0,0.6,1,0.5,0.4,0.3,0.4,0,0,0,0.4,0.5,1,
                      0.5,0.3,0.5,0,0,0,0.5,0.4,0.5,1,0.4,0.5,0.3,0.4,0.2,0.4,0.3,0.3,0.4,1,
                      0.4,0.3,0.3,0.1,0.2,0.4,0.5,0.5,0.4,1,0.5,0.3,0.1,0,0,0,0.3,0.3,0.5,1,0.6,
                      0.4,0,0,0,0.4,0.3,0.3,0.6,1,0.3,0,0,0,0.2,0.1,0.1,0.4,0.3,1), 3*J, 3*J)
SIGMA13 <- matrix(c(0.3, 0.2, 0.3, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2), 3, 3) 
SIGMA[1:3,7:9] <- SIGMA13
SIGMA[7:9,1:3] <- t(SIGMA13)
SIGMApop <- vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)
b2 <- c(0.5, 1.5, -1)
b3 <- c(1, 1, -1)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)
d2 <- c(1, 2, 2)
d3 <- c(1, 1.5, 1.8)

THETApop <- c(a1, a2, a3, b1, b2, b3, d1, d2, d3)
THETAhat <- PostResults$THETApost[,1,]
cbind(THETApop, rowMeans(THETAhat))
# SIGMAhat <- PostResults$SIGMApost[,1,]
SIGMAhat <- sapply(1:100, function(l){vech(diag(1/diag(xpnd(PostResults$SIGMApostNOst[,1,l]))^0.5)%*%xpnd(PostResults$SIGMApostNOst[,1,l])%*%diag(1/diag(xpnd(PostResults$SIGMApostNOst[,1,l]))^0.5))})
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars){
  RMSE <- (mean(sapply(1:length(pars), function(i){(pop-pars[i])^2})))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i){RMSEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
RMSESigma <- sapply(1:length(SIGMApop), function(i){RMSEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 
RMSESigmanew <- sapply(1:length(SIGMApop), function(i){RMSEfunct(SIGMApop[i],SIGMAhat[i,1])}) 
cbind(RMSESigma, RMSESigmanew)
plot(RMSESigma, type = "l")
lines(RMSESigmanew, col = "red") # It's better to standarized at the end!!!

MAPEfunct <- function(pop, pars){
  MAPE <- mean(sapply(1:length(pars), function(i){abs((pop-pars[i])/pop)}))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i){MAPEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i){MAPEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 
max(MAPESigma)
MAPESigmanew <- sapply(1:length(SIGMApop), function(i){MAPEfunct(SIGMApop[i],SIGMAhat[i,1])}) 
max(MAPESigmanew)

CovFunct <- function(pars){
  if(pars[1]<= pars[2] && pars[3] >= pars[2]){
  #if(pars[1]-0.001<= pars[2] && pars[3]+0.001 >= pars[2]){ # To avoid problems due to numerical approximation in cor. coef. = 1 in corr matrix
    Cov <- 1
  }else{
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],THETApop[i],PostResults$THETApost[i,3,]), 1, CovFunct))})
CoverageSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],SIGMApop[i],PostResults$SIGMApost[i,3,]), 1, CovFunct))})
SIGMAhatInf <- sapply(1:100, function(l){vech(diag(1/diag(xpnd(PostResults$SIGMApostNOst[,2,l]))^0.5)%*%xpnd(PostResults$SIGMApostNOst[,2,l])%*%diag(1/diag(xpnd(PostResults$SIGMApostNOst[,2,l]))^0.5))})
SIGMAhatSup <- sapply(1:100, function(l){vech(diag(1/diag(xpnd(PostResults$SIGMApostNOst[,3,l]))^0.5)%*%xpnd(PostResults$SIGMApostNOst[,3,l])%*%diag(1/diag(xpnd(PostResults$SIGMApostNOst[,3,l]))^0.5))})

CoverageSigmanew <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(SIGMAhatInf[i,],SIGMApop[i],SIGMAhatSup[i,]), 1, CovFunct))})

LengthFunct <- function(pars){
  leng <- abs(pars[2]-pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],PostResults$THETApost[i,3,]), 1, LengthFunct))})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],PostResults$SIGMApost[i,3,]), 1, LengthFunct))})
IntLengthSigmanew <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(SIGMAhatInf[i,],SIGMAhatSup[i,]), 1, LengthFunct))})

####################################################################
rm(list = ls())
load("PostResultsV2CnewV2.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
a2 <- c(0.8, -1.2)
a3 <- c(1.1, -0.7)
J <- 3
rho <- 1
SIGMA <- rho*matrix(c(1,0.6,0.4,0.5,0.4,0.2,0,0,0,0.6,1,0.5,0.4,0.3,0.4,0,0,0,0.4,0.5,1,
                      0.5,0.3,0.5,0,0,0,0.5,0.4,0.5,1,0.4,0.5,0.3,0.4,0.2,0.4,0.3,0.3,0.4,1,
                      0.4,0.3,0.3,0.1,0.2,0.4,0.5,0.5,0.4,1,0.5,0.3,0.1,0,0,0,0.3,0.3,0.5,1,0.6,
                      0.4,0,0,0,0.4,0.3,0.3,0.6,1,0.3,0,0,0,0.2,0.1,0.1,0.4,0.3,1), 3*J, 3*J)
SIGMA13 <- matrix(c(0.3, 0.2, 0.3, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2), 3, 3) 
SIGMA[1:3,7:9] <- SIGMA13
SIGMA[7:9,1:3] <- t(SIGMA13)
SIGMApop <- vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)
b2 <- c(0.5, 1.5, -1)
b3 <- c(1, 1, -1)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)
d2 <- c(1, 2, 2)
d3 <- c(1, 1.5, 1.8)

ALPHApop <- c(a1, a2, a3)
THETApop <- c(b1, b2, b3, d1, d2, d3)
ALPHAhat <- rowMeans(PostResults$ALPHApost[,1,])
THETAhat <- rowMeans(PostResults$THETApost[,1,])
cbind(ALPHApop, ALPHAhat)
cbind(THETApop, THETAhat)

RMSEfunct <- function(pop, pars){
  RMSE <- (mean(sapply(1:length(pars), function(i){(pop-pars[i])^2})))^0.5
  return(RMSE)
}

RMSEAlpha1 <- sapply(1:length(ALPHApop), function(i){RMSEfunct(ALPHApop[i],PostResults$ALPHApost[i,1,])})
RMSETheta1 <- sapply(1:length(THETApop), function(i){RMSEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
RMSE1 <- cbind(RMSETheta, c(RMSEAlpha1, RMSETheta1))

MAPEfunct <- function(pop, pars){
  MAPE <- mean(sapply(1:length(pars), function(i){abs((pop-pars[i])/pop)}))
  return(MAPE)
}
MAPETheta1 <- sapply(1:length(THETApop), function(i){MAPEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
MAPEAlpha1 <- sapply(1:length(ALPHApop), function(i){MAPEfunct(ALPHApop[i],PostResults$ALPHApost[i,1,])}) 
MAPE1 <- cbind(MAPETheta, c(MAPEAlpha1, MAPETheta1))

CovFunct <- function(pars){
  if(pars[1]<= pars[2] && pars[3] >= pars[2]){
    #if(pars[1]-0.001<= pars[2] && pars[3]+0.001 >= pars[2]){ # To avoid problems due to numerical approximation in cor. coef. = 1 in corr matrix
    Cov <- 1
  }else{
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta1 <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],THETApop[i],PostResults$THETApost[i,3,]), 1, CovFunct))})
CoverageAlpha1 <- sapply(1:length(ALPHApop), function(i){mean(apply(cbind(PostResults$ALPHApost[i,2,],ALPHApop[i],PostResults$ALPHApost[i,3,]), 1, CovFunct))})
Cov1 <- cbind(CoverageTheta, c(CoverageAlpha1, CoverageTheta1))

LengthFunct <- function(pars){
  leng <- abs(pars[2]-pars[1])
  return(leng)
}
IntLengthTheta1 <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],PostResults$THETApost[i,3,]), 1, LengthFunct))})
IntLengthAlpha1 <- sapply(1:length(ALPHApop), function(i){mean(apply(cbind(PostResults$ALPHApost[i,2,],PostResults$ALPHApost[i,3,]), 1, LengthFunct))})
RatioLength1 <- c(IntLengthTheta1, IntLengthAlpha1)/IntLengthTheta
cbind(RMSE1, MAPE1, Cov1, RatioLength1)

####################################################################
rm(list = ls())
load("PostResultsV1CnewV1.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho*matrix(c(1,0.5,0.0,0.5, 1, 0.3, 0.0, 0.3, 1), 3*J, 3*J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[,1,]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[,1,]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars){
  RMSE <- (mean(sapply(1:length(pars), function(i){(pop-pars[i])^2})))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i){RMSEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
RMSESigma <- sapply(1:length(SIGMApop), function(i){RMSEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 

MAPEfunct <- function(pop, pars){
  MAPE <- mean(sapply(1:length(pars), function(i){abs((pop-pars[i])/pop)}))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i){MAPEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i){MAPEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 
max(MAPESigma)

CovFunct <- function(pars){
  if(pars[1]<= pars[2] & pars[3] >= pars[2]){
    Cov <- 1
  }else{
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],THETApop[i],PostResults$THETApost[i,3,]), 1, CovFunct))})
CoverageSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],SIGMApop[i],PostResults$SIGMApost[i,3,]), 1, CovFunct))})

LengthFunct <- function(pars){
  leng <- abs(pars[2]-pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],PostResults$THETApost[i,3,]), 1, LengthFunct))})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],PostResults$SIGMApost[i,3,]), 1, LengthFunct))})

####################################################################
rm(list = ls())
load("PostResultsV1B.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho*matrix(c(1,0.6,0.4,0.6, 1, 0.7, 0.4, 0.7, 1), 3*J, 3*J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[,1,]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[,1,]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars){
  RMSE <- (mean(sapply(1:length(pars), function(i){(pop-pars[i])^2})))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i){RMSEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
RMSESigma <- sapply(1:length(SIGMApop), function(i){RMSEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 

MAPEfunct <- function(pop, pars){
  MAPE <- mean(sapply(1:length(pars), function(i){abs((pop-pars[i])/pop)}))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i){MAPEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i){MAPEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 
max(MAPESigma)

CovFunct <- function(pars){
  if(pars[1]<= pars[2] & pars[3] >= pars[2]){
    Cov <- 1
  }else{
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],THETApop[i],PostResults$THETApost[i,3,]), 1, CovFunct))})
CoverageSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],SIGMApop[i],PostResults$SIGMApost[i,3,]), 1, CovFunct))})

LengthFunct <- function(pars){
  leng <- abs(pars[2]-pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],PostResults$THETApost[i,3,]), 1, LengthFunct))})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],PostResults$SIGMApost[i,3,]), 1, LengthFunct))})


####################################################################
rm(list = ls())
load("PostResultsV1BNOexclusion.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho*matrix(c(1,0.6,0.4,0.6, 1, 0.7, 0.4, 0.7, 1), 3*J, 3*J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[,1,]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[,1,]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars){
  RMSE <- (mean(sapply(1:length(pars), function(i){(pop-pars[i])^2})))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i){RMSEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
RMSESigma <- sapply(1:length(SIGMApop), function(i){RMSEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 

MAPEfunct <- function(pop, pars){
  MAPE <- mean(sapply(1:length(pars), function(i){abs((pop-pars[i])/pop)}))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i){MAPEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i){MAPEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 
max(MAPESigma)

CovFunct <- function(pars){
  if(pars[1]<= pars[2] & pars[3] >= pars[2]){
    Cov <- 1
  }else{
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],THETApop[i],PostResults$THETApost[i,3,]), 1, CovFunct))})
CoverageSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],SIGMApop[i],PostResults$SIGMApost[i,3,]), 1, CovFunct))})

LengthFunct <- function(pars){
  leng <- abs(pars[2]-pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],PostResults$THETApost[i,3,]), 1, LengthFunct))})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],PostResults$SIGMApost[i,3,]), 1, LengthFunct))})

####################################################################
rm(list = ls())
load("PostResultsV1BNew.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho*matrix(c(1,0.6,0.4,0.6, 1, 0.7, 0.4, 0.7, 1), 3*J, 3*J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[,1,]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[,1,]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars){
  RMSE <- (mean(sapply(1:length(pars), function(i){(pop-pars[i])^2})))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i){RMSEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
RMSESigma <- sapply(1:length(SIGMApop), function(i){RMSEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 

MAPEfunct <- function(pop, pars){
  MAPE <- mean(sapply(1:length(pars), function(i){abs((pop-pars[i])/pop)}))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i){MAPEfunct(THETApop[i],PostResults$THETApost[i,1,])}) 
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i){MAPEfunct(SIGMApop[i],PostResults$SIGMApost[i,1,])}) 
max(MAPESigma)

CovFunct <- function(pars){
  if(pars[1]<= pars[2] & pars[3] >= pars[2]){
    Cov <- 1
  }else{
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],THETApop[i],PostResults$THETApost[i,3,]), 1, CovFunct))})
CoverageSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],SIGMApop[i],PostResults$SIGMApost[i,3,]), 1, CovFunct))})

LengthFunct <- function(pars){
  leng <- abs(pars[2]-pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i){mean(apply(cbind(PostResults$THETApost[i,2,],PostResults$THETApost[i,3,]), 1, LengthFunct))})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i){mean(apply(cbind(PostResults$SIGMApost[i,2,],PostResults$SIGMApost[i,3,]), 1, LengthFunct))})

