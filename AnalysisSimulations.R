### Analysis simulation ###
rm(list = ls())
library(MCMCpack)
library(matrixcalc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
library(latex2exp)


##### Uni-product 3 stages ####
#1) Baseline
load("PostResultsOneProduct3Stages.RData")
J <- 1
a1 <- c(1, -1, 1)
b1 <- c(1, -0.5, -1.2)
d1 <- c(1, -0.5, 1)

rho <- 1
SIGMA <- rho*matrix(c(1,0.7,0.6,0.7,1, 0.8, 0.6, 0.8, 1), 3*J, 3*J)
SIGMApop <- vech(SIGMA)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhat))
dist(rbind(THETApop, rowMeans(THETAhat)))
# SIGMAhat <- PostResults$SIGMApost[,1,]
SIGMAhat <- sapply(1:100, function(l) {
  vech(diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 1, l]))^0.5) %*% xpnd(PostResults$SIGMApostNOst[, 1, l]) %*% diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 1, l]))^0.5))
})
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSESigma <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
RMSESigmanew <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], SIGMAhat[i, 1])
})
cbind(RMSESigma, RMSESigmanew)
plot(RMSESigma, type = "l")
lines(RMSESigmanew, col = "red") # It's better to standarized at the end!!!

par(mar = c(1, 1, 1, 1))
par(mfrow = c(1, 2))
boxplot(RMSETheta[-c(1, 3, 5, 7, 10, 13, 16, 19, 22)])
boxplot(RMSETheta)

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
max(MAPESigma)
MAPESigmanew <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], SIGMAhat[i, 1])
})
max(MAPESigmanew)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] && pars[3] >= pars[2]) {
    # if(pars[1]-0.001<= pars[2] && pars[3]+0.001 >= pars[2]){ # To avoid problems due to numerical approximation in cor. coef. = 1 in corr matrix
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], SIGMApop[i], PostResults$SIGMApost[i, 3, ]), 1, CovFunct))
})
SIGMAhatInf <- sapply(1:100, function(l) {
  vech(diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 2, l]))^0.5) %*% xpnd(PostResults$SIGMApostNOst[, 2, l]) %*% diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 2, l]))^0.5))
})
SIGMAhatSup <- sapply(1:100, function(l) {
  vech(diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 3, l]))^0.5) %*% xpnd(PostResults$SIGMApostNOst[, 3, l]) %*% diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 3, l]))^0.5))
})

CoverageSigmanew <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(SIGMAhatInf[i, ], SIGMApop[i], SIGMAhatSup[i, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthThetaBP <- sapply(1:length(THETApop), function(i) {
  apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct)
})

IntLengthTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], PostResults$SIGMApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigmanew <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(SIGMAhatInf[i, ], SIGMAhatSup[i, ]), 1, LengthFunct))
})

#2) Independent stages (no access and no selection)
THETAhatIND <- PostResults$THETApostUnivar[, 1, ]
cbind(THETApop, rowMeans(THETAhatIND))
dist(rbind(THETApop, rowMeans(THETAhatIND)))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSEThetaIND <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApostUnivar[i, 1, ])
})


MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPEThetaIND <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApostUnivar[i, 1, ])
})

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] && pars[3] >= pars[2]) {
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageThetaIND <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApostUnivar[i, 2, ], THETApop[i], PostResults$THETApostUnivar[i, 3, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthThetaBPIND <- sapply(1:length(THETApop), function(i) {
  apply(cbind(PostResults$THETApostUnivar[i, 2, ], PostResults$THETApostUnivar[i, 3, ]), 1, LengthFunct)
})

IntLengthThetaIND <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApostUnivar[i, 2, ], PostResults$THETApostUnivar[i, 3, ]), 1, LengthFunct))
})

#3) Independent stages (no access and no selection)
load("PostResultsOneProduct3StagesNOexcrest.RData")
THETAhatNOexc <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhatNOexc))
dist(rbind(THETApop, rowMeans(THETAhatNOexc)))

RMSEThetaNOexc <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})


MAPEThetaNOexc <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})

CoverageThetaNOexc <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})

IntLengthThetaBPNOexc <- sapply(1:length(THETApop), function(i) {
  apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct)
})

IntLengthThetaNOexc <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})

#4) Exogenous access
load("PostResultsOneProduct3StagesExAccess.RData")
THETAhatExAc <- rbind(PostResults$ALPHApost[, 1, ], PostResults$THETApost[, 1, ])
cbind(THETApop, rowMeans(THETAhatExAc))
dist(rbind(THETApop, rowMeans(THETAhatExAc)))

RMSEThetaExAc <- sapply(1:length(THETApop[4:9]), function(i) {
  RMSEfunct(THETApop[3+i], PostResults$THETApost[i, 1, ])
})

MAPEAlphaExAc <- sapply(1:length(THETApop[1:3]), function(i) {
  MAPEfunct(THETApop[i], PostResults$ALPHApost[i, 1, ])
})

RMSEThetaExAc <- c(MAPEAlphaExAc, RMSEThetaExAc)

MAPEThetaExAc <- sapply(1:length(THETApop[4:9]), function(i) {
  MAPEfunct(THETApop[3+i], PostResults$THETApost[i, 1, ])
})

MAPEAlphaExAc <- sapply(1:length(THETApop[1:3]), function(i) {
  MAPEfunct(THETApop[i], PostResults$ALPHApost[i, 1, ])
})

MAPEThetaExAc <- c(MAPEAlphaExAc, MAPEThetaExAc)

CoverageThetaExAc <- sapply(1:length(THETApop[4:9]), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[3+i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})

CoverageAlphaExAc <- sapply(1:length(THETApop[1:3]), function(i) {
  mean(apply(cbind(PostResults$ALPHApost[i, 2, ], THETApop[i], PostResults$ALPHApost[i, 3, ]), 1, CovFunct))
})

CoverageThetaExAc <- c(CoverageAlphaExAc, CoverageThetaExAc)

IntLengthThetaBPExAc <- sapply(1:length(THETApop[4:9]), function(i) {
  apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct)
})

IntLengthAlphaBPExAc <- sapply(1:length(THETApop[1:3]), function(i) {
  apply(cbind(PostResults$ALPHApost[i, 2, ], PostResults$ALPHApost[i, 3, ]), 1, LengthFunct)
})

IntLengthThetaBPExAc <- cbind(IntLengthAlphaBPExAc, IntLengthThetaBPExAc)

IntLengthThetaExAc <- sapply(1:length(THETApop[4:9]), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})

IntLengthAlphaExAc <- sapply(1:length(THETApop[1:3]), function(i) {
  mean(apply(cbind(PostResults$ALPHApost[i, 2, ], PostResults$ALPHApost[i, 3, ]), 1, LengthFunct))
})

IntLengthThetaExAc <- c(IntLengthAlphaExAc, IntLengthThetaExAc)

#### Summary ######
cbind( RMSEThetaExAc/RMSETheta, RMSEThetaNOexc/RMSETheta, RMSEThetaIND/RMSETheta)
cbind(MAPEThetaExAc/MAPETheta, MAPEThetaNOexc/MAPETheta, MAPEThetaIND/MAPETheta)
cbind(CoverageTheta, CoverageThetaExAc, CoverageThetaNOexc, CoverageThetaIND)
cbind(IntLengthThetaExAc/IntLengthTheta, IntLengthThetaNOexc/IntLengthTheta, IntLengthThetaIND/IntLengthTheta)

################### Box plots: Coefficients Uni-product #################
S <- 100
a12Base <- THETAhat[2, ] # Baseline
a12ExAc <- THETAhatExAc[2, ] # ALPHAhat1[2, ] # Exogenous access
a12NoEx <- THETAhatNOexc[2, ] # No exclusion
a12Unv <- THETAhatIND[2, ] # Univariate

BoxPlota12 <- data.frame(c(a12Base, a12ExAc, a12NoEx, a12Unv))
BoxPlota12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota12) <- c("Coefficient", "Model")
BoxPlota12 <- BoxPlota12 %>%
  relocate(Model)

means <- BoxPlota12 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot1 <- ggplot(BoxPlota12, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a1[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab(expression(italic(alpha)[12]))
BoxPlot1

a13Base <- THETAhat[3, ] # Baseline
a13ExAc <- THETAhatExAc[3, ] # ALPHAhat1[2, ] # Exogenous access
a13NoEx <- THETAhatNOexc[3, ] # No exclusion
a13Unv <- THETAhatIND[3, ] # Univariate

BoxPlota13 <- data.frame(c(a13Base, a13ExAc, a13NoEx, a13Unv))
BoxPlota13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota13) <- c("Coefficient", "Model")
BoxPlota13 <- BoxPlota13 %>%
  relocate(Model)

means <- BoxPlota13 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot2 <- ggplot(BoxPlota13, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a1[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab(expression(italic(alpha)[13]))
BoxPlot2

b12Base <- THETAhat[5, ] # Baseline
b12ExAc <- THETAhatExAc[5, ] # ALPHAhat1[2, ] # Exogenous access
b12NoEx <- THETAhatNOexc[5, ] # No exclusion
b12Unv <- THETAhatIND[5, ] # Univariate

BoxPlotb12 <- data.frame(c(b12Base, b12ExAc, b12NoEx, b12Unv))
BoxPlotb12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb12) <- c("Coefficient", "Model")
BoxPlotb12 <- BoxPlotb12 %>%
  relocate(Model)

means <- BoxPlotb12 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot3 <- ggplot(BoxPlotb12, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b1[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab(expression(italic(delta)[12]))
BoxPlot3


b13Base <- THETAhat[6, ] # Baseline
b13ExAc <- THETAhatExAc[6, ] # ALPHAhat1[2, ] # Exogenous access
b13NoEx <- THETAhatNOexc[6, ] # No exclusion
b13Unv <- THETAhatIND[6, ] # Univariate

BoxPlotb13 <- data.frame(c(b13Base, b13ExAc, b13NoEx, b13Unv))
BoxPlotb13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb13) <- c("Coefficient", "Model")
BoxPlotb13 <- BoxPlotb13 %>%
  relocate(Model)

means <- BoxPlotb13 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot4 <- ggplot(BoxPlotb13, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b1[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab(expression(italic(delta)[13]))
BoxPlot4


d12Base <- THETAhat[8, ] # Baseline
d12ExAc <- THETAhatExAc[8, ] # ALPHAhat1[2, ] # Exogenous access
d12NoEx <- THETAhatNOexc[8, ] # No exclusion
d12Unv <- THETAhatIND[8, ] # Univariate

BoxPlotd12 <- data.frame(c(d12Base, d12ExAc, d12NoEx, d12Unv))
BoxPlotd12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd12) <- c("Coefficient", "Model")
BoxPlotd12 <- BoxPlotd12 %>%
  relocate(Model)

means <- BoxPlotd12 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot5 <- ggplot(BoxPlotd12, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d1[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab(expression(italic(beta)[12]))
BoxPlot5


d13Base <- THETAhat[9, ] # Baseline
d13ExAc <- THETAhatExAc[9, ] # ALPHAhat1[2, ] # Exogenous access
d13NoEx <- THETAhatNOexc[9, ] # No exclusion
d13Unv <- THETAhatIND[9, ] # Univariate

BoxPlotd13 <- data.frame(c(d13Base, d13ExAc, d13NoEx, d13Unv))
BoxPlotd13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd13) <- c("Coefficient", "Model")
BoxPlotd13 <- BoxPlotd13 %>%
  relocate(Model)

means <- BoxPlotd13 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot6 <- ggplot(BoxPlotd13, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d1[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab(expression(italic(beta)[13]))
BoxPlot6

ggarrange(BoxPlot1, BoxPlot2, BoxPlot3,
          BoxPlot4, BoxPlot5, BoxPlot6,
          labels = c(
            "A", "B", "C", "D", "E", "F"
          ),
          ncol = 2, nrow = 3,
          legend = "bottom",
          common.legend = TRUE
)

ggarrange(BoxPlot2, BoxPlot4, BoxPlot6,
          labels = c(
            "A", "B", "C"
          ),
          ncol = 3, nrow = 1,
          legend = "bottom",
          common.legend = TRUE
)


a12Base <- IntLengthThetaBP[, 2] # Baseline
a12ExAc <- IntLengthThetaBP[, 2] # ALPHAhat1[2, ] # Exogenous access
a12NoEx <- IntLengthThetaBPNOexc[, 2] # No exclusion
a12Unv <- IntLengthThetaBPIND[, 2] # Univariate

BoxPlota12 <- data.frame(c(a12Base, a12ExAc, a12NoEx, a12Unv))
BoxPlota12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota12) <- c("Interval_length", "Model")
BoxPlota12 <- BoxPlota12 %>%
  relocate(Model)

BoxPlot7 <- ggplot(BoxPlota12, aes(x = Model, y = Interval_length)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Interval length")
BoxPlot7

a13Base <- IntLengthThetaBP[, 3] # Baseline
a13ExAc <- IntLengthThetaBP[, 3] # ALPHAhat1[2, ] # Exogenous access
a13NoEx <- IntLengthThetaBPNOexc[, 3] # No exclusion
a13Unv <- IntLengthThetaBPIND[, 3] # Univariate

BoxPlota13 <- data.frame(c(a13Base, a13ExAc, a13NoEx, a13Unv))
BoxPlota13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota13) <- c("Interval_length", "Model")
BoxPlota13 <- BoxPlota13 %>%
  relocate(Model)

BoxPlot8 <- ggplot(BoxPlota13, aes(x = Model, y = Interval_length)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Interval length")
BoxPlot8

b12Base <- IntLengthThetaBP[, 4] # Baseline
b12ExAc <- IntLengthThetaBP[, 4] # ALPHAhat1[2, ] # Exogenous access
b12NoEx <- IntLengthThetaBPNOexc[, 4] # No exclusion
b12Unv <- IntLengthThetaBPIND[, 4] # Univariate

BoxPlotb12 <- data.frame(c(b12Base, b12ExAc, b12NoEx, b12Unv))
BoxPlotb12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb12) <- c("Interval_length", "Model")
BoxPlotb12 <- BoxPlotb12 %>%
  relocate(Model)

BoxPlot9 <- ggplot(BoxPlotb12, aes(x = Model, y = Interval_length)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Interval length")
BoxPlot9

b13Base <- IntLengthThetaBP[, 5] # Baseline
b13ExAc <- IntLengthThetaBP[, 5] # ALPHAhat1[2, ] # Exogenous access
b13NoEx <- IntLengthThetaBPNOexc[, 5] # No exclusion
b13Unv <- IntLengthThetaBPIND[, 5] # Univariate

BoxPlotb13 <- data.frame(c(b13Base, b13ExAc, b13NoEx, b13Unv))
BoxPlotb13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb13) <- c("Interval_length", "Model")
BoxPlotb13 <- BoxPlotb13 %>%
  relocate(Model)

BoxPlot10 <- ggplot(BoxPlotb13, aes(x = Model, y = Interval_length)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Interval length")
BoxPlot10

d12Base <- IntLengthThetaBP[, 8] # Baseline
d12ExAc <- IntLengthThetaBP[, 8] # ALPHAhat1[2, ] # Exogenous access
d12NoEx <- IntLengthThetaBPNOexc[, 8] # No exclusion
d12Unv <- IntLengthThetaBPIND[, 8] # Univariate

BoxPlotd12 <- data.frame(c(d12Base, d12ExAc, d12NoEx, d12Unv))
BoxPlotd12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd12) <- c("Interval_length", "Model")
BoxPlotd12 <- BoxPlotd12 %>%
  relocate(Model)

BoxPlot11 <- ggplot(BoxPlotd12, aes(x = Model, y = Interval_length)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Interval length")
BoxPlot11


d13Base <- IntLengthThetaBP[, 9] # Baseline
d13ExAc <- IntLengthThetaBP[, 9] # ALPHAhat1[2, ] # Exogenous access
d13NoEx <- IntLengthThetaBPNOexc[, 9] # No exclusion
d13Unv <- IntLengthThetaBPIND[, 9] # Univariate

BoxPlotd13 <- data.frame(c(d13Base, d13ExAc, d13NoEx, d13Unv))
BoxPlotd13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd13) <- c("Interval_length", "Model")
BoxPlotd13 <- BoxPlotd13 %>%
  relocate(Model)

BoxPlot12 <- ggplot(BoxPlotd13, aes(x = Model, y = Interval_length)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Interval length")
BoxPlot12

ggarrange(BoxPlot8, BoxPlot10, BoxPlot12,
          labels = c(
            "A", "B", "C"
          ),
          ncol = 3, nrow = 1,
          legend = "bottom",
          common.legend = TRUE
)

##### Multivariate probit: Access #####
load("PostResultsV3BnewV2.RData")
a1 <- c(1, -1)
a2 <- c(0.8, -1.2)
a3 <- c(1.1, -0.7)
J <- 3
rho <- 1
SIGMA <- rho * matrix(c(
  1, 0.6, 0.4, 0.5, 0.4, 0.2, 0, 0, 0, 0.6, 1, 0.5, 0.4, 0.3, 0.4, 0, 0, 0, 0.4, 0.5, 1,
  0.5, 0.3, 0.5, 0, 0, 0, 0.5, 0.4, 0.5, 1, 0.4, 0.5, 0.3, 0.4, 0.2, 0.4, 0.3, 0.3, 0.4, 1,
  0.4, 0.3, 0.3, 0.1, 0.2, 0.4, 0.5, 0.5, 0.4, 1, 0.5, 0.3, 0.1, 0, 0, 0, 0.3, 0.3, 0.5, 1, 0.6,
  0.4, 0, 0, 0, 0.4, 0.3, 0.3, 0.6, 1, 0.3, 0, 0, 0, 0.2, 0.1, 0.1, 0.4, 0.3, 1
), 3 * J, 3 * J)
SIGMA13 <- matrix(c(0.3, 0.2, 0.3, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2), 3, 3)
SIGMA[1:3, 7:9] <- SIGMA13
SIGMA[7:9, 1:3] <- t(SIGMA13)
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
THETAhat <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhat))
dist(rbind(THETApop, rowMeans(THETAhat)))
# SIGMAhat <- PostResults$SIGMApost[,1,]
SIGMAhat <- sapply(1:100, function(l) {
  vech(diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 1, l]))^0.5) %*% xpnd(PostResults$SIGMApostNOst[, 1, l]) %*% diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 1, l]))^0.5))
})
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSESigma <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
RMSESigmanew <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], SIGMAhat[i, 1])
})
cbind(RMSESigma, RMSESigmanew)
plot(RMSESigma, type = "l")
lines(RMSESigmanew, col = "red") # It's better to standarized at the end!!!

par(mar = c(1, 1, 1, 1))
par(mfrow = c(1, 2))
boxplot(RMSETheta[-c(1, 3, 5, 7, 10, 13, 16, 19, 22)])
boxplot(RMSETheta)

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
max(MAPESigma)
MAPESigmanew <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], SIGMAhat[i, 1])
})
max(MAPESigmanew)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] && pars[3] >= pars[2]) {
    # if(pars[1]-0.001<= pars[2] && pars[3]+0.001 >= pars[2]){ # To avoid problems due to numerical approximation in cor. coef. = 1 in corr matrix
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], SIGMApop[i], PostResults$SIGMApost[i, 3, ]), 1, CovFunct))
})
SIGMAhatInf <- sapply(1:100, function(l) {
  vech(diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 2, l]))^0.5) %*% xpnd(PostResults$SIGMApostNOst[, 2, l]) %*% diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 2, l]))^0.5))
})
SIGMAhatSup <- sapply(1:100, function(l) {
  vech(diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 3, l]))^0.5) %*% xpnd(PostResults$SIGMApostNOst[, 3, l]) %*% diag(1 / diag(xpnd(PostResults$SIGMApostNOst[, 3, l]))^0.5))
})

CoverageSigmanew <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(SIGMAhatInf[i, ], SIGMApop[i], SIGMAhatSup[i, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthThetaBP <- sapply(1:length(THETApop), function(i) {
  apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct)
})

IntLengthTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], PostResults$SIGMApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigmanew <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(SIGMAhatInf[i, ], SIGMAhatSup[i, ]), 1, LengthFunct))
})

####################################################################
# rm(list = ls())
load("PostResultsV2CnewV2.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
a2 <- c(0.8, -1.2)
a3 <- c(1.1, -0.7)
J <- 3
rho <- 1
SIGMA <- rho * matrix(c(
  1, 0.6, 0.4, 0.5, 0.4, 0.2, 0, 0, 0, 0.6, 1, 0.5, 0.4, 0.3, 0.4, 0, 0, 0, 0.4, 0.5, 1,
  0.5, 0.3, 0.5, 0, 0, 0, 0.5, 0.4, 0.5, 1, 0.4, 0.5, 0.3, 0.4, 0.2, 0.4, 0.3, 0.3, 0.4, 1,
  0.4, 0.3, 0.3, 0.1, 0.2, 0.4, 0.5, 0.5, 0.4, 1, 0.5, 0.3, 0.1, 0, 0, 0, 0.3, 0.3, 0.5, 1, 0.6,
  0.4, 0, 0, 0, 0.4, 0.3, 0.3, 0.6, 1, 0.3, 0, 0, 0, 0.2, 0.1, 0.1, 0.4, 0.3, 1
), 3 * J, 3 * J)
SIGMA13 <- matrix(c(0.3, 0.2, 0.3, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2), 3, 3)
SIGMA[1:3, 7:9] <- SIGMA13
SIGMA[7:9, 1:3] <- t(SIGMA13)
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
ALPHAhat1 <- PostResults$ALPHApost[, 1, ]
THETAhat1 <- PostResults$THETApost[, 1, ]
ALPHAhat1m <- rowMeans(PostResults$ALPHApost[, 1, ])
THETAhat1m <- rowMeans(PostResults$THETApost[, 1, ])
cbind(ALPHApop, ALPHAhat1m)
cbind(THETApop, THETAhat1m)
dist(rbind(c(ALPHApop, THETApop), c(ALPHAhat1m, THETAhat1m)))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}

RMSEAlpha1 <- sapply(1:length(ALPHApop), function(i) {
  RMSEfunct(ALPHApop[i], PostResults$ALPHApost[i, 1, ])
})
RMSETheta1 <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSE1 <- cbind(RMSETheta, c(RMSEAlpha1, RMSETheta1))
colMeans(RMSE1)
boxplot(c(RMSEAlpha1[-c(1, 3, 5)], RMSETheta1[-c(1, 4, 7, 10, 13, 16)]))
boxplot(c(RMSEAlpha1, RMSETheta1))
mean(c(RMSEAlpha1, RMSETheta1))

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta1 <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
MAPEAlpha1 <- sapply(1:length(ALPHApop), function(i) {
  MAPEfunct(ALPHApop[i], PostResults$ALPHApost[i, 1, ])
})
MAPE1 <- cbind(MAPETheta, c(MAPEAlpha1, MAPETheta1))
colMeans(MAPE1)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] && pars[3] >= pars[2]) {
    # if(pars[1]-0.001<= pars[2] && pars[3]+0.001 >= pars[2]){ # To avoid problems due to numerical approximation in cor. coef. = 1 in corr matrix
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta1 <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageAlpha1 <- sapply(1:length(ALPHApop), function(i) {
  mean(apply(cbind(PostResults$ALPHApost[i, 2, ], ALPHApop[i], PostResults$ALPHApost[i, 3, ]), 1, CovFunct))
})
Cov1 <- cbind(CoverageTheta, c(CoverageAlpha1, CoverageTheta1))
colMeans(Cov1)

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthTheta1BP <- sapply(1:length(THETApop), function(i) {
  apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct)
})
IntLengthAlpha1BP <- sapply(1:length(ALPHApop), function(i) {
  apply(cbind(PostResults$ALPHApost[i, 2, ], PostResults$ALPHApost[i, 3, ]), 1, LengthFunct)
})

IntLengthTheta1 <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthAlpha1 <- sapply(1:length(ALPHApop), function(i) {
  mean(apply(cbind(PostResults$ALPHApost[i, 2, ], PostResults$ALPHApost[i, 3, ]), 1, LengthFunct))
})
RatioLength1 <- c(IntLengthTheta1, IntLengthAlpha1) / IntLengthTheta
mean(RatioLength1)
cbind(RMSE1, MAPE1, Cov1, RatioLength1)

####################################################################
rm(list = ls())
load("PostResultsV1CnewV1.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho * matrix(c(1, 0.5, 0.0, 0.5, 1, 0.3, 0.0, 0.3, 1), 3 * J, 3 * J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[, 1, ]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSESigma <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
max(MAPESigma)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] & pars[3] >= pars[2]) {
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], SIGMApop[i], PostResults$SIGMApost[i, 3, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], PostResults$SIGMApost[i, 3, ]), 1, LengthFunct))
})

####################################################################
rm(list = ls())
load("PostResultsV1B.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho * matrix(c(1, 0.6, 0.4, 0.6, 1, 0.7, 0.4, 0.7, 1), 3 * J, 3 * J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[, 1, ]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSESigma <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
max(MAPESigma)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] & pars[3] >= pars[2]) {
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], SIGMApop[i], PostResults$SIGMApost[i, 3, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], PostResults$SIGMApost[i, 3, ]), 1, LengthFunct))
})


####################################################################
rm(list = ls())
load("PostResultsV1BNOexclusion.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho * matrix(c(1, 0.6, 0.4, 0.6, 1, 0.7, 0.4, 0.7, 1), 3 * J, 3 * J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[, 1, ]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSESigma <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
max(MAPESigma)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] & pars[3] >= pars[2]) {
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], SIGMApop[i], PostResults$SIGMApost[i, 3, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], PostResults$SIGMApost[i, 3, ]), 1, LengthFunct))
})
####################################################################
rm(list = ls())
load("PostResultsV1BNew.RData")
##### Multivariate probit: Access #####
a1 <- c(1, -1)
J <- 1
rho <- 1
SIGMA <- rho * matrix(c(1, 0.6, 0.4, 0.6, 1, 0.7, 0.4, 0.7, 1), 3 * J, 3 * J)
SIGMApop <- matrixcalc::vech(SIGMA)
##### Multivariate probit: Selection #####
b1 <- c(1, 1, -0.5)

##### SUR: Outcome #####
d1 <- c(1, 1.7, 1.5)

THETApop <- c(a1, b1, d1)
THETAhat <- PostResults$THETApost[, 1, ]
cbind(THETApop, rowMeans(THETAhat))
SIGMAhat <- PostResults$SIGMApost[, 1, ]
cbind(SIGMApop, rowMeans(SIGMAhat))

RMSEfunct <- function(pop, pars) {
  RMSE <- (mean(sapply(1:length(pars), function(i) {
    (pop - pars[i])^2
  })))^0.5
  return(RMSE)
}
RMSETheta <- sapply(1:length(THETApop), function(i) {
  RMSEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
RMSESigma <- sapply(1:length(SIGMApop), function(i) {
  RMSEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})

MAPEfunct <- function(pop, pars) {
  MAPE <- mean(sapply(1:length(pars), function(i) {
    abs((pop - pars[i]) / pop)
  }))
  return(MAPE)
}
MAPETheta <- sapply(1:length(THETApop), function(i) {
  MAPEfunct(THETApop[i], PostResults$THETApost[i, 1, ])
})
max(MAPETheta)
MAPESigma <- sapply(1:length(SIGMApop), function(i) {
  MAPEfunct(SIGMApop[i], PostResults$SIGMApost[i, 1, ])
})
max(MAPESigma)

CovFunct <- function(pars) {
  if (pars[1] <= pars[2] & pars[3] >= pars[2]) {
    Cov <- 1
  } else {
    Cov <- 0
  }
  return(Cov)
}
CoverageTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], THETApop[i], PostResults$THETApost[i, 3, ]), 1, CovFunct))
})
CoverageSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], SIGMApop[i], PostResults$SIGMApost[i, 3, ]), 1, CovFunct))
})

LengthFunct <- function(pars) {
  leng <- abs(pars[2] - pars[1])
  return(leng)
}
IntLengthTheta <- sapply(1:length(THETApop), function(i) {
  mean(apply(cbind(PostResults$THETApost[i, 2, ], PostResults$THETApost[i, 3, ]), 1, LengthFunct))
})
IntLengthSigma <- sapply(1:length(SIGMApop), function(i) {
  mean(apply(cbind(PostResults$SIGMApost[i, 2, ], PostResults$SIGMApost[i, 3, ]), 1, LengthFunct))
})

################### Box plots: Coefficients #################
S <- 100
a12Base <- THETAhat[2, ] # Baseline
a12ExAc <- ALPHAhat1[2, ] # Exogenous access
a12NoEx <- THETAhat[2, ] # No exclusion
a12Unv <- ALPHAhat1[2, ] # Univariate

BoxPlota12 <- data.frame(c(a12Base, a12ExAc, a12NoEx, a12Unv))
BoxPlota12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota12) <- c("Coefficient", "Model")
BoxPlota12 <- BoxPlota12 %>%
  relocate(Model)

means <- BoxPlota12 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot1 <- ggplot(BoxPlota12, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a1[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot1

a22Base <- THETAhat[4, ] # Baseline
a22ExAc <- ALPHAhat1[4, ] # Exogenous access
a22NoEx <- THETAhat[4, ] # No exclusion
a22Unv <- ALPHAhat1[4, ] # Univariate

BoxPlota22 <- data.frame(c(a22Base, a22ExAc, a22NoEx, a22Unv))
BoxPlota22$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota22) <- c("Coefficient", "Model")
BoxPlota22 <- BoxPlota22 %>%
  relocate(Model)

means <- BoxPlota22 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot2 <- ggplot(BoxPlota22, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a2[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot2

a32Base <- THETAhat[6, ] # Baseline
a32ExAc <- ALPHAhat1[6, ] # Exogenous access
a32NoEx <- THETAhat[6, ] # No exclusion
a32Unv <- ALPHAhat1[6, ] # Univariate

BoxPlota32 <- data.frame(c(a32Base, a32ExAc, a32NoEx, a32Unv))
BoxPlota32$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota32) <- c("Coefficient", "Model")
BoxPlota32 <- BoxPlota32 %>%
  relocate(Model)

means <- BoxPlota32 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot3 <- ggplot(BoxPlota32, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a3[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot3

##################
b12Base <- THETAhat[8, ] # Baseline
b12ExAc <- THETAhat1[2, ] # Exogenous access
b12NoEx <- THETAhat[8, ] # No exclusion
b12Unv <- THETAhat1[2, ] # Univariate

BoxPlotb12 <- data.frame(c(b12Base, b12ExAc, b12NoEx, b12Unv))
BoxPlotb12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb12) <- c("Coefficient", "Model")
BoxPlotb12 <- BoxPlotb12 %>%
  relocate(Model)

BoxPlot4 <- ggplot(BoxPlotb12, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b1[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot4

b22Base <- THETAhat[11, ] # Baseline
b22ExAc <- THETAhat1[5, ] # Exogenous access
b22NoEx <- THETAhat[11, ] # No exclusion
b22Unv <- THETAhat1[5, ] # Univariate

BoxPlotb22 <- data.frame(c(b22Base, b22ExAc, b22NoEx, b22Unv))
BoxPlotb22$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb22) <- c("Coefficient", "Model")
BoxPlotb22 <- BoxPlotb22 %>%
  relocate(Model)

BoxPlot5 <- ggplot(BoxPlotb22, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b2[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot5

b32Base <- THETAhat[14, ] # Baseline
b32ExAc <- THETAhat1[8, ] # Exogenous access
b32NoEx <- THETAhat[14, ] # No exclusion
b32Unv <- THETAhat1[8, ] # Univariate

BoxPlotb32 <- data.frame(c(b32Base, b32ExAc, b32NoEx, b32Unv))
BoxPlotb32$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb32) <- c("Coefficient", "Model")
BoxPlotb32 <- BoxPlotb32 %>%
  relocate(Model)

BoxPlot6 <- ggplot(BoxPlotb32, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b3[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot6

b13Base <- THETAhat[9, ] # Baseline
b13ExAc <- THETAhat1[3, ] # Exogenous access
b13NoEx <- THETAhat[9, ] # No exclusion
b13Unv <- THETAhat1[3, ] # Univariate

BoxPlotb13 <- data.frame(c(b13Base, b13ExAc, b13NoEx, b13Unv))
BoxPlotb13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb13) <- c("Coefficient", "Model")
BoxPlotb13 <- BoxPlotb13 %>%
  relocate(Model)

BoxPlot7 <- ggplot(BoxPlotb13, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b1[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot7

b23Base <- THETAhat[12, ] # Baseline
b23ExAc <- THETAhat1[6, ] # Exogenous access
b23NoEx <- THETAhat[12, ] # No exclusion
b23Unv <- THETAhat1[6, ] # Univariate

BoxPlotb23 <- data.frame(c(b23Base, b23ExAc, b23NoEx, b23Unv))
BoxPlotb23$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb23) <- c("Coefficient", "Model")
BoxPlotb23 <- BoxPlotb23 %>%
  relocate(Model)

BoxPlot8 <- ggplot(BoxPlotb23, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b2[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot8

b33Base <- THETAhat[15, ] # Baseline
b33ExAc <- THETAhat1[9, ] # Exogenous access
b33NoEx <- THETAhat[15, ] # No exclusion
b33Unv <- THETAhat1[9, ] # Univariate

BoxPlotb33 <- data.frame(c(b33Base, b33ExAc, b33NoEx, b33Unv))
BoxPlotb33$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb33) <- c("Coefficient", "Model")
BoxPlotb33 <- BoxPlotb33 %>%
  relocate(Model)

BoxPlot9 <- ggplot(BoxPlotb33, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b3[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot9

b13Base <- THETAhat[9, ] # Baseline
b13ExAc <- THETAhat1[3, ] # Exogenous access
b13NoEx <- THETAhat[9, ] # No exclusion
b13Unv <- THETAhat1[3, ] # Univariate

BoxPlotb13 <- data.frame(c(b13Base, b13ExAc, b13NoEx, b13Unv))
BoxPlotb13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb13) <- c("Coefficient", "Model")
BoxPlotb13 <- BoxPlotb13 %>%
  relocate(Model)

BoxPlot7 <- ggplot(BoxPlotb13, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b1[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot7

b23Base <- THETAhat[12, ] # Baseline
b23ExAc <- THETAhat1[6, ] # Exogenous access
b23NoEx <- THETAhat[12, ] # No exclusion
b23Unv <- THETAhat1[6, ] # Univariate

BoxPlotb23 <- data.frame(c(b23Base, b23ExAc, b23NoEx, b23Unv))
BoxPlotb23$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb23) <- c("Coefficient", "Model")
BoxPlotb23 <- BoxPlotb23 %>%
  relocate(Model)

BoxPlot8 <- ggplot(BoxPlotb23, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b2[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot8

d12Base <- THETAhat[17, ] # Baseline
d12ExAc <- THETAhat1[11, ] # Exogenous access
d12NoEx <- THETAhat[17, ] # No exclusion
d12Unv <- THETAhat1[11, ] # Univariate

BoxPlotd12 <- data.frame(c(d12Base, d12ExAc, d12NoEx, d12Unv))
BoxPlotd12$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd12) <- c("Coefficient", "Model")
BoxPlotd12 <- BoxPlotd12 %>%
  relocate(Model)

BoxPlot10 <- ggplot(BoxPlotd12, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d1[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot10

d22Base <- THETAhat[20, ] # Baseline
d22ExAc <- THETAhat1[14, ] # Exogenous access
d22NoEx <- THETAhat[20, ] # No exclusion
d22Unv <- THETAhat1[14, ] # Univariate

BoxPlotd22 <- data.frame(c(d22Base, d22ExAc, d22NoEx, d22Unv))
BoxPlotd22$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd22) <- c("Coefficient", "Model")
BoxPlotd22 <- BoxPlotd22 %>%
  relocate(Model)

BoxPlot11 <- ggplot(BoxPlotd22, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d2[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot11

d32Base <- THETAhat[23, ] # Baseline
d32ExAc <- THETAhat1[17, ] # Exogenous access
d32NoEx <- THETAhat[23, ] # No exclusion
d32Unv <- THETAhat1[17, ] # Univariate

BoxPlotd32 <- data.frame(c(d32Base, d32ExAc, d32NoEx, d32Unv))
BoxPlotd32$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd32) <- c("Coefficient", "Model")
BoxPlotd32 <- BoxPlotd32 %>%
  relocate(Model)

BoxPlot12 <- ggplot(BoxPlotd32, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d3[2], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot12

d13Base <- THETAhat[18, ] # Baseline
d13ExAc <- THETAhat1[11, ] # Exogenous access
d13NoEx <- THETAhat[18, ] # No exclusion
d13Unv <- THETAhat1[11, ] # Univariate

BoxPlotd13 <- data.frame(c(d13Base, d13ExAc, d13NoEx, d13Unv))
BoxPlotd13$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd13) <- c("Coefficient", "Model")
BoxPlotd13 <- BoxPlotd13 %>%
  relocate(Model)

BoxPlot13 <- ggplot(BoxPlotd13, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d1[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot13

d23Base <- THETAhat[21, ] # Baseline
d23ExAc <- THETAhat1[15, ] # Exogenous access
d23NoEx <- THETAhat[21, ] # No exclusion
d23Unv <- THETAhat1[15, ] # Univariate

BoxPlotd23 <- data.frame(c(d23Base, d23ExAc, d23NoEx, d23Unv))
BoxPlotd23$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd23) <- c("Coefficient", "Model")
BoxPlotd23 <- BoxPlotd23 %>%
  relocate(Model)

BoxPlot14 <- ggplot(BoxPlotd23, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d2[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot14

d33Base <- THETAhat[24, ] # Baseline
d33ExAc <- THETAhat1[18, ] # Exogenous access
d33NoEx <- THETAhat[24, ] # No exclusion
d33Unv <- THETAhat1[18, ] # Univariate

BoxPlotd33 <- data.frame(c(d33Base, d33ExAc, d33NoEx, d33Unv))
BoxPlotd33$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd33) <- c("Coefficient", "Model")
BoxPlotd33 <- BoxPlotd33 %>%
  relocate(Model)

BoxPlot15 <- ggplot(BoxPlotd33, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d3[3], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot15

ggarrange(BoxPlot1, BoxPlot2, BoxPlot3,
  BoxPlot4, BoxPlot5, BoxPlot6,
  BoxPlot7, BoxPlot8, BoxPlot9,
  BoxPlot10, BoxPlot11, BoxPlot12,
  BoxPlot13, BoxPlot14, BoxPlot15,
  labels = c(
    "A", "B", "C", "D", "E", "F",
    "G", "H", "I", "J", "K", "L",
    "M", "N", "O"
  ),
  ncol = 3, nrow = 5,
  legend = "bottom",
  common.legend = TRUE
)

ggarrange(BoxPlot1, BoxPlot2, BoxPlot3,
  labels = c("A", "B", "C"),
  ncol = 3, nrow = 1,
  legend = "bottom",
  common.legend = TRUE
)

ggarrange(BoxPlot4, BoxPlot5, BoxPlot6,
  BoxPlot7, BoxPlot8, BoxPlot9,
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3, nrow = 2,
  legend = "bottom",
  common.legend = TRUE
)

ggarrange(BoxPlot10, BoxPlot11, BoxPlot12,
  BoxPlot13, BoxPlot14, BoxPlot15,
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3, nrow = 2,
  legend = "bottom",
  common.legend = TRUE
)
################# Box plots: Lengths ###############

a12BaseL <- IntLengthThetaBP[, 2] # Baseline
a12ExAcL <- IntLengthAlpha1BP[, 2] # Exogenous access
a12NoExL <- IntLengthThetaBP[, 2] # No exclusion
a12UnvL <- IntLengthAlpha1BP[, 2] # Univariate

BoxPlota12L <- data.frame(c(a12BaseL, a12ExAcL, a12NoExL, a12UnvL))
BoxPlota12L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota12L) <- c("Coefficient", "Model")
BoxPlota12L <- BoxPlota12L %>%
  relocate(Model)

BoxPlot16 <- ggplot(BoxPlota12L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot16

a22BaseL <- IntLengthThetaBP[, 4] # Baseline
a22ExAcL <- IntLengthAlpha1BP[, 4] # Exogenous access
a22NoExL <- IntLengthThetaBP[, 4] # No exclusion
a22UnvL <- IntLengthAlpha1BP[, 4] # Univariate

BoxPlota22L <- data.frame(c(a22BaseL, a22ExAcL, a22NoExL, a22UnvL))
BoxPlota22L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota22L) <- c("Coefficient", "Model")
BoxPlota22L <- BoxPlota22L %>%
  relocate(Model)

BoxPlot17 <- ggplot(BoxPlota22L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot17

a32BaseL <- IntLengthThetaBP[, 6] # Baseline
a32ExAcL <- IntLengthAlpha1BP[, 6] # Exogenous access
a32NoExL <- IntLengthThetaBP[, 6] # No exclusion
a32UnvL <- IntLengthAlpha1BP[, 6] # Univariate

BoxPlota32L <- data.frame(c(a32BaseL, a32ExAcL, a32NoExL, a32UnvL))
BoxPlota32L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota32L) <- c("Coefficient", "Model")
BoxPlota32L <- BoxPlota32L %>%
  relocate(Model)

BoxPlot18 <- ggplot(BoxPlota32L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot18

ggarrange(BoxPlot16, BoxPlot17, BoxPlot18,
  labels = c(
    "A", "B", "C"
  ),
  ncol = 3, nrow = 1,
  legend = "bottom",
  common.legend = TRUE
)

b12BaseL <- IntLengthThetaBP[, 8] # Baseline
b12ExAcL <- IntLengthTheta1BP[, 2] # Exogenous access
b12NoExL <- IntLengthThetaBP[, 8] # No exclusion
b12UnvL <- IntLengthTheta1BP[, 2] # Univariate

BoxPlotb12L <- data.frame(c(b12BaseL, b12ExAcL, b12NoExL, b12UnvL))
BoxPlotb12L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb12L) <- c("Coefficient", "Model")
BoxPlotb12L <- BoxPlotb12L %>%
  relocate(Model)

BoxPlot19 <- ggplot(BoxPlotb12L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot19

b22BaseL <- IntLengthThetaBP[, 11] # Baseline
b22ExAcL <- IntLengthTheta1BP[, 5] # Exogenous access
b22NoExL <- IntLengthThetaBP[, 11] # No exclusion
b22UnvL <- IntLengthTheta1BP[, 5] # Univariate

BoxPlotb22L <- data.frame(c(b22BaseL, b22ExAcL, b22NoExL, b22UnvL))
BoxPlotb22L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb22L) <- c("Coefficient", "Model")
BoxPlotb22L <- BoxPlotb22L %>%
  relocate(Model)

BoxPlot20 <- ggplot(BoxPlotb22L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot20

b32BaseL <- IntLengthThetaBP[, 14] # Baseline
b32ExAcL <- IntLengthTheta1BP[, 8] # Exogenous access
b32NoExL <- IntLengthThetaBP[, 14] # No exclusion
b32UnvL <- IntLengthTheta1BP[, 8] # Univariate

BoxPlotb32L <- data.frame(c(b32BaseL, b32ExAcL, b32NoExL, b32UnvL))
BoxPlotb32L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb32L) <- c("Coefficient", "Model")
BoxPlotb32L <- BoxPlotb32L %>%
  relocate(Model)

BoxPlot21 <- ggplot(BoxPlotb32L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot21

b13BaseL <- IntLengthThetaBP[, 9] # Baseline
b13ExAcL <- IntLengthTheta1BP[, 3] # Exogenous access
b13NoExL <- IntLengthThetaBP[, 9] # No exclusion
b13UnvL <- IntLengthTheta1BP[, 3] # Univariate

BoxPlotb13L <- data.frame(c(b13BaseL, b13ExAcL, b13NoExL, b13UnvL))
BoxPlotb13L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb13L) <- c("Coefficient", "Model")
BoxPlotb13L <- BoxPlotb13L %>%
  relocate(Model)

BoxPlot22 <- ggplot(BoxPlotb13L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot22

b23BaseL <- IntLengthThetaBP[, 12] # Baseline
b23ExAcL <- IntLengthTheta1BP[, 6] # Exogenous access
b23NoExL <- IntLengthThetaBP[, 12] # No exclusion
b23UnvL <- IntLengthTheta1BP[, 6] # Univariate

BoxPlotb23L <- data.frame(c(b23BaseL, b23ExAcL, b23NoExL, b23UnvL))
BoxPlotb23L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb23L) <- c("Coefficient", "Model")
BoxPlotb23L <- BoxPlotb23L %>%
  relocate(Model)

BoxPlot23 <- ggplot(BoxPlotb23L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot23

b33BaseL <- IntLengthThetaBP[, 15] # Baseline
b33ExAcL <- IntLengthTheta1BP[, 9] # Exogenous access
b33NoExL <- IntLengthThetaBP[, 15] # No exclusion
b33UnvL <- IntLengthTheta1BP[, 9] # Univariate

BoxPlotb33L <- data.frame(c(b33BaseL, b33ExAcL, b33NoExL, b33UnvL))
BoxPlotb33L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb33L) <- c("Coefficient", "Model")
BoxPlotb33L <- BoxPlotb33L %>%
  relocate(Model)

BoxPlot24 <- ggplot(BoxPlotb33L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot24

ggarrange(BoxPlot19, BoxPlot20, BoxPlot21,
  BoxPlot22, BoxPlot23, BoxPlot24,
  labels = c(
    "A", "B", "C", "D", "E", "F"
  ),
  ncol = 3, nrow = 2,
  legend = "bottom",
  common.legend = TRUE
)

d12BaseL <- IntLengthThetaBP[, 17] # Baseline
d12ExAcL <- IntLengthTheta1BP[, 11] # Exogenous access
d12NoExL <- IntLengthThetaBP[, 17] # No exclusion
d12UnvL <- IntLengthTheta1BP[, 11] # Univariate

BoxPlotd12L <- data.frame(c(d12BaseL, d12ExAcL, d12NoExL, d12UnvL))
BoxPlotd12L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd12L) <- c("Coefficient", "Model")
BoxPlotd12L <- BoxPlotd12L %>%
  relocate(Model)

BoxPlot25 <- ggplot(BoxPlotd12L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot25

d22BaseL <- IntLengthThetaBP[, 20] # Baseline
d22ExAcL <- IntLengthTheta1BP[, 14] # Exogenous access
d22NoExL <- IntLengthThetaBP[, 20] # No exclusion
d22UnvL <- IntLengthTheta1BP[, 14] # Univariate

BoxPlotd22L <- data.frame(c(d22BaseL, d22ExAcL, d22NoExL, d22UnvL))
BoxPlotd22L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd22L) <- c("Coefficient", "Model")
BoxPlotd22L <- BoxPlotd22L %>%
  relocate(Model)

BoxPlot26 <- ggplot(BoxPlotd22L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot26

d32BaseL <- IntLengthThetaBP[, 23] # Baseline
d32ExAcL <- IntLengthTheta1BP[, 17] # Exogenous access
d32NoExL <- IntLengthThetaBP[, 23] # No exclusion
d32UnvL <- IntLengthTheta1BP[, 17] # Univariate

BoxPlotd32L <- data.frame(c(d32BaseL, d32ExAcL, d32NoExL, d32UnvL))
BoxPlotd32L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd32L) <- c("Coefficient", "Model")
BoxPlotd32L <- BoxPlotd32L %>%
  relocate(Model)

BoxPlot27 <- ggplot(BoxPlotd32L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot27

d13BaseL <- IntLengthThetaBP[, 18] # Baseline
d13ExAcL <- IntLengthTheta1BP[, 12] # Exogenous access
d13NoExL <- IntLengthThetaBP[, 18] # No exclusion
d13UnvL <- IntLengthTheta1BP[, 12] # Univariate

BoxPlotd13L <- data.frame(c(d13BaseL, d13ExAcL, d13NoExL, d13UnvL))
BoxPlotd13L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd13L) <- c("Coefficient", "Model")
BoxPlotd13L <- BoxPlotd13L %>%
  relocate(Model)

BoxPlot28 <- ggplot(BoxPlotd13L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot28

d23BaseL <- IntLengthThetaBP[, 21] # Baseline
d23ExAcL <- IntLengthTheta1BP[, 15] # Exogenous access
d23NoExL <- IntLengthThetaBP[, 21] # No exclusion
d23UnvL <- IntLengthTheta1BP[, 15] # Univariate

BoxPlotd23L <- data.frame(c(d23BaseL, d23ExAcL, d23NoExL, d23UnvL))
BoxPlotd23L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd23L) <- c("Coefficient", "Model")
BoxPlotd23L <- BoxPlotd23L %>%
  relocate(Model)

BoxPlot29 <- ggplot(BoxPlotd23L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot29

d33BaseL <- IntLengthThetaBP[, 24] # Baseline
d33ExAcL <- IntLengthTheta1BP[, 18] # Exogenous access
d33NoExL <- IntLengthThetaBP[, 24] # No exclusion
d33UnvL <- IntLengthTheta1BP[, 18] # Univariate

BoxPlotd33L <- data.frame(c(d33BaseL, d33ExAcL, d33NoExL, d33UnvL))
BoxPlotd33L$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd33L) <- c("Coefficient", "Model")
BoxPlotd33L <- BoxPlotd33L %>%
  relocate(Model)

BoxPlot30 <- ggplot(BoxPlotd33L, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot30

ggarrange(BoxPlot25, BoxPlot26, BoxPlot27,
  BoxPlot28, BoxPlot29, BoxPlot30,
  labels = c(
    "A", "B", "C", "D", "E", "F"
  ),
  ncol = 3, nrow = 2,
  legend = "bottom",
  common.legend = TRUE
)

########## Box plots: Intercepts #############
S <- 100
a11Base <- THETAhat[1, ] # Baseline
a11ExAc <- ALPHAhat1[1, ] # Exogenous access
a11NoEx <- THETAhat[1, ] # No exclusion
a11Unv <- ALPHAhat1[1, ] # Univariate

BoxPlota11 <- data.frame(c(a11Base, a11ExAc, a11NoEx, a11Unv))
BoxPlota11$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota11) <- c("Coefficient", "Model")
BoxPlota11 <- BoxPlota11 %>%
  relocate(Model)

means <- BoxPlota11 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot31 <- ggplot(BoxPlota11, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a1[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot31

a21Base <- THETAhat[3, ] # Baseline
a21ExAc <- ALPHAhat1[3, ] # Exogenous access
a21NoEx <- THETAhat[3, ] # No exclusion
a21Unv <- ALPHAhat1[3, ] # Univariate

BoxPlota21 <- data.frame(c(a21Base, a21ExAc, a21NoEx, a21Unv))
BoxPlota21$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota21) <- c("Coefficient", "Model")
BoxPlota21 <- BoxPlota21 %>%
  relocate(Model)

means <- BoxPlota21 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot32 <- ggplot(BoxPlota21, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a2[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot32

a31Base <- THETAhat[5, ] # Baseline
a31ExAc <- ALPHAhat1[5, ] # Exogenous access
a31NoEx <- THETAhat[5, ] # No exclusion
a31Unv <- ALPHAhat1[5, ] # Univariate

BoxPlota31 <- data.frame(c(a31Base, a31ExAc, a31NoEx, a31Unv))
BoxPlota31$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlota31) <- c("Coefficient", "Model")
BoxPlota31 <- BoxPlota31 %>%
  relocate(Model)

means <- BoxPlota31 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot33 <- ggplot(BoxPlota31, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = a3[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot33

b11Base <- THETAhat[7, ] # Baseline
b11ExAc <- THETAhat1[1, ] # Exogenous access
b11NoEx <- THETAhat[7, ] # No exclusion
b11Unv <- THETAhat1[1, ] # Univariate

BoxPlotb11 <- data.frame(c(b11Base, b11ExAc, b11NoEx, b11Unv))
BoxPlotb11$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb11) <- c("Coefficient", "Model")
BoxPlotb11 <- BoxPlotb11 %>%
  relocate(Model)

means <- BoxPlotb11 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot34 <- ggplot(BoxPlotb11, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b1[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot34

b21Base <- THETAhat[10, ] # Baseline
b21ExAc <- THETAhat1[4, ] # Exogenous access
b21NoEx <- THETAhat[10, ] # No exclusion
b21Unv <- THETAhat1[4, ] # Univariate

BoxPlotb21 <- data.frame(c(b21Base, b21ExAc, b21NoEx, b21Unv))
BoxPlotb21$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb21) <- c("Coefficient", "Model")
BoxPlotb21 <- BoxPlotb21 %>%
  relocate(Model)

means <- BoxPlotb21 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot35 <- ggplot(BoxPlotb21, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b2[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot35

b31Base <- THETAhat[13, ] # Baseline
b31ExAc <- THETAhat1[7, ] # Exogenous access
b31NoEx <- THETAhat[13, ] # No exclusion
b31Unv <- THETAhat1[7, ] # Univariate

BoxPlotb31 <- data.frame(c(b31Base, b31ExAc, b31NoEx, b31Unv))
BoxPlotb31$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotb31) <- c("Coefficient", "Model")
BoxPlotb31 <- BoxPlotb31 %>%
  relocate(Model)

means <- BoxPlotb31 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot36 <- ggplot(BoxPlotb31, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = b3[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot36

d11Base <- THETAhat[16, ] # Baseline
d11ExAc <- THETAhat1[10, ] # Exogenous access
d11NoEx <- THETAhat[16, ] # No exclusion
d11Unv <- THETAhat1[10, ] # Univariate

BoxPlotd11 <- data.frame(c(d11Base, d11ExAc, d11NoEx, d11Unv))
BoxPlotd11$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd11) <- c("Coefficient", "Model")
BoxPlotd11 <- BoxPlotd11 %>%
  relocate(Model)

means <- BoxPlotd11 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot37 <- ggplot(BoxPlotd11, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d1[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot37

d21Base <- THETAhat[19, ] # Baseline
d21ExAc <- THETAhat1[13, ] # Exogenous access
d21NoEx <- THETAhat[19, ] # No exclusion
d21Unv <- THETAhat1[13, ] # Univariate

BoxPlotd21 <- data.frame(c(d21Base, d21ExAc, d21NoEx, d21Unv))
BoxPlotd21$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd21) <- c("Coefficient", "Model")
BoxPlotd21 <- BoxPlotd21 %>%
  relocate(Model)

means <- BoxPlotd21 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot38 <- ggplot(BoxPlotd21, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d2[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot38

d31Base <- THETAhat[22, ] # Baseline
d31ExAc <- THETAhat1[16, ] # Exogenous access
d31NoEx <- THETAhat[22, ] # No exclusion
d31Unv <- THETAhat1[16, ] # Univariate

BoxPlotd31 <- data.frame(c(d31Base, d31ExAc, d31NoEx, d31Unv))
BoxPlotd31$Model <- c(
  rep("Baseline", S), rep("Exogenous Access", S),
  rep("No Exclusion", S), rep("Univariate", S)
)
colnames(BoxPlotd31) <- c("Coefficient", "Model")
BoxPlotd31 <- BoxPlotd31 %>%
  relocate(Model)

means <- BoxPlotd31 %>%
  group_by(Model) %>%
  summarise(Means = round(mean(Coefficient), 2))

BoxPlot39 <- ggplot(BoxPlotd31, aes(x = Model, y = Coefficient)) +
  geom_boxplot(aes(fill = Model)) +
  stat_summary(
    fun = mean, colour = "darkred", geom = "point",
    shape = 18, size = 3, show.legend = FALSE
  ) +
  geom_hline(yintercept = d3[1], colour = "blue") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
BoxPlot39