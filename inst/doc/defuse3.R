## ----echo=F--------------------------------------------------------------
### get knitr just the way we like it

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)

## ------------------------------------------------------------------------
library(ASSISTant)
##Fix randomization vector N, errors, eps
trialParameters <- list(N = c(200, 340, 476), type1Error = 0.025,
                        eps = 1/2, type2Error = 0.1)

## ------------------------------------------------------------------------
designParameters <- list(
    nul0 = list(prevalence = rep(1/6, 6), mean = matrix(0, 2, 6),
                sd = matrix(1, 2, 6)),
    alt1 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                       c(0.5, 0.4, 0.3, 0, 0, 0)),
                sd = matrix(1, 2, 6)),
    alt2 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                     c(0.5, 0.5, 0, 0, 0, 0)),
                sd = matrix(1,2, 6)),
    alt3 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.36, 6)),
                sd = matrix(1,2, 6)),
    alt4 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.30, 6)),
                sd = matrix(1,2, 6)),
    alt5 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                       c(0.4, 0.3, 0.2, 0, 0, 0)),
                sd = matrix(1,2, 6)),
    alt6 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                       c(0.5, 0.5, 0.3, 0.3, 0.1, 0.1)),
                sd = matrix(1,2, 6))
)

## ------------------------------------------------------------------------
defuse3 <- DEFUSE3Design$new(trialParameters = trialParameters,
                             numberOfSimulations = 500,
                             designParameters = designParameters$nul0,
                             showProgress = FALSE)
print(defuse3)

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 500,
                          showProgress = FALSE,
                          rngSeed = 28912)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result1 <- defuse3$explore(numberOfSimulations = 500,
                           trueParameters = designParameters$alt1,
                           showProgress = FALSE,
                           rngSeed = 737218)
analysis1 <- defuse3$analyze(result1)
print(defuse3$summary(analysis1))

## ------------------------------------------------------------------------
result2 <- defuse3$explore(numberOfSimulations = 500,
                           trueParameters = designParameters$alt2,
                           showProgress = FALSE,
                          rngSeed = 928812)
analysis2 <- defuse3$analyze(result2)
print(defuse3$summary(analysis2))

## ------------------------------------------------------------------------
null.uniform <- rep(1, 7L) ## uniform on 7 support points
hourglass <- c(1, 2, 2, 1, 2, 2, 1)
inverted.hourglass <- c(2, 1, 1, 2, 1, 1, 2)
bottom.heavy <- c(2, 2, 2, 1, 1, 1, 1)
bottom.heavier <- c(3, 3, 2, 2, 1, 1, 1)
bottom.loaded <- c(4, 4, 3, 3, 2, 1, 1)
top.heavy <- c(1, 1, 1, 1, 2, 2, 2)
top.heavier <- c(1, 1, 1, 2, 2, 3, 3)
top.loaded <- c(1, 1, 2, 3, 3, 4, 4)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
knitr::kable(
           sapply(list(null = null.uniform,
                       hourglass = hourglass,
                       inv.hourglass = inverted.hourglass,
                       bot.heavy = bottom.heavy,
                       bot.heavier = bottom.heavier,
                       bot.loaded = bottom.loaded,
                       top.heavy = top.heavy,
                       top.heavier = top.heavier,
                       top.loaded = top.loaded),
                  computeMeanAndSD)
       )

## ------------------------------------------------------------------------
designParameters <- list(
    nul0 = list(prevalence = rep(1, 2),
                ctlDist = null.uniform,
                trtDist = cbind(null.uniform,
                                null.uniform)),
    alt1 = list(prevalence = rep(1, 2), 
                ctlDist = null.uniform,
                trtDist = cbind(top.loaded,
                                null.uniform)),
    alt2 = list(prevalence = rep(1, 2), 
                ctlDist = null.uniform,
                trtDist = cbind(null.uniform,
                                top.loaded))
)

## ------------------------------------------------------------------------
discDefuse3 <- DEFUSE3Design$new(trialParameters = trialParameters,
                                 numberOfSimulations = 5000,
                                 discreteData = TRUE,
                                 designParameters = designParameters$nul0,
                                 showProgress = FALSE)
print(discDefuse3)

## ------------------------------------------------------------------------
result <- discDefuse3$explore(numberOfSimulations = 50,
                              showProgress = FALSE,
                              rngSeed = 3783)
analysis <- discDefuse3$analyze(result)
print(discDefuse3$summary(analysis))

## ------------------------------------------------------------------------
result1 <- discDefuse3$explore(numberOfSimulations = 50,
                               trueParameters = designParameters$alt1,
                               showProgress = FALSE,
                               rngSeed = 28912)
analysis1 <- discDefuse3$analyze(result1)
print(discDefuse3$summary(analysis1))

## ------------------------------------------------------------------------
result2 <- discDefuse3$explore(numberOfSimulations = 50,
                               trueParameters = designParameters$alt2,
                               showProgress = FALSE,
                               rngSeed = 931)
analysis2 <- discDefuse3$analyze(result2)
print(discDefuse3$summary(analysis2))

