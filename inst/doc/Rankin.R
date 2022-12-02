## ----echo=F-------------------------------------------------------------------
### get knitr just the way we like it

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)

## -----------------------------------------------------------------------------
library(ASSISTant)
null.uniform <- rep(1, 7L) ## uniform on 7 support points
hourglass <- c(1, 2, 2, 1, 2, 2, 1)
inverted.hourglass <- c(2, 1, 1, 2, 1, 1, 2)
bottom.heavy <- c(2, 2, 2, 1, 1, 1, 1)
bottom.heavier <- c(3, 3, 2, 2, 1, 1, 1)
top.heavy <- c(1, 1, 1, 1, 2, 2, 2)
top.heavier <- c(1, 1, 1, 2, 2, 3, 3)

## -----------------------------------------------------------------------------
ctlDist <- null.uniform
trtDist <- cbind(null.uniform, null.uniform, null.uniform,
                 hourglass, hourglass, hourglass)

##d <- generateDiscreteRankinScores(rep(1, 6), 10, ctlDist, trtDist)

## -----------------------------------------------------------------------------
data(LLL.SETTINGS)
designParameters <- list(prevalence = rep(1/6, 6),
                         ctlDist = ctlDist,
                         trtDist = trtDist)

designA <- ASSISTDesign$new(trialParameters = LLL.SETTINGS$trialParameters,
                            designParameters = designParameters, discreteData = TRUE)
print(designA)

## -----------------------------------------------------------------------------
result <- designA$explore(numberOfSimulations = 5000, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

