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
## Various settings
settings <- list(setting1 = list(N = c(250, 400, 550), type1Error = 0.025,
                                 eps = 1/2, type2Error = 0.1),
                 setting2 = list(N = c(250, 400, 550), type1Error = 0.05,
                                 eps = 1/2, type2Error = 0.1),
                 setting3 =  list(N = c(250, 400, 550), type1Error = 0.1,
                                  eps = 1/2, type2Error = 0.2),
                 setting4 = list(N = c(250, 400, 550), type1Error = 0.2,
                                 eps = 1/2, type2Error = 0.3))



## ------------------------------------------------------------------------
scenarios <- list(
    scenario0 = list(prevalence = rep(1/6, 6), mean = matrix(0, 2, 6),
                     sd = matrix(1, 2, 6)),
    scenario1 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                            c(0.5, 0.4, 0.3, 0, 0, 0)),
                     sd = matrix(1, 2, 6)),
    scenario2 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                            c(0.3, 0.3, 0, 0, 0, 0)),
                     sd = matrix(1, 2, 6)),
    scenario3 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.3, 6)),
                     sd = matrix(1, 2, 6)),
    scenario4 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                            c(0.4, 0.3, 0.2, 0, 0, 0)),
                     sd = matrix(1, 2, 6)),
    scenario5 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                            c(0.5, 0.5, 0.3, 0.3, 0.1, 0.1)),
                     sd = matrix(1, 2, 6)),
    scenario6 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                            c(0.6, 0.6, -0.3, -0.3, -0.3, -0.3)),
                     sd = matrix(1, 2, 6)),
    scenario7 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.01, 6)),
                     sd = matrix(1, 2, 6)), ## very small effect
    scenario8 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.3, 6)),
                     sd = matrix(1, 2, 6)), ## moderate negative effect
    scenario9 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                            c(0.9, 0.3, 0, -0.1, -0.4, -0.7)),
                     sd = matrix(1, 2, 6)), ## single strong effect with negatives thrown in
    scenario10 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(-0.01, 6)),
                      sd = matrix(1, 2, 6)) ## very small negative effect
)

## ------------------------------------------------------------------------
rngSeed <- 2128783
set.seed(rngSeed)
for (setting in names(settings)) {
    trialParameters <- settings[[setting]]
    for (scenario in names(scenarios)) {
        designParameters <- scenarios[[scenario]]
        cat("##############################\n")
        print(sprintf("%s/%s", setting, scenario))
        cat("##############################\n")
        designA <- ASSISTDesign$new(trialParameters = trialParameters,
                                    designParameters = designParameters)
        print(designA)
        result <- designA$explore(numberOfSimulations = 5000,
                                  rngSeed = rngSeed,
                                  showProgress = FALSE)
        analysis <- designA$analyze(result)
        print(designA$summary(analysis))
        rngSeed <- floor(runif(100000 * runif(1)))
    }
}

