#' A class to encapsulate the adaptive clinical trial design of Lai, Lavori and Liao
#'
#' @description `ASSISTDesign` objects are used to design, simulate and analyze
#' adaptive group sequential clinical trial with three stages. For details refer to the paper
#' _Adaptive Choice of Patient Subgroup for Comparing Two Treatments_
#' by Tze Leung Lai and Philip W. Lavori and Olivia Yueh-Wen Liao. Contemporary Clinical Trials,
#' Vol. 39, No. 2, pp 191-200 (2014). 
#'
#' @seealso `LLL.SETTINGS` for an explanation of trial parameters
#' @importFrom R6 R6Class
#' @importFrom dplyr mutate summarize filter group_by ungroup n select
#' @importFrom magrittr %>%
#' @importFrom knitr kable
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @export
#' @examples
#' \dontrun{
#' data(LLL.SETTINGS)
#' prevalence <- LLL.SETTINGS$prevalences$table1
#' scenario <- LLL.SETTINGS$scenarios$S0
#' designParameters <- list(prevalence = prevalence,
#'                          mean = scenario$mean,
#'                          sd = scenario$sd)
#' designA <- ASSISTDesign$new(trialParameters = LLL.SETTINGS$trialParameters,
#'                               designParameters = designParameters)
#' print(designA)
## A realistic design uses 5000 simulations or more!
#' result <- designA$explore(showProgress = interactive())
#' analysis <- designA$analyze(result)
#' designA$summary(analysis)
#' }
## For full examples, try:
## browseURL(system.file("full_doc/ASSISTant.html", package="ASSISTant"))
ASSISTDesign <-
  R6Class(classname = "ASSISTDesign",
          private = list(
            designParameters = NA,
            trialParameters = NA,
            boundaries = NA,
            discreteData = FALSE,
            checkParameters = function(designParameters, trialParameters, discreteData) {
              if (length(trialParameters$N) != NUM_STAGES) {
                stop(sprintf("Improper sample size vector; this design assumes %d stages", NUM_STAGES))
              }
              if (!integerInRange(trialParameters$N, low = 1)) {
                stop("Improper values for sample sizes")
              }
              if (!identical(order(trialParameters$N), seq_along(trialParameters$N))) {
                stop("Improper values for sample sizes; need increasing sequence")
              }
              prevalence <- designParameters$prevalence
              J <- length(prevalence)
              if (!scalarInRange(J, low = 2, high = 10)) {
                stop("Improper number of subgroups; need at least 2; max 10")
              }
              if (!scalarInRange(trialParameters$type1Error, low=0.0001, 0.2)) {
                stop("Improper type 1 error")
              }
              if (!scalarInRange(trialParameters$type2Error, low=0.0001, 0.3)) {
                stop("Improper type 2 error")
              }
              if (!scalarInRange(trialParameters$eps, low=1e-5, high = 1 - 1e-5)) {
                stop("Improper epsilon specified")
              }
              if (any(prevalence <= 0)) {
                stop("Improper prevalence specified")
              }
              if (discreteData) {
                support <- designParameters$distSupport
                ## Assume Rankin is 0:6 unless specified in designParameters
                if (is.null(support)) {
                  support <- 0L:6L
                }
                K <- length(support)
                ctlDist <- designParameters$ctlDist
                if (!is.matrix(ctlDist)) {
                  if ((length(ctlDist) != K) || any(ctlDist < 0)) {
                    stop(sprintf("Improper ctlDist; need a %d-length probability vector for Rankin scores", K))
                  }
                } else {
                  if ((nrow(ctlDist) != K) || (ncol(ctlDist) != J) || any(ctlDist < 0)) {
                    stop(sprintf("Improper ctlDist; need a %d x %d matrix, with each column a probability vector", K, J))
                  }
                }
                trtDist <- designParameters$trtDist
                if ((nrow(trtDist) != K) || (ncol(trtDist) != J) || any(trtDist < 0)) {
                  stop(sprintf("Improper trtDist; need a %d x %d matrix, with each column a probability vector", K, J))
                }
              } else {
                if (!all.equal(dim(designParameters$mean), c(2, J))) {
                  stop("Mean dimension does not match number of groups")
                }
                if (!all.equal(dim(designParameters$sd), c(2, J))) {
                  stop("Mean dimension does not match number of groups")
                }
                if (!all(designParameters$sd > 0)) {
                  stop("SDs are not all positive")
                }
              }
              TRUE
            },
            selectSubgroup = function (data) {
              data <- data[order(data$subGroup), ]
              ## ASSUME EACH GROUP is represented!!
              ## FIX needed if you use names for groups instead of 0, 1, 2, .., J
              counts <- cumsum(table(data$subGroup))
              wcx <- sapply(counts[-length(counts)],
                            function(n) {
                              d <- data[seq_len(n), ]
                              score <- split(d$score, d$trt)
                              wilcoxon(score$`1`, score$`0`)
                            })
              which.max(wcx)
            },
            doInterimLook = function (data, stage, recordStats = FALSE) {
              d <- split(data$score, data$trt)
              nc <- length(d$`0`) ## control
              nt <- length(d$`1`) ## treatment
              Nl <- sqrt(nt * nc / (nt +  nc))
              wcx <- wilcoxon(d$`1`, d$`0`)
              wcx.fut <- NA
              bdy <- private$boundaries
              if (stage < NUM_STAGES) { ## all but last stage
                if (wcx >= bdy["b"]) { ## Reject
                  decision <- 1
                } else {
                  effectSize <- private$trialParameters$effectSize
                  wcx.fut <- wilcoxon(d$`1`, d$`0`, theta = effectSize)
                  if (wcx.fut < bdy["btilde"]) { ## Futility, so accept
                    decision <- -1
                  } else {
                    decision <- 0 ## continue
                  }
                }
              } else { ## the last stage
                if (wcx >= bdy["c"]) { ## Final boundary
                  decision <- 1 ## reject
                } else {
                  decision <- -1 ## accept
                }
              }
              if (recordStats) {
                stats <- private$getStat(data, stage)
              } else {
                stats <- NA
              }

              list(decision = decision, wcx = wcx, wcx.fut = wcx.fut,
                   Nl = Nl, stats = stats)
            },
            getStat = function(d, stage, sigma = 1) {
              ## THIS HAS TO BE ARTICULATED WITH the trial History columns
              ## IN explore()!!
              ## number of columns for each stage is
              ## decision, wilcoxon, wilcoxon.futility, Nl (= 4L)
              ## means, sds and N for control and treatment, wilcoxon (= 7L) for each group
              ## decision, wilcoxon, wilcoxon.futility NL (= 4L more) for selected subgroup IHat
              ## trialHistory <- matrix(NA, nrow = numberOfSimulations,
              ##                       ncol = NUM_STAGES * (4L + 7L * J) + 4L)
              ##
              ## d is the data so far in the trial
              ## stage is the stage at which the trial stopped (1, 2, 3)
              N <- private$trialParameters$N
              J <- private$designParameters$J
              subGroup <- max(d$subGroup) ## Restrict to the groups we see
              result <- unlist(lapply(seq_len(J),
                                      function(j) {
                                        if (j <= subGroup) {
                                          data <- subset(d, subGroup <= j)
                                          splitData <- split(data$score, data$trt)
                                          c(wilcoxon(splitData$`1`, splitData$`0`),
                                            length(splitData$`0`),
                                            length(splitData$`1`),
                                            mean(splitData$`0`),
                                            mean(splitData$`1`),
                                            sd(splitData$`0`),
                                            sd(splitData$`1`))
                                        } else {
                                          rep(NA, 7L)
                                        }
                                      }))
              ## Drop the first three columns because of above
              names(result) <- colNamesForStage(stage, J)[-(1:4)]
              result
            }
          ),
          public = list(
            #' @description
            #' Create a new `ASSISTDesign` instance using the parameters specified.
            #' @param designParameters parameters of the experimental design. Must contain apropriate distributions to sample from, if `discreteData = TRUE`
            #' @param trialParameters the trial parameters, such as sample size etc.
            #' @param discreteData a flag indicating that a discrete distribution is to be used for the Rankin scores
            #' @param boundaries decision boundaries to use for interim looks, a named vector of `btilde`, `b` and `c` values
            #' @return a new `AssistDesign` object
            initialize = function(designParameters, trialParameters, discreteData = FALSE, boundaries) {
              ## Check parameters
              private$checkParameters(designParameters, trialParameters, discreteData)
              private$discreteData <- discreteData
              ## Conform parameters
              private$designParameters  <- conformParameters(designParameters, discreteData)
              trialParameters$effectSize <- (qnorm(1 - trialParameters$type1Error) +
                                               qnorm(1 - trialParameters$type2Error)) /
                sqrt(3 * trialParameters$N[3])

              private$trialParameters <- trialParameters

              if (missing(boundaries)) {
                private$boundaries <- self$computeCriticalValues()
              } else {
                self$setBoundaries(boundaries)
              }
            },
            #' @description
            #' return the designParameters field
            getDesignParameters = function() private$designParameters,
            #' @description
            #' return the trialParameters field
            getTrialParameters = function() private$trialParameters,
            #' @description
            #' return the boundaries field
            getBoundaries  = function() private$boundaries,
            #' @description
            #' Set the boundaries field
            #' @param value a named vector of `btilde`, `b` and `c` values
            setBoundaries  = function(value) {
              if (!all(names(value) == c("btilde", "b", "c"))) {
                stop("setBoundaries: Need names 'btilde', 'b', and 'c' in order")
              }
              private$boundaries <- value
            },
            #' @description
            #' Print details of the design to console
            print = function() {
              designParameters <- private$designParameters
              cat("Design Parameters:\n")
              cat(sprintf(" Number of Groups: %d\n", designParameters$J))
              cat(" Prevalence:")
              prevalence <- matrix(designParameters$prevalence, nrow = 1)
              colnames(prevalence) <- names(designParameters$prevalence)
              print(knitr::kable(prevalence))
              cat(sprintf("\n Using Discrete Rankin scores? %s\n\n", private$discreteData))
              if (private$discreteData) {
                support <- designParameters$distSupport
                cat(" Null Rankin Distribution:")
                print(knitr::kable(designParameters$ctlDist))
                cat(" Null Mean and SD")
                print(knitr::kable(apply(designParameters$ctlDist, 2,
                                         computeMeanAndSD, support = support)))
                cat(" Alternative Rankin Distribution:\n")
                print(knitr::kable(designParameters$trtDist))
                cat(" Alternative Mean and SD")
                print(knitr::kable(apply(designParameters$trtDist, 2,
                                         computeMeanAndSD, support = support)))
              } else {
                cat(" Normal Rankin Distribution means (null row, alt. row):\n")
                print(knitr::kable(designParameters$mean))
                cat("\n Normal Rankin Distribution SDs (null row, alt. row):\n")
                print(knitr::kable(designParameters$sd))
              }
              cat("\nTrial Parameters:\n")
              str(private$trialParameters)
              
              cat("\nBoundaries:\n")
              boundaries <- matrix(private$boundaries, nrow = 1)
              colnames(boundaries) <- names(private$boundaries)
              print(knitr::kable(boundaries))
            },
            #' @description
            #' Compute the critical boundary values \eqn{\tilde{b}}, \eqn{b} and \eqn{c} for futility, efficacy and final efficacy decisions. This is time consuming so cache where possible.
            #' @return a named vector of critical values with names `btilde`, `b`, and `c` as in the paper
            computeCriticalValues = function() {
              trialParameters <- private$trialParameters
              computeMHPBoundaries(prevalence = private$designParameters$prevalence,
                                   N = trialParameters$N,
                                   alpha = trialParameters$type1Error,
                                   beta = trialParameters$type2Error,
                                   eps = trialParameters$eps)
            },
            #' @description
            #' Explore the design using the specified number of simulations and random number seed and other parameters.
            #' @param numberOfSimulations default number of simulations is 5000
            #' @param rngSeed default seed is 12345
            #' @param trueParameters the state of nature, by default the value of `self$getDesignParameters()` as would be the case for a Type I error calculation. If changed, would yield power.
            #' @param recordStats a boolean flag (default `TRUE`) to record statistics
            #' @param showProgress a boolean flag to show progress, default `TRUE`
            #' @param fixedSampleSize a bollean flag indicating that patients lost after a futile overall look are not made up, default `FALSE`.
            #' @param saveRawData a flag (default `FALSE`) to indicate if raw data has to be saved
            #' @return a list of results
            explore = function (numberOfSimulations = 5000, rngSeed = 12345,
                                trueParameters = self$getDesignParameters(),
                                recordStats = TRUE,
                                showProgress = TRUE,
                                fixedSampleSize = FALSE,
                                saveRawData = FALSE) {
              ##browser()
              ## Save rng state
              oldRngState <- if (exists(".Random.seed", envir = .GlobalEnv)) {
                get(x = ".Random.seed", envir=.GlobalEnv)
              } else {
                NULL
              }
              ## set our seed
              set.seed(seed = rngSeed, normal.kind = NULL)
              
              ## SOME CHECKS needed here when trueParameters is provided
              ## for conformity
              trialParameters <- private$trialParameters
              trueParameters  <- conformParameters(trueParameters, private$discreteData)
              J  <- trueParameters$J
              support <- trueParameters$distSupport
              glrBoundary <- private$boundaries
              prevalence <- trueParameters$prevalence
              
              ## We record the entire trial history
              ## number of columns for each stage is
              ## decision, wilcoxon, Nl (= 3L)
              ## means, sds and N for control and treatment, wilcoxon (= 7L) for each group
              ## decision, wilcoxon, NL (= 3L more) for selected subgroup IHat
              ## + 1 for confidence interval bounds
              ##
              trialHistoryColumnNames <- c(unlist(lapply(seq_len(NUM_STAGES),
                                                         colNamesForStage, J)),
                                           IHAT_COL_NAMES,
                                           CI_COL_NAME,
                                           STAGE_COL_NAME)
              trialHistory <- matrix(NA, nrow = numberOfSimulations,
                                     ncol = length(trialHistoryColumnNames))
              colnames(trialHistory) <- trialHistoryColumnNames
              
              if (showProgress) {
                pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
              }
              
              if (saveRawData) {
                rawData <- vector(mode = "list", length = numberOfSimulations)
              }
              discreteData <- private$discreteData
              
              
              for (i in seq_len(numberOfSimulations)) {
                ## Generate Empty dataset
                if (discreteData) {
                  thisTrialData <-
                    trialData <- generateDiscreteData(prevalence = prevalence,
                                                      N = 0,
                                                      support = support,
                                                      ctlDist = trueParameters$ctlDist,
                                                      trtDist = trueParameters$trtDist)
                  
                } else {
                  thisTrialData <-
                    trialData <- generateNormalData(prevalence = prevalence,
                                                    N = 0,
                                                    mean = trueParameters$mean,
                                                    sd = trueParameters$sd)
                }
                ## H_J is tested first
                subGroup <- J
                N <- c(0, trialParameters$N)
                ## decisions always follow: 0 = continue, 1 = reject, -1 = accept
                for (stage in seq_len(NUM_STAGES)) {
                  ## Generate data for this stage
                  groupIndices <- seq_len(subGroup)
                  if (fixedSampleSize) {
                    sampleSizeForThisStage <- N[stage + 1] - N[stage]
                  } else {
                    sampleSizeForThisStage <- N[stage + 1] - nrow(trialData)
                  }
                  
                  if (discreteData) {
                    thisStageData <- generateDiscreteData(prevalence = prevalence[groupIndices],
                                                          N = sampleSizeForThisStage,
                                                          support = support,
                                                          ctlDist = trueParameters$ctlDist,
                                                          trtDist = trueParameters$trtDist[, groupIndices,
                                                                                           drop = FALSE])
                  } else {
                    thisStageData <- generateNormalData(prevalence = prevalence[groupIndices],
                                                        N = sampleSizeForThisStage,
                                                        mean = trueParameters$mean[, groupIndices,
                                                                                   drop = FALSE],
                                                        sd = trueParameters$sd[, groupIndices,
                                                                               drop = FALSE])
                  }
                  ## Combine it with previous data
                  trialData <- rbind(trialData, thisStageData)
                  if (saveRawData) {
                    thisTrialData <- rbind(thisTrialData, thisStageData)
                  }
                  
                  ## doInterimLook is guaranteed to return a decision 1 or -1 at stage 3
                  interimResult <- private$doInterimLook(data = trialData,
                                                         stage = stage,
                                                         recordStats = recordStats)
                  resultNames <- colNamesForStage(stage, J)
                  if (recordStats) {
                    trialHistory[i, resultNames[-(1:4)]] <- interimResult$stats
                  }
                  trialHistory[i, resultNames[1:4] ] <- unlist(interimResult[1:4], use.names = FALSE)
                  
                  if (interimResult$decision == 1L) {
                    ## H_{subGroup} was rejected
                    ## so trial stops
                    break
                  } else if (interimResult$decision == -1L) {
                    if (subGroup == J) {
                      ## Select a subgroup and perform an interim look using that subgroup
                      subGroup <- private$selectSubgroup(trialData)
                      ## Restrict the data from now on to those in the subGroup
                      ## So we lose some patients when we restrict to subGroup!
                      prevN <- nrow(trialData)
                      trialData <- trialData[trialData$subGroup <= subGroup, ]
                      interimResult <- private$doInterimLook(data = trialData,
                                                             stage = stage,
                                                             recordStats = FALSE)
                      ## Append this sub-result to the interim result already obtained
                      ## so that both the overall and the subgroup results are retained
                      trialHistory[i, IHAT_COL_NAMES] <- c(unlist(interimResult[1:4], use.names = FALSE),
                                                           subGroup, stage, prevN - nrow(trialData))
                      if (interimResult$decision != 0L) {
                        ## H_{\hat{I}} was accepted or rejected
                        ## so trial stops
                        break
                      }
                    } else {
                      ## Trial was futile for H_{\hat{I}}
                      ## So stop
                      break
                    }
                  } else {
                    ## Trial continues to next stage
                    ##
                  }
                }
                ## Record stage at which trial stopped
                trialHistory[i, STAGE_COL_NAME] <- stage
                ##
                ## sigmahat is fixed at 1 for now
                ##
                sigmahat <- 1
                if (stage == 3) {
                  cut <- glrBoundary["c"] # c
                } else {
                  cut <- glrBoundary["b"]
                }
                trialHistory[i, CI_COL_NAME] <- (interimResult$wcx - cut) * sigmahat / interimResult$Nl
                
                if (showProgress) {
                  setTxtProgressBar(pb, i)
                }
                if (saveRawData) {
                  rawData[[i]] <- thisTrialData
                }
              }
              if (showProgress) {
                close(pb)
              }
              ## Restore rng state
              if (is.null(oldRngState)) {
                rm(".Random.seed", envir = .GlobalEnv)
              } else {
                assign(x = ".Random.seed", value = oldRngState, envir = .GlobalEnv)
              }
              if (saveRawData) {
                list(trialHistory = trialHistory, trueParameters = trueParameters,
                     rawData = rawData)
              } else {
                list(trialHistory = trialHistory, trueParameters = trueParameters)
              }
            },
            #' @description
            #' Perform an interim look on trial data
            #' @param trialData trial data frame
            #' @param stage the trial stage
            #' @param recordStats a boolean flag to record all statistics
            #' @param fixedSampleSize a flag to use a fixed sample size to account for loss to follow up
            #' @return the trial history
            performInterimLook = function (trialData, stage, recordStats = FALSE, fixedSampleSize = FALSE) {
              ## Functionally equivalent to private$doInterimLook, but
              ## more error checking is done on data and stage
              expectedDFNames <- c("subGroup", "trt", "score")
              if (any(is.na(match(expectedDFNames, names(trialData))))) {
                stop("Data is missing one or more of columns 'subGroup', 'trt', 'score'")
              }
              if (!(scalarIntegerInRange(stage, low = 1, high = 3))) {
                stop("Stage has to be between 1 and 3 (inclusive)")
              }
              n <- nrow(trialData)
              N <- c(0, private$trialParameters$N)
              ## if (nrow(trialData) != N[stage]) {
              ##     stop("Data size does not match design sample size!")
              ## }

              ## We have to perform interim looks at all previous stages as well..
              ## Remember that the data is cumulative, includes _ALL_ subjects
              subGroup <- J <- private$designParameters$J

              trialHistoryColumnNames <- c(unlist(lapply(seq_len(stage),
                                                         colNamesForStage, J)),
                                           IHAT_COL_NAMES,
                                           CI_COL_NAME,
                                           STAGE_COL_NAME)
              trialHistory <- matrix(NA, nrow = 1,
                                     ncol = length(trialHistoryColumnNames))
              colnames(trialHistory) <- trialHistoryColumnNames

              ##
              ## sigmahat is fixed at 1 for now
              ##
              sigmahat <- 1

              ## Initialize
              stageData <- trialData[0, ]
              dataPtr <- 0
              for (st in seq_len(stage)) {
                if (fixedSampleSize) {
                  sampleSizeForThisStage <- N[stage + 1] - N[stage]
                } else {
                  sampleSizeForThisStage <- N[stage + 1] - nrow(stageData)
                }

                stageData <- rbind(stageData ,
                                   trialData[(dataPtr + 1):(dataPtr + sampleSizeForThisStage), ])
                dataPtr <- dataPtr + sampleSizeForThisStage
                interimResult <- private$doInterimLook(stageData, stage, recordStats)
                resultNames <- colNamesForStage(st, J)
                if (recordStats) {
                  trialHistory[1, resultNames[-(1:4)]] <- interimResult$stats
                }
                trialHistory[1, resultNames[1:4] ] <- unlist(interimResult[1:4], use.names = FALSE)

                if (interimResult$decision == 1L) {
                  ## H_{subGroup} was rejected
                  ## so trial stops
                  trialHistory[1, STAGE_COL_NAME] <- st
                  if (st == 3) {
                    cut <- glrBoundary["c"] # c
                  } else {
                    cut <- glrBoundary["b"]
                  }
                  trialHistory[1, CI_COL_NAME] <- (interimResult$wcx - cut) * sigmahat / interimResult$Nl
                  break
                } else if (interimResult$decision == -1L) {
                  if (subGroup == J) {
                    ## Select a subgroup and perform an interim look using that subgroup
                    subGroup <- private$selectSubgroup(stageData)
                    ## Restrict the data from now on to those in the subGroup
                    ## So we lose some patients when we restrict to subGroup!
                    prevN <- nrow(stageData)
                    stageData <- stageData[stageData$subGroup <= subGroup, ]
                    interimResult <- private$doInterimLook(data = stageData,
                                                           stage = st,
                                                           recordStats = FALSE)
                    ## Append this sub-result to the interim result already obtained
                    ## so that both the overall and the subgroup results are retained
                    trialHistory[1, IHAT_COL_NAMES] <- c(unlist(interimResult[1:4], use.names = FALSE),
                                                         subGroup, st, prevN - nrow(stageData))
                    if (interimResult$decision != 0L) {
                      ## H_{\hat{I}} was accepted or rejected
                      ## so trial stops
                      trialHistory[1, STAGE_COL_NAME] <- st
                      if (st == 3) {
                        cut <- glrBoundary["c"] # c
                      } else {
                        cut <- glrBoundary["b"]
                      }
                      trialHistory[1, CI_COL_NAME] <- (interimResult$wcx - cut) * sigmahat / interimResult$Nl
                      break
                    }
                  } else {
                    ## Trial was futile for H_{\hat{I}}
                    ## So stop
                    trialHistory[1, STAGE_COL_NAME] <- st
                    if (st == 3) {
                      cut <- glrBoundary["c"] # c
                    } else {
                      cut <- glrBoundary["b"]
                    }
                    trialHistory[1, CI_COL_NAME] <- (interimResult$wcx - cut) * sigmahat / interimResult$Nl
                    break
                  }
                }
              }
              trialHistory
            },
            #' @description
            #' Analyze the exploration data from trial
            #' @param trialExploration the result of a call to `explore()` to simulate the design
            #' @return Return a list of summary quantities
            analyze = function (trialExploration) {
              trialHistory <- as.data.frame(trialExploration$trialHistory)
              trueParameters <- trialExploration$trueParameters
              numberOfSimulations <- nrow(trialHistory)
              trialParameters <- private$trialParameters
              designParameters <- private$designParameters
              J <- designParameters$J
              if (private$discreteData) {
                mu <- computeMeanAndSD(probVec = designParameters$ctlDist,
                                       support = designParameters$distSupport)["mean"]
                trueTheta <- cumsum(designParameters$prevalence * mu)
              } else {
                ## These two lines were in versions prior to 1.3-15. They seem wrong!
                ## trueTheta <- cumsum(designParameters$mean[2, ]) / seq_len(J)
                ## trueDelta <- designParameters$mean[2, ]
                trueTheta <- cumsum(designParameters$prevalence * designParameters$mean[2, ])
              }

              trialHistory %>%
                dplyr::mutate(
                  ## Compute reject.ITT
                  ## i.e. no subgroup is chosen and a decision is only made on H_J!
                  reject.ITT = (is.na(Ihat) &
                                  (decision_1 == 1 | decision_2 == 1 | decision_3 == 1)),

                  ## Compute reject.subgp
                  ## i.e. A subgroup is chosen and a rejection is made on the subgroup
                  reject.subgp = !is.na(Ihat) &
                    ((stage_Ihat == 1) &
                       (decision_Ihat == 1 | decision_2 == 1 | decision_3 == 1)) |
                    ((stage_Ihat == 2) & (decision_Ihat == 1 | decision_3 == 1)) |
                    ((stage_Ihat == 3) & (decision_Ihat == 1)),
                  ## Fix up the NAs
                  reject.subgp = ifelse(is.na(reject.subgp), FALSE, reject.subgp),
                  ## Did the trial stop before the last stage?
                  earlyStop = (exitStage < 3),
                  ## Fix up lost, which is 0 if NA
                  lost = ifelse(is.na(lost), 0, lost),
                  ## Did the trial reject the overall or subgroup null?
                  reject = (reject.ITT | reject.subgp),
                  ## Did the trial stop before the last stage for efficacy?
                  earlyStopEff = (reject & earlyStop),
                  ## Did the trial stop before the last stage for futility?
                  earlyStopFut = (!reject & earlyStop),
                  ## If H_J is tested, Ihat is NA, so create a group variable
                  group = ifelse(is.na(Ihat), J, Ihat),
                  ## for each group, we have the relevant true theta
                  theta = trueTheta[group]
                ) -> result

              result %>%
                dplyr::summarize(Rej_H0_ITT = mean(reject.ITT),
                                 Rej_H0_subgp = mean(reject.subgp),
                                 Rej_H0 = mean(reject)) ->
                rejectStats

              result %>%
                dplyr::summarize(earlyStopEff = mean(earlyStopEff),
                                 earlyStopFut = mean(earlyStopFut)) ->
                earlyStopStats

              ## Proportion of rejections by subgroup
              result %>%
                dplyr::filter(reject) %>%
                dplyr::group_by(group) %>%
                dplyr::summarize(count = n()) %>%
                dplyr::mutate(proportion = count / numberOfSimulations) ->
                popReject

              ## Sample size at trial exit
              exitRandSS <- trialParameters$N[result$exitStage]
              ## Sample size at exit taking loss into account
              exitAnalyzeSS <- exitRandSS - result$lost

              ##  Table of exit Stage and proportion of occurrence
              result %>%
                dplyr::group_by(exitStage) %>%
                dplyr::summarize(count = n()) %>%
                dplyr::mutate(proportion = count / numberOfSimulations) %>%
                dplyr::select(exitStage, proportion) ->
                stageAtExitProportion

              ## Table of futility by stage
              result %>%
                dplyr::group_by(stage_Ihat, Ihat) %>%
                dplyr::summarize(count = n()) %>%
                dplyr::filter(!is.na(stage_Ihat)) %>%
                dplyr::select(stage_Ihat, Ihat, count) %>%
                as.matrix -> temp

              futilityTable <- matrix(0L, nrow = 3, ncol = J - 1)
              for (i in seq_len(nrow(temp))) {
                futilityTable[temp[i, 1], temp[i, 2]] <- temp[i, 3]
              }
              futilityTable <- cbind(seq_len(3), futilityTable)
              colnames(futilityTable) <- c("FutilityStage", paste0("G", seq_len(J-1)))

              ## Table of loss statistics by stage (mean and sd)
              result %>%
                dplyr::group_by(stage_Ihat, Ihat) %>%
                dplyr::filter(!is.na(stage_Ihat)) %>%
                dplyr::summarize(mean = mean(lost), sd = sd(lost)) %>%
                dplyr::rename(FutilityStage = stage_Ihat, selectedGroup = Ihat) ->
                lossTable

              ## CI Report
              result %>%
                dplyr::group_by(group) %>%
                dplyr::summarize(coverage = mean(bounds <= theta),
                                 selectedCount = n(),
                                 rejectedCount = sum(reject)) %>%
                dplyr::select(coverage, selectedCount, rejectedCount) ->
                coverage

              result %>%
                summarize(overall = mean(bounds <= theta ),
                          rejection = mean( bounds[reject] <= theta[reject] , na.rm = TRUE)) ->
                overallAndRejectionCoverage

              list(numberOfSimulations = numberOfSimulations,
                   reject = rejectStats,
                   earlyStopStats = earlyStopStats,
                   popReject = popReject,
                   lost_Stats= list(mean = mean(result$lost), sd = sd(result$lost)),
                   exitRandSS_Stats = list(mean = mean(exitRandSS), sd = sd(exitRandSS)),
                   exitAnalyzeSS_Stats = list(mean = mean(exitAnalyzeSS), sd = sd(exitAnalyzeSS)),
                   futilityTable = futilityTable,
                   lossTable = lossTable,
                   stageAtExitProportion = stageAtExitProportion,
                   coverage = coverage,
                   overallAndRejectionCoverage = overallAndRejectionCoverage)
            },
            
            #' @description
            #' Print the operating characteristics of the design using the analysis data
            #' @param analysis the analysis result from the `analyze()` call
            summary = function(analysis) {
              with(analysis, {
                cat(sprintf("P(Reject H0_ITT) = %f; P(Reject H0_subgp) = %f; P(Reject H0) = %f\n",
                            reject$Rej_H0_ITT, reject$Rej_H0_subgp, reject$Rej_H0))
                cat(sprintf("P(Early stop for efficacy [futility]) = %f [%f]\n",
                            earlyStopStats$earlyStopEff, earlyStopStats$earlyStopFut))
                cat(sprintf("Mean [SD] Randomized N = %f [%f]\n",
                            exitRandSS_Stats$mean, exitRandSS_Stats$sd))
                cat("\nStage at exit (proportion)\n")
                print(knitr::kable(stageAtExitProportion))
                cat(sprintf("\nMean [SD] Lost N = %f [%f]\n",
                            lost_Stats$mean, lost_Stats$sd))
                cat(sprintf("Mean [SD] Analyzed N = %f [%f]\n",
                            exitAnalyzeSS_Stats$mean, exitAnalyzeSS_Stats$sd))
                cat("\nMean loss by futility stage and subgroup\n")
                print(knitr::kable(lossTable))
                cat("\nChance of each subpopulation rejected\n")
                print(knitr::kable(popReject))
                cat("\nCounts by futility stage and subgroup choice\n")
                print(knitr::kable(futilityTable))
                cat("\nCI Statistics:")
                cat('\nOverall coverage and coverage for rejections:')
                print(knitr::kable(overallAndRejectionCoverage))
                cat('\nP(theta_test is in the confidence interval)\n')
                print(knitr::kable(coverage))
                invisible()
              })
            }
          ))


#' A fixed sample design to compare against the adaptive clinical
#' trial design
#'
#' @description `ASSISTDesignB` objects are used to design a trial
#'   with certain characteristics provided in the object instantiation
#'   method. This design differs from `ASSISTDesign` in only how it
#'   computes the critical boundaries, how it performs the interim
#'   look, and what quantities are computed in a trial run.
#'
#' @seealso `ASSISTDesign` which is a superclass of this object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @export
#' @examples
#' \dontrun{
#' data(LLL.SETTINGS)
#' prevalence <- LLL.SETTINGS$prevalences$table1
#' scenario <- LLL.SETTINGS$scenarios$S0
#' designParameters <- list(prevalence = prevalence,
#'                        mean = scenario$mean,
#'                        sd = scenario$sd)
#' designB <- ASSISTDesignB$new(trialParameters = LLL.SETTINGS$trialParameters,
#'                             designParameters = designParameters)
#' print(designB)
#' ## A realistic design uses 5000 simulations or more!
#' result <- designB$explore(showProgress = interactive())
#' analysis <- designB$analyze(result)
#' designB$summary(analysis)
#' }
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/ASSISTant.html", package="ASSISTant"))
#'
ASSISTDesignB <-
  R6Class("ASSISTDesignB",
          inherit = ASSISTDesign,
          private = list(
            doInterimLook = function (data) {
              d <- split(data$score, data$trt)
              wcx <- wilcoxon(d$`1`, d$`0`)
              bdy <- private$boundaries
              if (wcx >= bdy["cAlpha"]) { ## Final boundary
                decision <- 1 ## reject
              } else {
                decision <- -1 ## accept
              }

              list(decision = decision, wcx = wcx)
            },
            ## Function for computing futility boundary btilde
            mHP.ITT = function (mu.prime, Sigma.prime, alpha) {
              J <- private$designParameters$J
              ## Derive interim eff boundary b.I for subgp
              crossingProb <- function(c) {
                f <- function(i) { #i=sub-population selected
                  integrate(function(x) {
                    sapply(x, function(x)
                      private$den.vs(x, i, mu.prime, Sigma.prime, c))}, c, Inf)$value
                }
                sum(sapply(seq_len(J - 1), function(i) f(i))) +
                  pnorm(c, lower.tail = FALSE) - alpha
              }
              uniroot(f = crossingProb, lower = 1, upper = 4, maxiter = 20)$root
            }
          ),
          public = list(
            #' @description
            #' Compute the critical boundary value \eqn{c_\alpha}
            #' @return a named vector of a single value containing the value for `c`
            computeCriticalValues = function() {
              trialParameters <- private$trialParameters
              designParameters <- private$designParameters
              computeMHPBoundaryITT(prevalence = private$designParameters$prevalence,
                                    alpha = private$trialParameters$type1Error)
            },
            #' @description
            #' Explore the design using the specified number of simulations, random number seed, and
            #' further parameters.
            #' @param numberOfSimulations default number of simulations is 100
            #' @param rngSeed default seed is 12345
            #' @param trueParameters the state of nature, by default the value of `self$getDesignParameters()` as would be the case for a Type I error calculation. If changed, would yield power.
            #' @param showProgress a boolean flag to show progress, default `TRUE`
            #' @param saveRawData a flag (default `FALSE`) to indicate if raw data has to be saved
            #' @return a list of results
            explore = function (numberOfSimulations = 100, rngSeed = 12345,
                                trueParameters = self$getDesignParameters(),
                                showProgress = TRUE,
                                saveRawData = FALSE) {
              ## Save rng state
              oldRngState <- if (exists(".Random.seed", envir = .GlobalEnv)) {
                get(x = ".Random.seed", envir=.GlobalEnv)
              } else {
                NULL
              }
              ## set our seed
              set.seed(seed = rngSeed, normal.kind = NULL)
              
              trialParameters <- private$trialParameters
              
              ## SOME CHECKS needed here when trueParameters is provided
              ## for conformity
              trueParameters  <- conformParameters(trueParameters, private$discreteData)
              support <- trueParameters$distSupport
              J <- trueParameters$J
              
              glrBoundary <- private$boundaries
              support <- trueParameters$distSupport
              
              naVec <- rep(NA, numberOfSimulations)
              zeroVec <- integer(numberOfSimulations)
              trialHistory <- data.frame(decision = naVec, select = naVec,
                                         statistic = naVec,
                                         matrix(0, numberOfSimulations, J))
              
              if (showProgress) {
                pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
              }
              
              if (saveRawData) {
                rawData <- data.frame(simId = integer(0),
                                      subGroup = integer(0),
                                      trt = integer(0),
                                      score = numeric(0))
              }
              
              for (i in seq_len(numberOfSimulations)) {
                if (private$discreteData) {
                  dataSoFar <- generateDiscreteData(prevalence = trueParameters$prevalence,
                                                    N = trialParameters$N[3],
                                                    support = support,
                                                    ctlDist = trueParameters$ctlDist,
                                                    trtDist = trueParameters$trtDist)
                  
                } else {
                  dataSoFar <- generateNormalData(prevalence = trueParameters$prevalence,
                                                  N = trialParameters$N[3],
                                                  mean = trueParameters$mean,
                                                  sd = trueParameters$sd)
                }
                
                if (saveRawData) {
                  rawData <- rbind(rawData,
                                   data.frame(simId = i, trialData))
                }
                
                interim <- private$doInterimLook(dataSoFar)
                subGroup <- J ## Last group
                if (interim$decision == -1) { ## continue
                  subGroup <- private$selectSubgroup(dataSoFar)
                  interim <- private$doInterimLook(dataSoFar[dataSoFar$subGroup <= subGroup, ])
                }
                trialHistory[i, ] <- c(decision = interim$decision,
                                       select = subGroup,
                                       statistic = interim$wcx,
                                       table(dataSoFar$subGroup))
                if (showProgress) {
                  setTxtProgressBar(pb, i)
                }
              }
              if (showProgress) {
                close(pb)
              }
              ## Restore rng state
              if (is.null(oldRngState)) {
                rm(".Random.seed", envir = .GlobalEnv)
              } else {
                assign(x = ".Random.seed", value = oldRngState, envir = .GlobalEnv)
              }
              names(trialHistory) <- c("decision", "select", "statistic",
                                       sapply(seq_len(J), function(i) paste0("G", i)))
              if (saveRawData) {
                list(trialHistory = trialHistory, trueParameters = trueParameters,
                     rawData = rawData)
              } else {
                list(trialHistory = trialHistory, trueParameters = trueParameters)
              }
            },
            #' @description
            #' Analyze the exploration data from trial
            #' @param trialExploration the result of a call to `explore()` to simulate the design
            #' @return Return a list of summary quantities
            analyze = function (trialExploration) {
              J <- private$designParameters$J
              trialHistory <- trialExploration$trialHistory
              numberOfSimulations <- nrow(trialHistory)
              reject <- (trialHistory$decision == 1)
              rejectGroupTable <- table(trialHistory$select[reject])

              list(reject = reject, rejectGroupTable = rejectGroupTable,
                   rejectSubgroup = sum(rejectGroupTable[-J]) / numberOfSimulations)
            },
            #' @description
            #' Print the operating characteristics of the design using the analysis data
            #' @param analysis the analysis result from the `analyze()` call
           summary = function(analysis) {
              numberOfSimulations <- length(analysis$reject)
              cat(sprintf("P(Reject H0) = %f\n",
                          mean(analysis$reject)))
              cat(sprintf("P(Reject H0_ITT) = %f\n",
                          mean(analysis$reject) - analysis$rejectSubgroup))
              cat(sprintf("P(Reject H0_subgp) = %f\n",
                          analysis$rejectSubgroup))
              cat("\nChance of each subpopulation rejected\n")
              print(analysis$rejectGroupTable / numberOfSimulations)
            }
          ))

#' A fixed sample RCT design to compare against the adaptive clinical
#' trial design of Lai, Lavori and Liao.
#'
#' @description `ASSISTDesignC` objects are used to design a trial
#'   with certain characteristics provided in the object instantiation
#'   method. This design differs from `ASSISTDesign` in only how it
#'   computes the critical boundaries, how it performs the interim
#'   look, and what quantities are computed in a trial run.
#'
#' @seealso `ASSISTDesignB` which is a superclass of this object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @export
#' @examples
#' data(LLL.SETTINGS)
#' prevalence <- LLL.SETTINGS$prevalences$table1
#' scenario <- LLL.SETTINGS$scenarios$S0
#' designParameters <- list(prevalence = prevalence,
#'                        mean = scenario$mean,
#'                        sd = scenario$sd)
#' ## A realistic design uses 5000 simulations or more!
#' designC <- ASSISTDesignC$new(trialParameters = LLL.SETTINGS$trialParameters,
#'                             designParameters = designParameters)
#' print(designC)
#' result <- designC$explore(numberOfSimulations = 100, showProgress = interactive())
#' analysis <- designC$analyze(result)
#' designC$summary(analysis)
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/ASSISTant.html", package="ASSISTant"))
#'
ASSISTDesignC <-
  R6Class("ASSISTDesignC",
          inherit = ASSISTDesignB,
          public = list(
            #' @description
            #' Compute the critical boundary values \eqn{\tilde{b}}, \eqn{b} and \eqn{c} for futility, efficacy and final efficacy decisions. This is time consuming so cache where possible.
            #' @return a named list containing the critical value `cAlpha`
            computeCriticalValues = function() {
              list(cAlpha = qnorm(1 - private$trialParameters$type1Error))
            },
            #' @description
            #' Explore the design using the specified number of simulations and random number seed and other parameters.
            #' @param numberOfSimulations default number of simulations is 5000
            #' @param rngSeed default seed is 12345
            #' @param trueParameters the state of nature, by default the value of `self$getDesignParameters()` as would be the case for a Type I error calculation. If changed, would yield power.
            #' @param showProgress a boolean flag to show progress, default `TRUE`
            #' @param saveRawData a flag (default `FALSE`) to indicate if raw data has to be saved
            #' @return a list of results
            explore = function (numberOfSimulations = 5000, rngSeed = 12345,
                                trueParameters = self$getDesignParameters(),
                                showProgress = TRUE,
                                saveRawData = FALSE) {
              ## Save rng state
              oldRngState <- if (exists(".Random.seed", envir = .GlobalEnv)) {
                get(x = ".Random.seed", envir=.GlobalEnv)
              } else {
                NULL
              }
              ## set our seed
              set.seed(seed = rngSeed, normal.kind = NULL)
              trialParameters <- private$trialParameters
              trueParameters  <- conformParameters(trueParameters, private$discreteData)
              J <- trueParameters$J

              glrBoundary <- private$boundaries
              support <- trueParameters$distSupport

              naVec <- rep(NA, numberOfSimulations)
              zeroVec <- integer(numberOfSimulations)
              trialHistory <- data.frame(decision = naVec,
                                         statistic = naVec,
                                         matrix(0, numberOfSimulations, J))

              if (showProgress) {
                pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
              }

              if (saveRawData) {
                rawData <- data.frame(simId = integer(0),
                                      subGroup = integer(0),
                                      trt = integer(0),
                                      score = numeric(0))
              }
              for (i in seq_len(numberOfSimulations)) {
                if (private$discreteData) {
                  dataSoFar <- generateDiscreteData(prevalence = trueParameters$prevalence,
                                                    N = trialParameters$N[3],
                                                    support = support,
                                                    ctlDist = trueParameters$ctlDist,
                                                    trtDist = trueParameters$trtDist)

                } else {
                  dataSoFar <- generateNormalData(prevalence = trueParameters$prevalence,
                                                  N = trialParameters$N[3],
                                                  mean = trueParameters$mean,
                                                  sd = trueParameters$sd)
                }
                if (saveRawData) {
                  rawData <- rbind(rawData,
                                   data.frame(simId = i, trialData))
                }
                interim <- private$doInterimLook(dataSoFar)
                trialHistory[i, ] <- c(decision = interim$decision,
                                       statistic = interim$wcx,
                                       table(dataSoFar$subGroup))
                if (showProgress) {
                  setTxtProgressBar(pb, i)
                }
              }
              if (showProgress) {
                close(pb)
              }
              ## Restore rng state
              if (is.null(oldRngState)) {
                rm(".Random.seed", envir = .GlobalEnv)
              } else {
                assign(x = ".Random.seed", value = oldRngState, envir = .GlobalEnv)
              }
              names(trialHistory) <- c("decision", "statistic",
                                       sapply(seq_len(J), function(i) paste0("G", i)))
              if (saveRawData) {
                list(trialHistory = trialHistory, trueParameters = trueParameters,
                     rawData = rawData)
              } else {
                list(trialHistory = trialHistory, trueParameters = trueParameters)
              }
            },
            #' @description
            #' Analyze the design given the `trialExploration` data
            #' @param trialExploration the results from a call to `explore()` to simulate the design
            #' @return a named list of rejections
            analyze = function (trialExploration) {
              J <- private$designParameters$J
              trialHistory = trialExploration$trialHistory
              numberOfSimulations <- nrow(trialHistory)
              reject <- (trialHistory$decision == 1)
              list(reject = reject)
            },
            #' @description
            #' Print the operating characteristics of the design using the analysis data
            #' @param analysis the analysis result from the `analyze()` call
            #' @return no value, just print
            summary = function(analysis) {
              numberOfSimulations <- length(analysis$reject)
              cat(sprintf("P(Reject H0) = %f\n",
                          mean(analysis$reject)))
            }
          ))

#' The DEFUSE3 design
#'
#'
#' @description `DEFUSE3Design` is a slight variant of the the adaptive
#' clinical trial design of Lai, Lavori and Liao. Simulation is used to compute
#' the expected maximum sample size and the boundary for early futility is adjusted to
#' account as well.
#'
#' @seealso `ASSISTDesign` which is a superclass of this object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#'
#' @export
#' @examples
#' trialParameters <- list(N = c(200, 340, 476), type1Error = 0.025,
#'                         eps = 1/2, type2Error = 0.1)
#' designParameters <- list(
#'    nul0 = list(prevalence = rep(1/6, 6), mean = matrix(0, 2, 6),
#'                sd = matrix(1, 2, 6)),
#'    alt1 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.5, 0.4, 0.3, 0, 0, 0)),
#'                sd = matrix(1, 2, 6)),
#'    alt2 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.5, 0.5, 0, 0, 0, 0)),
#'                sd = matrix(1,2, 6)),
#'    alt3 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.36, 6)),
#'                sd = matrix(1,2, 6)),
#'    alt4 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.30, 6)),
#'                sd = matrix(1,2, 6)),
#'    alt5 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.4, 0.3, 0.2, 0, 0, 0)),
#'                sd = matrix(1,2, 6)),
#'    alt6 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.5, 0.5, 0.3, 0.3, 0.1, 0.1)),
#'                sd = matrix(1,2, 6)))
#'
#'\dontrun{
#' ## A realistic design uses 5000 simulations or more!
#' defuse3 <- DEFUSE3Design$new(trialParameters = trialParameters,
#'                              numberOfSimulations = 25,
#'                              designParameters = designParameters$nul0,
#'                              showProgress = FALSE)
#' print(defuse3)
#' result <- defuse3$explore(showProgress = interactive())
#' analysis <- defuse3$analyze(result)
#' print(defuse3$summary(analysis))
#' }
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/defuse3.html", package="ASSISTant"))
#'
DEFUSE3Design <-
  R6Class("DEFUSE3Design",
          inherit = ASSISTDesign,
          private = list(
            originalBoundaries = NA
          ),
          public = list(

            #' @description
            #' Return the original boundaries for the design
            #' @return a named vector of values for `b`, `btilde` and `c`
            getOriginalBoundaries = function() private$originalBoundaries,

            #' @description
            #' Create a `DEFUSE3Design` object
            #' @param designParameters parameters of the experimental design. Must contain apropriate distributions to sample from, if `discreteData = TRUE`
            #' @param trialParameters the trial parameters, such as sample size etc.
            #' @param discreteData a flag indicating that a discrete distribution is to be used for the Rankin scores
            #' @param numberOfSimulations the number of simulations to use, default 5000
            #' @param rngSeed the random number generator seed
            #' @param showProgress a boolean flag to show progress (default `TRUE`)
            #' @param trueParameters a list of true parameter values reflecting the state of nature
            #' @param boundaries decision boundaries to use for interim looks, a named vector of `btilde`, `b` and `c` values
            #' @return a new `AssistDesign` object
            initialize = function(designParameters, trialParameters, discreteData = FALSE,
                                  numberOfSimulations = 5000, rngSeed = 54321,
                                  showProgress = TRUE,
                                  trueParameters = NULL,
                                  boundaries) {
              super$initialize(designParameters, trialParameters, discreteData, boundaries)
              ## Save original Effect sizes
              ##browser()
              private$originalBoundaries <- private$boundaries
              private$trialParameters$originalEffectSize <- private$trialParameters$effectSize
              if (missing(boundaries)) {
                self$adjustCriticalValues(numberOfSimulations, rngSeed, showProgress)
              }
            },
            
            #' @description
            #' Adjust critical values to account for sample size loss due to futility
            #' @param numberOfSimulations the number of simulations to use
            #' @param rngSeed the random number generator seed
            #' @param showProgress a boolean flag for showing progress
            #' @return the adjusted boundaries
            adjustCriticalValues = function(numberOfSimulations, rngSeed, showProgress) {
              designParameters <- private$designParameters
              trialParameters <- private$trialParameters
              ## Run simulations to estimate expect max sample sizes
              result <- as.data.frame(
                self$explore(numberOfSimulations = numberOfSimulations,
                             rngSeed = rngSeed,
                             recordStats = FALSE,
                             showProgress = showProgress)$trialHistory
              )
              q <- cumsum(designParameters$prevalence)
              N <- trialParameters$N
              simDN <- matrix(NA, nrow = numberOfSimulations, ncol = 3L)

              for (i in seq_len(numberOfSimulations)) {
                j <- result$exitStage[i]  ## the stage at which the trial ended
                seq_j <- seq_len(j)
                ## incremental recruitment per stage per norm
                simDN[i, seq_j] <- diff(c(0, N[seq_j]))
                jfut <- result$stage_Ihat[i] ## Stage at which we had futility
                if (!is.na(jfut)) {
                  seq_jfut <- seq_len(jfut)
                  ## sample size adjustment for loss
                  simDN[i, seq_jfut] <- simDN[i, seq_jfut] * q[result$Ihat[i]]
                }
              }

              ## End Simulation

              DEM <- apply(simDN, 2, mean, na.rm = TRUE)
              EM <- floor(cumsum(DEM)) ## Expected N actually
              J <- designParameters$J
              expectedEffectSize <- (qnorm(1 - trialParameters$type1Error) +
                                       qnorm(1 - trialParameters$type2Error)) /
                sqrt(3 * EM[3])
              ## Update the effect size
              private$trialParameters$effectSize <- expectedEffectSize
              private$boundaries <- computeMHPBoundaries(prevalence = designParameters$prevalence,
                                                         N = EM,
                                                         alpha = trialParameters$type1Error,
                                                         beta = trialParameters$type2Error,
                                                         eps = trialParameters$eps)
            },

            #' @description
            #' Explore the design using the specified number of simulations and random number seed and other parameters.
            #' @param numberOfSimulations default number of simulations is 5000
            #' @param rngSeed default seed is 12345
            #' @param trueParameters the state of nature, by default the value of `self$getDesignParameters()` as would be the case for a Type I error calculation. If changed, would yield power.
            #' @param showProgress a boolean flag to show progress, default `TRUE`
            #' @param recordStats a boolean flag (default `TRUE`) to record statistics
            #' @param saveRawData a flag (default `FALSE`) to indicate if raw data has to be saved
            #' @return a list of results
            explore = function (numberOfSimulations = 5000, rngSeed = 12345,
                                trueParameters = self$getDesignParameters(),
                                recordStats = TRUE,
                                showProgress = TRUE,
                                saveRawData = FALSE) {
              super$explore(numberOfSimulations, rngSeed, trueParameters,
                            recordStats, showProgress, fixedSampleSize = TRUE,
                            saveRawData = saveRawData)
            },
            
            #' @description
            #' Perform an interim look for futility 
            #' @param trialData trial data frame
            #' @param stage the trial stage
            #' @param recordStats a boolean flag to record all statistics
            #' @return the trial history
            performInterimLook = function (trialData, stage, recordStats = FALSE) {
              super$performInterimLook(trialData, stage, recordStats, fixedSampleSize = TRUE)
            }

          ))

