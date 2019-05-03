## Global constant, the number of stages for this design
NUM_STAGES <- 3L

## Global constant, the names of columns in trial history relating to Ihat, the subgroup chosen
IHAT_COL_NAMES = c("decision_Ihat", "wcx_Ihat", "wcx.fut_Ihat", "Nl_Ihat", "Ihat", "stage_Ihat", "lost")

## Global constant, the name of the CI column
CI_COL_NAME = c("bounds")

## Global constant, the name of the column for stage at which the trial exits
STAGE_COL_NAME = "exitStage"


#' Compute the standardized Wilcoxon test statistic for two samples
#'
#' We compute the standardized Wilcoxon test statistic with mean 0 and
#' and standard deviation 1 for samples \eqn{x} and \eqn{y}.  The R function
#' [stats::wilcox.test()] returns the statistic
#'
#' \deqn{
#' U = \sum_i R_i - \frac{m(m + 1)}{2}
#' }{
#' U = (sum over i) R_i - m(m + 1) / 2
#' }
#'
#' where \eqn{R_i} are the ranks of the first sample \eqn{x} of size
#' \eqn{m}. We compute
#'
#' \deqn{
#' \frac{(U - mn(1/2 + \theta))}{\sqrt{mn(m + n + 1) / 12}}
#' }{
#' (U - mn(1/2 + theta)) / (mn(m + n + 1) / 12)^(1/2)
#' }
#'
#' where \eqn{\theta} is the alternative hypothesis shift on the
#'     probability scale, i.e. \eqn{P(X > Y) = 1/2 + \theta}.
#'
#' @param x a sample numeric vector
#' @param y a sample numeric vector
#' @param theta a value > 0 but < 1/2.
#' @return the standardized Wilcoxon statistic
#'
#' @importFrom stats wilcox.test
#' @export
#' @md
wilcoxon <- function(x, y, theta = 0) {
    r <- rank(c(x, y))
    n.x <- as.double(length(x))
    n.y <- as.double(length(y))
    STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x +
                                                     1)/2)
    TIES <- (length(r) != length(unique(r)))
    NTIES <- table(r)
    z <- STATISTIC - n.x * n.y * (1/2 + theta)
    SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) -
                                    sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y -
                                                                         1))))
    ##CORRECTION <- sign(z) * 0.5
    z / SIGMA
}


#' Compute the sample size for any group at a stage assuming a nested
#' structure as in the paper.
#'
#' In the three stage design under consideration, the groups are
#' nested with assumed prevalences and fixed total sample size at each
#' stage. This function returns the sample size for a specified group
#' at a given stage, where the futility stage for the overall group
#' test may be specified along with the chosen subgroup.
#'
#' @param prevalence the vector of prevalence, will be normalized if
#'     not already so. The length of this vector implicitly indicates
#'     the number of groups J.
#' @param N an integer vector of length 3 indicating total sample size
#'     at each of the three stages
#' @param stage the stage of the trial
#' @param group the group whose sample size is desired
#' @param HJFutileAtStage is the stage at which overall futility
#'     occured. Default `NA` indicating it did not occur. Also
#'     ignored if stage is 1.
#' @param chosenGroup the selected group if HJFutilityAtStage is not
#'     `NA`. Ignored if stage is 1.
#' @return the sample size for group
#'
#' @export
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
groupSampleSize <- function(prevalence, N, stage, group, HJFutileAtStage = NA, chosenGroup = NA) {
    if (!integerInRange(N, low = 1) || length(N) != NUM_STAGES) {
        stop("Improper values for sample size N")
    }
    if (!identical(order(N), seq_along(N))) {
        stop("Sample size vector N is not monotone increasing sequence")
    }

    J <- length(prevalences)

    if (!scalarInRange(J, low = 2, high = 10)) {
        stop("Improper number of subgroups; need at least 2; max 10")
    }

    if (any(prevalence <= 0)) {
        stop("Improper prevalence specified")
    }
    prevalences <- prevalence / sum(prevalence)
    q <- cumsum(prevalence)

    if (stage == 1 || is.na(HJFutileAtStage) || stage == HJFutileAtStage) {
        ## catches stage = 1 and all cases where stage == HJFutileAtStage
        N[stage] * q[group]
    } else {
        ## stage > 1 && stage > HJFutileAtStage
        stopifnot(stage <= 3 && 1 <= HJFutileAtStage && HJFutileAtStage < stage)
        if (stage == 2) {
            ## HJFutileAtStage = 1 for sure
            if (group <= chosenGroup) {
                qq <- prevalence[seq_len(chosenGroup)]
                qq <- cumsum(qq / sum(qq))
                s1 <- N[1] * q[chosenGroup]
                (N[2] - s1) * qq[group] + N[1] * q[group]
            } else {
                N[1] * q[group] + (N[2] - N[1] * q[chosenGroup])
            }
        } else {
            ## stage == 3 here
            if (HJFutileAtStage == 2) {
                if (group <= chosenGroup) {
                    qq <- prevalence[seq_len(chosenGroup)]
                    qq <- cumsum(qq / sum(qq))
                    s2 <- N[2] * q[chosenGroup]
                    (N[3] - s2) * qq[group] + N[2] * q[group]
                } else {
                    N[2] * q[group] + (N[3] - N[2]* q[chosenGroup])
                }
            } else {
                ## HJFutileAtStage = 1
                if (group <= chosenGroup) {
                    qq <- prevalence[seq_len(chosenGroup)]
                    qq <- cumsum(qq / sum(qq))
                    s1 <- N[1] * q[chosenGroup]
                    (N[3] - s1) * qq[group] + N[1] * q[group]
                } else {
                    N[1] * q[group] + (N[3] - N[1] * q[chosenGroup])
                }
            }
        }
    }
}

#' Conditional probability of \eqn{i}-th subgroup statistic being
#' chosen given the appropriate mean, covariance matrix and futility
#' boundary \eqn{\tilde{b}}{btilde} at \eqn{v}.
#'
#' The computation involves a \eqn{J-1} multivariate normal integral
#' of the conditional density of the \eqn{i}-th subgroup statistic
#' given that it was maximal among all subgroups:
#'
#'\deqn{
#' \phi_i(v)(\int_0^v\int_0^v\ldots
#' \int_0^{\tilde{b}} \phi_v(z_{-i})dz_{-i})
#' }
#'
#' where \eqn{z_{-i}} denotes all subgroups other than \eqn{i}.
#'
#' @param v the value of the statistic
#' @param i the subgroup
#' @param mu.prime the conditional mean vector of the distribution of
#'     length \eqn{J - 1}; needs to be multipled by the conditional
#'     value, the parameter `v`.
#' @param Sigma.prime the conditional covariance matrix of dimension
#'     \eqn{J-1} by \eqn{J-1}
#' @param fut the futility boundary, which is \eqn{\tilde{b}}{btilde}
#'     for stages 1 and 2, but \eqn{c} for stage 3
#' @return the conditional probability
#'
#' @rdname ASSISTant-internal
#'
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats dnorm
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
den.vs <- function(v, i, mu.prime, Sigma.prime, fut ) {
    ## Density function used in integration
    mu.prime <- mu.prime * v
    mvtnorm::pmvnorm(upper = c(rep(v, nrow(mu.prime) - 1), fut), mean = mu.prime[, i],
                     sigma = Sigma.prime[[i]], algorithm = mvtnorm::Miwa()) * stats::dnorm(v)
}

#' Compute the futility boundary (modified Haybittle-Peto) for the
#' first two stages
#'
#' The futility boundary \eqn{\tilde{b}}{btilde} is computed by
#' solving (under the alternative)
#'
#' \deqn{
#' P(\tilde{Z}_J^1\le\tilde{b} or \tilde{Z}_J^2\le\tilde{b}) = \epsilon\beta }
#'
#' where the superscripts denote the stage and \eqn{\epsilon} is the
#' fraction of the type I error (\eqn{\alpha}) spent and \eqn{\beta}
#' is the type II error. We make use of the joint normal density of
#' \eqn{Z_{J}} (the overall group) at each of the three stages and the
#' fact that the \eqn{\tilde{Z_J}} is merely a translation of
#' \eqn{Z_J}. So here the calculation is based on a mean of zero and
#' has to be translated during use!
#'
#' @param beta the type II error
#' @param cov.J the 3 x 3 covariance matrix
#'
#' @importFrom stats uniroot qnorm
#' @importFrom mvtnorm pmvnorm Miwa
#' @export
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
mHP.btilde <- function (beta, cov.J) {
    sigma <- cov.J[-NUM_STAGES, -NUM_STAGES]
    btilde <- stats::uniroot(f = function(btilde) {
        1 - mvtnorm::pmvnorm(lower = rep(btilde, NUM_STAGES - 1),
                             upper = rep(Inf, NUM_STAGES - 1),
                             sigma = sigma,
                             algorithm = Miwa()) -
            beta },
        lower = stats::qnorm(beta) - 1,
        upper = stats::qnorm(beta^(1 / (NUM_STAGES - 1))) + 1)
    btilde$root
}

#' Compute the efficacy boundary (modified Haybittle-Peto) for the
#' first two stages
#'
#' @param prevalence the vector of prevalences between 0 and 1 summing
#'     to 1. \eqn{J}, the number of groups, is implicitly the length
#'     of this vector and should be at least 2.
#' @param N a three-vector of total sample size at each stage
#' @param cov.J the 3 x 3 covariance matrix for Z_J at each of the
#'     three stages
#' @param mu.prime a list of \eqn{J} mean vectors, each of length
#'     \eqn{J-1} representing the conditional means of all the other
#'     \eqn{Z_j} given \eqn{Z_i}. This mean does not account for the
#'     conditioned value of \eqn{Z_i} and so has to be multiplied by
#'     that during use!
#' @param Sigma.prime a list of \eqn{J} covariance matrices, each
#'     \eqn{J-1} by \eqn{J-1} representing the conditional covariances
#'     all the other \eqn{Z_j} given \eqn{Z_i}
#' @param alpha the amount of type I error to spend
#' @param btilde the futility boundary
#' @param theta the effect size on the probability scale
#' @importFrom stats integrate pnorm uniroot
#' @importFrom mvtnorm pmvnorm Miwa
#' @export
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
mHP.b <- function (prevalence, N, cov.J, mu.prime, Sigma.prime, alpha, btilde, theta) {
    J <- length(prevalence)
    q <- cumsum(prevalence / sum(prevalence))

    crossingProb <- function(b) {
        ##
        ## Function to compute conditional probability of rejecting
        ## the subgroup hypothesis for group i at stage (should be 1
        ## or 2), given that the trial was futile at stage.accept.J.
        ##
        f <- function(stage, stage.accept.J, i) {
            ## Translate btilde appropriately from the theta
            ## (probability) scale to the standard scale; see writeup.
            btilde <- btilde + theta * sqrt(3 * N[stage.accept.J])
            ## Adjust the sample size to account for the loss, once
            ## HJ is accepted
            ssi <- replace(N, stage.accept.J, N[stage.accept.J] * q[i])
            if (stage == stage.accept.J) {
                ## Rejection of subgroup at same stage as futility stage
                ## So this is just an integral of the conditional joint
                ## distribution
                ##
                stats::integrate(
                           function(x) {
                               sapply(x, function(x) den.vs(x, i, mu.prime,
                                                            Sigma.prime, fut = btilde))
                           },
                           lower = b, upper = Inf)$value
            } else {
                ## Rejection of subgroup at the subsequent stage,
                ## that is, stage === stage.accept.J + 1
                ## For efficiency, we don't do explicit checking of such
                ## conditions, but the invocation below is implicitly expected
                ## to respect this fact.
                ##
                ## So here we have to integrate the product of the
                ## conditional joint distribution and the probability
                ## of the i-th group statistic exceeding the at the
                ## next stage.
                ##
                sigma <- sqrt(ssi[stage - 1] / ssi[stage])
                integrand <- function(u) {
                    den.vs(u, i, mu.prime, Sigma.prime, fut = btilde) *
                        stats::pnorm(b, mean = u * sigma, sd = sqrt(1 - sigma^2),
                              lower.tail = FALSE)
                }
                stats::integrate(
                           function(u) sapply(u, integrand),
                           lower = -Inf, upper = b)$value
            }
        }
        ## Type I error at interim stages 1 and 2 =
        ## P(accept H_J at stage 1, reject H_I at stage 1) +
        ## P(accept H_J at stage 1, reject H_I at stage 2) +
        ## P(accept H_J at stage 2, reject H_I at stage 2)
        ##          for I in 1:(J-1) +
        ## P(reject H_J at stage 1 or 2)
        ##
        sum(sapply(
            seq_len(J - 1), function(i)  f(1, 1, i) + f(2, 1, i) + f(2, 2, i))) +
            ##
            1 - mvtnorm::pmvnorm(lower = -Inf, upper = b, mean = rep(0, NUM_STAGES - 1),
                                 sigma = cov.J[-NUM_STAGES, -NUM_STAGES],
                                 algorithm = mvtnorm::Miwa()) -
            ##
            alpha
    }
    stats::uniroot(f = crossingProb, lower = 1.0, upper = 4.0)$root
}


#' Compute the efficacy boundary (modified Haybittle-Peto) for the
#' final (third) stage
#'
#' @param prevalence the vector of prevalences between 0 and 1 summing
#'     to 1. \eqn{J}, the number of groups, is implicitly the length
#'     of this vector and should be at least 2.
#' @param N a three-vector of total sample size at each stage
#' @param cov.J the 3 x 3 covariance matrix for Z_J at each of the
#'     three stages
#' @param mu.prime a list of \eqn{J} mean vectors, each of length
#'     \eqn{J-1} representing the conditional means of all the other
#'     \eqn{Z_j} given \eqn{Z_i}. This mean does not account for the
#'     conditioned value of \eqn{Z_i} and so has to be multiplied by
#'     that during use!
#' @param Sigma.prime a list of \eqn{J} covariance matrices, each
#'     \eqn{J-1} by \eqn{J-1} representing the conditional covariances
#'     all the other \eqn{Z_j} given \eqn{Z_i}
#' @param alpha the amount of type I error to spend
#' @param btilde the futility boundary
#' @param b the efficacy boundary for the first two stages
#' @param theta the effect size on the probability scale
#'
#' @importFrom stats uniroot qnorm integrate
#' @importFrom mvtnorm pmvnorm Miwa
#' @export
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
mHP.c <- function (prevalence, N, cov.J, mu.prime, Sigma.prime, alpha, btilde, b, theta) {
    ## Function for computing final boundary c
    J <- length(prevalence)
    q <- cumsum(prevalence / sum(prevalence))

    crossingProb <- function(c) {
        ##
        ## Function for bounding the probability of rejecting either
        ## H_J or H_I at the third stage.
        ##
        f <- function(stage.accept.J, i) {
            ## i is sub-population selected
            ##
            ## Translate btilde appropriately from the theta
            ## (probability) scale to the standard scale; see writeup.
            ##
            btilde <- btilde + theta * sqrt(3 * N[stage.accept.J])
            ## Adjust the sample size to account for the loss, once
            ## HJ is accepted
            ssi <- replace(N, stage.accept.J, N[stage.accept.J] * q[i])

            if (stage.accept.J == 3 ) {
                ## Rejection of subgroup at same stage as futility
                ## stage, that is, stage = 3.  So this is just an
                ## integral of the conditional joint distribution,
                ## except that at the third stage, the critical
                ## boundary is c.
                ##
                stats::integrate(
                           function(x) {
                               sapply(x, function(x)
                                   den.vs(x, i, mu.prime, Sigma.prime, fut = c))
                           },
                           lower = c,
                           upper = Inf)$value
            } else if (stage.accept.J == 2 ) {
                ## Rejection of subgroup at the subsequent stage,
                ## that is, H_J was futile at stage 2, and H_I is rejected
                ## at stage 3.
                ##
                ## So here we have to integrate the product of the
                ## conditional joint distribution at stage 2 and the
                ## probability of the i-th group statistic exceeding
                ## the critical value stage 3. The latter is a
                ## one-dimensional integral in this case, with
                ## appropriate standard deviation. Note that the
                ## critical boundary is c in the last stage!
                ##
                sigma <- sqrt(ssi[2] / ssi[3])
                integrand <- function(u) {
                    den.vs(u, i, mu.prime, Sigma.prime, btilde) *
                        stats::pnorm(c, mean = u * sigma, sd = sqrt(1 - sigma^2), lower.tail = FALSE)
                }
                stats::integrate(function(u) sapply(u, integrand),
                                 lower = -Inf,
                                 upper = b)$value
            } else {
                ## Rejection of subgroup at the third stage, while
                ## H_J was futile at stage 1, and H_I is rejected
                ## at stage 3.
                ##
                ## So here we have to integrate the product of the
                ## conditional joint distribution at stage 1 and the
                ## probability of the i-th group statistic exceeding
                ## the critical value at stage 3, but not stage 2. The
                ## latter probability is a 2-dimensional integral with
                ## an appropriate covariance structure. Once again,
                ## note that the critical boundary at stage 3 is c.
                ##
                v23 <- sqrt(c(ssi[1] / ssi[2], ssi[1] / ssi[3]))
                sigma23 <- matrix(sqrt(ssi[2] / ssi[3]), 2, 2)
                diag(sigma23) <- 1
                sigma <- sigma23 - v23 %*% t(v23)
                integrand <- function(u) {
                    den.vs(u, i, mu.prime, Sigma.prime, btilde) *
                        mvtnorm::pmvnorm(lower = c(-Inf, c),
                                         upper = c(b, Inf),
                                         mean = u * v23,
                                         sigma = sigma,
                                         algorithm = mvtnorm::Miwa())
                }
                stats::integrate(function(u) sapply(u, integrand),
                                 lower = -Inf,
                                 upper = b)$value
            }
        }
        ##
        ## Type I error at final stage =
        ## P(accept H_J at stage 1, reject H_I at stage 3) +
        ## P(accept H_J at stage 2, reject H_I at stage 3) +
        ## P(accept H_J at stage 3, reject H_I at stage 3) or
        ##             for I in 1:(J-1) +
        ## P(reject H_J at stage 3)
        ##
        sum(sapply(seq_len(J - 1), function(i) f(1, i) + f(2, i) + f(3, i))) +
            ##
            mvtnorm::pmvnorm(lower = c(rep(-Inf, NUM_STAGES - 1), c),
                             upper = c(rep(b, NUM_STAGES - 1), Inf), sigma = cov.J,
                             mean = rep(0, NUM_STAGES),
                             algorithm = mvtnorm::Miwa()) -
            ##
            alpha
    }
    stats::uniroot(f = crossingProb,
                   lower = min(0.0, b - 2.0),
                   upper = max(b + 2.0, 4.0))$root
}

#' Compute the three modified Haybittle-Peto boundaries
#'
#' @param prevalence the vector of prevalences between 0 and 1 summing
#'     to 1. \eqn{J}, the number of groups, is implicitly the length
#'     of this vector and should be at least 2.
#' @param N a three-vector of total sample size at each stage
#' @param alpha the type I error
#' @param beta the type II error
#' @param eps the fraction (between 0 and 1) of the type 1 error to
#'     spend in the interim stages 1 and 2
#' @param futilityOnly a logical value indicating only the futility
#'     boundary is to be computed; default `FALSE`
#' @return a named vector of three values containing
#'     \eqn{\tilde{b}}{btilde}, b, c
#'
#' @importFrom stats integrate pnorm uniroot
#' @export
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
computeMHPBoundaries <- function(prevalence, N, alpha, beta, eps, futilityOnly = FALSE) {
    J <- length(prevalence)
    q <- cumsum(prevalence / sum(prevalence))

    theta <- (qnorm(1 - alpha) + qnorm(1 - beta)) / sqrt(3 * N[3])

    ## Sigma = covariance matrix between subgroup,
    ## which is roughly stage independent

    Sigma <- matrix(0, J, J)
    for (i in seq_len(J - 1)) {
        for (j in (i + 1):J) {
            Sigma[i, j] <- sqrt(q[i] / q[j])
        }
    }
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- 1

    mu.prime <- matrix(0, (J - 1), J)
    Sigma.prime <- vector("list", J)
    for (i in seq_len(J)) {
        mu.prime[, i] <- Sigma[-i, i]
        Sigma.prime[[i]] <- Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
    }


    cov.J <- matrix(0, NUM_STAGES, NUM_STAGES)
    for (i in seq_len(NUM_STAGES - 1)) {
        for (j in (i + 1):NUM_STAGES) {
            cov.J[i, j] <- sqrt((N[i]^2 * (N[j] + 1)) / (N[j]^2 * (N[i] + 1)))
        }
    }
    cov.J <- cov.J + t(cov.J)
    diag(cov.J) <- 1

    btilde <- mHP.btilde(beta = beta * eps, cov.J = cov.J)
    if (futilityOnly) {
        b <- c <- NA
    } else {
        b <- mHP.b(prevalence = prevalence,
                   N = N,
                   cov.J = cov.J,
                   mu.prime = mu.prime,
                   Sigma.prime = Sigma.prime,
                   alpha = alpha * eps,
                   btilde = btilde,
                   theta = theta)
        c <- mHP.c(prevalence = prevalence,
                   N = N,
                   cov.J = cov.J,
                   mu.prime = mu.prime,
                   Sigma.prime = Sigma.prime,
                   alpha = alpha * (1 - eps),
                   btilde = btilde,
                   b = b,
                   theta = theta)
    }
    c(btilde = btilde, b = b, c = c)
}

#' Compute the three modified Haybittle-Peto boundaries and effect size
#'
#' @param prevalence the vector of prevalences between 0 and 1 summing
#'     to 1. \eqn{J}, the number of groups, is implicitly the length
#'     of this vector and should be at least 2.
#' @param alpha the type I error
#' @return a named vector of a single value containing the value for `c`
#' @export
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'     Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'     Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'     191-200
#'     (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @md
computeMHPBoundaryITT <- function(prevalence, alpha) {
    J <- length(prevalence)
    q <- cumsum(prevalence / sum(prevalence))

    ## Sigma = covariance matrix between subgroup,
    ## which is roughly stage independent

    Sigma <- matrix(0, J, J)
    for (i in seq_len(J - 1)) {
        for (j in (i + 1):J) {
            Sigma[i, j] <- sqrt(q[i] / q[j])
        }
    }
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- 1

    mu.prime <- matrix(0, (J - 1), J)
    Sigma.prime <- vector("list", J)
    for (i in seq_len(J)) {
        mu.prime[, i] <- Sigma[-i, i]
        Sigma.prime[[i]] <- Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
    }

    ## Derive interim eff boundary b.I for subgp
    crossingProb <- function(c) {
        f <- function(i) {
            ##i=sub-population selected
            stats::integrate(
                       function(x) {
                           sapply(x, function(x)
                               den.vs(x, i, mu.prime, Sigma.prime, fut = c))
                       },
                       lower = c,
                       upper = Inf)$value
        }
        sum(sapply(seq_len(J - 1), function(i) f(i))) +
            stats::pnorm(c, lower.tail = FALSE) - alpha
    }
    c(cAlpha = stats::uniroot(f = crossingProb, lower = 1, upper = 4)$root)
}


#' Return a vector of column names for statistics for a given stage
#'
#' @param stage the trial stage (1 to 3 inclusive).
#' @param J the number of subgroups
#' @return a character vector of the column names
#' @export
#' @md
colNamesForStage <- function(stage, J) {
    seqJ <- seq_len(J)
    c(paste(c("decision", "wcx", "wcx.fut", "Nl"), stage, sep = "_"),
      sapply(seqJ, function(group) {
          c(sprintf("wcx_%d_%d", stage, group),
            sprintf("nc_%d_%d", stage, group),
            sprintf("nt_%d_%d", stage, group),
            sprintf("muc_%d_%d", stage, group),
            sprintf("mut_%d_%d", stage, group),
            sprintf("sdc_%d_%d", stage, group),
            sprintf("sdt_%d_%d", stage, group))
      }))
}

#' A data generation function using a discrete distribution for Rankin
#' score rather than a normal distribution
#'
#' @param prevalence a vector of group prevalences (length denoted by J below)
#' @param N the sample size to generate
#' @param support the support values of the discrete distribution (length K), default 0:6
#' @param ctlDist a probability vector of length K denoting the Rankin score distribution for control.
#' @param trtDist an K x J probability matrix with each column is the Rankin distribution for the associated group
#' @return a three-column data frame of `subGroup`, `trt` (0 or 1), and `score`
#' @examples
#' # Simulate data from a discrete distribution for the Rankin scores,
#' # which are typically ordinal integers from 0 to 6 in the following
#' # simulations. So we define a few scenarios.
#' library(ASSISTant)
#' null.uniform <- rep(1, 7L) ## uniform on 7 support points
#' hourglass <- c(1, 2, 2, 1, 2, 2, 1)
#' inverted.hourglass <- c(2, 1, 1, 2, 1, 1, 2)
#' bottom.heavy <- c(2, 2, 2, 1, 1, 1, 1)
#' bottom.heavier <- c(3, 3, 2, 2, 1, 1, 1)
#' top.heavy <- c(1, 1, 1, 1, 2, 2, 2)
#' top.heavier <- c(1, 1, 1, 2, 2, 3, 3)
#' ctlDist <- null.uniform
#' trtDist <- cbind(null.uniform, null.uniform, hourglass, hourglass) ## 4 groups
#' generateDiscreteData(prevalence = rep(1, 4), N = 10, ctlDist = ctlDist,
#'                      trtDist = trtDist) ## default support is 0:6
#' trtDist <- cbind(bottom.heavy, bottom.heavy, top.heavy, top.heavy)
#' generateDiscreteData(prevalence = rep(1, 4), N = 10, ctlDist = ctlDist,
#'                      trtDist = trtDist)
#' support <- c(-2, -1, 0, 1, 2) ## Support of distribution
#' top.loaded <- c(1, 1, 1, 3, 3) ## Top is heavier
#' ctl.dist <- c(1, 1, 1, 1, 1) ## null on 5 support points
#' trt.dist <- cbind(ctl.dist, ctl.dist, top.loaded) ## 3 groups
#' generateDiscreteData(prevalence = rep(1, 3), N = 10, support = support,
#'                      ctlDist = ctl.dist, trtDist = trt.dist)
#' @export
#' @md
generateDiscreteData <- function(prevalence, N, support = 0L:6L, ctlDist, trtDist) {
    nR <- length(support)
    stopifnot(length(ctlDist) == nR && nrow(trtDist) == nR)
    J <- length(prevalence)
    stopifnot(J == ncol(trtDist))
    null <- sapply(seq_len(J), function(x) ctlDist)
    dists <- cbind(null, trtDist)

    if (N == 0) {
        data.frame(subGroup = integer(0), trt = integer(0),
                   score = integer(0))
    } else {
        subGroup <- sample.int(n = J, size = N, replace = TRUE,
                               prob = prevalence)
        ## Scale trt to 0:1 from 1:2.
        trt <- sample.int(n = 2L, size = N, replace = TRUE) - 1L

        rankin <- sapply(J * trt + subGroup,
                         function(k) sample(support, size = 1, prob = dists[, k] ))
        data.frame(subGroup = subGroup, trt = trt, score = rankin)
    }
}

#' A data generation function along the lines of what was used in the Lai, Lavori, Liao paper.
#' score rather than a normal distribution
#'
#' @param prevalence a vector of group prevalences (length denoted by J below)
#' @param N the sample size to generate
#' @param mean a 2 x J matrix of means under the null (first row) and alternative for each group
#' @param sd a 2 x J matrix of standard deviations under the null (first row) and alternative for each group
#' @return a three-column data frame of `subGroup`, `trt` (0 or 1), and `score`
#' @export
#' @md
generateNormalData <- function(prevalence, N, mean, sd) {
    J <- length(prevalence)
    stopifnot((J == ncol(mean)) && (J == ncol(sd)))
    if (N == 0) {
        data.frame(subGroup = integer(0), trt = integer(0),
                   score = numeric(0))
    } else {
        subGroup <- sample.int(n = J, size = N, replace = TRUE,
                               prob = prevalence)
        trt <- sample.int(n = 2L, size = N, replace = TRUE) - 1L
        rankin <- unlist(
            Map(function(i, j)
                rnorm(n = 1, mean = mean[i, j], sd = sd[i, j]),
                trt + 1, subGroup))
        data.frame(subGroup = subGroup, trt = trt, score = rankin)
    }
}

#' Compute the mean and sd of a discrete Rankin distribution
#' @param probVec a probability vector of length equal to length of support,
#'        default is uniform
#' @param support a vector of support values (default 0:6 for Rankin Scores)
#' @return a named vector of `mean` and `sd`
#' @export
#' @md
computeMeanAndSD <- function(probVec = rep(1, 7L), support = 0L:6L) {
    stopifnot(all(probVec >= 0))
    probVec <- probVec / sum(probVec)
    mean <- sum(support * probVec)
    sd <- sqrt(sum(probVec * support^2) - mean^2)
    c(mean = mean, sd = sd)
}
