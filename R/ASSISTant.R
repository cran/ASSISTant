#' Three stage group sequential adaptive design with subgroup selection
#'
#' \code{ASSISTant} is a package that implements a three-stage
#' adaptive clinical trial design with provision for subgroup
#' selection where the treatment may be effective; see Lai, Lavori and
#' Liao (\doi{10.1016/j.cct.2014.09.001}). The main design object is
#' an \code{R6} class that can be instantiated and manipulated to
#' obtain the operating characteristics. A vignette is provided
#' showing the use of this package for designing the DEFUSE-3 trial,
#' described in the paper by Lai, Lavori and Liao. The package
#' contains everything necessary to reproduce the results of the
#' paper.
#' @docType package
#' @name ASSISTant
#' @references Adaptive Choice of Patient Subgroup for Comparing Two
#'   Treatments by Tze Leung Lai and Philip W. Lavori and Olivia
#'   Yueh-Wen Liao. Contemporary Clinical Trials, Vol. 39, No. 2, pp
#'   191-200 (2014, \doi{10.1016/j.cct.2014.09.001}).
#' @references Adaptive design of confirmatory trials: Advances and
#'   challenges by Tze Leung Lai and Philip W. Lavori and Ka Wai
#'   Tsang. Contemporary Clinical Trials, Vol. 45, Part A, pp 93-102
#'   (2015, \doi{10.1016/j.cct.2015.06.007}).
#'
NULL

#' Is a scalar quantity is a specified range?
#' @description Check if the argument is a scalar within specified range
#'
#' @param low the lower bound, default \code{-Inf}
#' @param high the upper bound, default \code{Inf}
#' @return \code{TRUE} or \code{FALSE}
#'
#' @rdname ASSISTant-internal
#'
#' @keywords internal
scalarInRange <- function(x, low = -Inf, high = Inf) {
  ## #'
  ## #' @examples
  ## #' ASSISTant:::scalarInRange(x = 10, low = 2)  ## TRUE
  ## #' ASSISTant:::scalarInRange(x = 10, high = 2) ## FALSE
  (length(x) == 1) && (x >= low) && (x <= high)
}

#' Is the numeric quantity is a specified range?
#' @description Check if the argument is within specified range
#'
#' @param low the lower bound, default \code{-Inf}
#' @param high the upper bound, default \code{Inf}
#' @return \code{TRUE} or \code{FALSE}
#'
#' @rdname ASSISTant-internal
#'
#' @keywords internal
numberInRange <- function(x, low = -Inf, high = Inf) {
  ## #'
  ## #' @examples
  ## #' ASSISTant:::numberInRange(x = 2:10, low = 2)  ## TRUE
  ## #' ASSISTant:::numberInRange(x = 10:15, high = 2) ## FALSE
  all((x >= low) & (x <= high))
}

#' Is a scalar quantity is an integer in specified range?
#' @description Check if the argument is a scalar integer within specified range
#'
#' @param low the lower bound, default \code{-Inf}
#' @param high the upper bound, default \code{Inf}
#' @return \code{TRUE} or \code{FALSE}
#'
#' @rdname ASSISTant-internal
#'
#' @keywords internal
scalarIntegerInRange <- function(x, low = -Inf, high = Inf) {
  ## #'
  ## #' @examples
  ## #' ASSISTant:::integerInRange(x = 10, low = 2)  ## TRUE
  ## #' ASSISTant:::integerInRange(x = 10.5) ## FALSE
  (length(x) == 1) && (x == trunc(x)) && (x >= low) && (x <= high)
}

#' Is the numeric quantity is an integer vector in specified range?
#' @description Check if the argument is an integer vector within specified range
#'
#' @param low the lower bound, default \code{-Inf}
#' @param high the upper bound, default \code{Inf}
#' @return \code{TRUE} or \code{FALSE}
#'
#' @rdname ASSISTant-internal
#'
#' @keywords internal
integerInRange <- function(x, low = -Inf, high = Inf) {
  ## #' @examples
  ## #' ASSISTant:::integerInRange(x = 2:10, low = 2)  ## TRUE
  ## #' ASSISTant:::integerInRange(x = 10:15, high = 2) ## FALSE
  ## #' ASSISTant:::integerInRange(x = c(0.5, 1.5), high = 2) ## FALSE
  ## #'
  all((x == trunc(x)) & (x >= low) & (x <= high))
}

