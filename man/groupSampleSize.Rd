% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{groupSampleSize}
\alias{groupSampleSize}
\title{Compute the sample size for any group at a stage assuming a nested
structure as in the paper.}
\usage{
groupSampleSize(
  prevalence,
  N,
  stage,
  group,
  HJFutileAtStage = NA,
  chosenGroup = NA
)
}
\arguments{
\item{prevalence}{the vector of prevalence, will be normalized if
not already so. The length of this vector implicitly indicates
the number of groups J.}

\item{N}{an integer vector of length 3 indicating total sample size
at each of the three stages}

\item{stage}{the stage of the trial}

\item{group}{the group whose sample size is desired}

\item{HJFutileAtStage}{is the stage at which overall futility
occured. Default \code{NA} indicating it did not occur. Also
ignored if stage is 1.}

\item{chosenGroup}{the selected group if HJFutilityAtStage is not
\code{NA}. Ignored if stage is 1.}
}
\value{
the sample size for group
}
\description{
In the three stage design under consideration, the groups are
nested with assumed prevalences and fixed total sample size at each
stage. This function returns the sample size for a specified group
at a given stage, where the futility stage for the overall group
test may be specified along with the chosen subgroup.
}
