% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/riley.r
\name{summary.riley}
\alias{summary.riley}
\title{Parameter summaries
Provides the summary estimates of the alternative model for bivariate random-effects meta-analysis by Riley et al. 
(2008) with their corresponding standard errors (derived from the inverse Hessian). For confidence intervals,
asymptotic normality is assumed.}
\usage{
\method{summary}{riley}(object, ...)
}
\arguments{
\item{object}{A \code{riley} object}

\item{\dots}{Arguments to be passed on to other functions (currently ignored)}
}
\value{
array with confidence intervals for the estimated model parameters. For diagnostic test accuracy data, 
the resulting summary sensitivity and false positive rate are included.
}
\description{
Parameter summaries
Provides the summary estimates of the alternative model for bivariate random-effects meta-analysis by Riley et al. 
(2008) with their corresponding standard errors (derived from the inverse Hessian). For confidence intervals,
asymptotic normality is assumed.
}
\details{
For meta-analysis of diagnostic test accuracy data, \code{beta1} equals the logit sensitivity (Sens) and 
\code{beta2} equals the logit false positive rate (FPR).
}
\note{
For the overall correlation (\code{rho}) confidence intervals are derived using the transformation 
\code{logit((rho+1)/2)}. Similarly, the logit transformation is used to derive confidence intervals for the summary
sensitivity and false positive rate.
}
\references{
Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the 
within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}
\author{
Thomas Debray <thomas.debray@gmail.com>
}
