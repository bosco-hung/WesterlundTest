\name{westerlund_test_reg}
\alias{westerlund_test_reg}
\title{Formatted Coefficient Table for Westerlund Test Reporting}
\description{
Builds a Stata-style coefficient table using either Z inference or t inference.
Mean-group reporting uses Z inference, while individual OLS ECM output uses t
inference with residual degrees of freedom.
}
\usage{
westerlund_test_reg(
  b,
  V,
  verbose = FALSE,
  distribution = c("z", "t"),
  df = NULL
)
}
\arguments{
\item{b}{Named numeric coefficient vector.}
\item{V}{Variance-covariance matrix corresponding to \code{b}.}
\item{verbose}{Logical. Print the coefficient table.}
\item{distribution}{Either \code{"z"} or \code{"t"}.}
\item{df}{Positive residual degrees of freedom, required for \code{distribution="t"}.}
}
\details{
Standard errors are the square roots of \code{diag(V)}. The reported statistic is
coefficient divided by standard error. Two-sided p-values and 95 percent
confidence intervals are computed from either the standard normal or Student t
distribution according to \code{distribution}.
}
\value{A numeric matrix containing coefficients, standard errors, the selected test statistic, two-sided p-values, and 95 percent confidence intervals.}
\examples{
b <- c(alpha = -0.2, beta = 1.05)
V <- diag(c(0.03^2, 0.10^2))
names(b) <- rownames(V) <- colnames(V) <- c("alpha", "beta")
westerlund_test_reg(b, V)
westerlund_test_reg(b, V, distribution = "t", df = 25)
}
\seealso{\code{\link{westerlund_test_mg}}, \code{\link{westerlund_test}}}
