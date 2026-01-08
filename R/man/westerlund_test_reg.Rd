\name{westerlund_test_reg}
\alias{westerlund_test_reg}
\title{Formatted Coefficient Table for Westerlund Test Reporting}
\usage{
westerlund_test_reg(b, V)
}
\arguments{
\item{b}{A named numeric vector of coefficients.}

\item{V}{A numeric variance-covariance matrix corresponding to \code{b}.}
}
\value{
This function is called for its side effect (printing).
It invisibly returns \code{NULL}.
}
\description{
Internal helper that calculates standard errors, z-statistics, p-values, and
confidence intervals for a given coefficient vector and covariance matrix,
then prints a Stata-style regression table.
}
\details{
\strong{Calculation Logic.}
The function replicates the behavior of Stata's \code{ereturn post} when no
degrees of freedom are specified, using the standard normal distribution (Z)
for all inference:

\itemize{
  \item \strong{Standard Errors:} Derived as the square root of the diagonal of \code{V}.
  \item \strong{z-statistics:} Calculated as \eqn{z = \hat{\beta} / SE(\hat{\beta})}.
  \item \strong{P-values:} Two-tailed probabilities from the standard normal distribution.
  \item \strong{Confidence Intervals:} Calculated as \eqn{\hat{\beta} \pm 1.96 \times SE(\hat{\beta})}.
}

\strong{Formatting.}
The table is printed using \code{stats::printCoefmat()} to ensure clean
alignment and decimal consistency. It includes columns for the Coefficient,
Standard Error, z-statistic, P-value, and the 95\% Confidence Interval.

\strong{Intended use.}
This is an internal utility called by \code{\link{westerlund_test_mg}}. It
provides a standardized way to report results across different parts of the
Westerlund cointegration test output.
}
\section{Reporting Style}{
This section describes the alignment with econometric software reporting standards.

\subsection{Stata Compatibility}{
The output format (column headers and statistical assumptions) is designed to
match the output seen in Stata's \code{xtcpmg} or \code{xtwest} routines. This
ensures that users transitioning from or comparing results with Stata find the
output familiar.
}
}
\examples{
## Define simple coefficient vector and VCV
b <- c(alpha = -0.20, beta = 1.05)
V <- diag(c(0.03^2, 0.10^2))
names(b) <- rownames(V) <- colnames(V) <- c("alpha", "beta")

## Print the formatted table
westerlund_test_reg(b, V)
}
\seealso{
\code{\link{westerlund_test_mg}},
\code{\link{westerlund_test}}
}