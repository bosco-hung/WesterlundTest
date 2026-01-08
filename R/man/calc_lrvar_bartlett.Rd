\name{calc_lrvar_bartlett}
\alias{calc_lrvar_bartlett}
\title{Long-Run Variance Estimation with Bartlett Kernel}
\description{
Computes a Bartlett-kernel (Newey--West style) long-run variance estimate for a
univariate series using Stata-like conventions: missing values are removed
(\code{na.omit}), autocovariances are scaled by \eqn{1/n} (not \eqn{1/(n-j)}),
and optional centering is controlled by \code{nodemean}.
}
\usage{
calc_lrvar_bartlett(x, maxlag, nodemean = FALSE)
}
\arguments{
  \item{x}{A numeric vector. Missing values are removed prior to computation.}

  \item{maxlag}{Non-negative integer. Maximum lag order \eqn{m} used in the Bartlett kernel.
  If \code{maxlag=0}, the function returns the sample variance estimate \eqn{\gamma_0}.}

  \item{nodemean}{Logical. If \code{FALSE} (default), the series is centered by subtracting its mean
  before computing autocovariances. If \code{TRUE}, the function assumes \code{x} is already centered
  and does not de-mean it.}
}
\details{
Let \eqn{x_t} be the input series after removing missing values. If
\code{nodemean=FALSE}, the function replaces \eqn{x_t} with
\eqn{x_t - \bar{x}}.

Define the lag-\eqn{j} autocovariance using the Stata-style scaling:
\deqn{
\gamma_j = \frac{1}{n}\sum_{t=j+1}^{n} x_t x_{t-j}, \qquad j=0,1,\dots,m,
}
where \eqn{n} is the length of the cleaned series and \eqn{m=\code{maxlag}}.
Note that the scaling uses \eqn{1/n} for all \eqn{j} (rather than \eqn{1/(n-j)}).

The Bartlett kernel weight at lag \eqn{j} is:
\deqn{
w_j = 1 - \frac{j}{m+1}.
}
The long-run variance estimate is then:
\deqn{
\widehat{\Omega} = \gamma_0 + 2\sum_{j=1}^{m} w_j \gamma_j.
}

If the cleaned series has zero length, the function returns \code{NA}.
}
\value{
A single numeric value:
\itemize{
  \item The Bartlett-kernel long-run variance estimate \eqn{\widehat{\Omega}} (scalar),
  \item or \code{NA} if \code{x} contains no non-missing values.
}
}
\section{Technical Notes}{
This section illustrates how \code{calc_lrvar_bartlett()} computes a Bartlett-kernel
long-run variance (LRV) estimate and how its options map to common time-series
preprocessing choices.

\subsection{Stata-like conventions}{
This helper function is designed to match Stata-style calculations:
\itemize{
  \item \strong{Missing values are dropped}: \code{na.omit(x)} is applied first.
  \item \strong{Scaling by \eqn{1/n}}: autocovariances at all lags divide by \eqn{n}, not by \eqn{n-j}.
  \item \strong{Optional de-meaning}: controlled by \code{nodemean}.
}
}

\subsection{Centering}{
By default (\code{nodemean=FALSE}), the series is centered: \eqn{x_t \leftarrow x_t - \bar{x}}.
Set \code{nodemean=TRUE} when you have already centered the series elsewhere.
}

\subsection{Choosing maxlag}{
\code{maxlag} sets the truncation point \eqn{m}. Larger \eqn{m} captures more
serial correlation but increases estimation noise.
}
}
\examples{
## Example 1: Basic usage
x <- rnorm(200)
calc_lrvar_bartlett(x, maxlag = 4)

## Example 2: maxlag = 0 returns gamma_0
calc_lrvar_bartlett(x, maxlag = 0)

## Example 3: Handle missing values (they are removed)
x_na <- x
x_na[c(5, 10, 50)] <- NA
calc_lrvar_bartlett(x_na, maxlag = 4)

## Example 4: Compare centering choices
## Default: de-mean internally
lr1 <- calc_lrvar_bartlett(x, maxlag = 4, nodemean = FALSE)

## Pre-center and skip de-meaning
x_centered <- x - mean(x)
lr2 <- calc_lrvar_bartlett(x_centered, maxlag = 4, nodemean = TRUE)
}
\seealso{
\code{\link{westerlund_test}}
}
