\name{get_diff}
\alias{get_diff}
\title{First Difference with Strict Time Indexing}
\description{
Computes the first difference of a time-indexed series using strict time-based
lagging. The function respects gaps in the time index and returns \code{NA} when
the previous time period does not exist, mirroring Stata’s \code{D.} operator.
}
\usage{
get_diff(vec, tvec)
}
\arguments{
  \item{vec}{A numeric (or atomic) vector of observations.}

  \item{tvec}{A vector of time indices corresponding one-to-one with \code{vec}.
  Each value must uniquely identify a time period within the series.}
}
\details{
This helper function computes first differences as:
\deqn{
\Delta x_t = x_t - x_{t-1},
}
where the lagged value \eqn{x_{t-1}} is obtained using \code{\link{get_lag}},
which performs strict time-based lookup.

Internally, the function calls:
\preformatted{
val_t_minus_1 <- get_lag(vec, tvec, 1)
}
and then subtracts this lagged vector from \code{vec}. If the time index contains
gaps, or if the previous time period does not exist for a given observation,
the lagged value is \code{NA} and the corresponding difference is also
\code{NA}.

No interpolation or implicit shifting is performed; missing time periods
propagate as missing differences.
}
\value{
A vector of the same length as \code{vec}, containing the first differences
aligned by the time index. Elements are \code{NA} when the previous time period
does not exist.
}
\section{Time Indexing Logic}{
This section explains how \code{get_diff()} computes first differences and why
strict time indexing matters in the presence of gaps.

\subsection{Relation to Stata’s D. operator}{
The function replicates the behaviour of Stata’s first-difference operator
\code{D.x}. When time periods are missing, Stata returns missing values rather
than differencing across gaps. Because \code{get_diff()} relies on
\code{\link{get_lag}}, it follows the same rule.
}

\subsection{Why not use diff()?}{
The base R function \code{diff()} computes differences based on vector positions.
This implicitly assumes a complete and regularly spaced time index. When time
periods are missing, \code{diff()} can produce misleading results by differencing
across gaps. \code{get_diff()} avoids this by differencing only when the
previous time period exists.
}
}
\examples{
## Example 1: Regular time series
t <- 1:5
x <- c(10, 20, 30, 40, 50)

get_diff(x, t)
# [1] NA 10 10 10 10

## Example 2: Time series with a gap
t_gap <- c(1, 2, 4, 5)
x_gap <- c(10, 20, 40, 50)

get_diff(x_gap, t_gap)
# [1] NA 10 NA 10

## Explanation:
## At t = 4, the previous period t-1 = 3 does not exist, so the difference is NA.

## Example 3: Comparison with diff()
diff(x_gap)
# [1] 10 20 10
}
\seealso{
\code{\link{get_lag}},
\code{\link{get_ts_val}},
\code{\link{westerlund_test}}
}
