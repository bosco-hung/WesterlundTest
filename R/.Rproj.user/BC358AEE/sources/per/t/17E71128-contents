\name{get_lag}
\alias{get_lag}
\title{Time-Indexed Lag Extraction with Explicit Gaps}
\description{
Extracts lagged (or led) values from a time-indexed vector using a strict
time-based lookup. The function respects gaps in the time index and returns
\code{NA} when the requested lagged time does not exist, mirroring the behaviour
of Stata’s time-series lag/lead operators.
}
\usage{
get_lag(vec, tvec, k)
}
\arguments{
  \item{vec}{A numeric (or atomic) vector of observations.}

  \item{tvec}{A vector of time indices corresponding one-to-one with \code{vec}.
  Each value must uniquely identify a time period within the series.}

  \item{k}{Integer. The lag order. Positive values correspond to lags
  (e.g., \code{k = 1} for \eqn{t-1}), while negative values correspond to leads
  (e.g., \code{k = -1} for \eqn{t+1}).}
}
\details{
This helper function performs strict time-based lagging rather than
position-based shifting. Internally, it constructs a mapping from time indices
to observed values:
\preformatted{
val_map <- setNames(vec, tvec)
}
For each observation at time \eqn{t}, the function retrieves the value associated
with time \eqn{t - k}. If that time value is not present in \code{tvec}, the
result is \code{NA}.

This behaviour is particularly important when working with irregular time
series or panel data with gaps. In such cases, simple vector shifting can
incorrectly carry values across missing time periods, while \code{get_lag()}
preserves the correct alignment.

No interpolation, padding, or reordering is performed; missing time periods
propagate as \code{NA} in the lagged (or led) series.
}
\value{
A vector of the same length as \code{vec}, containing the lagged (or led) values
aligned by the time index. Elements are \code{NA} when the requested time period
does not exist.
}
\section{Technical Background}{
This section illustrates the logic of \code{get_lag()} and explains how
it differs from standard position-based lagging.

\subsection{Why not simple shifting?}{
In many applications, lags are computed by shifting vectors by positions.
This implicitly assumes a complete and regular time index. When time periods
are missing, such shifting produces incorrect lag values. \code{get_lag()}
avoids this by explicitly matching on time values.
}

\subsection{Relation to Stata operators}{
The function mirrors the behaviour of Stata’s time-series operators:
\itemize{
  \item \code{L.x}: lagged value at \eqn{t-1},
  \item \code{L2.x}: lagged value at \eqn{t-2},
  \item \code{F.x}: lead value at \eqn{t+1}.
}
As in Stata, gaps in the time index result in missing values in the lagged
series.
}
}
\examples{
## Example 1: Regular time series
t <- 1:5
x <- c(10, 20, 30, 40, 50)

get_lag(x, t, k = 1)
# [1] NA 10 20 30 40

get_lag(x, t, k = -1)
# [1] 20 30 40 50 NA

## Example 2: Time series with a gap
t_gap <- c(1, 2, 4, 5)
x_gap <- c(10, 20, 40, 50)

get_lag(x_gap, t_gap, k = 1)
# [1] NA 10 NA 40

## Explanation:
## At t = 4, the requested time t-1 = 3 does not exist, so NA is returned.

## Example 3: Higher-order lags
get_lag(x_gap, t_gap, k = 2)
# [1] NA NA NA 20
}
\seealso{
\code{\link{get_ts_val}},
\code{\link{westerlund_test}},
\code{\link{calc_lrvar_bartlett}}
}
