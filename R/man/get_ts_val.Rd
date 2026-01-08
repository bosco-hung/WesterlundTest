\name{get_ts_val}
\alias{get_ts_val}
\title{Strict Time-Series Mapping with Explicit Gaps}
\description{
Retrieves lagged (or led) values from a time-indexed vector using a strict
time-based mapping.
}
\usage{
get_ts_val(vec, tvec, lag)
}
\arguments{
  \item{vec}{A numeric (or atomic) vector of observations.}

  \item{tvec}{A vector of time indices corresponding one-to-one with \code{vec}.
  Values must uniquely identify time periods within the series.}

  \item{lag}{Integer. The lag order. Positive values correspond to lags
  (e.g., \code{lag = 1} for \eqn{t-1}), while negative values correspond to leads
  (e.g., \code{lag = -1} for \eqn{t+1}).}
}
\details{
This helper function implements a strict time-based lookup rather than
position-based indexing. Internally, it constructs a named mapping from time
indices to observed values:
\preformatted{
val_map <- setNames(vec, tvec)
}
For each observation at time \eqn{t}, the function computes the target time
\eqn{t - \text{lag}} and retrieves the corresponding value from the map.

If the target time does not exist in \code{tvec}, the function returns
\code{NA} for that observation. This behaviour exactly mirrors Stata’s
lag/lead operators in the presence of gaps.

Importantly, the function does \strong{not} interpolate or shift values when
time periods are missing. A gap in the time index propagates as \code{NA} in the
lagged (or led) series.
}
\value{
A vector of the same length as \code{vec}, containing the lagged (or led)
values aligned by time index. Elements are \code{NA} when the requested
time period does not exist.
}
\section{Implementation Details}{
This section explains why \code{get_ts_val()} is useful and how it differs from
standard lagging approaches based on vector positions.

\subsection{Why strict time-based lagging matters}{
In panel and time-series econometrics, lagged variables should be defined with
respect to the time index, not the row position. When data contain gaps in time,
simple shifting (e.g., \code{c(NA, x[-length(x)])}) can produce incorrect values.
\code{get_ts_val()} avoids this by explicitly matching on time values.
}

\subsection{Relation to Stata operators}{
This function replicates the behaviour of Stata’s time-series operators:
\itemize{
  \item \code{L.x}: lagged value at \eqn{t-1},
  \item \code{L2.x}: lagged value at \eqn{t-2},
  \item \code{F.x}: lead value at \eqn{t+1}.
}
When time periods are missing, Stata returns missing values rather than shifting
across gaps.
}
}
\examples{
## Example 1: Regular time series
t <- 1:5
x <- c(10, 20, 30, 40, 50)

## Lag by one period
get_ts_val(x, t, lag = 1)
# [1] NA 10 20 30 40

## Lead by one period
get_ts_val(x, t, lag = -1)
# [1] 20 30 40 50 NA

## Example 2: Time series with a gap
t_gap <- c(1, 2, 4, 5)
x_gap <- c(10, 20, 40, 50)

## Lag by one period: note the NA at time 4
get_ts_val(x_gap, t_gap, lag = 1)
# [1] NA 10 NA 40

## Explanation:
## - At t = 4, t - 1 = 3 does not exist in t_gap, so NA is returned.

## Example 3: Higher-order lags
get_ts_val(x_gap, t_gap, lag = 2)
# [1] NA NA NA 20
}
\seealso{
\code{\link{westerlund_test}},
\code{\link{calc_lrvar_bartlett}}
}
