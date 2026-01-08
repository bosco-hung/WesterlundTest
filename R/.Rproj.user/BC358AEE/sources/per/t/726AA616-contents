\name{WesterlundPlain}
\alias{WesterlundPlain}
\title{Compute Raw Westerlund ECM Panel Cointegration Statistics (Plain Routine)}
\description{
Internal plain (non-bootstrap) routine for computing the four Westerlund (2007)
ECM-based panel cointegration test statistics \eqn{G_t}, \eqn{G_a}, \eqn{P_t},
and \eqn{P_a}. The function estimates unit-specific ECM regressions to form the
mean-group statistics and then constructs pooled (panel) statistics using
cross-unit aggregation and partialling-out steps. Time indexing is handled
strictly via gap-aware lag/difference helpers.
}
\usage{
WesterlundPlain(
  data,
  touse,
  idvar,
  timevar,
  yvar,
  xvars,
  constant = FALSE,
  trend = FALSE,
  lags,
  leads = NULL,
  lrwindow = 2,
  westerlund = FALSE,
  bootno = FALSE,
  noisily = FALSE,
  debug = FALSE,
  mg = FALSE,
  detrend = FALSE
)
}
\arguments{
  \item{data}{A \code{data.frame} containing panel data.}

  \item{touse}{Logical vector of length \code{nrow(data)} indicating rows eligible for estimation. Rows are further filtered to remove missing \code{yvar} and \code{xvars}.}

  \item{idvar}{String. Column identifying cross-sectional units.}

  \item{timevar}{String. Column identifying time.}

  \item{yvar}{String. Name of the dependent variable (levels).}

  \item{xvars}{Character vector. Names of regressors in the long-run relationship (levels).}

  \item{constant}{Logical. If \code{TRUE}, includes a constant term in the ECM design matrix.}

  \item{trend}{Logical. If \code{TRUE}, includes a linear time trend in the ECM design matrix.}

  \item{lags}{Integer or length-2 integer vector. Fixed lag order or range \code{c(min,max)} for short-run dynamics. If a range is supplied, the routine performs an information-criterion search over candidate lag/lead combinations.}

  \item{leads}{Integer or length-2 integer vector, or \code{NULL}. Fixed lead order or range \code{c(min,max)}. If \code{NULL}, defaults to 0.}

  \item{lrwindow}{Integer. Bartlett kernel window (maximum lag) used in long-run variance calculations via \code{\link{calc_lrvar_bartlett}}.}

  \item{westerlund}{Logical. If \code{TRUE}, uses a Westerlund-specific information criterion and trimming logic for variance estimation.}

  \item{bootno}{Logical. If \code{TRUE}, prints a short header and progress dots (intended for higher-level routines).}

  \item{noisily}{Logical. If \code{TRUE}, allows more console output during computation.}

  \item{debug}{Logical. If \code{TRUE}, suppresses progress dots in some branches and can be used for debugging prints.}

  \item{mg}{Logical. Reserved for future use / compatibility; currently not used.}

  \item{detrend}{Logical. Reserved for future use / compatibility; currently not used.}
}
\details{
\strong{Purpose and status.}
\code{WesterlundPlain()} is typically called internally by \code{\link{westerlund_test}}.
It returns the four \emph{raw} test statistics and lag/lead diagnostics needed
for printing and standardization.

\strong{Workflow overview.}
The routine proceeds in two main stages:

\enumerate{
  \item \strong{Unit-specific ECM regressions (Loop 1):} For each cross-sectional unit, it constructs an ECM with
  \eqn{\Delta y_t} as the dependent variable and includes deterministic terms (optional), \eqn{y_{t-1}},
  \eqn{x_{t-1}}, lagged \eqn{\Delta y_t}, and leads/lags of \eqn{\Delta x_t}. Lags and leads are computed using
  strict time-indexed helpers (\code{\link{get_lag}}, \code{\link{get_diff}}), which respect gaps in the time index.
  If \code{lags} and/or \code{leads} are provided as ranges, an information-criterion search selects the
  lag/lead orders for each unit. The routine stores the unit-level error-correction estimate \eqn{\hat{\alpha}_i}
  and its standard error.
  \item \strong{Pooled (panel) aggregation (Loop 2):} Using the mean of selected lag/lead orders across units,
  the routine constructs pooled quantities needed for \eqn{P_t} and \eqn{P_a} via partialling-out regressions
  and cross-unit aggregation of residual products.
}

\strong{Long-run variance calculations.}
Long-run variances are computed using \code{\link{calc_lrvar_bartlett}} with
\code{maxlag = lrwindow}. In \code{westerlund=TRUE} mode, the routine applies
Stata-like trimming at the start/end of the differenced series based on selected
lags/leads prior to long-run variance estimation.

\strong{Returned statistics.}
Let \eqn{\hat{\alpha}_i} denote the unit-specific error-correction coefficient
on \eqn{y_{t-1}} (as constructed in the ECM), with standard error \eqn{\widehat{\mathrm{se}}(\hat{\alpha}_i)}.
The routine computes:
\itemize{
  \item \eqn{G_t}: the mean of the individual t-ratios \eqn{\hat{\alpha}_i/\widehat{\mathrm{se}}(\hat{\alpha}_i)},
  \item \eqn{G_a}: a scaled mean-group statistic using a unit-specific normalization factor derived from long-run variances,
  \item \eqn{P_t}: a pooled t-type statistic based on a pooled \eqn{\hat{\alpha}} and its pooled standard error,
  \item \eqn{P_a}: a pooled scaled statistic using an average effective time dimension.
}

\strong{Notes on deterministics and trend handling.}
The function constructs \code{cons} and \code{tren} columns in the working data.
When \code{constant=TRUE} and \code{trend=TRUE} in non-\code{westerlund} mode,
the routine includes a detrending/centering adjustment in the long-run variance
preprocessing steps (as implemented in the code).

\strong{Intended use.}
For user-facing workflows, call \code{\link{westerlund_test}}; it handles
input checks, time-continuity checks, minimum-observation checks, optional
bootstrapping, and formatting of results.
}
\value{
A list containing:
\itemize{
  \item \code{Gt}, \code{Ga}, \code{Pt}, \code{Pa}: raw Westerlund test statistics,
  \item \code{meanlag}, \code{meanlead}: integer means of selected unit-level lags/leads,
  \item \code{realmeanlag}, \code{realmeanlead}: non-integer means of selected unit-level lags/leads,
  \item \code{auto}: logical indicating whether lag/lead ranges were supplied (automatic selection used).
}
}
\section{Internal Logic}{
This section provides a high-level guide to what \code{WesterlundPlain()} does
internally.

\subsection{Two-stage structure}{
\strong{Loop 1 (mean-group)} estimates unit-specific ECMs. Each unit produces an
estimated error-correction coefficient on \eqn{y_{t-1}} and an associated standard
error. These are aggregated into \eqn{G_t} and \eqn{G_a}.

\strong{Loop 2 (pooled)} fixes a common short-run structure based on the average
selected lag/lead orders and constructs pooled residual products to obtain \eqn{P_t} and \eqn{P_a}.
}

\subsection{Strict time indexing and gaps}{
All lags and differences are computed using strict time-based helpers
(\code{\link{get_lag}}, \code{\link{get_diff}}). This ensures that gaps in the
time index propagate as missing values rather than shifting across gaps.
}

\subsection{Automatic lag/lead selection}{
If \code{lags} and/or \code{leads} are supplied as ranges, the routine performs a grid search
over candidate pairs and selects the combination minimizing the information criterion.
}
}
\examples{
\dontrun{
## Suppose `touse` has been constructed and inputs validated upstream.
plain_res <- WesterlundPlain(
  data       = df,
  touse      = touse,
  idvar      = "id",
  timevar    = "t",
  yvar       = "y",
  xvars      = c("x1","x2"),
  constant   = TRUE,
  trend      = FALSE,
  lags       = 1,
  leads      = 0,
  lrwindow   = 2,
  westerlund = FALSE,
  bootno     = TRUE,
  noisily    = FALSE
)

c(Gt = plain_res$Gt, Ga = plain_res$Ga, Pt = plain_res$Pt, Pa = plain_res$Pa)
}
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{westerlund_test}},
\code{\link{WesterlundBootstrap}},
\code{\link{get_lag}},
\code{\link{get_diff}},
\code{\link{calc_lrvar_bartlett}}
}
