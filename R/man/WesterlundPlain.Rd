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
  aic = TRUE,
  bootno = FALSE,
  indiv.ecm = FALSE,
  verbose = FALSE
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

  \item{aic}{Logical. If \code{TRUE}, uses AIC for lag/lead selection when ranges. If \code{FALSE}, uses BIC.}

  \item{bootno}{Logical. If \code{TRUE}, prints a short header and progress dots (intended for higher-level routines).}

  \item{indiv.ecm}{Logical. If \code{TRUE}, gets output of individual ECM regressions.}

  \item{verbose}{Logical. If \code{TRUE}, prints additional output.}
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
}
\value{
A nested list containing:
\itemize{
  \item \code{stats}: A list of the four raw Westerlund test statistics:
    \itemize{
      \item \code{Gt}: Mean-group tau statistic.
      \item \code{Ga}: Mean-group alpha statistic.
      \item \code{Pt}: Pooled tau statistic.
      \item \code{Pa}: Pooled alpha statistic.
    }
  \item \code{indiv_data}: A named list where each element corresponds to a cross-sectional unit (ID), containing:
     \itemize{
       \item \code{ai}: The estimated speed of adjustment (alpha).
       \item \code{seai}: The standard error of alpha (adjusted for degrees of freedom).
       \item \code{betai}: Vector of long-run coefficients (\eqn{\beta = -\lambda / \alpha}).
       \item \code{blag, blead}: The lags and leads selected for that specific unit.
       \item \code{ti}: Raw observation count for the unit.
       \item \code{tnorm}: Degrees of freedom used for normalization.
       \item \code{reg_coef}: If \code{indiv.ecm = TRUE}, the full coefficient matrix from \code{westerlund_test_reg}.
     }
  \item \code{results_df}: A summary \code{data.frame} containing all unit-level results in vectorized format.
  \item \code{settings}: A list of routine metadata:
    \itemize{
      \item \code{meanlag}, \code{meanlead}: Integer averages of the selected unit lags/leads.
      \item \code{realmeanlag}, \code{realmeanlead}: Numeric averages of the selected unit lags/leads.
      \item \code{auto}: Logical; \code{TRUE} if automatic selection (ranges) was used.
    }
}
}
\section{Internal Logic}{
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
}
\examples{
\donttest{
set.seed(123)
N <- 5
T <- 20
df <- data.frame(
  id = rep(1:N, each = T),
  t  = rep(1:T, N),
  y  = rnorm(N * T),
  x1 = rnorm(N * T),
  x2 = rnorm(N * T)
)

touse <- rep(TRUE, nrow(df))

plain_res <- WesterlundPlain(
  data       = df,
  touse      = touse,
  idvar      = "id",
  timevar    = "t",
  yvar       = "y",
  xvars      = c("x1","x2"),
  lags       = 1,
  leads      = 0
)

# Accessing results from the nested structure:
stats <- plain_res$stats
print(c(Gt = stats$Gt, Ga = stats$Ga, Pt = stats$Pt, Pa = stats$Pa))

# Checking unit-specific coefficients for ID '101'
unit_101 <- plain_res$indiv_data[["101"]]
print(unit_101$ai)
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
