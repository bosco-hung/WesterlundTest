\name{WesterlundBootstrap}
\alias{WesterlundBootstrap}
\title{Bootstrap Routine for Westerlund ECM Panel Cointegration Tests}
\description{
Internal bootstrap routine for generating a bootstrap distribution of the
Westerlund (2007) ECM-based panel cointegration statistics under the null
hypothesis of no error correction. The routine estimates short-run dynamics for
each unit, constructs residual-based bootstrap samples in a Stata-aligned manner,
and re-computes \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a} by calling
\code{\link{WesterlundPlain}} on each bootstrap sample.
}
\usage{
WesterlundBootstrap(
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
  westerlund = FALSE,
  aic = TRUE,
  bootstrap = 100,
  indiv.ecm = FALSE,
  lrwindow = 2,
  verbose = FALSE
)
}
\arguments{
  \item{data}{A \code{data.frame} containing the original panel data.}

  \item{touse}{Logical vector of length \code{nrow(data)} indicating rows eligible for estimation. The routine further removes rows with missing \code{yvar} or \code{xvars}.}

  \item{idvar}{String. Column identifying cross-sectional units.}

  \item{timevar}{String. Column identifying time. Within-unit time gaps are handled by Stata-style differencing in the bootstrap setup.}

  \item{yvar}{String. Name of the dependent variable (levels) in the original data.}

  \item{xvars}{Character vector. Names of regressors (levels) in the original data.}

  \item{constant}{Logical. Include a constant term in the short-run regression used to obtain bootstrap residuals and coefficients.}

  \item{trend}{Logical. Include a trend term in the short-run regression used to obtain bootstrap residuals and coefficients.}

  \item{lags}{Integer or length-2 integer vector. Fixed lag order or range \code{c(min,max)} used for selecting short-run dynamics in the bootstrap setup.}

  \item{leads}{Integer or length-2 integer vector, or \code{NULL}. Fixed lead order or range \code{c(min,max)} used for selecting short-run dynamics. If \code{NULL}, defaults to 0.}

  \item{westerlund}{Logical. If \code{TRUE}, uses a Westerlund-style information criterion and trimming logic in the bootstrap setup, matching the \code{westerlund} mode used elsewhere.}

  \item{aic}{Logical. If \code{TRUE}, uses AIC for lag/lead selection when ranges. If \code{FALSE}, uses BIC.}

  \item{bootstrap}{Integer. Number of bootstrap replications.}

  \item{indiv.ecm}{Logical. If \code{TRUE}, gets output of individual ECM regressions..}

  \item{lrwindow}{Integer. Bartlett kernel window (maximum lag) forwarded to \code{\link{WesterlundPlain}} during bootstrap re-estimation.}

  \item{verbose}{Logical. If \code{TRUE}, prints additional output.}
}
\details{
\strong{Purpose and status.}
\code{WesterlundBootstrap()} is an internal engine typically called by
\code{\link{westerlund_test}} when bootstrap inference is requested.
It returns a matrix of bootstrap test statistics that can be used to compute
bootstrap p-values for the \emph{raw} statistics.

\strong{High-level algorithm.}
For each bootstrap replication, the routine:
\enumerate{
  \item Filters the original data using \code{touse} and removes missing \code{yvar}/\code{xvars}.
  \item For each unit, estimates a short-run model for \eqn{\Delta y_t} using Stata-style differencing that respects time gaps:
  if \eqn{t} is not consecutive, the corresponding difference is set to missing.
  \item Stores the estimated coefficients for lagged \eqn{\Delta y_t} and leads/lags/current \eqn{\Delta x_t}, and constructs residuals \eqn{e}.
  \item Demeans residuals \eqn{e} within each unit and constructs within-unit centered differenced regressors on the residual-support rows.
  \item Resamples residual-support clusters to form bootstrap innovations, recursively generates bootstrap \eqn{\Delta y_t} with an AR recursion implied by the estimated \eqn{\Delta y} coefficients, integrates to levels to obtain bootstrap \code{booty} and \code{bootx} series, and enforces Stata-style truncation rules.
  \item Calls \code{\link{WesterlundPlain}} on the constructed bootstrap panel to compute \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a}.
}
The resulting statistics are stored row-by-row in \code{BOOTSTATS}.

\strong{Time indexing in bootstrap setup.}
Unlike the main estimation routines (which use strict time-indexed helpers), this bootstrap routine uses a local differencing helper \code{diff_ts()} that mimics Stata's \code{D.} operator under gaps: differences are marked \code{NA} when \code{diff(timevec) != 1}.

\strong{Lag/lead selection in the bootstrap setup.}
If \code{lags} and/or \code{leads} are supplied as ranges, the routine performs an information-criterion search over candidate short-run specifications for each unit. In \code{westerlund=TRUE} mode, a Westerlund-style criterion is used; otherwise, a standard AIC-based criterion is applied.

\strong{Random number generation.}
The routine sets the RNG kind to Mersenne-Twister and uses a deterministic sequence given the current RNG state. Users should set \code{set.seed()} upstream if reproducibility is desired.
}
\value{
A list containing:
\itemize{
  \item \code{BOOTSTATS}: a numeric matrix with \code{bootstrap} rows and 4 columns,
  containing bootstrap replications of \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a}
  (in that order).
}
}
\section{Vignette}{
This section describes how \code{WesterlundBootstrap()} constructs bootstrap
samples and how it is typically used in the overall testing workflow.

\subsection{How does this relate to \code{westerlund_test()}?}{
When the user-facing \code{\link{westerlund_test}} is called with \code{bootstrap > 0}, it typically:
\enumerate{
  \item obtains a bootstrap distribution via \code{WesterlundBootstrap()},
  \item computes the observed statistics via \code{\link{WesterlundPlain}},
  \item derives bootstrap p-values by comparing observed raw statistics to the bootstrap distribution.
}
}

\subsection{Residual-based bootstrap under the null}{
The routine first estimates a short-run model for \eqn{\Delta y_t} within each
unit and extracts residuals. It then resamples these residuals (clustered on the
residual-support rows) and uses the stored short-run coefficients to propagate
bootstrap dynamics.
}

\subsection{Handling time gaps}{
To mirror Stata semantics in differencing, the bootstrap setup uses a local
difference rule: if two adjacent observations are not exactly one time unit apart,
the difference is set to missing for that location.
}

\subsection{Example}{
\preformatted{
## Upstream: set a seed for reproducibility
set.seed(123)

boot_res <- WesterlundBootstrap(
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
  westerlund = FALSE,
  bootstrap  = 399,
  indiv.ecm  = TRUE,
  lrwindow   = 2,
  verbose    = FALSE
)

## Example left-tail p-value for Gt (assuming 'obs' exists)
# mean(boot_res$BOOTSTATS[,1] <= obs$Gt, na.rm = TRUE)
}
}
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{westerlund_test}},
\code{\link{WesterlundPlain}}
}
