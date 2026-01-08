\name{westerlund_test}
\alias{westerlund_test}
\title{Westerlund ECM-Based Panel Cointegration Test}
\description{
Implements the error-correction-based panel cointegration tests proposed by
Westerlund (2007). The function computes both mean-group (\eqn{G_t}, \eqn{G_a})
and pooled (\eqn{P_t}, \eqn{P_a}) statistics based on unit-specific ECMs, with
Stata-like input validation, time-continuity checks, and optional bootstrap
p-values for the raw statistics.
}
\usage{
westerlund_test(
  data,
  yvar,
  xvars,
  idvar,
  timevar,
  constant = FALSE,
  trend = FALSE,
  lags,
  leads = NULL,
  westerlund = FALSE,
  bootstrap = -1,
  noisily = FALSE,
  lrwindow = 2
)
}
\arguments{
  \item{data}{A \code{data.frame} containing the panel data. Panels may be unbalanced, but each unit's time index must be continuous (no holes) after removing missing observations on model variables.}

  \item{yvar}{String. Name of the dependent variable.}

  \item{xvars}{Character vector. Names of regressors entering the long-run relationship. A maximum of 6 regressors is allowed. If \code{westerlund=TRUE}, at most one regressor may be included.}

  \item{idvar}{String. Column identifying cross-sectional units. Missing values are not allowed.}

  \item{timevar}{String. Column identifying time. Within each unit, the time index must be strictly increasing and continuous (differences of 1).}

  \item{constant}{Logical. Include an intercept in the cointegrating relationship.}

  \item{trend}{Logical. Include a linear trend in the cointegrating relationship. Setting \code{trend=TRUE} requires \code{constant=TRUE}.}

  \item{lags}{Integer or length-2 integer vector. Fixed lag order or range \code{c(min,max)} for the short-run lag length \eqn{p}. This option must be provided.}

  \item{leads}{Integer or length-2 integer vector, or \code{NULL}. Fixed lead order or range \code{c(min,max)} for the short-run lead length \eqn{q}. If \code{NULL}, defaults to 0.}

  \item{westerlund}{Logical. If \code{TRUE}, enforces additional restrictions: at least a constant must be included (constant and/or trend) and at most one regressor may be specified.}

  \item{bootstrap}{Integer. If \code{bootstrap > 0}, runs an internal bootstrap routine with \code{bootstrap} replications and returns bootstrap p-values for the raw statistics. If \code{bootstrap <= 0} (default \code{-1}), no bootstrap is performed.}

  \item{noisily}{Logical. If \code{TRUE}, prints progress output from internal routines.}

  \item{lrwindow}{Integer. Window parameter used for long-run variance estimation in internal routines.}
}
\details{
The test is based on estimating unit-specific error-correction models (ECMs) and
testing the null hypothesis of no error correction for all cross-sectional units.
Four statistics are reported:
\itemize{
  \item \eqn{G_t}: mean of individual t-ratios of the error-correction coefficient,
  \item \eqn{G_a}: mean of individually scaled error-correction coefficients,
  \item \eqn{P_t}: pooled t-type statistic based on a common error-correction coefficient,
  \item \eqn{P_a}: pooled statistic based on the scaled common coefficient.
}

\strong{Input validation and data checks.}
The function enforces several Stata-like guardrails:
\itemize{
  \item A maximum of 6 regressors can be specified in \code{xvars}.
  \item If \code{trend=TRUE}, then \code{constant=TRUE} is required.
  \item If \code{westerlund=TRUE}, at least a constant must be included and \code{xvars} must contain at most one regressor.
  \item The panel must be declared by specifying \code{idvar} and \code{timevar}, and \code{idvar} may not contain missing values.
  \item Rows with missing values are excluded, and continuity checks are applied to the resulting usable sample.
}

\strong{Inference.}
If \code{bootstrap <= 0}, the function returns standardized Z-scores and
asymptotic (left-tail) p-values. If \code{bootstrap > 0},
it additionally produces bootstrap p-values for the raw statistics.
}
\value{
A list containing the following components:
\itemize{
  \item \code{test_stats}: A numeric vector containing:
    \itemize{
      \item \code{gt, ga, pt, pa}: Raw test statistics.
      \item \code{gt_z, ga_z, pt_z, pa_z}: Standardized Z-scores.
      \item \code{gt_pval, ga_pval, pt_pval, pa_pval}: Asymptotic left-tail p-values.
    }
  \item \code{boot_pvals}: A list containing \code{gt_pvalboot, ga_pvalboot, pt_pvalboot, pa_pvalboot}
    (only populated if \code{bootstrap > 0}).
  \item \code{bootstrap_distributions}: A matrix of the bootstrap replicates for each statistic.
  \item \code{unit_data}: A \code{data.frame} with unit-specific ECM results (alpha and beta estimates).
  \item \code{mean_group}: A list containing the Mean Group estimates (\code{mg_alpha}, \code{mg_betas})
    and their respective standard errors.
  \item \code{settings}: A list of the lag/lead specifications and internal parameters used.
}
}
\section{Vignette}{
This section demonstrates how to use \code{westerlund_test()} and interprets common outputs.

\subsection{What the test does}{
The test evaluates the null hypothesis of \strong{no error correction for all cross-sectional units}.
Rejection suggests that at least some units exhibit error-correcting behaviour.
}

\subsection{Data requirements: unbalanced panels and time continuity}{
Unbalanced panels are permitted as long as each unit's time index is continuous
after removing rows with missing values. If any unit has gaps (e.g., time jumps from 2001 to 2003),
the function stops and reports the unit identifier.
}

\subsection{Examples}{
\preformatted{
## Plain (asymptotic) test
res_plain <- westerlund_test(
  data       = df,
  yvar       = "y",
  xvars      = c("x1","x2"),
  idvar      = "id",
  timevar    = "t",
  constant   = TRUE,
  trend      = FALSE,
  lags       = 1,
  leads      = 0,
  westerlund = FALSE,
  bootstrap  = -1,
  noisily    = FALSE,
  lrwindow   = 2
)
}
}

\subsection{Interpreting common errors}{
\itemize{
  \item \strong{Continuous time-series are required}: at least one unit has holes
  in the time index. Fix by ensuring the usable sample is continuous.
  \item \strong{At least ... observations are required}: units are too short
  for the requested lag/lead orders.
  \item \strong{If a trend is included, a constant must be included}:
  set \code{constant=TRUE} when \code{trend=TRUE}.
}
}
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{WesterlundPlain}},
\code{\link{WesterlundBootstrap}},
\code{\link{DisplayWesterlund}}
}
