\name{DisplayWesterlund}
\alias{DisplayWesterlund}
\title{Standardize and Display Westerlund ECM Panel Cointegration Test Results}
\description{
Formats, standardizes, and prints the Westerlund (2007) ECM-based panel
cointegration test results for the null hypothesis of no cointegration. Given
raw statistics \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a}, the function
computes standardized Z-statistics using tabulated (or hard-coded) asymptotic
moments and reports left-tail asymptotic p-values. If a bootstrap distribution
is supplied, it also computes bootstrap (robust) p-values for the raw statistics.
}
\usage{
DisplayWesterlund(
  stats,
  bootstats = NULL,
  nobs,
  nox,
  constant = FALSE,
  trend = FALSE,
  meanlag = -1,
  meanlead = -1,
  realmeanlag = -1,
  realmeanlead = -1,
  auto = 0,
  westerlund = FALSE,
  mg = FALSE
)
}
\arguments{
  \item{stats}{A numeric vector of length 4 or a \code{1x4} numeric matrix containing the raw test statistics in the order \code{c(Gt, Ga, Pt, Pa)}.}

  \item{bootstats}{Optional. A numeric matrix with bootstrap replications of the raw statistics, with 4 columns in the order \code{[Gt, Ga, Pt, Pa]}. If provided, robust (bootstrap) p-values are computed and displayed.}

  \item{nobs}{Integer. Number of valid cross-sectional units (series) used in the test.}

  \item{nox}{Integer. Number of covariates in the long-run relationship (length of \code{xvars} in upstream calls). Used to select asymptotic moments for standardization in the non-\code{westerlund} case.}

  \item{constant}{Logical. Indicates whether a constant was included in the cointegrating relationship. Affects the deterministic-case row used for asymptotic moment lookup.}

  \item{trend}{Logical. Indicates whether a trend was included. If \code{TRUE}, \code{constant} is assumed to be \code{TRUE} and the deterministic-case row for moment lookup changes accordingly.}

  \item{meanlag}{Integer. Mean selected lag length (rounded down to integer) reported in the header when \code{auto} is non-zero.}

  \item{meanlead}{Integer. Mean selected lead length (rounded down to integer) reported in the header when \code{auto} is non-zero.}

  \item{realmeanlag}{Numeric. Unrounded mean selected lag length (for display when \code{auto} is non-zero).}

  \item{realmeanlead}{Numeric. Unrounded mean selected lead length (for display when \code{auto} is non-zero).}

  \item{auto}{Logical/integer. If non-zero, the header prints the average AIC-selected lag and lead lengths based on \code{realmeanlag} and \code{realmeanlead}.}

  \item{westerlund}{Logical. If \code{TRUE}, uses the hard-coded Westerlund-specific moment constants for standardization (separate constants for trend vs no-trend). If \code{FALSE}, uses lookup tables indexed by deterministic case and number of covariates.}

  \item{mg}{Logical. Reserved for future use / compatibility; currently not used in the computation.}
}
\details{
\strong{What this function does.}
\code{DisplayWesterlund()} takes raw Westerlund ECM test statistics and produces:
\itemize{
  \item standardized Z-statistics for \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, \eqn{P_a},
  \item left-tail asymptotic p-values using \code{pnorm()},
  \item optionally, bootstrap (robust) p-values when \code{bootstats} is provided,
  \item a console table summarizing the results.
}

\strong{Standardization and deterministic cases.}
The function selects a deterministic-case index:
\itemize{
  \item Row 1: no constant, no trend,
  \item Row 2: constant, no trend,
  \item Row 3: constant and trend,
}
implemented as \code{row_idx = (as.integer(constant) + as.integer(trend)) + 1}.
In the non-\code{westerlund} case, asymptotic means and variances are taken from
hard-coded lookup matrices indexed by \code{row_idx} and \code{nox}. Z-statistics are computed using:
\deqn{
Z = \frac{\sqrt{N}S - \sqrt{N}\mu}{\sqrt{\sigma^2}}
}
for mean-group statistics and analogous formulas for pooled statistics as
implemented in the code, where \eqn{N =} \code{nobs} and \eqn{S} is the raw statistic.

In the \code{westerlund=TRUE} case, the function uses a separate set of
hard-coded mean/variance constants, with different values depending on whether
\code{trend} is included.

\strong{Asymptotic p-values.}
All asymptotic p-values are computed as left-tail probabilities:
\code{pnorm(Z)}.

\strong{Bootstrap p-values (robust p-values).}
If \code{bootstats} is supplied, bootstrap p-values are computed for each raw
statistic using a left-tail empirical rule:
\deqn{
\hat{p} = \frac{r + 1}{B + 1}, \qquad r = \sum_{b=1}^{B} \mathbb{I}(S_b \le S_{\text{obs}}),
}
after dropping non-finite bootstrap draws. This is a common finite-sample
correction used in Stata-style bootstrap code.

\strong{Return values.}
In addition to printing, the function returns a list containing the raw
statistics, Z-statistics, asymptotic p-values, and (if applicable) bootstrap
p-values.
}
\value{
A list containing at minimum:
\itemize{
  \item \code{gt, ga, pt, pa}: raw test statistics,
  \item \code{gt_z, ga_z, pt_z, pa_z}: standardized Z-statistics,
  \item \code{gt_pval, ga_pval, pt_pval, pa_pval}: left-tail asymptotic p-values.
}
If \code{bootstats} is provided, the list also includes:
\itemize{
  \item \code{gt_pvalboot, ga_pvalboot, pt_pvalboot, pa_pvalboot}: bootstrap (robust) p-values for the raw statistics.
}
}
\section{Vignette}{
This section explains how to use \code{DisplayWesterlund()} and how its output
connects to the broader testing workflow.

\subsection{Where does \code{DisplayWesterlund()} fit?}{
Typically, the workflow is:
\enumerate{
  \item Compute observed raw statistics via \code{\link{WesterlundPlain}}.
  \item Optionally compute a bootstrap distribution via \code{\link{WesterlundBootstrap}}.
  \item Call \code{DisplayWesterlund()} to standardize, print, and return p-values.
}
The user-facing \code{\link{westerlund_test}} wraps these steps and collects the
returned scalars.
}

\subsection{Input format for \code{stats}}{
\code{stats} can be either a numeric vector \code{c(Gt, Ga, Pt, Pa)} or a \code{1x4} matrix with \code{[1,]} corresponding to \code{Gt, Ga, Pt, Pa}.
}

\subsection{Asymptotic standardization and left-tail p-values}{
The function converts raw statistics to Z-statistics using hard-coded asymptotic
means and variances. P-values are computed as \code{pnorm(Z)}, i.e., a left-tail test.
}

\subsection{Bootstrap (robust) p-values}{
If \code{bootstats} is provided, robust p-values are computed by comparing the
observed raw statistic to its bootstrap distribution, using a finite-sample
correction:
\eqn{(r+1)/(B+1)} where \eqn{r} is the number of bootstrap draws less than or
equal to the observed statistic.
}

\subsection{Examples}{
\preformatted{
## Example 1: Asymptotic-only display (no bootstrap)
stats <- c(Gt = -2.1, Ga = -9.5, Pt = -1.8, Pa = -6.2)

res1 <- DisplayWesterlund(
  stats      = stats,
  bootstats = NULL,
  nobs       = 18,
  nox        = 2,
  constant   = TRUE,
  trend      = FALSE,
  auto       = 0,
  westerlund = FALSE
)

## Example 2: With bootstrap distribution (robust p-values)
set.seed(123)
bootstats <- cbind(
  rnorm(399, mean = -1.8, sd = 0.9),
  rnorm(399, mean = -8.0, sd = 3.0),
  rnorm(399, mean = -1.5, sd = 1.0),
  rnorm(399, mean = -5.0, sd = 3.5)
)

res2 <- DisplayWesterlund(
  stats      = stats,
  bootstats = bootstats,
  nobs       = 18,
  nox        = 2,
  constant   = TRUE,
  trend      = FALSE,
  auto       = 1,
  realmeanlag  = 1.35,
  realmeanlead = 0.40,
  westerlund = FALSE
)
}
}
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{westerlund_test}},
\code{\link{WesterlundPlain}},
\code{\link{WesterlundBootstrap}}
}
