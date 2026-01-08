\name{plot_westerlund_bootstrap}
\alias{plot_westerlund_bootstrap}
\title{Plot Bootstrap Distributions for Westerlund ECM Panel Cointegration Tests}
\description{
Plots kernel density estimates of the bootstrap distributions for the four
Westerlund (2007) ECM-based panel cointegration statistics (\eqn{G_t}, \eqn{G_a},
\eqn{P_t}, \eqn{P_a}) in a 2x2 panel. Each subplot shows the bootstrap density,
the observed test statistic, and a left-tail critical value at the specified
confidence level.
}
\usage{
plot_westerlund_bootstrap(
  results,
  title = "Westerlund Test: Bootstrap Distributions",
  conf_level = 0.05,
  colors = list(
    obs = "#D55E00",
    crit = "#0072B2",
    fill = rgb(0.8, 0.8, 0.8, 0.5),
    density = "grey"
  ),
  lwd = list(obs = 3, crit = 2, density = 2),
  show_grid = TRUE,
  ...
)
}
\arguments{
  \item{results}{A list-like object containing bootstrap output. It must include:
  \code{results$bootstrap_distributions} (a numeric matrix with 4 columns for \code{Gt, Ga, Pt, Pa})
  and \code{results$test_stats} (a named list or vector containing the observed statistics).}

  \item{title}{Character. Overall figure title added to the outer margin.}

  \item{conf_level}{Numeric in (0,1). Left-tail quantile used as the critical value. For example,
  \code{conf_level = 0.05} plots the 5\% left-tail critical value of each bootstrap distribution.}

  \item{colors}{A named list of plotting colors. Expected elements are:
  \code{obs} (observed-stat line), \code{crit} (critical-value line),
  \code{fill} (density fill), and \code{density} (density outline).}

  \item{lwd}{A named list of line widths. Expected elements are:
  \code{obs}, \code{crit}, and \code{density}.}

  \item{show_grid}{Logical. If \code{TRUE}, adds a background grid to each subplot.}

  \item{\dots}{Additional arguments passed to \code{\link[stats]{plot.density}} (e.g., \code{xlim}, \code{ylim}).}
}
\details{
\strong{Expected structure of \code{results}.}
This function assumes a results object with two key components:
\itemize{
  \item \code{bootstrap_distributions}: a numeric matrix with 4 columns corresponding to
  \code{Gt}, \code{Ga}, \code{Pt}, and \code{Pa} (in that order),
  \item \code{test_stats}: the observed statistics, accessible as \code{results$test_stats[["Gt"]]}
  etc. (a named vector or list is acceptable).
}
If \code{results$bootstrap_distributions} is missing, the function stops with an error.

\strong{What is plotted.}
For each statistic:
\enumerate{
  \item Non-finite bootstrap draws are removed.
  \item A kernel density estimate is computed via \code{\link[stats]{density}}.
  \item The observed statistic is drawn as a solid vertical line.
  \item The left-tail critical value \code{quantile(data, conf_level)} is drawn as a dashed vertical line.
  \item A legend reports the observed statistic and the critical value.
}
The four plots are arranged in a 2x2 grid using \code{\link[graphics]{par}} with
\code{mfrow = c(2,2)}. The original graphics parameters are restored on exit.

\strong{Interpretation.}
The dashed line is a left-tail critical value from the bootstrap distribution.
If the observed statistic lies to the left of this line (smaller than the
critical value), it indicates rejection of the null at approximately
\code{conf_level} under a left-tail decision rule.

\strong{Dependencies.}
The legend label uses \code{scales::percent()} for formatting. Ensure the
\pkg{scales} package is available if you keep this call in your function.
}
\value{
This function is called for its side effect (plotting).
It invisibly returns \code{NULL}.
}
\section{Visualization Notes}{
This section explains how \code{plot_westerlund_bootstrap()} visualizes bootstrap
evidence for the Westerlund ECM test statistics.

\subsection{When is this plot useful?}{
Bootstrap p-values are often preferred for Westerlund-type tests in finite
samples or under cross-sectional dependence. Visualizing the bootstrap
distribution alongside the observed statistic helps diagnose rejection regions
and potential outliers.
}

\subsection{Required inputs}{
You need a results object that contains a \code{bootstrap_distributions} matrix
with four columns and a \code{test_stats} container holding the observed
statistics. These are typically produced by the main \code{westerlund_test} function.
}
}
\examples{
## Example: minimal mock results object
set.seed(123)
mock <- list(
  bootstrap_distributions = cbind(
    Gt = rnorm(399, mean = -2, sd = 1),
    Ga = rnorm(399, mean = -8, sd = 3),
    Pt = rnorm(399, mean = -2, sd = 1),
    Pa = rnorm(399, mean = -6, sd = 3)
  ),
  test_stats = list(Gt = -3.1, Ga = -10.2, Pt = -2.7, Pa = -7.8)
)

## Basic plot (5\% left-tail critical values)
plot_westerlund_bootstrap(mock)

## Customize critical level and x-limits passed via ...
plot_westerlund_bootstrap(mock, conf_level = 0.10, xlim = c(-6, 2))
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{westerlund_test}},
\code{\link{WesterlundBootstrap}},
\code{\link{DisplayWesterlund}}
}
