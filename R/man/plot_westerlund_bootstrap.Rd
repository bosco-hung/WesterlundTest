\name{plot_westerlund_bootstrap}
\alias{plot_westerlund_bootstrap}
\title{Plot Bootstrap Distributions for Westerlund ECM Panel Cointegration Tests}
\description{
Creates a faceted \code{ggplot2} visualization of the bootstrap distributions for
the four Westerlund (2007) panel cointegration statistics (\eqn{G_t}, \eqn{G_a},
\eqn{P_t}, \eqn{P_a}). The plot displays the kernel density of bootstrap
replications, the observed test statistic as a solid line, and the left-tail
bootstrap critical value as a dashed line.
}
\usage{
plot_westerlund_bootstrap(
  results,
  title = "Westerlund Test: Bootstrap Distributions",
  conf_level = 0.05,
  colors = list(
    obs = "#D55E00",
    crit = "#0072B2",
    fill = "grey80",
    density = "grey30"
  ),
  lwd = list(obs = 1, crit = 0.8, density = 0.5),
  show_grid = TRUE
)
}
\arguments{
  \item{results}{A list-like object containing bootstrap output. It must include:
  \code{results$bootstrap_distributions} (a numeric matrix with 4 columns for \code{Gt, Ga, Pt, Pa})
  and \code{results$test_stats} (a named list or vector containing the observed statistics).}

  \item{title}{Character. The main title of the plot.}

  \item{conf_level}{Numeric in (0,1). Left-tail quantile used as the critical value.
  For example, \code{conf_level = 0.05} identifies the 5\% left-tail critical value.}

  \item{colors}{A named list of colors for plot elements. Expected elements:
  \code{obs} (observed line), \code{crit} (critical value line),
  \code{fill} (density area fill), and \code{density} (density curve outline).}

  \item{lwd}{A named list of line widths (or sizes) for \code{obs}, \code{crit}, and \code{density}.}

  \item{show_grid}{Logical. If \code{TRUE}, displays major panel grids.}

}
\details{
\strong{Implementation Details.}
The function converts the bootstrap distribution matrix into a long-format
data frame to utilize \code{ggplot2::facet_wrap()}.



Unlike the base R version, this
implementation:
\itemize{
  \item Uses \code{ggplot2::geom_density()} for empirical distribution estimation.
  \item Arranges the four subplots into a 2x2 grid with independent x-axis scales
  using \code{facet_wrap(~Statistic, scales = "free")}.
  \item Adds text annotations to each facet displaying the exact numerical
  values of the Observed Stat and the Critical Value.
}

\strong{Interpretation.}
The dashed line represents the bootstrap critical value for a left-tailed test.
If the observed statistic (solid line) is to the left of the dashed line, it
indicates a rejection of the null hypothesis (\eqn{H_0}: no cointegration)
at the \code{conf_level} significance level.

\strong{Dependencies.}
The function requires the \pkg{ggplot2}, \pkg{tidyr}, and \pkg{scales} packages.
}
\value{
Returns a \code{ggplot} object. This allows the user to apply further
customizations, such as changing themes (\code{+ theme_bw()}) or adding
additional layers.
}
\section{Inference and Visualization}{
\subsection{Why use bootstrap distributions?}{
Westerlund's pooled statistics can be sensitive to nuisance parameters and
cross-sectional dependence in finite samples. Visualizing the bootstrap
density provides a more robust reference than asymptotic normal curves.
}

\subsection{Customization}{
Because the function returns a \code{ggplot} object, users can modify the
output easily. For example, to change the theme:
\code{plot_westerlund_bootstrap(res) + ggplot2::theme_dark()}.
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

## Generate the faceted ggplot
p <- plot_westerlund_bootstrap(mock)

## Display plot
if (interactive()) print(p)

## Add ggplot2 layers
# p + ggplot2::theme_classic()
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
