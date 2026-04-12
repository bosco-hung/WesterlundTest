\name{plot.westerlund_test}
\alias{plot.westerlund_test}
\title{Plot Bootstrap Distributions for Westerlund ECM Panel Cointegration Tests}
\description{
This function provides an S3 method for visualizing objects of class \code{westerlund_test}.
It creates a faceted \code{ggplot2} visualization of the bootstrap distributions for
the four Westerlund (2007) panel cointegration statistics: \eqn{G_t}, \eqn{G_a},
\eqn{P_t}, and \eqn{P_a}. The plot displays the kernel density of bootstrap
replications, the observed test statistic as a solid line, and the left-tail
bootstrap critical value as a dashed line.
}
\usage{
\method{plot}{westerlund_test}(
  x,
  title = "Westerlund Test: Bootstrap Distributions",
  conf_level = 0.05,
  colors = list(
    obs = "#D55E00",
    crit = "#0072B2",
    fill = "grey80",
    density = "grey30"
  ),
  lwd = list(obs = 1, crit = 0.8, density = 0.5),
  show_grid = TRUE,
  ...
)
}
\arguments{
  \item{x}{An object of class \code{westerlund_test}. This list-like object must
  include \code{x$bootstrap_distributions} (a numeric matrix with 4 columns for
  \code{Gt, Ga, Pt, Pa}) and \code{x$test_stats} (a named list or vector
  containing the observed statistics).}

  \item{title}{Character. The main title of the plot.}

  \item{conf_level}{Numeric in (0,1). Left-tail quantile used as the critical value.
  For example, \code{conf_level = 0.05} identifies the 5\% left-tail critical value.}

  \item{colors}{A named list of colors for plot elements. Expected elements include
  \code{obs} for the observed line, \code{crit} for the critical value line,
  \code{fill} for the density area fill, and \code{density} for the density curve outline.}

  \item{lwd}{A named list of line widths or sizes for the \code{obs}, \code{crit},
  and \code{density} elements.}

  \item{show_grid}{Logical. If \code{TRUE}, displays major panel grids.}

  \item{...}{Additional arguments passed to or from other methods.}
}
\details{
The function converts the bootstrap distribution matrix into a long-format data
frame to utilize the faceting capabilities of the \pkg{ggplot2} package. This
implementation uses \code{ggplot2::geom_density()} for empirical distribution
estimation and arranges the four subplots into a 2x2 grid. By using
\code{facet_wrap(~Statistic, scales = "free")}, the plots maintain independent
x-axis scales to better visualize the specific range of each statistic.
Text annotations are automatically added to each facet to display the exact
numerical values of the Observed Stat and the Critical Value for quick reference.

The dashed line represents the bootstrap critical value for a left-tailed test.
If the observed statistic, represented by the solid line, is located to the left
of the dashed line, it indicates a rejection of the null hypothesis of no
cointegration at the specified \code{conf_level} significance level. This
visual approach allows for an immediate assessment of the test results relative
to the empirically derived distribution.

The function requires the \pkg{ggplot2}, \pkg{tidyr}, and \pkg{scales} packages
to be installed. Because it is an S3 method, it is typically invoked by calling
\code{plot()} directly on the result object of a Westerlund test.
}
\value{
Returns a \code{ggplot} object. This allows the user to apply further
customizations, such as changing themes with \code{+ theme_bw()} or adding
additional layers and labels.
}
\section{Inference and Visualization}{
\subsection{Why use bootstrap distributions?}{
Westerlund's pooled statistics can be sensitive to nuisance parameters and
cross-sectional dependence in finite samples. Visualizing the bootstrap
density provides a more robust reference than asymptotic normal curves,
offering a clearer picture of the distribution under the null hypothesis.
}

\subsection{Customization}{
Because the function returns a \code{ggplot} object, users can modify the
output easily. For example, the theme can be changed by adding a theme
layer, such as \code{plot(res) + ggplot2::theme_dark()}.
}
}
\examples{
## Example: generating a mock westerlund_test object
set.seed(123)
mock_res <- list(
  bootstrap_distributions = cbind(
    Gt = rnorm(399, mean = -2, sd = 1),
    Ga = rnorm(399, mean = -8, sd = 3),
    Pt = rnorm(399, mean = -2, sd = 1),
    Pa = rnorm(399, mean = -6, sd = 3)
  ),
  test_stats = list(Gt = -3.1, Ga = -10.2, Pt = -2.7, Pa = -7.8)
)
class(mock_res) <- "westerlund_test"

## Use the S3 plot method
p <- plot(mock_res)

## Display plot
if (interactive()) print(p)
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
