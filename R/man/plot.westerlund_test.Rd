\name{plot.westerlund_test}
\alias{plot.westerlund_test}
\title{Plot Bootstrap Distributions for Westerlund ECM Panel Cointegration Tests}
\description{
Creates a 2-by-2 bootstrap density display for \code{Gt}, \code{Ga}, \code{Pt},
and \code{Pa}, including the observed statistic, a configurable lower-tail
bootstrap critical value, and optional robust p-value annotations.
}
\usage{
\method{plot}{westerlund_test}(
  x,
  title = "Westerlund Test: Bootstrap Distributions",
  conf_level = 0.05,
  save_path = NULL,
  dpi = 300,
  figsize = c(12, 10),
  colors = list(obs = "#D55E00", crit = "#0072B2", fill = "grey80", density = "grey30"),
  lwd = list(obs = 1, crit = 0.8, density = 0.5),
  alpha = 0.5,
  show_grid = TRUE,
  show_robust_p = TRUE,
  ...
)
}
\arguments{
\item{x}{A \code{westerlund_test} object containing bootstrap distributions.}
\item{title}{Main plot title.}
\item{conf_level}{Lower-tail bootstrap critical probability in (0,1).}
\item{save_path}{Optional output path. If non-\code{NULL}, the plot is saved with \code{ggsave}.}
\item{dpi}{Saved image resolution.}
\item{figsize}{Numeric width and height in inches.}
\item{colors}{Named colors for observed line, critical line, density fill, and density outline.}
\item{lwd}{Named line widths for observed, critical, and density lines.}
\item{alpha}{Density fill transparency.}
\item{show_grid}{Logical; display major grid lines.}
\item{show_robust_p}{Logical; annotate bootstrap p-values in the facets.}
\item{...}{Reserved for S3 compatibility.}
}
\value{A \code{ggplot} object. If \code{save_path} is supplied, the same plot is also saved to disk.}
\details{
The observed statistic is shown as a solid line and the empirical
\code{conf_level} quantile as a dashed line. Since the Westerlund tests are
lower-tail tests, rejection occurs when the observed statistic lies to the left
of the critical value.
}
\seealso{\code{\link{westerlund_test}}}
