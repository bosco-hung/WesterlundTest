\name{summary.westerlund_test}
\alias{summary.westerlund_test}
\alias{print.summary.westerlund_test}
\title{Summary Method for Westerlund ECM Panel Cointegration Tests}
\description{
Provides a detailed summary of the Westerlund (2007) panel cointegration test
results. This method extends the basic \code{print} output by including
standardized Z-scores, asymptotic p-values, bootstrap p-values (if applicable),
and the Mean Group (MG) parameter estimates.
}
\usage{
\method{summary}{westerlund_test}(object, ...)

\method{print}{summary.westerlund_test}(x, ...)
}
\arguments{
  \item{object}{An object of class \code{westerlund_test}.}
  \item{x}{An object of class \code{summary.westerlund_test}.}
  \item{...}{Additional arguments passed to or from other methods.}
}
\details{
The \code{summary} method acts as a bridge between the raw test results and
a human-readable report. While the standard \code{print} method focuses on
the raw statistics, \code{summary} organizes the data into tables that
compare the observed statistics against their standardized counterparts and
associated p-values.

Specifically, the output includes a summary of the model settings (such as
the inclusion of constants or trends), the main test statistics table, and the
Mean Group results which highlight the average error correction speed (\eqn{\alpha})
and long-run relationships (\eqn{\beta}). If a bootstrap was performed, the
bootstrap p-values are appended to the main statistics table for direct
comparison with asymptotic values.
}
\value{
The \code{summary} method returns a list of class \code{summary.westerlund_test}
containing:
\itemize{
  \item \code{stats_table}: A data frame combining raw stats, Z-scores, and p-values.
  \item \code{settings}: A list of the parameters used in the test.
  \item \code{mg_results}: A data frame of the Mean Group estimation results.
  \item \code{n_units}: The number of cross-sectional units in the panel.
}
}
\examples{
## Assuming 'res' is the output from westerlund_test()
# res <- westerlund_test(data, ...)
# summary(res)
}
\seealso{
\code{\link{westerlund_test}}, \code{\link{plot.westerlund_test}}
}
