\name{summary.westerlund_test}
\alias{summary.westerlund_test}
\alias{print.summary.westerlund_test}
\title{Summary Method for Westerlund ECM Panel Cointegration Tests}
\description{
Returns and prints a compact report containing the raw statistics, standardized
Z-scores, asymptotic p-values, optional bootstrap p-values, settings, and
mean-group estimates.
}
\usage{
\method{summary}{westerlund_test}(object, ...)

\method{print}{summary.westerlund_test}(x, ...)
}
\arguments{
\item{object}{A \code{westerlund_test} object.}
\item{x}{A \code{summary.westerlund_test} object.}
\item{...}{Additional or unused arguments.}
}
\value{
\code{summary.westerlund_test()} returns a list of class
\code{summary.westerlund_test} containing \code{stats_table}, \code{settings},
\code{mean_group}, \code{mg_results}, \code{mg_tables}, and \code{n_units}. The print method returns
its input invisibly.
}
\seealso{\code{\link{westerlund_test}}, \code{\link{plot.westerlund_test}}}
