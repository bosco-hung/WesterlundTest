\name{print.westerlund_test}
\alias{print.westerlund_test}
\title{Print Method for Westerlund ECM Panel Cointegration Tests}
\description{
This function provides a clean and concise console output for objects of class
\code{westerlund_test}. Instead of displaying the raw list structure, it
summarizes the key results of the Westerlund (2007) panel cointegration test,
including the four observed test statistics and information regarding the
bootstrap replications if they were performed.
}
\usage{
\method{print}{westerlund_test}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{westerlund_test}. This is typically the output
  returned by a call to \code{westerlund_test} or a manually constructed list
  assigned that class.}

  \item{...}{Further arguments passed to or from other methods. Currently unused
  but maintained for S3 generic compatibility.}
}
\details{
The print method is designed to enhance readability by presenting the test
results in a formatted table-like structure. It first identifies the test type
and then displays the observed values for the group-mean statistics (\eqn{G_t}
and \eqn{G_a}) and the pooled panel statistics (\eqn{P_t} and \eqn{P_a}). This
summary provides an immediate overview of the test's findings without
requiring the user to navigate the internal list structure of the object.

If the test was conducted with a bootstrap distribution, the method also
reports the total number of successful replications. This is particularly
useful for verifying that the requested number of bootstrap draws was
attained after accounting for any non-finite results. Finally, the output
includes a brief instructional note suggesting the use of the \code{plot()}
method for visual inference, which is the recommended way to interpret the
bootstrap results relative to the observed statistics.
}
\value{
The function returns the input object \code{x} invisibly (via \code{invisible(x)}).
Its primary purpose is the side effect of printing a formatted summary to the
R console.
}
\examples{
## Example: generating a mock westerlund_test object
set.seed(123)
mock_res <- list(
  test_stats = list(Gt = -3.14, Ga = -10.22, Pt = -2.71, Pa = -7.89),
  bootstrap_distributions = matrix(rnorm(400), ncol = 4)
)
class(mock_res) <- "westerlund_test"

## Calling the object directly invokes print.westerlund_test
mock_res

## Or call it explicitly
print(mock_res)
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{westerlund_test}},
\code{\link{plot.westerlund_test}}
}
