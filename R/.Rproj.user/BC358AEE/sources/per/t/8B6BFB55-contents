\name{westerlund_test_mg}
\alias{westerlund_test_mg}
\title{Print Mean-Group ECM Output and Long-Run Relationship Tables}
\usage{
westerlund_test_mg(b, V, b2, V2, auto)
}
\arguments{
\item{b}{A numeric vector of coefficients for the long-run relationship and short-run adjustment terms (as reported in the second table).}

\item{V}{A numeric variance-covariance matrix corresponding to \code{b}.}

\item{b2}{A numeric vector of coefficients for the mean-group error-correction model (as reported in the first table).}

\item{V2}{A numeric variance-covariance matrix corresponding to \code{b2}.}

\item{auto}{Logical/integer. If non-zero, prints a note indicating that short-run coefficients (apart from the error-correction term) are omitted because selected lag/lead lengths may differ across units.}
}
\value{
This function is called for its side effect (printing).
It invisibly returns \code{NULL}.
}
\description{
Internal reporting helper that prints two Stata-style coefficient tables:
(1) the mean-group error-correction model (short-run ECM output) and
(2) the estimated long-run relationship with short-run adjustment.
This function delegates all table computation and formatting to
\code{\link{westerlund_test_reg}}.
}
\details{
\strong{What it prints.}
The function prints two sections:

\enumerate{
  \item \strong{Mean-group error-correction model:} This corresponds to reporting the mean-group ECM coefficients (often denoted \code{b2} with variance \code{V2}). If \code{auto} is non-zero, a note is printed to remind that short-run coefficients may be omitted due to unit-specific lag/lead selection.
  \item \strong{Estimated long-run relationship and short run adjustment:} This corresponds to reporting the long-run relationship and the adjustment dynamics (\code{b} with variance \code{V}).
}

In each section, \code{\link{westerlund_test_reg}} is called to compute standard
errors, z-statistics, two-sided p-values, and 95\% confidence intervals, and to
print the formatted coefficient table using \code{stats::printCoefmat()}.

\strong{Intended use.}
This is an internal display helper used in Stata-aligned reporting flows. It is
not required for computing the Westerlund test statistics themselves.
}
\section{Vignette}{
This section shows the expected calling pattern for \code{westerlund_test_mg()}
and explains the role of the \code{auto} flag.

\subsection{Two-table reporting}{
In mean-group implementations, it is common to show:
\itemize{
  \item a short-run ECM table (possibly with omitted coefficients when lag/lead
  lengths differ by unit),
  \item a long-run relationship table with the short-run adjustment term.
}
\code{westerlund_test_mg()} prints these tables sequentially.
}

\subsection{Examples}{
\preformatted{
## Long-run relationship table inputs
b <- c(alpha = -0.20, beta = 1.05)
V <- diag(c(0.03^2, 0.10^2))
rownames(V) <- colnames(V) <- names(b)

## Mean-group ECM table inputs
b2 <- c(ec = -0.18)
V2 <- matrix(0.04^2, nrow = 1, ncol = 1)
rownames(V2) <- colnames(V2) <- names(b2)

## Print both sections (auto = TRUE prints the omission note)
westerlund_test_mg(b, V, b2, V2, auto = TRUE)
}
}
}
\seealso{
\code{\link{westerlund_test_reg}},
\code{\link{westerlund_test}}
}
