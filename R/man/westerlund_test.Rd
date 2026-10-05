\name{westerlund_test}
\alias{westerlund_test}
\title{Westerlund ECM-Based Panel Cointegration Test}
\description{
Computes the Westerlund (2007) ECM-based panel cointegration statistics
\eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a}, with optional automatic lag/lead
selection, bootstrap inference, unit-level ECM output, and mean-group summaries.
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
  lags = 1,
  leads = NULL,
  westerlund = FALSE,
  aic = TRUE,
  bootstrap = -1,
  indiv.ecm = FALSE,
  lrwindow = 2,
  seed = NULL,
  verbose = FALSE
)
}
\arguments{
\item{data}{A \code{data.frame} containing the panel data.}
\item{yvar}{String naming the dependent variable.}
\item{xvars}{Character vector naming the long-run regressors. At most six are supported; with \code{westerlund=TRUE}, at most one is allowed.}
\item{idvar}{String naming the panel identifier.}
\item{timevar}{String naming the time variable.}
\item{constant}{Logical. Include a constant.}
\item{trend}{Logical. Include a trend. Requires \code{constant=TRUE}.}
\item{lags}{Integer or length-2 integer vector giving a fixed lag or selection range. Defaults to 1.}
\item{leads}{Integer or length-2 integer vector, or \code{NULL}. Defaults to 0.}
\item{westerlund}{Logical. Use the Westerlund-specific selection/variance conventions.}
\item{aic}{Logical. Use AIC for automatic selection; if \code{FALSE}, use BIC.}
\item{bootstrap}{Integer number of bootstrap replications; values \eqn{<=0} disable bootstrapping.}
\item{indiv.ecm}{Logical. Store individual ECM coefficient tables.}
\item{lrwindow}{Non-negative integer Bartlett long-run variance window.}
\item{seed}{Optional integer bootstrap seed. If \code{NULL}, the current R RNG state is used without resetting it.}
\item{verbose}{Logical. Print detailed output.}
}
\details{
The function permits unbalanced panels, but the usable time series within each
unit must be continuous after removing missing model observations. Automatic
lag/lead selection is performed unit by unit when ranges are supplied.

Bootstrap p-values use the finite-sample correction \eqn{(r+1)/(B+1)}, where
\eqn{r} is the number of finite bootstrap statistics less than or equal to the
observed statistic. This convention is intentionally retained in both the R and
Python implementations.
}
\value{
An object of class \code{westerlund_test}, a list containing:
\itemize{
\item \code{test_stats}: raw \code{Gt}, \code{Ga}, \code{Pt}, and \code{Pa}.
\item \code{z_scores}: standardized Z-scores.
\item \code{p_values}: asymptotic left-tail p-values.
\item \code{boot_pvals}: finite-sample-corrected bootstrap p-values when requested.
\item \code{bootstrap_distributions}: bootstrap draws of the four raw statistics.
\item \code{unit_data}: unit-level alpha estimates, standard errors, beta vectors,
selected lags/leads, normalization terms, long-run variance components, and
observation counts. Cross-language aliases are included for compatibility.
\item \code{indiv_data}: detailed per-unit intermediate results.
\item \code{mean_group}: mean-group alpha/beta estimates and their standard errors.
\item \code{settings}: model and selection settings, including the selection method,
bandwidth, bootstrap count, seed, and average panel length.
\item \code{metadata}: compact metadata aligned with the Python implementation.
\item \code{mg_results}: mean-group reporting tables.
\item \code{mg_tables}: alias of the mean-group reporting tables for cross-language consistency.
\item \code{indiv_reg}: stored individual ECM tables when \code{indiv.ecm=TRUE}.
}
}
\examples{
\donttest{
# 1. Generate a small synthetic panel dataset
set.seed(123)
N <- 10
T <- 30

df <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  y = rnorm(N * T),
  x1 = rnorm(N * T)
)

# 2. Run Asymptotic (Plain) Test
res_plain <- westerlund_test(
  data = df,
  yvar = "y",
  xvars = "x1",
  idvar = "id",
  timevar = "time",
  constant = TRUE,
  lags = 1,
  leads = 0,
  verbose = FALSE
)

print(res_plain$test_stats)

# 3. Run Bootstrap Test with automatic lag selection
# Note: bootstrap replications are kept low for example purposes
res_boot <- westerlund_test(
  data = df,
  yvar = "y",
  xvars = "x1",
  idvar = "id",
  timevar = "time",
  constant = TRUE,
  lags = c(0, 1),
  bootstrap = 50,
  verbose = TRUE
)

summary(res_boot)
}
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{
\code{\link{WesterlundPlain}}, \code{\link{WesterlundBootstrap}},
\code{\link{DisplayWesterlund}}, \code{\link{plot.westerlund_test}}
}
