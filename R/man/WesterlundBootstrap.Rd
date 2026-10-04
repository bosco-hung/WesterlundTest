\name{WesterlundBootstrap}
\alias{WesterlundBootstrap}
\title{Bootstrap Routine for Westerlund ECM Panel Cointegration Tests}
\description{
Generates the null bootstrap distribution of the Westerlund ECM panel
cointegration statistics using Stata-aligned restricted regressions, time-cluster
resampling, and bootstrap reconstruction of the panel.
}
\usage{
WesterlundBootstrap(
  data,
  touse,
  idvar,
  timevar,
  yvar,
  xvars,
  constant = FALSE,
  trend = FALSE,
  lags = 1,
  leads = NULL,
  westerlund = FALSE,
  aic = TRUE,
  bootstrap = 100,
  indiv.ecm = FALSE,
  lrwindow = 2,
  seed = NULL,
  verbose = FALSE
)
}
\arguments{
\item{data}{Original panel \code{data.frame}.}
\item{touse}{Logical vector identifying usable rows.}
\item{idvar}{Panel identifier column name.}
\item{timevar}{Time column name.}
\item{yvar}{Dependent variable name.}
\item{xvars}{Regressor names.}
\item{constant}{Include a constant in the restricted bootstrap regression.}
\item{trend}{Include a trend.}
\item{lags}{Fixed lag or length-2 lag range. Defaults to 1.}
\item{leads}{Fixed lead or length-2 lead range; \code{NULL} means 0.}
\item{westerlund}{Use the Westerlund-specific restricted-model information criterion.}
\item{aic}{Use AIC rather than BIC outside Westerlund-specific mode.}
\item{bootstrap}{Number of bootstrap replications.}
\item{indiv.ecm}{Retained for interface compatibility.}
\item{lrwindow}{Bartlett long-run variance window passed to re-estimation.}
\item{seed}{Optional integer seed. If \code{NULL}, the current R RNG state is used.}
\item{verbose}{Print progress information.}
}
\details{
For each unit the routine estimates the restricted short-run model under the
null of no error correction. Residuals are demeaned within unit and
\eqn{\Delta x} is centered over all usable observations. The bootstrap then
reproduces the Stata-style \code{expandcl 2}, actual-time \code{cluster(t)}
sampling, and \code{newttt}/\code{tussent}/\code{newtt} temporal alignment while
retaining unit-specific panel lengths.

When ranges are supplied, lag/lead selection is performed in the restricted
model. In \code{westerlund=TRUE} mode the criterion is
\deqn{\log(RSS/(T_i-p-q-1)) + 2(p+q+c+r+1)/(T_i-p_{max}-q_{max}).}
For literal compatibility with the 2010 \code{xtwest} source, automatic lead
reconstruction uses the final processed unit's selected \code{currlead} as the
global lead-loop upper bound.

The routine sets the RNG kind to Mersenne-Twister. Supplying \code{seed}
restarts that generator reproducibly; with \code{seed=NULL}, the existing RNG
state is used.
}
\value{
A list with \code{BOOTSTATS}, a \code{bootstrap} by 4 numeric matrix containing
\code{Gt}, \code{Ga}, \code{Pt}, and \code{Pa} draws.
}
\references{
Westerlund, J. (2007). Testing for error correction in panel data.
\emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
}
\seealso{\code{\link{westerlund_test}}, \code{\link{WesterlundPlain}}}
