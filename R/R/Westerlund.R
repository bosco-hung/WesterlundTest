# ==============================================================================
# Main Function: westerlund_test
# ==============================================================================
#' Westerlund ECM-Based Panel Cointegration Test
#'
#' Implements the error-correction-based panel cointegration tests proposed by
#' Westerlund (2007). The function computes both mean-group (\eqn{G_t}, \eqn{G_a})
#' and pooled (\eqn{P_t}, \eqn{P_a}) statistics based on unit-specific ECMs, with
#' Stata-like input validation, time-continuity checks, and optional bootstrap
#' p-values for the raw statistics.
#'
#' @param data A \code{data.frame} containing the panel data. Panels may be unbalanced,
#'   but each unit's time index must be continuous (no holes) after removing missing
#'   observations on model variables.
#' @param yvar String. Name of the dependent variable.
#' @param xvars Character vector. Names of regressors entering the long-run relationship.
#'   A maximum of 6 regressors is allowed. If \code{westerlund=TRUE}, at most one
#'   regressor may be included.
#' @param idvar String. Column identifying cross-sectional units. Missing values
#'   are not allowed.
#' @param timevar String. Column identifying time. Within each unit, the time index
#'   must be strictly increasing and continuous (differences of 1).
#' @param constant Logical. Include an intercept in the cointegrating relationship.
#' @param trend Logical. Include a linear trend in the cointegrating relationship.
#'   Setting \code{trend=TRUE} requires \code{constant=TRUE}.
#' @param lags Integer or length-2 integer vector. Fixed lag order or range
#'   \code{c(min,max)} for the short-run lag length \eqn{p}. Defaults to 1.
#' @param leads Integer or length-2 integer vector, or \code{NULL}. Fixed lead order
#'   or range \code{c(min,max)} for the short-run lead length \eqn{q}. If \code{NULL},
#'   defaults to 0.
#' @param westerlund Logical. If \code{TRUE}, enforces additional restrictions: at
#'   least a constant must be included (constant and/or trend) and at most one
#'   regressor may be specified.
#' @param aic Logical. If \code{TRUE}, uses AIC for lag/lead selection when ranges.
#'   If \code{FALSE}, uses BIC.
#' @param bootstrap Integer. If \code{bootstrap > 0}, runs an internal residual-based
#'   bootstrap with \code{bootstrap} replications. The bootstrap resamples actual
#'   time clusters across units to preserve contemporaneous cross-sectional dependence,
#'   retains unit-specific panel lengths, and returns bootstrap p-values for the raw
#'   statistics. If \code{bootstrap <= 0} (default \code{-1}), no bootstrap is performed.
#' @param indiv.ecm Logical. If \code{TRUE}, gets output of individual ECM regressions.
#' @param lrwindow Integer. Window parameter used for long-run variance estimation
#'   in internal routines. Default is 2.
#' @param seed Optional integer random seed used for bootstrap reproducibility. If
#'   \code{NULL}, the current R random-number state is used without resetting it.
#' @param verbose Logical. If \code{TRUE}, prints additional information.
#'
#' @details
#' The test evaluates the null hypothesis of no error correction for all
#' cross-sectional units. Rejection suggests that at least some units exhibit
#' error-correcting behavior.
#'
#' \strong{Statistics:}
#' \itemize{
#'   \item \eqn{G_t}: mean of individual t-ratios of the error-correction coefficient.
#'   \item \eqn{G_a}: mean of individually scaled error-correction coefficients.
#'   \item \eqn{P_t}: pooled t-type statistic based on a common error-correction coefficient.
#'   \item \eqn{P_a}: pooled statistic based on the scaled common coefficient.
#' }
#'
#' \strong{Input validation and data checks:}
#' The function enforces several guardrails:
#' \itemize{
#'   \item \strong{Regressor Limits:} A maximum of 6 regressors in \code{xvars}.
#'   \item \strong{Deterministic Logic:} If \code{trend=TRUE}, then \code{constant=TRUE} is required.
#'   \item \strong{Westerlund Restrictions:} If \code{westerlund=TRUE}, at least a constant
#'     must be included and \code{xvars} must contain at most one regressor.
#'   \item \strong{Time Continuity:} The panel must be continuous within each unit. If
#'     any unit contains time-series holes (time difference > 1), the function stops.
#'   \item \strong{Observation Sufficiency:} A minimum number of observations per unit
#'     is required based on: \code{maxlag + maxlead + 1 + (constant + trend) + 1 + nox +
#'     maxlag + nox * (maxlag + maxlead + 1) + 1}.
#' }
#'
#' \strong{Inference:}
#' If \code{bootstrap <= 0}, standardized Z-scores and asymptotic (left-tail)
#' p-values are returned. If \code{bootstrap > 0}, the function additionally
#' constructs a residual-based null distribution using Stata-aligned time-cluster
#' resampling and reports bootstrap p-values for the raw statistics. The bootstrap
#' p-values in this R implementation use the finite-sample correction
#' \eqn{(r+1)/(B+1)}.
#'
#' @return An object of class \code{westerlund_test}. This is a list containing
#' the following components:
#' \itemize{
#'   \item \code{test_stats}: Named raw statistics \code{Gt}, \code{Ga}, \code{Pt}, and \code{Pa}.
#'   \item \code{z_scores}: Named standardized Z-scores for the four statistics.
#'   \item \code{p_values}: Named asymptotic left-tail p-values for the four statistics.
#'   \item \code{boot_pvals}: A named list containing bootstrap p-values under the
#'     names \code{Gt}, \code{Ga}, \code{Pt}, and \code{Pa}, calculated only if
#'     \code{bootstrap > 0}.
#'   \item \code{bootstrap_distributions}: A matrix of dimension \code{[bootstrap x 4]}
#'     containing the bootstrap distributions of \code{Gt}, \code{Ga}, \code{Pt}, and
#'     \code{Pa} across replications.
#'   \item \code{unit_data}: A \code{data.frame} containing unit-specific alpha estimates,
#'     standard errors, beta vectors, selected lags/leads, normalization terms,
#'     long-run variance components, and observation counts. 
#'   \item \code{indiv_data}: A list of length equal to the number of cross-sectional units storing unit-specific results.
#'   \item \code{mean_group}: A list containing:
#'     \itemize{
#'       \item \code{mg_alpha}: The Mean Group estimate of the error correction coefficient.
#'       \item \code{se_mg_alpha}: Standard error of the MG alpha.
#'       \item \code{mg_betas}: Vector of Mean Group estimates for the long-run regressors.
#'       \item \code{se_mg_betas}: Standard errors for the MG betas.
#'     }
#'   \item \code{settings}: A list of internal parameters, lag/lead selections,
#'     selection method, bandwidth, bootstrap count, seed, and average panel length.
#'   \item \code{metadata}: Compact metadata.
#'   \item \code{mg_results}: Mean-group reporting tables.
#'   \item \code{mg_tables}: Alias of the two mean-group reporting tables for
#'     cross-language consistency.
#'   \item \code{indiv_reg}: Stored individual ECM tables when requested.
#' }
#'
#' @references
#' Westerlund, J. (2007). Testing for error correction in panel data.
#' \emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
#'
#' @seealso
#' \code{\link{WesterlundPlain}}, \code{\link{WesterlundBootstrap}}, \code{\link{DisplayWesterlund}}
#'
#' @examples
#' \donttest{
#' # 1. Generate a small synthetic panel dataset
#' set.seed(123)
#' N <- 10; T <- 30
#' df <- data.frame(
#'   id = rep(1:N, each = T),
#'   time = rep(1:T, N),
#'   y = rnorm(N * T),
#'   x1 = rnorm(N * T)
#' )
#'
#' # 2. Run Asymptotic (Plain) Test
#' res_plain <- westerlund_test(
#'   data = df, yvar = "y", xvars = "x1",
#'   idvar = "id", timevar = "time",
#'   constant = TRUE, lags = 1, leads = 0,
#'   verbose = FALSE
#' )
#' print(res_plain$test_stats)
#'
#' # 3. Run Bootstrap Test with automatic lag selection
#' # Note: bootstrap replications are kept low for example purposes
#' res_boot <- westerlund_test(
#'   data = df, yvar = "y", xvars = "x1",
#'   idvar = "id", timevar = "time",
#'   constant = TRUE, lags = c(0, 1),
#'   bootstrap = 50, verbose = TRUE
#' )
#' }
#'
#' @export
westerlund_test <- function(data,
                   yvar,
                   xvars,
                   idvar,
                   timevar,
                   constant = FALSE,
                   trend = FALSE,
                   lags = 1,          # Can be a single integer or a vector c(min, max)
                   leads = NULL,      # Can be a single integer or a vector c(min, max)
                   westerlund = FALSE,
                   aic = TRUE,
                   bootstrap = -1,    # Default -1
                   indiv.ecm = FALSE,
                   lrwindow = 2,
                   seed = NULL,
                   verbose = FALSE) {

  # ----------------------------------------------------------------------------
  # 1. Input Parsing and Validation
  # ----------------------------------------------------------------------------
  nox <- length(xvars)
  if (nox > 6) stop("Error: No more than 6 covariates can be specified. (Exit 459)")

  withconstant <- constant
  withtrend <- trend
  withwesterlund <- westerlund

  if (withwesterlund) {
    if (!withconstant && !withtrend) stop("Error: if westerlund is specified, at least a constant must be added.")
    if (nox > 1) stop("Error: if westerlund is specified at most one x-variable may be included.")
  }

  if (withtrend && !withconstant) stop("Error: If a trend is included, a constant must be included as well.")
  if (length(lrwindow) != 1 || !is.finite(lrwindow) || lrwindow < 0 || lrwindow != as.integer(lrwindow))
    stop("Error: lrwindow must be a non-negative integer.")
  if (length(bootstrap) != 1 || !is.finite(bootstrap) || bootstrap != as.integer(bootstrap))
    stop("Error: bootstrap must be an integer.")
  if (!is.null(seed) && (length(seed) != 1 || !is.finite(seed) || seed != as.integer(seed)))
    stop("Error: seed must be NULL or a single integer.")

  # ----------------------------------------------------------------------------
  # 2. Data Preparation
  # ----------------------------------------------------------------------------
  if (!(idvar %in% names(data)) || !(timevar %in% names(data))) stop("Error: Specify idvar and timevar.")
  data <- data[order(data[[idvar]], data[[timevar]]), ]
  if (any(is.na(data[[idvar]]))) stop("Error: missing values in idvar not allowed.")

  touse <- rep(TRUE, nrow(data))
  touse[is.na(data[[idvar]]) | is.na(data[[yvar]])] <- FALSE
  for (x in xvars) touse[is.na(data[[x]])] <- FALSE

  # ----------------------------------------------------------------------------
  # 3. Time Series Continuity Check
  # ----------------------------------------------------------------------------
  data_check <- data[touse, ]
  ids <- unique(data_check[[idvar]])
  for (i in ids) {
    sub_t <- data_check[[timevar]][data_check[[idvar]] == i]
    if (length(sub_t) > 1 && any(diff(sub_t) > 1)) stop(paste("Error: Holes in time series for ID", i))
  }

  # ----------------------------------------------------------------------------
  # 4. Lags and Leads Parsing
  # ----------------------------------------------------------------------------
  unique_levels <- unique(data[[idvar]][touse])
  nobs <- length(unique_levels)

  minlag <- if(length(lags) > 1) min(lags) else lags[1]
  maxlag <- if(length(lags) > 1) max(lags) else lags[1]

  minlead <- 0; maxlead <- 0
  if (!is.null(leads)) {
    minlead <- if(length(leads) > 1) min(leads) else leads[1]
    maxlead <- if(length(leads) > 1) max(leads) else leads[1]
  }
  auto <- (minlead != maxlead) || (minlag != maxlag)
  if (minlag < 0 || maxlag < 0 || minlead < 0 || maxlead < 0 ||
      any(c(minlag, maxlag, minlead, maxlead) != as.integer(c(minlag, maxlag, minlead, maxlead))))
    stop("Error: lags and leads must be non-negative integers.")

  # ----------------------------------------------------------------------------
  # 5. Observation Sufficiency Check
  # ----------------------------------------------------------------------------
  numbobs_table <- table(data_check[[idvar]])
  minobs <- maxlag + maxlead + 1 + (as.integer(constant) + as.integer(trend)) + 1 + nox + maxlag + nox * (maxlag + maxlead + 1) + 1
  if (min(numbobs_table) < minobs) stop("Error: Insufficient observations for selected lags/leads.")

  # ----------------------------------------------------------------------------
  # 6. Branching Logic
  # ----------------------------------------------------------------------------
  boot_dist <- NULL
  boot_pvals <- list()

  if (bootstrap > 0) {
    if(verbose)cat("\nCalculating Westerlund ECM panel cointegration tests\n\n")
    boot_res <- WesterlundBootstrap(data, touse, idvar, timevar, yvar, xvars, constant, trend, lags, leads, westerlund, aic, bootstrap, indiv.ecm = FALSE, lrwindow = lrwindow, seed = seed, verbose = verbose)
    boot_dist <- boot_res$BOOTSTATS
  }

  plain_res <- WesterlundPlain(data, touse, idvar, timevar, yvar, xvars, constant, trend, lags, leads, lrwindow, westerlund, aic, bootno = TRUE, indiv.ecm = indiv.ecm, verbose = verbose)

  # ----------------------------------------------------------------------------
  # 7. Mean Group (MG) Estimation and Display
  # ----------------------------------------------------------------------------
  # Pass 1.5: Aggregating individual results stored in plain_res$indiv_data

  # A. Error Correction term (Alpha)
  all_alphas <- sapply(plain_res$indiv_data, function(x) x$ai)
  mg_alpha <- mean(all_alphas)
  se_mg_alpha <- sd(all_alphas) / sqrt(length(all_alphas))

  b_ec <- c("ec (alpha)" = mg_alpha)
  V_ec <- matrix(se_mg_alpha^2, 1, 1, dimnames = list("ec (alpha)", "ec (alpha)"))

  # B. Long-Run relationship (Betas)
  # Extract beta vectors and convert to matrix [N x nox]
  beta_matrix <- do.call(rbind, lapply(plain_res$indiv_data, function(x) x$betai))
  mg_betas <- colMeans(beta_matrix)

  # Calculate Variance of the Mean Group Beta (Group Variation / N)
  # Using diag() ensures the matrix is square for the regdisplay checks
  if (nox > 1) {
    se_mg_betas <- apply(beta_matrix, 2, sd) / sqrt(nobs)
    V_lr <- diag(se_mg_betas^2)
  } else {
    se_mg_betas <- sd(beta_matrix) / sqrt(nobs)
    V_lr <- matrix(se_mg_betas^2, 1, 1)
  }
  names(mg_betas) <- xvars
  dimnames(V_lr) <- list(xvars, xvars)

  # Display the MG Results (Stata: eret disp)
  mg_results = westerlund_test_mg(b = mg_betas, V = V_lr, b2 = b_ec, V2 = V_ec, auto = auto, verbose = verbose)

  # ----------------------------------------------------------------------------
  # 8. Test Statistics Display
  # ----------------------------------------------------------------------------
  if (bootstrap > 0) {
    for (j in 1:4) {
      stat_name <- c("Gt", "Ga", "Pt", "Pa")[j]
      obs <- plain_res$stats[[stat_name]]
      boots <- boot_dist[, j]; boots <- boots[is.finite(boots)]
      boot_pvals[[stat_name]] <- (sum(boots <= obs) + 1) / (length(boots) + 1)
    }
  }

  display_res <- DisplayWesterlund(
    stats = plain_res$stats,
    bootstats = boot_dist,
    nobs = nobs,
    nox = nox,
    constant = constant,
    trend = trend,
    westerlund = westerlund,
    aic = aic,
    meanlag = plain_res$settings$meanlag,
    meanlead = plain_res$settings$meanlead,
    realmeanlag = plain_res$settings$realmeanlag,
    realmeanlead = plain_res$settings$realmeanlead,
    auto = plain_res$settings$auto,
    verbose = verbose
  )

  # ----------------------------------------------------------------------------
  # 9. Construct Final Bundle
  # ----------------------------------------------------------------------------
  z_scores <- c(Gt = display_res$gt_z, Ga = display_res$ga_z,
                Pt = display_res$pt_z, Pa = display_res$pa_z)
  p_values <- c(Gt = display_res$gt_pval, Ga = display_res$ga_pval,
                Pt = display_res$pt_pval, Pa = display_res$pa_pval)
  selection_method <- if (westerlund) "Westerlund IC" else if (aic) "AIC" else "BIC"
  avg_obs <- mean(as.numeric(numbobs_table))
  model_type <- if (trend) "trend" else if (constant) "constant" else "none"
  settings <- c(plain_res$settings, list(
    minlag = minlag, maxlag = maxlag, minlead = minlead, maxlead = maxlead,
    selection_method = selection_method, lrwindow = lrwindow,
    bootstrap = bootstrap, seed = seed, T = avg_obs
  ))
  metadata <- list(
    bandwidth = lrwindow, n_groups = nobs, avg_obs = avg_obs,
    model_type = model_type, selection_method = selection_method, seed = seed
  )
  indiv_reg <- lapply(plain_res$indiv_data, function(z) z$reg_coef)
  unit_data <- plain_res$results_df
  unit_data$alpha_i <- unit_data$ai
  unit_data$se_alpha_i <- unit_data$seai
  unit_data$lag <- unit_data$lags
  unit_data$lead <- unit_data$leads
  unit_data$obs <- vapply(as.character(unit_data$id), function(uid) plain_res$indiv_data[[uid]]$ti, numeric(1))
  unit_data$ti <- unit_data$obs
  unit_data$beta_i <- I(lapply(as.character(unit_data$id), function(uid) plain_res$indiv_data[[uid]]$betai))

  results_bundle <- list(
    test_stats = plain_res$stats,
    z_scores = z_scores,
    p_values = p_values,
    boot_pvals = boot_pvals,
    bootstrap_distributions = boot_dist,
    unit_data = unit_data,
    indiv_data = plain_res$indiv_data,
    mean_group = list(
      mg_alpha = mg_alpha,
      se_mg_alpha = se_mg_alpha,
      mg_betas = mg_betas,
      se_mg_betas = se_mg_betas
    ),
    settings = settings,
    metadata = metadata,
    mg_results = mg_results,
    mg_tables = mg_results,
    indiv_reg = indiv_reg
  )

  class(results_bundle) <- "westerlund_test"

  return(results_bundle)

}

# ------------------------------------------------------------------------------
# HELPER 1: Bartlett Kernel Long-Run Variance
# ------------------------------------------------------------------------------
#' Long-Run Variance Estimation with Bartlett Kernel
#'
#' Computes a Bartlett-kernel (Newey-West style) long-run variance estimate for a
#' univariate series. This function follows Stata-like conventions: missing values
#' are removed, autocovariances are scaled by \eqn{1/n} (not \eqn{1/(n-j)}), and
#' optional centering is provided.
#'
#' @param x A numeric vector. Missing values are removed prior to computation
#'   using \code{na.omit}.
#' @param maxlag Non-negative integer. Maximum lag order \eqn{m} used in the
#'   Bartlett kernel. If \code{maxlag = 0}, the function returns the sample
#'   variance estimate \eqn{\gamma_0}.
#' @param nodemean Logical. If \code{FALSE} (default), the series is centered
#'   by subtracting its mean before computing autocovariances. If \code{TRUE},
#'   the function assumes \code{x} is already centered.
#'
#' @details
#' The Bartlett kernel provides a consistent estimator of the long-run variance
#' (the spectral density at frequency zero) by calculating a weighted sum of
#' autocovariances.
#'
#' \strong{Mathematical Formulation:}
#' Let \eqn{x_t} be the cleaned series of length \eqn{n}. If \code{nodemean = FALSE},
#' \eqn{x_t} is transformed to \eqn{x_t - \bar{x}}.
#'
#' The lag-\eqn{j} autocovariance is calculated as:
#' \deqn{\gamma_j = \frac{1}{n}\sum_{t=j+1}^{n} x_t x_{t-j}}
#'
#' The Bartlett kernel weight for lag \eqn{j} given maximum lag \eqn{m} is:
#' \deqn{w_j = 1 - \frac{j}{m+1}}
#'
#' Finally, the long-run variance estimate is:
#' \deqn{\widehat{\Omega} = \gamma_0 + 2\sum_{j=1}^{m} w_j \gamma_j}
#'
#' \strong{Comparison with other implementations:}
#' This function specifically replicates Stata's behavior. Most notably, it uses
#' a fixed \eqn{1/n} denominator for all autocovariances rather than the
#' \eqn{1/(n-j)} degrees-of-freedom correction often found in other R
#' implementations.
#'
#' @return A single numeric value representing the Bartlett-kernel long-run
#'   variance estimate. Returns \code{NA} if the input vector contains only
#'   missing values.
#'
#' @examples
#' ## Basic usage with a white noise series
#' x <- rnorm(200)
#' calc_lrvar_bartlett(x, maxlag = 4)
#'
#' ## Setting maxlag = 0 returns the sample variance (gamma_0)
#' calc_lrvar_bartlett(x, maxlag = 0)
#'
#' ## Example with missing values (NA values are removed before calculation)
#' x_na <- x
#' x_na[c(5, 10, 50)] <- NA
#' calc_lrvar_bartlett(x_na, maxlag = 4)
#'
#' ## Pre-centered data using nodemean = TRUE
#' x_centered <- x - mean(x)
#' calc_lrvar_bartlett(x_centered, maxlag = 4, nodemean = TRUE)
#'
#' @seealso \code{\link{westerlund_test}}
#'
#' @export
calc_lrvar_bartlett <- function(x, maxlag, nodemean = FALSE) {
  # Remove NAs strictly
  x_clean <- na.omit(x)
  n <- length(x_clean)
  if (n == 0) return(NA)

  # Stata 'nodemean': If FALSE, center the variable. If TRUE, assume already centered.
  if (!nodemean) {
    x_clean <- x_clean - mean(x_clean)
  }

  # Variance (Gamma 0)
  gamma_0 <- sum(x_clean^2) / n

  if (maxlag == 0) return(gamma_0)

  # Weighted Sum of Autocovariances
  sum_weighted <- 0
  for (j in 1:maxlag) {
    # Gamma_j: sum(x_t * x_{t-j}) / n
    # Note: Stata divides by N, not N-k
    gamma_j <- sum(x_clean[(j + 1):n] * x_clean[1:(n - j)]) / n

    # Bartlett Kernel Weight: 1 - j / (m + 1)
    weight <- 1 - (j / (maxlag + 1))
    sum_weighted <- sum_weighted + (weight * gamma_j)
  }

  # LRV = Gamma0 + 2 * Sum(Weighted Gamma_j)
  lrvar <- gamma_0 + 2 * sum_weighted
  return(lrvar)
}

#' Strict Time-Series Mapping with Explicit Gaps
#'
#' Retrieves lagged (or led) values from a time-indexed vector using a strict
#' time-based mapping. This function replicates the behavior of Stata’s \code{L.}
#' and \code{F.} operators by respecting gaps in the time index and returning
#' \code{NA} when the target time does not exist.
#'
#' @param vec A numeric (or atomic) vector of observations.
#' @param tvec A vector of time indices corresponding one-to-one with \code{vec}.
#'   Values must uniquely identify time periods within the series.
#' @param lag Integer. The lag order. Positive values correspond to lags
#'   (e.g., \code{lag = 1} for \eqn{t-1}), while negative values correspond to
#'   leads (e.g., \code{lag = -1} for \eqn{t+1}).
#'
#' @details
#' In panel and time-series econometrics, lagging by position can lead to
#' incorrect results if the time series is not perfectly continuous.
#' \code{get_ts_val} implements a strict lookup based on the values in
#' \code{tvec}.
#'
#' For each observation at time \eqn{t}, the function computes the target time
#' \eqn{t^* = t - \text{lag}}. It then retrieves the value associated with
#' \eqn{t^*} from the original series. If \eqn{t^*} is not present in
#' \code{tvec}, the function returns \code{NA} for that observation.
#'
#' \strong{Comparison with Shifting:}
#' Unlike simple position-based shifts (e.g., \code{dplyr::lag}), this function
#' ensures that if a series jumps from \eqn{t=2} to \eqn{t=4}, the lag-1 value
#' for \eqn{t=4} is \code{NA} (because \eqn{t=3} is missing) rather than
#' erroneously returning the value from \eqn{t=2}.
#'
#' @return A vector of the same length as \code{vec}, containing the lagged
#'   (or led) values aligned by the time index. Elements are \code{NA} when
#'   the requested time period does not exist in \code{tvec}.
#'
#' @examples
#' ## Example 1: Regular time series
#' t <- 1:5
#' x <- c(10, 20, 30, 40, 50)
#'
#' # Lag by one period: NA 10 20 30 40
#' get_ts_val(x, t, lag = 1)
#'
#' ## Example 2: Time series with a gap (missing t=3)
#' t_gap <- c(1, 2, 4, 5)
#' x_gap <- c(10, 20, 40, 50)
#'
#' # Lag by one period: note the NA at time 4 because t=3 is missing
#' # Result: NA 10 NA 40
#' get_ts_val(x_gap, t_gap, lag = 1)
#'
#' # Lead by one period: 20 NA 50 NA
#' get_ts_val(x_gap, t_gap, lag = -1)
#'
#' @seealso \code{\link{westerlund_test}}, \code{\link{calc_lrvar_bartlett}}
#'
#' @export
# Strict Time-Series Mapping (Replicates Stata L. and D. operators with gaps)
get_ts_val <- function(vec, tvec, lag) {
  val_map <- setNames(vec, tvec)
  target_time <- tvec - lag
  unname(val_map[as.character(target_time)])
}

# ------------------------------------------------------------------------------
# HELPER 2: Stata-style Time Series Operators (Handles Gaps)
# ------------------------------------------------------------------------------
#' Time-Indexed Lag Extraction with Explicit Gaps
#'
#' Extracts lagged (or led) values from a time-indexed vector using a strict
#' time-based lookup. The function respects gaps in the time index and returns
#' \code{NA} when the requested lagged time does not exist, mirroring the behavior
#' of Stata’s time-series lag/lead operators.
#'
#' @param vec A numeric (or atomic) vector of observations.
#' @param tvec A vector of time indices corresponding one-to-one with \code{vec}.
#'   Each value must uniquely identify a time period within the series.
#' @param k Integer. The lag order. Positive values correspond to lags
#'   (e.g., \code{k = 1} for \eqn{t-1}), while negative values correspond to leads
#'   (e.g., \code{k = -1} for \eqn{t+1}).
#'
#' @details
#' This helper function performs strict time-based lagging rather than
#' position-based shifting. It ensures that variables are aligned according to the
#' actual time index, which is critical for irregular time series or panels with
#' missing periods.
#'
#' \strong{Logic:}
#' For each observation at time \eqn{t}, the function seeks the value associated
#' with time \eqn{t - k}. If that specific time value is not found in \code{tvec},
#' the function returns \code{NA}.
#'
#' \strong{Comparison with Position Shifting:}
#' In standard vector shifting (like \code{x[i-1]}), a gap in time (e.g., jumping
#' from 2005 to 2007) would result in the 2005 value being incorrectly assigned
#' as the "lag" for 2007. \code{get_lag} correctly identifies that 2006 is
#' missing and returns \code{NA} for the 2007 lag-1 observation.
#'
#' @return A vector of the same length as \code{vec}, containing the lagged (or led)
#'   values aligned by the time index. Elements are \code{NA} when the requested
#'   time period does not exist.
#'
#' @examples
#' ## Example 1: Regular time series
#' t <- 1:5
#' x <- c(10, 20, 30, 40, 50)
#'
#' # Lag by one period: NA 10 20 30 40
#' get_lag(x, t, k = 1)
#'
#' # Lead by one period: 20 30 40 50 NA
#' get_lag(x, t, k = -1)
#'
#' ## Example 2: Time series with a gap (t=3 is missing)
#' t_gap <- c(1, 2, 4, 5)
#' x_gap <- c(10, 20, 40, 50)
#'
#' # Lag by one period: note the NA at t=4 because t=3 is missing
#' # Result: NA 10 NA 40
#' get_lag(x_gap, t_gap, k = 1)
#'
#' ## Example 3: Higher-order lags
#' # Result: NA NA NA 20
#' get_lag(x_gap, t_gap, k = 2)
#'
#' @seealso \code{\link{get_ts_val}}, \code{\link{westerlund_test}},
#'   \code{\link{calc_lrvar_bartlett}}
#'
#' @export
get_lag <- function(vec, tvec, k) {
  val_map <- setNames(vec, tvec)
  unname(val_map[as.character(tvec - k)])
}

#' First Difference with Strict Time Indexing
#'
#' Computes the first difference of a time-indexed series using strict time-based
#' lagging. The function respects gaps in the time index and returns \code{NA} when
#' the previous time period does not exist, mirroring the behavior of Stata’s
#' \code{D.} operator.
#'
#' @param vec A numeric (or atomic) vector of observations.
#' @param tvec A vector of time indices corresponding one-to-one with \code{vec}.
#'   Each value must uniquely identify a time period within the series.
#'
#' @details
#' This helper function computes first differences as:
#' \deqn{\Delta x_t = x_t - x_{t-1}}
#' where the lagged value \eqn{x_{t-1}} is obtained using \code{\link{get_lag}}.
#'
#' \strong{Time-Based Logic vs. Position-Based Logic:}
#' Standard R functions like \code{diff()} compute differences based on the
#' position of elements in a vector. This assumes that observations are
#' contiguous in time. If your data has a gap (e.g., observations for 2010 and
#' 2012 but not 2011), \code{diff()} will subtract the 2010 value from the
#' 2012 value.
#'
#' In contrast, \code{get_diff()} identifies that the observation for \eqn{t-1}
#' is missing and returns \code{NA} for that period. This ensures that
#' short-run dynamics are only estimated between adjacent time periods.
#'
#' No interpolation or implicit shifting is performed; missing time periods
#' propagate as \code{NA} values in the resulting vector.
#'
#' @return A vector of the same length as \code{vec}, containing the first
#'   differences aligned by the time index. Elements are \code{NA} when the
#'   previous time period (\eqn{t-1}) does not exist in \code{tvec}.
#'
#' @examples
#' ## Example 1: Regular time series
#' t <- 1:5
#' x <- c(10, 20, 30, 40, 50)
#'
#' # Result: NA 10 10 10 10
#' get_diff(x, t)
#'
#' ## Example 2: Time series with a gap (t=3 is missing)
#' t_gap <- c(1, 2, 4, 5)
#' x_gap <- c(10, 20, 40, 50)
#'
#' # Result: NA 10 NA 10
#' # Note the NA at t=4 because the difference requires t=3.
#' get_diff(x_gap, t_gap)
#'
#' ## Comparison with base R diff()
#' # diff(x_gap) would return: 10 20 10 (losing the time context and vector length)
#'
#' @seealso \code{\link{get_lag}}, \code{\link{get_ts_val}}, \code{\link{westerlund_test}}
#'
#' @export
get_diff <- function(vec, tvec) {
  val_t_minus_1 <- get_lag(vec, tvec, 1)
  vec - val_t_minus_1
}

# ==============================================================================
# Main Function: WesterlundPlain
# ==============================================================================
#' Compute Raw Westerlund ECM Panel Cointegration Statistics (Plain Routine)
#'
#' Internal plain (non-bootstrap) routine for computing the four Westerlund (2007)
#' ECM-based panel cointegration test statistics: \eqn{G_t}, \eqn{G_a}, \eqn{P_t},
#' and \eqn{P_a}. This function estimates unit-specific Error Correction Models (ECMs)
#' to form mean-group statistics and then constructs pooled statistics using
#' cross-unit aggregation and partialling-out steps.
#'
#' @param data A \code{data.frame} containing the panel data.
#' @param touse Logical vector of length \code{nrow(data)} indicating rows eligible
#'   for estimation.
#' @param idvar String. Column identifying cross-sectional units.
#' @param timevar String. Column identifying time.
#' @param yvar String. Name of the dependent variable (levels).
#' @param xvars Character vector. Names of regressors in the long-run relationship (levels).
#' @param constant Logical. If \code{TRUE}, includes a constant term in the ECM design matrix.
#' @param trend Logical. If \code{TRUE}, includes a linear time trend in the ECM design matrix.
#' @param lags Integer or length-2 integer vector. Fixed lag order or range
#'   \code{c(min,max)} for short-run dynamics. If a range is supplied, the routine
#'   performs an information-criterion search.
#' @param leads Integer or length-2 integer vector, or \code{NULL}. Fixed lead order
#'   or range \code{c(min,max)}. Defaults to 0 if \code{NULL}.
#' @param lrwindow Integer. Bartlett kernel window (maximum lag) used in long-run
#'   variance calculations via \code{\link{calc_lrvar_bartlett}}.
#' @param westerlund Logical. If \code{TRUE}, uses a Westerlund-specific
#'   information criterion and trimming logic for variance estimation.
#' @param aic Logical. If \code{TRUE}, uses AIC for lag/lead selection when ranges.
#'   If \code{FALSE}, uses BIC.
#' @param bootno Logical. If \code{TRUE}, suppresses certain header outputs,
#'   typically used when called from a bootstrap routine.
#' @param indiv.ecm Logical. If \code{TRUE}, gets output of individual ECM regressions.
#' @param verbose Logical. If \code{TRUE}, prints progress messages to the console.
#'
#' @details
#' \code{WesterlundPlain} is the computational core for the Westerlund test. It
#' processes the data in two major stages:
#'
#' \strong{Stage 1: Unit-Specific ECMs (Loop 1)}
#' For each unit \eqn{i}, the following ECM is estimated:
#' \deqn{\Delta y_{it} = \delta'_i d_t + \alpha_i y_{i,t-1} + \lambda'_i x_{i,t-1} + \sum_{j=1}^{p_i} \phi_{ij} \Delta y_{i,t-j} + \sum_{j=-q_i}^{q_i} \gamma_{ij} \Delta x_{i,t-j} + e_{it}}
#' where \eqn{d_t} contains the deterministic components (constant/trend). The
#' routine selects \eqn{p_i} (lags) and \eqn{q_i} (leads) per unit if ranges are
#' provided. The unit-level statistics are:
#' \itemize{
#'   \item \eqn{G_t = \frac{1}{N} \sum \frac{\hat{\alpha}_i}{SE(\hat{\alpha}_i)}}
#'   \item \eqn{G_a = \frac{1}{N} \sum \frac{T \hat{\alpha}_i}{\hat{\alpha}_i(1)}}
#' }
#'
#' \strong{Stage 2: Pooled Aggregation (Loop 2)}
#' The routine performs partialling-out regressions to remove the effects of
#' short-run dynamics and deterministics. It then aggregates residual products
#' across all units to compute the pooled statistics \eqn{P_t} and \eqn{P_a}.
#'
#'
#'
#' \strong{Time-Series Handling:}
#' All transformations (differences and lags) are performed using \code{\link{get_lag}}
#' and \code{\link{get_diff}}, ensuring that time-series gaps correctly result
#' in \code{NA} values rather than shifts across missing periods.
#'
#' @return A list containing three main components:
#' \itemize{
#'   \item \code{stats}: A list of the four raw Westerlund test statistics:
#'     \itemize{
#'       \item \code{Gt}: Mean-group statistic (tau).
#'       \item \code{Ga}: Mean-group statistic (alpha).
#'       \item \code{Pt}: Pooled statistic (tau).
#'       \item \code{Pa}: Pooled statistic (alpha).
#'     }
#'   \item \code{indiv_data}: A named list where each element corresponds to a cross-sectional unit (ID), containing:
#'     \itemize{
#'       \item \code{ai}: The estimated speed of adjustment (alpha).
#'       \item \code{seai}: The standard error of alpha (adjusted for degrees of freedom).
#'       \item \code{betai}: Vector of long-run coefficients (\eqn{\beta = -\lambda / \alpha}).
#'       \item \code{blag, blead}: The lags and leads selected for that specific unit.
#'       \item \code{ti}: Raw observation count for the unit.
#'       \item \code{tnorm}: Degrees of freedom used for normalization.
#'       \item \code{reg_coef}: If \code{indiv.ecm = TRUE}, the full coefficient matrix from \code{westerlund_test_reg}.
#'     }
#'   \item \code{results_df}: A summary \code{data.frame} containing all unit-level results in vectorized format.
#'   \item \code{settings}: A list of parameters used in the routine:
#'     \itemize{
#'       \item \code{constant}: Logical indicating if a constant was included.
#'       \item \code{trend}: Logical indicating if a trend was included.
#'       \item \code{meanlag, meanlead}: Integer averages of the selected lags/leads.
#'       \item \code{realmeanlag, realmeanlead}: Precise (floating-point) averages of lags/leads.
#'       \item \code{auto}: Logical indicating if automatic lag/lead selection was active.
#'     }
#' }
#'
#' @references
#' Westerlund, J. (2007). Testing for error correction in panel data.
#' \emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
#'
#' @seealso \code{\link{westerlund_test}}, \code{\link{WesterlundBootstrap}},
#'   \code{\link{calc_lrvar_bartlett}}, \code{\link{get_lag}}
#'
#' @importFrom stats BIC
#' @export
WesterlundPlain <- function(data, touse, idvar, timevar, yvar, xvars,
                            constant = FALSE, trend = FALSE, lags, leads = NULL,
                            lrwindow = 2, westerlund = FALSE, aic = TRUE,
                            bootno = FALSE, indiv.ecm = FALSE, verbose = FALSE) {

  if (bootno & verbose) cat("\nCalculating Westerlund ECM panel cointegration tests\n\n")

  # 1. Data Cleaning and Preparation
  # Sort strictly to ensure loop order matches Stata
  data <- data[order(data[[idvar]], data[[timevar]]), ]

  # Apply 'touse' and check valid Y/X
  valid_rows <- touse & !is.na(data[[yvar]])
  for (x in xvars) valid_rows <- valid_rows & !is.na(data[[x]])
  work_data <- data[valid_rows, ]

  # Calculate T_i (Average time periods per group)
  ti_counts <- table(work_data[[idvar]])
  meanbigT <- mean(as.numeric(ti_counts))

  nox <- length(xvars)

  # 2. Lag/Lead Handling (Stata Logic)
  numlags <- length(lags)
  minlag <- lags[1]; maxlag <- lags[1]

  if (numlags > 1) {
    minlag <- lags[1]; maxlag <- lags[2]

    # if min > max in syntax, swaps them
    if (minlag > maxlag) { temp <- maxlag; maxlag <- minlag; minlag <- temp }
  }

  minlead <- 0; maxlead <- 0

  if (!is.null(leads)) {
    if (length(leads) > 1) {
      minlead <- leads[1]; maxlead <- leads[2]
      if (minlead > maxlead) {t<-maxlead; maxlead<-minlead; minlead<-t}
    } else {
      minlead <- leads[1]; maxlead <- leads[1]
    }
  }

  auto <- (minlead != maxlead) | (minlag != maxlag)
  indiv_data <- list()
  ids <- unique(work_data[[idvar]])

  # Prepare Results Storage
  results_df <- data.frame(id=ids, ai=NA, seai=NA, aonesemi=NA, lags=NA, leads=NA, wysq=NA, wusq=NA, tnorm=NA)
  work_data$cons <- if(constant) 1 else 0
  work_data$tren <- 0
  if (trend) {
    work_data$tren <- ave(work_data[[timevar]], work_data[[idvar]],
                          FUN = function(z) seq_along(z))
  }
  counter <- 0; dots <- 0

  # ============================================================================
  # LOOP 1: Individual Regressions (G Stats)
  # ============================================================================
  for (i in seq_along(ids)) {
    uid <- ids[i]; counter <- counter + 1
    idx_rows <- which(work_data[[idvar]] == uid)
    sub_df <- work_data[idx_rows, ]

    # Extract time vector for gap-aware lag/diff
    tv <- sub_df[[timevar]]
    thisti <- length(tv)

    # Generate Time Series Vectors using Helper (Handles Gaps)
    dy <- get_diff(sub_df[[yvar]], tv)
    ly <- get_lag(sub_df[[yvar]], tv, 1)
    lx_list <- list()
    for (xv in xvars) lx_list[[xv]] <- get_lag(sub_df[[xv]], tv, 1)

    currlag <- maxlag; currlead <- maxlead

    # --- AIC/BIC Selection Loop ---
    if (auto) {
      curroptic <- Inf; curroptlag <- NA; curroptlead <- NA

      for (l in maxlag:minlag) {
        for (ld in maxlead:minlead) {
          # Build Regressor Matrix
          X_mat <- matrix(nrow=thisti, ncol=0)
          if (constant) X_mat <- cbind(X_mat, sub_df$cons)
          if (trend)    X_mat <- cbind(X_mat, sub_df$tren)
          X_mat <- cbind(X_mat, ly)
          for(xv in xvars) X_mat <- cbind(X_mat, lx_list[[xv]])

          # Lagged dy: L(1/l)D.y
          if (l > 0) {
            for (k in 1:l) {
              lag_dy <- get_lag(dy, tv, k)
              X_mat <- cbind(X_mat, lag_dy)
            }
          }

          # Dx Leads/Lags: L(-lead/lag)D.x
          for (xv in xvars) {
            dx <- get_diff(sub_df[[xv]], tv)

            # Leads (F.dx)
            if (ld > 0) {
              for (k in ld:1) {
                # Re-map value at time + k
                lead_dx <- sapply(tv, function(t) {
                  val_t_k <- dx[which(tv == t + k)]
                  if(length(val_t_k)==0) return(NA) else return(val_t_k)
                })
                X_mat <- cbind(X_mat, lead_dx)
              }
            }

            # Contemporary dx (L0)
            X_mat <- cbind(X_mat, dx)

            # Lags (L.dx)
            if (l > 0) {
              for (k in 1:l) {
                lag_dx <- get_lag(dx, tv, k)
                X_mat <- cbind(X_mat, lag_dx)
              }
            }
          }

          valid <- complete.cases(dy, X_mat)

          # Check df
          if (sum(valid) > ncol(X_mat) + 2) {
            # Run Regression (no intercept, as cons is manual)
            mod <- lm(dy[valid] ~ 0 + X_mat[valid, ])
            rss <- sum(resid(mod)^2)
            n_obs <- sum(valid)

            # --- AIC/BIC Calculation ---
            if(westerlund) {
              # Penalty: 2*(currlag + currlead + cons + tren + 1)
              nparams_stata <- l + ld + (if(constant) 1 else 0) + (if(trend) 1 else 0) + 1
              ic <- log(rss / (thisti - l - ld - 1)) + 2 * (nparams_stata) / (thisti - maxlag - maxlead)
            } else {
              # Stata Standard AIC (estat ic)
              # Formula: N * ln(RSS/N) + 2*k
              # Note: Stata's 'estat ic' typically uses -2*LL + 2*k.
              # For OLS, -2LL = N + N*log(2*pi) + N*log(RSS/N).
              # We use R's AIC() which follows -2LL + 2k.
              ic <- if(aic) AIC(mod) else BIC(mod)
            }
            if (ic < curroptic) {
              curroptic <- ic; curroptlag <- l; curroptlead <- ld
            }
          }
        }
      }
      currlag <- curroptlag; currlead <- curroptlead
    }

    # --- Recalculate Long Run Variance (wysq) ---
    dytmp <- dy

    # Stata trims based on selected lags (Line 214)
    # Using row indices relative to the FULL sub_df is risky with gaps.
    # However, Stata's code sets indices to missing based on absolute position.
    trim_start <- currlag + 1
    trim_end <- if (currlead > 0) currlead else 0

    if (westerlund) {
      if (trim_start > 0) dytmp[1:min(length(dytmp), trim_start)] <- NA
      if (trim_end > 0 && trim_end < length(dytmp)) {
        dytmp[(length(dytmp)-trim_end+1):length(dytmp)] <- NA
      }
    }

    # Detrending (Non-Westerlund only)
    if (!westerlund && constant && trend) {
      v <- !is.na(dytmp)
      if(sum(v) > 1) {
        # Regress on constant only (center)
        dytmp[v] <- resid(lm(dytmp[v] ~ 0 + sub_df$cons[v]))
      }
    }

    wysq <- calc_lrvar_bartlett(dytmp, maxlag=lrwindow, nodemean=TRUE)

    # --- Final Individual Regression ---

    # Rebuild X_mat with optimal lags
    X_mat <- matrix(nrow=thisti, ncol=0); cn <- c()

    if (constant) { X_mat<-cbind(X_mat, sub_df$cons); cn<-c(cn,"cons") }
    if (trend)    { X_mat<-cbind(X_mat, sub_df$tren); cn<-c(cn,"tren") }

    X_mat <- cbind(X_mat, ly); cn<-c(cn, paste0("l.", yvar))
    for(xv in xvars) { X_mat <- cbind(X_mat, lx_list[[xv]]); cn<-c(cn, paste0("l.", xv)) }

    if (currlag > 0) {
      for (k in 1:currlag) {
        X_mat <- cbind(X_mat, get_lag(dy, tv, k))
        cn <- c(cn, paste0("L",k,".dy"))
      }
    }

    for (xv in xvars) {
      dx <- get_diff(sub_df[[xv]], tv)
      # Leads
      if (currlead > 0) {
        for (k in currlead:1) {
          lead_dx <- sapply(tv, function(t) {
            val <- dx[which(tv == t + k)]
            if(length(val)==0) return(NA) else return(val)
          })
          X_mat <- cbind(X_mat, lead_dx)
          cn <- c(cn, paste0("F",k,".d.",xv))
        }
      }
      # Current
      X_mat <- cbind(X_mat, dx)
      cn <- c(cn, paste0("L0.d.",xv))
      # Lags
      if (currlag > 0) {
        for (k in 1:currlag) {
          X_mat <- cbind(X_mat, get_lag(dx, tv, k))
          cn <- c(cn, paste0("L",k,".d.",xv))
        }
      }
    }
    colnames(X_mat) <- cn

    valid <- complete.cases(dy, X_mat)
    if (sum(valid) > ncol(X_mat)) {
      reg_df <- as.data.frame(X_mat[valid, , drop=FALSE])
      reg_df$dep_var_dy <- dy[valid]

      # 1. Run Individual Unit Regression
      mod <- lm(dep_var_dy ~ 0 + ., data=reg_df)

      reg_coef <- NULL
      if (indiv.ecm) {
          if(verbose)cat(sprintf("\nIndividual ECM Regression for Unit: %s\n", uid))
          # Pass the estimated coefficients and the VCE matrix
          reg_coef <- westerlund_test_reg(b = coef(mod), V = vcov(mod), verbose = verbose,
                                           distribution = "t", df = df.residual(mod))
      }

      sc <- summary(mod)$coefficients
      cfs <- coef(mod)

      # 2. Dynamic Coefficient Extraction
      target_y <- paste0("l.", yvar)
      alpha_idx <- match(target_y, rownames(sc))

      if (!is.na(alpha_idx)) {
        # Speed of adjustment (Alpha)
        alpha_i <- sc[alpha_idx, "Estimate"]
        se_alpha_i <- sc[alpha_idx, "Std. Error"]

        # Long-run relationship (Beta = -gamma / alpha)
        beta_is <- numeric(length(xvars))
        names(beta_is) <- xvars
        for (xv in xvars) {
          target_x <- paste0("l.", xv)
          gamma_idx <- match(target_x, rownames(sc))
          if (!is.na(gamma_idx)) {
            beta_is[xv] <- -cfs[target_x] / alpha_i
          } else {
            beta_is[xv] <- 0
          }
        }

        # 3. Store Metadata for results_bundle
        indiv_data[[as.character(uid)]] <- list(
          gid = uid,
          ai = alpha_i,
          seai = se_alpha_i,
          betai = beta_is,
          blag = currlag,
          blead = currlead,
          ti = thisti,
          tnorm = 0, # Calculated below
          reg_coef = reg_coef
        )

        # 4. Extract results for Gt/Ga calculation (results_df)
        results_df$ai[i] <- alpha_i
        results_df$seai[i] <- se_alpha_i

        # --- Calculate u (Residuals for Sieve Bootstrap) ---
        u <- dy
        if (constant && "cons" %in% names(cfs)) u <- u - cfs["cons"] * sub_df$cons
        if (trend && "tren" %in% names(cfs)) u <- u - cfs["tren"] * sub_df$tren
        u <- u - cfs[target_y] * ly

        for (x in xvars) {
          nm <- paste0("l.", x)
          if (nm %in% names(cfs)) u <- u - cfs[nm] * lx_list[[x]]
        }

        if (currlag > 0) {
          for (v in 1:currlag) {
            nm <- paste0("L", v, ".dy")
            if (nm %in% names(cfs)) u <- u - cfs[nm] * get_lag(dy, tv, v)
          }
        }

        tmpu <- u
        # Trimming logic matches Stata's _n constraints
        trim_start <- currlag + 1
        trim_end <- if (currlead > 0) currlead else 0

        if (westerlund) {
          if (trim_start > 0) tmpu[1:min(length(tmpu), trim_start)] <- NA
          if (trim_end > 0 && trim_end < length(tmpu)) {
            tmpu[(length(tmpu) - trim_end + 1):length(tmpu)] <- NA
          }
        }

        # 5. Long-Run Variances (wysq, wusq)
        wusq <- calc_lrvar_bartlett(tmpu, maxlag = lrwindow, nodemean = TRUE)

        # Recalculate wysq with same trimming for consistency
        dytmp <- dy
        if (westerlund) {
          if (trim_start > 0) dytmp[1:min(length(dytmp), trim_start)] <- NA
          if (trim_end > 0 && trim_end < length(dytmp)) {
            dytmp[(length(dytmp) - trim_end + 1):length(dytmp)] <- NA
          }
        } else if (constant && trend) {
          v <- !is.na(dytmp)
          if(sum(v) > 1) dytmp[v] <- dytmp[v] - mean(dytmp[v], na.rm = TRUE)
        }
        wysq <- calc_lrvar_bartlett(dytmp, maxlag = lrwindow, nodemean = TRUE)

        results_df$aonesemi[i] <- sqrt(wusq / wysq)
        results_df$wysq[i] <- wysq
        results_df$wusq[i] <- wusq
        results_df$lags[i] <- currlag
        results_df$leads[i] <- currlead
      }
    }

    if (bootno && !indiv.ecm && counter > (dots * length(ids) / 10) && verbose) {
      cat(".")
      dots <- dots + 1
    }
  }
  if (bootno & verbose) cat("\n")

  # ============================================================================
  # LOOP 2: Pooled Regressions (P Stats)
  # ============================================================================

  meanlag <- as.integer(mean(results_df$lags, na.rm=TRUE))
  realmeanlag <- mean(results_df$lags, na.rm=TRUE)
  meanlead <- as.integer(mean(results_df$leads, na.rm=TRUE))
  realmeanlead <- mean(results_df$leads, na.rm=TRUE)

  pooledalphatop <- 0; pooledalphabottom <- 0; sum_sisq <- 0

  for (i in seq_along(ids)) {
    if (is.na(results_df$ai[i])) next
    uid <- ids[i]; idx <- which(work_data[[idvar]] == uid)
    sub_df <- work_data[idx, ]
    tv <- sub_df[[timevar]]; thisti <- length(tv)

    dy <- get_diff(sub_df[[yvar]], tv)
    ly <- get_lag(sub_df[[yvar]], tv, 1)

    # Prep dytmp
    dytmp <- dy
    trim_start <- meanlag + 1
    trim_end <- if (meanlead > 0) meanlead else 0

    if (westerlund) {
      if (trim_start > 0) dytmp[1:min(length(dytmp), trim_start)] <- NA
      if (trim_end > 0 && trim_end < length(dytmp)) {
        dytmp[(length(dytmp)-trim_end+1):length(dytmp)] <- NA
      }
    }

    # Detrending (matches Stata 'egen mean = mean(), by(id)')
    if (!westerlund && constant && trend) {
      dytmp <- dytmp - mean(dytmp, na.rm=TRUE)
    }

    wysq <- calc_lrvar_bartlett(dytmp, maxlag=lrwindow, nodemean=TRUE)

    # Build Z_mat (Common structure based on meanlag/meanlead)
    Z_mat <- matrix(nrow=thisti, ncol=0)

    # Stata: reg d.y cons tren l.y l.xvars L(1/meanlag)D.y L(-meanlead/meanlag)D.x

    if (constant) Z_mat <- cbind(Z_mat, sub_df$cons)
    if (trend)    Z_mat <- cbind(Z_mat, sub_df$tren)
    for (x in xvars) Z_mat <- cbind(Z_mat, get_lag(sub_df[[x]], tv, 1)) # l.x

    if (meanlag > 0) {
      for (k in 1:meanlag) Z_mat <- cbind(Z_mat, get_lag(dy, tv, k)) # L.dy
    }

    # Dx: Leads, Current, Lags
    for (x in xvars) {
      dx <- get_diff(sub_df[[x]], tv)
      if (meanlead > 0) {
        for (k in meanlead:1) {
          lead_dx <- sapply(tv, function(t) {
            val <- dx[which(tv == t + k)]
            if(length(val)==0) return(NA) else return(val)
          })
          Z_mat <- cbind(Z_mat, lead_dx)
        }
      }
      Z_mat <- cbind(Z_mat, dx)
      if (meanlag > 0) {
        for (k in 1:meanlag) Z_mat <- cbind(Z_mat, get_lag(dx, tv, k))
      }
    }

    # Auxiliary regressions (Partialling out Z_mat)
    # Stata: reg d.y Z_mat (excluding l.y from Z_mat)
    v_dy <- complete.cases(dy, Z_mat)
    v_ly <- complete.cases(ly, Z_mat)

    if (sum(v_dy) <= ncol(Z_mat) || sum(v_ly) <= ncol(Z_mat)) next

    dyresid <- rep(NA, thisti)
    dyresid[v_dy] <- resid(lm(dy[v_dy] ~ 0 + Z_mat[v_dy, ]))

    yresid <- rep(NA, thisti)
    yresid[v_ly] <- resid(lm(ly[v_ly] ~ 0 + Z_mat[v_ly, ]))

    # Full ECM Regression for Sigma/u
    full_X <- cbind(ly, Z_mat)
    v_f <- complete.cases(dy, full_X)
    if (sum(v_f) <= ncol(full_X)) next

    mod_full <- lm(dy[v_f] ~ 0 + full_X[v_f, ])
    cfs <- coef(mod_full)

    # Calculate Pooled u
    # Stata subtracts: cons, tren, l.x, L(meanlag)D.y
    # Crucially, it relies on column order of variables in the regression
    # Stata reg: ly, cons, tren, lx, ldy, leads/lags_dx

    # Define count of "Long Run" parameters to subtract
    count_sub <- 1 + (if(constant) 1 else 0) + (if(trend) 1 else 0) + nox + meanlag

    # Note: full_X = [ly, cons, tren, l.x, l.dy, sr.dx]
    # The first 'count_sub' columns of full_X correspond exactly to Stata's subtracted variables

    if (ncol(full_X) >= count_sub) {
      u <- dy - (full_X[, 1:count_sub, drop=FALSE] %*% cfs[1:count_sub])[,1]
    } else {
      # Fallback if structure is weird, though shouldn't happen with strict counts
      u <- dy - resid(mod_full)
    }

    tmpu <- u
    if (westerlund) {
      if (trim_start > 0) tmpu[1:min(length(tmpu), trim_start)] <- NA
      if (trim_end > 0 && trim_end < length(tmpu)) {
        tmpu[(length(tmpu)-trim_end+1):length(tmpu)] <- NA
      }
    }

    wusq_pool <- calc_lrvar_bartlett(tmpu, maxlag=lrwindow, nodemean=TRUE)
    aonesemipool <- sqrt(wusq_pool / wysq)

    vp <- complete.cases(yresid, dyresid)
    if (sum(vp) > 0) {
      pooledalphatop <- pooledalphatop + sum((1/aonesemipool) * yresid[vp] * dyresid[vp])
      pooledalphabottom <- pooledalphabottom + sum(yresid[vp]^2)
    }

    sigmasqi <- sum(resid(mod_full)^2)

   if (westerlund) {
      tnorm <- meanbigT - meanlag - meanlead - 1
    } else {
      kp <- (if(constant) 1 else 0) + (if(trend) 1 else 0) + nox + meanlag + nox*(meanlag+meanlead+1)
      tnorm <- meanbigT - meanlag - meanlead - 1 - kp - 1
    }

    if (tnorm > 0) {
      val <- (sqrt(sigmasqi/tnorm) / aonesemipool)^2
      if (!is.na(val)) sum_sisq <- sum_sisq + val
    }
  }

  # ============================================================================
  # FINAL STATISTICS
  # ============================================================================

  valid_res <- !is.na(results_df$ai)
  n_groups <- sum(valid_res)

  pooledalpha <- pooledalphatop / pooledalphabottom
  se_pooledalpha <- sqrt(sum_sisq / n_groups) / sqrt(pooledalphabottom)

  for (i in 1:nrow(results_df)) {
    if (!valid_res[i]) next

    ti <- ti_counts[as.character(results_df$id[i])]
    l <- results_df$lags[i]; ld <- results_df$leads[i]
    kp <- (if(constant) 1 else 0) + (if(trend) 1 else 0) + nox + l + nox*(l+ld+1)

    if (westerlund) {
      tn <- ti - l - ld - 1
      alt <- tn - (kp + 1)
      results_df$seai[i] <- results_df$seai[i] * sqrt(alt) / sqrt(tn)
      results_df$tnorm[i] <- tn
    } else {
      results_df$tnorm[i] <- ti - l - ld - 1 - kp - 1
    }
  }

  # Synchronize indiv_data with final adjusted results_df
  for (i in seq_along(ids)) {
    uid_str <- as.character(ids[i])
    if (!is.null(indiv_data[[uid_str]])) {
      indiv_data[[uid_str]]$seai <- results_df$seai[i]
      indiv_data[[uid_str]]$tnorm <- results_df$tnorm[i]
      indiv_data[[uid_str]]$aonesemi <- results_df$aonesemi[i]
    }
  }

  Gt <- mean(results_df$ai / results_df$seai, na.rm=TRUE)
  Ga <- mean(results_df$tnorm * results_df$ai / results_df$aonesemi, na.rm=TRUE)
  Pt <- pooledalpha / se_pooledalpha

  kp <- (if(constant) 1 else 0) + (if(trend) 1 else 0) + nox + meanlag + nox*(meanlag+meanlead+1)
  tnorm_pool <- if (westerlund) meanbigT - meanlag - meanlead - 1 else meanbigT - meanlag - meanlead - 1 - kp - 1
  Pa <- tnorm_pool * pooledalpha

  return(list(
    stats = list(Gt=Gt, Ga=Ga, Pt=Pt, Pa=Pa),
    indiv_data = indiv_data,
    results_df = results_df,
    settings = list(constant=constant, trend=trend, meanlag=meanlag, meanlead=meanlead, realmeanlag=realmeanlag, realmeanlead=realmeanlead,auto=auto)
  ))
}

# ==============================================================================
# Main Function: WesterlundBootstrap
# ==============================================================================
#' Bootstrap Routine for Westerlund ECM Panel Cointegration Tests
#'
#' Internal bootstrap routine for generating a bootstrap distribution of the
#' Westerlund (2007) ECM-based panel cointegration statistics under the null
#' hypothesis of no error correction. The routine estimates short-run dynamics for
#' each unit, constructs residual-based bootstrap samples in a Stata-aligned manner,
#' and re-computes \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a} by calling
#' \code{\link{WesterlundPlain}} on each bootstrap sample.
#'
#' @param data A \code{data.frame} containing the original panel data.
#' @param touse Logical vector of length \code{nrow(data)} indicating rows eligible
#'   for estimation.
#' @param idvar String. Column identifying cross-sectional units.
#' @param timevar String. Column identifying time. Time gaps are handled via
#'   Stata-style differencing.
#' @param yvar String. Name of the dependent variable (levels).
#' @param xvars Character vector. Names of regressors (levels).
#' @param constant Logical. Include a constant in the short-run bootstrap regression.
#' @param trend Logical. Include a trend in the short-run bootstrap regression.
#' @param lags Integer or length-2 integer vector. Fixed lag order or range
#'   \code{c(min,max)} for selecting short-run dynamics in the bootstrap setup.
#' @param leads Integer or length-2 integer vector, or \code{NULL}. Fixed lead order
#'   or range. Defaults to 0 if \code{NULL}.
#' @param westerlund Logical. If \code{TRUE}, uses the Westerlund-specific
#'   lag/lead selection criterion in the restricted bootstrap regression:
#'   \eqn{\log(RSS/(T_i-p-q-1)) + 2(p+q+c+r+1)/(T_i-p_{\max}-q_{\max})},
#'   where \eqn{c} and \eqn{r} indicate the constant and trend terms.
#' @param aic Logical. If \code{TRUE}, uses AIC for lag/lead selection when ranges.
#'   If \code{FALSE}, uses BIC.
#' @param bootstrap Integer. Number of bootstrap replications.
#' @param indiv.ecm Logical. If \code{TRUE}, gets output of individual ECM regressions.
#' @param lrwindow Integer. Bartlett kernel window forwarded to
#'   \code{\link{WesterlundPlain}} during re-estimation.
#' @param seed Optional integer random seed. If \code{NULL}, the current R RNG state
#'   is used.
#' @param verbose Logical. If \code{TRUE}, prints additional output.
#'
#' @details
#' \code{WesterlundBootstrap} is an internal engine used by \code{\link{westerlund_test}}
#' to provide inference when asymptotic distributions may be unreliable.
#'
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item \strong{Residual Extraction:} Estimates a restricted short-run model for
#'     \eqn{\Delta y_{it}} under the null of no error correction and saves residuals
#'     \eqn{e_{it}} and short-run coefficients.
#'   \item \strong{Centering:} Demeans residuals within each unit and centers each
#'     \eqn{\Delta x_{it}} over all usable (\code{touse}) observations for that unit,
#'     matching the original \code{xtwest} bootstrap setup.
#'   \item \strong{Resampling and alignment:} Reproduces the Stata-style
#'     \code{expandcl 2}, actual-time \code{cluster(t)} sampling, and
#'     \code{newttt}/\code{tussent}/\code{newtt} alignment so that sampled time
#'     clusters are shared across units while unit-specific panel lengths are retained.
#'   \item \strong{Recursion:} Constructs bootstrap innovations and recursively
#'     generates \eqn{\Delta y} using the stored AR and regressor-difference coefficients.
#'   \item \strong{Integration:} Integrates bootstrap differences to levels and applies
#'     unit-specific Stata truncation rules.
#'   \item \strong{Re-Testing:} Calls \code{\link{WesterlundPlain}} on the simulated
#'     panel to store \eqn{G_t}, \eqn{G_a}, \eqn{P_t}, and \eqn{P_a}.
#' }
#'
#' \strong{Time Indexing:}
#' The bootstrap setup uses a local \code{diff_ts()} helper that mimics Stata's
#' \code{D.} operator under gaps: differences are marked \code{NA} when adjacent
#' time values are not exactly one unit apart.
#'
#' \strong{Lag/lead selection:}
#' When ranges are supplied, selection is performed separately for each unit using
#' the restricted bootstrap regression, without lagged level terms. In
#' \code{westerlund=TRUE} mode the criterion is
#' \deqn{IC(p,q)=\log\left(\frac{RSS}{T_i-p-q-1}\right)
#' + \frac{2(p+q+c+r+1)}{T_i-p_{\max}-q_{\max}}.}
#' Otherwise AIC or BIC is used according to \code{aic}.
#'
#' \strong{Compatibility note for automatic leads:}
#' For literal compatibility with the 2010 \code{xtwest} source, the bootstrap
#' reconstruction uses the final processed unit's selected \code{currlead} as the
#' global upper bound for lead terms. This reproduces an implementation detail of
#' that program rather than a separate theoretical requirement.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{BOOTSTATS}: A numeric matrix with \code{bootstrap} rows and 4 columns
#'     containing the simulated replications of (\eqn{G_t, G_a, P_t, P_a}).
#' }
#'
#' @references
#' Westerlund, J. (2007). Testing for error correction in panel data.
#' \emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
#'
#' @seealso \code{\link{westerlund_test}}, \code{\link{WesterlundPlain}}
#'
#' @export
WesterlundBootstrap <- function(data, touse, idvar, timevar, yvar, xvars,
                                constant = FALSE, trend = FALSE, lags = 1, leads = NULL,
                                westerlund = FALSE, aic = TRUE, bootstrap = 100,
                                indiv.ecm = FALSE, lrwindow = 2, seed = NULL, verbose = FALSE) {

  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  if (!is.null(seed)) set.seed(as.integer(seed))
  if (verbose) cat("\nBootstrapping critical values under H0\n")

  # --- Helper: Lag/Lead Vector ---
  L_vec <- function(x, k) {
    if (k == 0) return(x)
    n <- length(x)
    out <- rep(NA, n)
    if (k > 0) {
      if(n > k) out[(k+1):n] <- x[1:(n-k)]
    } else {
      k <- abs(k)
      if(n > k) out[1:(n-k)] <- x[(k+1):n]
    }
    return(out)
  }

  # --- Helper: Stata-style Difference (Respects Time Gaps) ---
  # Replicates D.var: returns NA if previous row is not exactly time-1
  diff_ts <- function(vec, timevec) {
    d <- diff(vec)
    dt <- diff(timevec)
    # If time gap != 1, the difference is invalid (NA)
    d[dt != 1] <- NA
    return(c(NA, d))
  }

  # 1. Setup
  data <- data[order(data[[idvar]], data[[timevar]]), ]

  # Stata marksample: exclude if any var is NA
  touse[is.na(data[[yvar]])] <- FALSE
  for (x in xvars) touse[is.na(data[[x]])] <- FALSE

  ti_counts <- table(data[[idvar]][touse])

  # Lags/Leads Logic
  numlags <- length(lags)
  if (numlags > 1) {
    minlag <- lags[1]; maxlag <- lags[2]
    if (minlag > maxlag) {t<-maxlag; maxlag<-minlag; minlag<-t}
  } else {
    minlag <- lags[1]; maxlag <- lags[1]
  }

  minlead <- 0; maxlead <- 0
  if (!is.null(leads)) {
    if(length(leads)>1) {
      minlead<-leads[1]; maxlead<-leads[2]
      if(minlead>maxlead) {t<-maxlead; maxlead<-minlead; minlead<-t}
    } else {
      minlead<-leads[1]; maxlead<-leads[1]
    }
  }

  auto <- (minlead != maxlead) | (minlag != maxlag)

  work_data <- data
  work_data$e <- NA
  valid_ids <- unique(work_data[[idvar]][touse])
  coeff_store <- list()

  # xtwest quirk: after the per-ID loop, currlead retains the value selected
  # for the final panel unit and controls the global bootstrap lead loop.
  stata_currlead <- maxlead

  # 2. Get Coefficients and Residuals
  for (uid in valid_ids) {
    idx <- which(work_data[[idvar]] == uid & touse)
    if (length(idx) == 0) next

    sub_df <- work_data[idx, ]
    sub_t  <- sub_df[[timevar]]

    # Use diff_ts to handle time gaps
    dy <- diff_ts(sub_df[[yvar]], sub_t)

    thisti <- length(idx)
    tren <- if (trend) 1:length(idx) else rep(0, length(idx))
    cons <- if (constant) rep(1, length(idx)) else rep(0, length(idx))

    # Optimization Loop
    currlag <- maxlag; currlead <- maxlead

    if (auto) {
      curroptic <- Inf; curroptlag <- NA; curroptlead <- NA
      # Stata loop: maxlag down to minlag
      for (l in maxlag:minlag) {
        for (ld in maxlead:minlead) {
          X_mat <- matrix(nrow=length(dy), ncol=0)

          if (constant) X_mat <- cbind(X_mat, cons)
          if (trend) X_mat <- cbind(X_mat, tren)
          if (l > 0) for (k in 1:l) X_mat <- cbind(X_mat, L_vec(dy, k))

          for (x in xvars) {
            dx <- diff_ts(sub_df[[x]], sub_t)
            if (ld > 0) for (k in ld:1) X_mat <- cbind(X_mat, L_vec(dx, -k))
            X_mat <- cbind(X_mat, L_vec(dx, 0))
            if (l > 0) for (k in 1:l) X_mat <- cbind(X_mat, L_vec(dx, k))
          }

          valid_rows <- complete.cases(X_mat, dy)
          if (sum(valid_rows) > (ncol(X_mat)+1)) {
            mod <- lm(dy[valid_rows] ~ 0 + X_mat[valid_rows, ])
            rss <- sum(resid(mod)^2)
            ic <- if(westerlund) {
              # Literal xtwest bootstrap IC:
              # log(RSS/(Ti-p-q-1)) + 2*(p+q+cons+trend+1)/(Ti-pmax-qmax)
              nparams_stata <- l + ld + as.integer(constant) + as.integer(trend) + 1
              log(rss/(thisti-l-ld-1)) + 2*nparams_stata/(thisti-maxlag-maxlead)
            } else {
              # AIC: Stata includes constant and sigma in k.
              # R's logLik includes them too.
              ic <- if(aic) AIC(mod) else BIC(mod)
            }

            if (!is.na(ic) && ic < curroptic) {
              curroptic <- ic; curroptlag <- l; curroptlead <- ld
            }
          }
        }
      }
      currlag <- curroptlag
      currlead <- curroptlead
    }

    # Preserve the literal 2010 xtwest behavior: this is overwritten for each
    # unit, leaving the last unit's selected lead in scope after the loop.
    stata_currlead <- currlead

    # Final Regression with Optimal Lags
    X_mat <- matrix(nrow=length(dy), ncol=0); cn <- c()
    if (constant) {X_mat<-cbind(X_mat, cons); cn<-c(cn, "cons")}
    if (trend) {X_mat<-cbind(X_mat, tren); cn<-c(cn, "tren")}
    if (currlag > 0) for (k in 1:currlag) {X_mat<-cbind(X_mat, L_vec(dy, k)); cn<-c(cn, paste0("L",k,".dy"))}
    for (x in xvars) {
      dx <- diff_ts(sub_df[[x]], sub_t)
      if (currlead > 0) for (k in currlead:1) {X_mat<-cbind(X_mat, L_vec(dx, -k)); cn<-c(cn, paste0("F",k,".d.",x))}
      X_mat<-cbind(X_mat, L_vec(dx, 0)); cn<-c(cn, paste0("L0.d.",x))
      if (currlag > 0) for (k in 1:currlag) {X_mat<-cbind(X_mat, L_vec(dx, k)); cn<-c(cn, paste0("L",k,".d.",x))}
    }

    colnames(X_mat) <- cn
    valid_rows <- complete.cases(X_mat, dy)

    if (sum(valid_rows) > 0) {
      mod <- lm(dy[valid_rows] ~ 0 + X_mat[valid_rows, ])
      coeffs <- coef(mod)
      coeffs[is.na(coeffs)] <- 0

      # Replicate Stata predict (score)
      fitted_vals <- X_mat %*% coeffs

      # Residuals
      work_data$e[idx] <- dy - fitted_vals

      # Store Coefficients
      cvals <- coef(mod)
      id_coeffs <- list(xvars=list(), yvars=numeric(maxlag))

      get_c <- function(n) {
        i <- which(cn == n)
        if(length(i)>0) {
          val <- cvals[i]
          if(is.na(val)) 0 else val
        } else 0
      }

      for (x in xvars) {
        xc <- list(L=numeric(maxlag+1), F=numeric(maxlead+1))
        for (l in 0:currlag) xc$L[l+1] <- get_c(paste0("L",l,".d.",x))
        if (currlead>0) for (l in 1:currlead) xc$F[l+1] <- get_c(paste0("F",l,".d.",x))
        id_coeffs$xvars[[x]] <- xc
      }
      if (currlag > 0) {
        for (l in 1:currlag) id_coeffs$yvars[l] <- get_c(paste0("L",l,".dy"))
      }
      coeff_store[[as.character(uid)]] <- id_coeffs
    }
  }

  # 3. Demeaning

  # Demean e within id over its nonmissing residual support
  means_e <- tapply(work_data$e, work_data[[idvar]], mean, na.rm = TRUE)
  work_data$e <- work_data$e - means_e[as.character(work_data[[idvar]])]

  centerd_names <- character(0)

  for (x in xvars) {
    dx_name <- paste0("d", x)
    cd_name <- paste0("centerd", x)

    # compute D.x within id
    work_data[[dx_name]] <- NA_real_

    for (uid in valid_ids) {
      id_rows <- which(work_data[[idvar]] == uid & touse)
      if (length(id_rows) == 0) next
      sub_t <- work_data[[timevar]][id_rows]
      sub_x <- work_data[[x]][id_rows]
      work_data[[dx_name]][id_rows] <- diff_ts(sub_x, sub_t)
    }

    # Stata: by id: egen meandx = mean(dx) if touse
    # Center D.x over all usable observations for the unit, not only the
    # narrower residual-support rows. Missing D.x values are ignored by mean().
    work_data[[cd_name]] <- NA_real_
    for (uid in valid_ids) {
      id_rows <- which(work_data[[idvar]] == uid & touse)
      if (length(id_rows) == 0) next
      dx_sub <- work_data[[dx_name]][id_rows]
      mu <- mean(dx_sub, na.rm = TRUE)
      work_data[[cd_name]][id_rows] <- dx_sub - mu
    }

    centerd_names <- c(centerd_names, cd_name)
  }


  BOOTSTATS <- matrix(NA, nrow = bootstrap, ncol = 4)

  # 4. Bootstrap Loop

  # Build the eligible pool exactly on the support used by Stata's
  #     bsample if e < ., cluster(t)
  # while retaining the ACTUAL time labels and each unit's original Ti.
  eligible_pool <- work_data[touse & !is.na(work_data$e),
                             c(idvar, timevar, "e", centerd_names),
                             drop = FALSE]
  if (nrow(eligible_pool) == 0) stop("No valid observations are available for bootstrap resampling.")

  ti_lookup <- as.integer(ti_counts)
  names(ti_lookup) <- names(ti_counts)
  eligible_pool$.__Ti <- ti_lookup[as.character(eligible_pool[[idvar]])]
  if (any(is.na(eligible_pool$.__Ti))) stop("Could not recover unit-specific Ti for bootstrap pool.")

  # Stata: expandcl 2, cluster(id).  The generated expansion marker is not
  # subsequently used by xtwest; two copies are sufficient for the next step.
  expanded_pool <- rbind(
    transform(eligible_pool, .__expanded = 1L),
    transform(eligible_pool, .__expanded = 2L)
  )
  rownames(expanded_pool) <- NULL

  # Stata's cluster(t) uses actual time values, not within-unit row positions.
  time_keys <- unique(eligible_pool[[timevar]])
  if (length(time_keys) == 0) stop("No eligible time clusters are available for bootstrap resampling.")
  expanded_pool$.__time_cluster <- match(expanded_pool[[timevar]], time_keys)
  rows_by_time <- split(expanded_pool, expanded_pool$.__time_cluster, drop = TRUE)

  for (b in seq_len(bootstrap)) {

    # Equivalent sampling unit to: bsample if e < ., cluster(t)
    draw_pos <- sample(seq_along(time_keys), size = length(time_keys), replace = TRUE)
    sampled_parts <- lapply(draw_pos, function(ii) rows_by_time[[as.character(ii)]])
    sampled <- do.call(rbind, sampled_parts)
    rownames(sampled) <- NULL

    # Stata: gen newttt = ceil(uniform()*(ti-maxlag-maxlead-1))
    upper <- as.integer(sampled$.__Ti) - maxlag - maxlead - 1L
    if (any(upper <= 0L)) stop("Insufficient observations for xtwest bootstrap shuffle.")
    sampled$.__newttt <- ceiling(runif(nrow(sampled)) * upper)

    # sort id, stable; by id: gen tussent = _n
    sampled$.__stable0 <- seq_len(nrow(sampled))
    sampled <- sampled[order(sampled[[idvar]], sampled$.__stable0), , drop = FALSE]
    rownames(sampled) <- NULL
    sampled$.__tussent <- ave(seq_len(nrow(sampled)), sampled[[idvar]], FUN = seq_along)

    # sort tussent, stable; by tussent: egen newtt = mean(newttt)
    sampled$.__stable1 <- seq_len(nrow(sampled))
    sampled <- sampled[order(sampled$.__tussent, sampled$.__stable1), , drop = FALSE]
    rownames(sampled) <- NULL
    sampled$.__newtt <- ave(sampled$.__newttt, sampled$.__tussent, FUN = mean)

    # sort id newtt, stable
    sampled$.__stable2 <- seq_len(nrow(sampled))
    sampled <- sampled[order(sampled[[idvar]], sampled$.__newtt, sampled$.__stable2), , drop = FALSE]
    rownames(sampled) <- NULL

    boot_data_list <- vector("list", length(valid_ids))
    names(boot_data_list) <- as.character(valid_ids)

    for (uid in valid_ids) {
      Ti_i <- as.integer(ti_lookup[as.character(uid)])
      if (is.na(Ti_i)) next

      # Stata: by id: keep if _n <= ti + maxlag + maxlead + 2
      sg <- sampled[sampled[[idvar]] == uid, , drop = FALSE]
      required_len <- Ti_i + maxlag + maxlead + 2L
      if (nrow(sg) > required_len) sg <- sg[seq_len(required_len), , drop = FALSE]
      if (nrow(sg) == 0) next

      Text_i <- nrow(sg)
      e_b <- as.numeric(sg$e)
      cdX_b <- as.matrix(sg[, centerd_names, drop = FALSE])
      storage.mode(cdX_b) <- "double"

      cfs_id <- coeff_store[[as.character(uid)]]
      if (is.null(cfs_id)) next

      # -------------------------- e -> u --------------------------
      u <- e_b
      if (maxlag > 0) u[seq_len(min(maxlag, Text_i))] <- NA_real_
      u[is.na(u)] <- 0

      # Stata loops over 0..maxlag globally; coefficients outside an ID's
      # selected lag are zero-filled. Do not skip zero coefficients: Stata's
      # missing * 0 remains missing and is only zeroed when dy is initialized.
      for (j in seq_along(xvars)) {
        x <- xvars[j]
        cfs <- cfs_id$xvars[[x]]
        if (is.null(cfs)) next
        v_sub <- cdX_b[, j]

        for (L in 0:maxlag) {
          beta <- cfs$L[L + 1L]
          if (is.na(beta)) beta <- 0
          u <- u + shiftNA(v_sub, L) * beta
        }

        # Literal xtwest currlead quirk: the final unit's selected lead controls
        # the upper limit for every unit; unused unit coefficients are zero.
        if (stata_currlead > 0) {
          for (F in seq_len(stata_currlead)) {
            beta <- cfs$F[F + 1L]
            if (is.na(beta)) beta <- 0
            u <- u + shiftNA(v_sub, -F) * beta
          }
        }
      }

      # ------------------------- u -> dy --------------------------
      # Stata: gen dy=u; replace dy=0 if dy>=.; then AR recursion.
      u_for_dy <- u
      u_for_dy[is.na(u_for_dy)] <- 0
      dy <- u_for_dy

      phi <- cfs_id$yvars
      if (length(phi) < maxlag) phi <- c(phi, rep(0, maxlag - length(phi)))

      if (maxlag > 0 && Text_i > maxlag) {
        for (tt in (maxlag + 1L):Text_i) {
          ar_term <- 0
          for (L in seq_len(maxlag)) {
            bL <- phi[L]
            if (is.na(bL)) bL <- 0
            ar_term <- ar_term + dy[tt - L] * bL
          }
          dy[tt] <- dy[tt] + ar_term
        }
      }
      if (maxlag > 0) dy[seq_len(min(maxlag, Text_i))] <- NA_real_

      # -------------------- integrate to levels -------------------
      booty <- cumsum(replace(dy, is.na(dy), 0))
      if (maxlag > 0) booty[seq_len(min(maxlag, Text_i))] <- NA_real_

      bootX <- matrix(NA_real_, nrow = Text_i, ncol = length(xvars))
      for (j in seq_along(xvars)) {
        bootX[, j] <- cumsum(replace(cdX_b[, j], is.na(cdX_b[, j]), 0))
      }
      if (maxlag > 0) bootX[seq_len(min(maxlag, Text_i)), ] <- NA_real_

      # Stata: replace bootstrap levels = . if _n > ti + maxlag
      upper_keep <- Ti_i + maxlag
      if (upper_keep < Text_i) {
        bad <- (upper_keep + 1L):Text_i
        booty[bad] <- NA_real_
        bootX[bad, ] <- NA_real_
      }

      # Like tsset id newt in Stata.  The generated time index is unit-specific
      # and starts at one after the common newtt sorting step.
      df_b <- data.frame(
        id = rep(as.character(uid), Text_i),
        t = as.integer(seq_len(Text_i)),
        booty = booty,
        stringsAsFactors = FALSE
      )
      for (j in seq_along(xvars)) df_b[[paste0("boot", xvars[j])]] <- bootX[, j]
      boot_data_list[[as.character(uid)]] <- df_b
    }

    boot_data_list <- Filter(Negate(is.null), boot_data_list)
    if (length(boot_data_list) == 0) stop("Bootstrap replication produced no panel observations.")
    boot_data <- do.call(rbind, boot_data_list)
    rownames(boot_data) <- NULL

    # WesterlundPlain's marksample excludes these boundary rows. Dropping them
    # here leaves the same usable support and preserves the generated newt labels.
    keep <- !is.na(boot_data$booty)
    for (j in seq_along(xvars)) keep <- keep & !is.na(boot_data[[paste0("boot", xvars[j])]])
    boot_data <- boot_data[keep, , drop = FALSE]

    touse_boot <- rep(TRUE, nrow(boot_data))

    res <- WesterlundPlain(
      boot_data, touse_boot,
      idvar = "id", timevar = "t",
      yvar = "booty", xvars = paste0("boot", xvars),
      constant = constant, trend = trend,
      lags = lags, leads = leads,
      lrwindow = lrwindow, westerlund = westerlund, aic = aic,
      bootno = FALSE, verbose = verbose
    )

    BOOTSTATS[b, ] <- c(res$stats$Gt, res$stats$Ga, res$stats$Pt, res$stats$Pa)
  }

  return(list(BOOTSTATS = BOOTSTATS))
}

#' NA-Padded Lag and Lead Operator
#'
#' Shifts a vector forward or backward by a specified number of positions and
#' fills out-of-range values with \code{NA}. This helper is designed for
#' position-based transformations where missing values should be explicitly
#' propagated as \code{NA}.
#'
#' @param v A vector (typically numeric, but can be any type).
#' @param k Integer. Shift order. Positive values correspond to lags
#' (shifting the series downward, i.e., $t-k$), while negative values
#' correspond to leads (shifting the series upward, i.e., $t+k$).
#'
#' @details
#' This function performs a position-based shift of the input vector \code{v}.
#' Unlike strict time-indexed helpers (such as \code{\link{get_lag}}), \code{shiftNA()}
#' does not rely on an explicit time index and does not propagate gaps based on
#' timestamps. Instead, it performs a simple index shift, padding the
#' resulting empty slots with \code{NA}.
#'
#' Let $n$ denote the length of \code{v}. The behavior is:
#' \itemize{
#'   \item \code{k = 0}: return \code{v} unchanged,
#'   \item \code{k > 0}: lag by \code{k} positions; the first \code{k} elements are \code{NA},
#'   \item \code{k < 0}: lead by \eqn{|k|} positions; the last \eqn{|k|} elements are \code{NA}.
#' }
#'
#' Formally, for $k > 0$:
#' \deqn{
#' \text{out}_t =
#' \begin{cases}
#' \text{NA}, & t \le k, \\
#' v_{t-k}, & t > k,
#' \end{cases}
#' }
#' and for $k < 0$:
#' \deqn{
#' \text{out}_t =
#' \begin{cases}
#' v_{t+|k|}, & t \le n-|k|, \\
#' \text{NA}, & t > n-|k|.
#' \end{cases}
#' }
#'
#' @return A vector of the same length as \code{v}, containing the shifted values
#' with \code{NA} inserted where the shift would otherwise go out of range.
#'
#' @section Usage and Propagation:
#' This section explains the implications of using \code{NA} padding in
#' recursive or difference-based calculations.
#'
#' \bold{Why NA-padding?}
#' \cr
#' Using \code{NA} is the standard behavior in R for out-of-bounds operations.
#' It ensures that subsequent calculations (like \code{v - shiftNA(v, 1)})
#' correctly result in \code{NA} for the boundary cases, preventing the
#' accidental use of arbitrary values (like zero) in statistical estimations
#' unless explicitly intended.
#'
#' \bold{Difference from get_lag()}
#' \cr
#' Unlike \code{\link{get_lag}}, which may handle time-series objects or
#' specific index matching, \code{shiftNA()} is a "dumb" position shifter.
#' It is faster for internal loops where the user guarantees that the
#' vector is already sorted and contiguous.
#'
#' @examples
#' v <- c(1, 2, 3, 4, 5)
#'
#' ## Lag by one (k = 1)
#' shiftNA(v, 1)
#' # [1] NA  1  2  3  4
#'
#' ## Lead by one (k = -1)
#' shiftNA(v, -1)
#' # [1]  2  3  4  5 NA
#'
#' ## Larger shifts
#' shiftNA(v, 2)
#' # [1] NA NA  1  2  3
#'
#' shiftNA(v, -2)
#' # [1]  3  4  5 NA NA
#'
#' @seealso \code{\link{WesterlundBootstrap}}, \code{\link{get_lag}}, \code{\link{get_diff}}
#' @export
shiftNA <- function(v, k) {
  n <- length(v)
  out <- rep(NA, n) # Start with NA
  if (k == 0) return(v)
  if (k > 0) {
    if (n > k) out[(k + 1):n] <- v[1:(n - k)]
  } else {
    kk <- abs(k)
    if (n > kk) out[1:(n - kk)] <- v[(kk + 1):n]
  }
  return(out)
}

# ============================================================================
# Main Function: DisplayWesterlund
# ==============================================================================
#' Standardize and Display Westerlund ECM Panel Cointegration Test Results
#'
#' Formats, standardizes, and prints the Westerlund (2007) ECM-based panel
#' cointegration test results. This function computes standardized Z-statistics
#' using asymptotic moments and reports left-tail asymptotic p-values. If a
#' bootstrap distribution is provided, it also computes robust bootstrap p-values.
#'
#' @param stats A numeric vector of length 4 or a \code{1x4} matrix containing
#'   the raw test statistics in the order: \code{Gt, Ga, Pt, Pa}.
#' @param bootstats Optional. A numeric matrix with bootstrap replications of the
#'   raw statistics (4 columns: \code{Gt, Ga, Pt, Pa}).
#' @param nobs Integer. Number of valid cross-sectional units (N) in the panel.
#' @param nox Integer. Number of regressors in the long-run relationship.
#' @param constant Logical. If \code{TRUE}, a constant was included in the model.
#' @param trend Logical. If \code{TRUE}, a linear trend was included.
#' @param meanlag Integer. Average selected lag length (rounded) for display.
#' @param meanlead Integer. Average selected lead length (rounded) for display.
#' @param realmeanlag Numeric. Exact average selected lag length.
#' @param realmeanlead Numeric. Exact average selected lead length.
#' @param auto Logical/Integer. If non-zero, the output header reports
#'   AIC/BIC-selected lag/lead averages.
#' @param westerlund Logical. If \code{TRUE}, uses hard-coded Westerlund-specific
#'   moment constants for standardization.
#' @param aic Logical. If \code{TRUE}, AIC was used for lag/lead selection; If
#'   \code{FALSE}, BIC was used.
#' @param verbose Logical. If \code{TRUE}, prints additional output.
#'
#' @details
#' \code{DisplayWesterlund} converts raw statistics into a standardized normal
#' distribution. The choice of asymptotic moments depends on the deterministic
#' specification (none, constant, or constant + trend) and the number of
#' regressors.
#'
#' \strong{Standardization Logic:}
#' The Z-statistics are calculated as:
#' \deqn{Z = \frac{\sqrt{N}(S - \mu)}{\sqrt{\sigma^2}}}
#' where \eqn{S} is the raw statistic and \eqn{\mu, \sigma^2} are the asymptotic
#' mean and variance.
#'
#'
#'
#' \strong{Bootstrap Inference:}
#' If \code{bootstats} is provided, the function calculates "robust" p-values.
#' These are empirical p-values adjusted for finite samples:
#' \deqn{p^* = \frac{r + 1}{B + 1}}
#' where \eqn{r} is the count of bootstrap replicates less than or equal to the
#' observed statistic and \eqn{B} is the total number of valid bootstrap draws.
#' This correction is intentionally retained in this R implementation; the
#' original 2010 \code{xtwest} display code uses the uncorrected proportion
#' \eqn{r/B}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{gt, ga, pt, pa}: Raw test statistics.
#'   \item \code{gt_z, ga_z, pt_z, pa_z}: Standardized Z-statistics.
#'   \item \code{gt_pval, ga_pval, pt_pval, pa_pval}: Asymptotic left-tail p-values.
#'   \item \code{gt_pvalboot, ga_pvalboot, pt_pvalboot, pa_pvalboot}: Bootstrap
#'     p-values (returned only if \code{bootstats} is provided).
#' }
#'
#' @references
#' Westerlund, J. (2007). Testing for error correction in panel data.
#' \emph{Oxford Bulletin of Economics and Statistics}, 69(6), 709--748.
#'
#' @seealso \code{\link{westerlund_test}}, \code{\link{WesterlundPlain}},
#'   \code{\link{WesterlundBootstrap}}
#'
#' @examples
#' ## Asymptotic-only display
#' stats <- c(Gt = -2.1, Ga = -9.5, Pt = -1.8, Pa = -6.2)
#' res <- DisplayWesterlund(stats, nobs = 20, nox = 1, constant = TRUE)
#'
#' ## Accessing calculated Z-scores
#' res$gt_z
#'
#' @export
DisplayWesterlund <- function(stats,
                              bootstats = NULL,
                              nobs,
                              nox,
                              constant = FALSE,
                              trend = FALSE,
                              meanlag = -1,
                              meanlead = -1,
                              realmeanlag = -1,
                              realmeanlead = -1,
                              auto = 0,
                              westerlund = FALSE,
                              aic = TRUE,
                              verbose = FALSE) {

  # Extract statistics
  gt <- stats$Gt
  ga <- stats$Ga
  pt <- stats$Pt
  pa <- stats$Pa

  # ----------------------------------------------------------------------------
  # 1. Header Display
  # ----------------------------------------------------------------------------
  if (verbose) {
    cat("\n\n")
    cat("Results for H0: no cointegration\n")

    if (nox > 1) {
      cat(sprintf("With %d series and %d covariates\n", nobs, nox))
    } else {
      cat(sprintf("With %d series and 1 covariate\n", nobs))
    }
  }


  if (auto) {
    # 1. Determine the correct label based on user settings
    # Uses AIC if either aic or westerlund is TRUE; otherwise defaults to BIC
    ic_label <- if (aic || westerlund) "AIC" else "BIC"

    # 2. Process Lag Length
    roundedrealmeanlag <- round(realmeanlag, 2)
    if (verbose) {
      cat(sprintf(
        "Average %s selected lag length: %s\n",
        ic_label,
        format(roundedrealmeanlag, nsmall = 2)
      ))
    }

    # 3. Process Lead Length
    roundedrealmeanlead <- round(realmeanlead, 2)
    if (verbose) {
      cat(sprintf(
        "Average %s selected lead length: %s\n",
        ic_label,
        format(roundedrealmeanlead, nsmall = 2)
      ))
    }
  }

  if (verbose) cat("\n")

  # ----------------------------------------------------------------------------
  # 2. Calculating Z-stats (Normalization)
  # ----------------------------------------------------------------------------

  gtnorm <- 0
  ganorm <- 0
  ptnorm <- 0
  panorm <- 0

  # Logic to determine row index for lookup matrices
  # Stata: constant + trend + 1
  # Since trend implies constant (checked in main), options are:
  # 1. No const, No trend (Row 1)
  # 2. Const, No trend (Row 2)
  # 3. Const, Trend (Row 3)
  row_idx <- (as.integer(constant) + as.integer(trend)) + 1

  if (westerlund) {
    if (!trend) {
      gtnorm <- (sqrt(nobs) * gt - sqrt(nobs) * (-1.793)) / sqrt(0.7904)
      ganorm <- (sqrt(nobs) * ga - sqrt(nobs) * (-7.2014)) / sqrt(29.3677)
      ptnorm <- (pt - sqrt(nobs) * (-1.4746)) / sqrt(1.0262)
      panorm <- (sqrt(nobs) * pa - sqrt(nobs) * (-4.3559)) / sqrt(21.0535)
    } else {
      gtnorm <- (sqrt(nobs) * gt - sqrt(nobs) * (-2.356)) / sqrt(0.6450)
      ganorm <- (sqrt(nobs) * ga - sqrt(nobs) * (-11.8978)) / sqrt(44.2471)
      ptnorm <- (pt - sqrt(nobs) * (-2.1128)) / sqrt(0.7371)
      panorm <- (sqrt(nobs) * pa - sqrt(nobs) * (-8.9536)) / sqrt(35.6802)
    }
  } else {
    # ---------------------------------------------------------
    # Lookup Matrices (3 rows, 6 columns)
    # ---------------------------------------------------------

    # GT Mean
    gtmean <- matrix(c(
      -0.9763, -1.3816, -1.7093, -1.9789, -2.1985, -2.4262,
      -1.7776, -2.0349, -2.2332, -2.4453, -2.6462, -2.8358,
      -2.3664, -2.5284, -2.7040, -2.8639, -3.0146, -3.1710
    ), nrow = 3, byrow = TRUE)

    # GA Mean
    gamean <- matrix(c(
      -3.8022, -5.8239, -7.8108, -9.8791, -11.7239, -13.8581,
      -7.1423, -9.1249, -10.9667, -12.9561, -14.9752, -17.0673,
      -12.0116, -13.6324, -15.5262, -17.3648, -19.2533, -21.2479
    ), nrow = 3, byrow = TRUE)

    # PT Mean
    ptmean <- matrix(c(
      -0.5105, -0.9370, -1.3169, -1.6167, -1.8815, -2.1256,
      -1.4476, -1.7131, -1.9206, -2.1484, -2.3730, -2.5765,
      -2.1124, -2.2876, -2.4633, -2.6275, -2.7858, -2.9537
    ), nrow = 3, byrow = TRUE)

    # PA Mean
    pamean <- matrix(c(
      -1.0263, -2.4988, -4.2699, -6.1141, -8.0317, -10.0074,
      -4.2303, -5.8650, -7.4599, -9.3057, -11.3152, -13.3180,
      -8.9326, -10.4874, -12.1672, -13.8889, -15.6815, -17.6515
    ), nrow = 3, byrow = TRUE)

    # GT Var
    gtvar <- matrix(c(
      1.0823, 1.0981, 1.0489, 1.0576, 1.0351, 1.0409,
      0.8071, 0.8481, 0.8886, 0.9119, 0.9083, 0.9236,
      0.6603, 0.7070, 0.7586, 0.8228, 0.8477, 0.8599
    ), nrow = 3, byrow = TRUE)

    # GA Var
    gavar <- matrix(c(
      20.6868, 29.9016, 39.0109, 50.5741, 58.9595, 69.5967,
      29.6336, 39.3428, 49.4880, 58.7035, 67.9499, 79.1093,
      46.2420, 53.7428, 64.5591, 74.7403, 84.7990, 94.0024
    ), nrow = 3, byrow = TRUE)

    # PT Var
    ptvar <- matrix(c(
      1.3624, 1.7657, 1.7177, 1.6051, 1.4935, 1.4244,
      0.9885, 1.0663, 1.1168, 1.1735, 1.1684, 1.1589,
      0.7649, 0.8137, 0.8857, 0.9985, 0.9918, 0.9898
    ), nrow = 3, byrow = TRUE)

    # PA Var
    pavar <- matrix(c(
      8.3827, 24.0223, 39.8827, 53.4518, 63.2406, 76.6757,
      19.7090, 31.2637, 42.9975, 57.4844, 69.4374, 81.0384,
      37.5948, 45.6890, 57.9985, 74.1258, 81.3934, 91.2392
    ), nrow = 3, byrow = TRUE)

    # Access Matrices (using nox as column index)
    # Stata: matrix[row, col]

    gtnorm <- (sqrt(nobs) * gt - sqrt(nobs) * gtmean[row_idx, nox]) / sqrt(gtvar[row_idx, nox])
    ganorm <- (sqrt(nobs) * ga - sqrt(nobs) * gamean[row_idx, nox]) / sqrt(gavar[row_idx, nox])
    ptnorm <- (pt - sqrt(nobs) * ptmean[row_idx, nox]) / sqrt(ptvar[row_idx, nox])
    panorm <- (sqrt(nobs) * pa - sqrt(nobs) * pamean[row_idx, nox]) / sqrt(pavar[row_idx, nox])
  }

  # ----------------------------------------------------------------------------
  # 3. Display Logic
  # ----------------------------------------------------------------------------

  # Return object initialization
  ret_list <- list(
    gt = gt, ga = ga, pt = pt, pa = pa,
    gt_z = gtnorm, ga_z = ganorm, pt_z = ptnorm, pa_z = panorm,
    gt_pval = pnorm(gtnorm),
    ga_pval = pnorm(ganorm),
    pt_pval = pnorm(ptnorm),
    pa_pval = pnorm(panorm)
  )

  # Check if Bootstrapping was done
  if (is.null(bootstats)) {
    if (verbose){
      # Standard Table

      cat("----------------------------------------------------\n")
      cat(" Statistic |   Value   |  Z-value  |  P-value  |\n")
      cat("----------------------------------------------------\n")

      # Helper to print row
      print_row <- function(name, val, z, p) {
        cat(sprintf("     %s    | %8.3f  | %8.3f  | %8.3f  |\n",
                    name, val, z, round(p, 4)))
      }

      print_row("Gt", gt, gtnorm, pnorm(gtnorm))
      print_row("Ga", ga, ganorm, pnorm(ganorm))
      print_row("Pt", pt, ptnorm, pnorm(ptnorm))
      print_row("Pa", pa, panorm, pnorm(panorm))

      cat("----------------------------------------------------\n")
    }

  } else {
    # Bootstrapped Table

    # bootstats is expected to be a matrix (rows=reps, cols=4 [Gt, Ga, Pt, Pa])
    rows <- nrow(bootstats)

    # Extract columns
    GTSTATS <- bootstats[, 1]
    GASTATS <- bootstats[, 2]
    PTSTATS <- bootstats[, 3]
    PASTATS <- bootstats[, 4]

    # Sorting and P-value calculation
    # Logic: pvalue = (count of boot_stats <= obs_stat) / rows
    # Stata implementation: Sort bootstats ascending. Loop i=1..rows.
    # if (boot[i] - obs > 0) break. position = i-1. pval = position/rows.
    # This implies pvalue = proportion of boot stats strictly less than or equal to observed.

    calc_boot_pval <- function(boot_vec, obs_val) {
      boot_vec <- boot_vec[is.finite(boot_vec)]   # drop NA/NaN/Inf like Stata implicitly does
      B <- length(boot_vec)
      if (B == 0) return(NA_real_)

      # Stata left-tail: count of boot <= observed
      r <- sum(boot_vec <= obs_val)

      # Finite-sample correction
      (r + 1) / (B + 1)
    }

    pvaluegtboot <- calc_boot_pval(GTSTATS, gt)
    pvaluegaboot <- calc_boot_pval(GASTATS, ga)
    pvalueptboot <- calc_boot_pval(PTSTATS, pt)
    pvaluepaboot <- calc_boot_pval(PASTATS, pa)

    # Display Table with Robust P-value
    if (verbose) {
      cat("--------------------------------------------------------------------------\n")
      cat(" Statistic |   Value   |  Z-value  |  P-value  | Robust P-value |\n")
      cat("--------------------------------------------------------------------------\n")

      print_boot_row <- function(name, val, z, p, bp) {
        cat(sprintf("     %s    | %8.3f  | %8.3f  | %8.3f  |      %8.3f    |\n",
                    name, val, z, round(p, 4), round(bp, 4)))
      }

      print_boot_row("Gt", gt, gtnorm, pnorm(gtnorm), pvaluegtboot)
      print_boot_row("Ga", ga, ganorm, pnorm(ganorm), pvaluegaboot)
      print_boot_row("Pt", pt, ptnorm, pnorm(ptnorm), pvalueptboot)
      print_boot_row("Pa", pa, panorm, pnorm(panorm), pvaluepaboot)

      cat("--------------------------------------------------------------------------\n")
    }

    # Add boot p-values to return list
    ret_list$gt_pvalboot <- pvaluegtboot
    ret_list$ga_pvalboot <- pvaluegaboot
    ret_list$pt_pvalboot <- pvalueptboot
    ret_list$pa_pvalboot <- pvaluepaboot
  }

  # Return the list of scalars
  return(ret_list)
}

# ==============================================================================
# Main Function: westerlund_test_reg
# ==============================================================================
#' Formatted Coefficient Table for Westerlund Test Reporting
#'
#' Internal helper that calculates standard errors, z-statistics, p-values, and
#' confidence intervals for a given coefficient vector and covariance matrix,
#' then prints a Stata-style regression table.
#'
#' @param b A named numeric vector of coefficients.
#' @param V A numeric variance-covariance matrix corresponding to \code{b}.
#' @param verbose Logical. If \code{TRUE}, prints additional output.
#' @param distribution Character, either \code{"z"} or \code{"t"}. Mean-group
#'   reporting uses Z inference; individual OLS ECM output uses t inference.
#' @param df Residual degrees of freedom, required when \code{distribution="t"}.
#'
#' @details
#' \strong{Calculation Logic.}
#' The function supports the two reporting conventions used by the package: Z inference
#' for mean-group \code{ereturn post}-style tables and t inference for individual OLS ECMs.
#' Standard errors are derived from \code{V}; two-sided p-values and 95\% confidence
#' intervals use the selected reference distribution.
#'
#' \strong{Formatting.}
#' The table is printed using \code{stats::printCoefmat()} to ensure clean
#' alignment and decimal consistency. It includes columns for the Coefficient,
#' Standard Error, selected test statistic, P-value, and the 95\% Confidence Interval.
#'
#' \strong{Intended use.}
#' This is an internal utility called by \code{\link{westerlund_test_mg}}. It
#' provides a standardized way to report results across different parts of the
#' Westerlund cointegration test output.
#'
#' @return A numeric matrix with coefficients, standard errors, the selected
#' test statistic, p-values, and 95\% confidence interval bounds.
#'
#' @section Reporting Style:
#' This section describes the alignment with econometric software reporting standards.
#'
#' \subsection{Stata Compatibility}{
#' The output format (column headers and statistical assumptions) is designed to
#' match the output seen in Stata's \code{xtcpmg} or \code{xtwest} routines. This
#' ensures that users transitioning from or comparing results with Stata find the
#' output familiar.
#' }
#'
#' @examples
#' ## Define simple coefficient vector and VCV
#' b <- c(alpha = -0.20, beta = 1.05)
#' V <- diag(c(0.03^2, 0.10^2))
#' names(b) <- rownames(V) <- colnames(V) <- c("alpha", "beta")
#'
#' ## Print the formatted table
#' westerlund_test_reg(b, V, verbose = TRUE)
#'
#' @seealso \code{\link{westerlund_test_mg}}, \code{\link{westerlund_test}}
#' @export
westerlund_test_reg <- function(b, V, verbose = FALSE,
                                distribution = c("z", "t"), df = NULL) {
  distribution <- match.arg(distribution)
  if (length(b) != nrow(V) || length(b) != ncol(V))
    stop("Dimensions of coefficient vector and variance matrix do not match.")
  if (distribution == "t" && (is.null(df) || length(df) != 1 || !is.finite(df) || df <= 0))
    stop("A positive residual df must be supplied when distribution='t'.")

  coefs <- as.numeric(b); names(coefs) <- names(b)
  se <- sqrt(diag(V))
  stat_vals <- coefs / se
  if (distribution == "t") {
    p_vals <- 2 * stats::pt(-abs(stat_vals), df = df)
    crit <- stats::qt(0.975, df = df)
    stat_name <- "t"; p_name <- "P>|t|"
  } else {
    p_vals <- 2 * stats::pnorm(-abs(stat_vals))
    crit <- stats::qnorm(0.975)
    stat_name <- "z"; p_name <- "P>|z|"
  }
  ci_lower <- coefs - crit * se
  ci_upper <- coefs + crit * se
  res_table <- cbind("Coef." = coefs, "Std. Err." = se, stat_vals, p_vals,
                     "[95% Conf." = ci_lower, "Interval]" = ci_upper)
  colnames(res_table)[3:4] <- c(stat_name, p_name)
  rownames(res_table) <- names(b)
  if (verbose) {
    cat("\n")
    printCoefmat(res_table, digits = 4, signif.stars = FALSE, na.print = "NA",
                 cs.ind = 1:2, tst.ind = 3, zap.ind = 4,
                 P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  }
  return(res_table)
}

# ==============================================================================
# Main Function: westerlund_test_mg
# ==============================================================================
#' Print Mean-Group ECM Output and Long-Run Relationship Tables
#'
#' Internal reporting helper that prints two Stata-style coefficient tables:
#' (1) the mean-group error-correction model (short-run ECM output) and
#' (2) the estimated long-run relationship with short-run adjustment.
#' This function delegates all table computation and formatting to
#' \code{\link{westerlund_test_reg}}.
#'
#' @param b A numeric vector of coefficients for the long-run relationship and
#'   short-run adjustment terms (as reported in the second table).
#' @param V A numeric variance-covariance matrix corresponding to \code{b}.
#' @param b2 A numeric vector of coefficients for the mean-group error-correction
#'   model (as reported in the first table).
#' @param V2 A numeric variance-covariance matrix corresponding to \code{b2}.
#' @param auto Logical/integer. If non-zero, prints a note indicating that
#'   short-run coefficients (apart from the error-correction term) are omitted
#'   because selected lag/lead lengths may differ across units.
#' @param verbose Logical. If \code{TRUE}, prints additional output.
#'
#' @details
#' This function is designed to replicate the reporting structure of mean-group
#' estimators in Stata. It organizes the output into two distinct sections:
#'
#' \strong{1. Mean-group error-correction model:}
#' This table displays the coefficients related to the short-run dynamics and the
#' error-correction term. When \code{auto} is enabled (suggesting unit-specific
#' lag/lead selection), the function informs the user that non-ec short-run
#' coefficients are not displayed because they are not uniform across the panel.
#'
#' \strong{2. Estimated long-run relationship and short-run adjustment:}
#' This table displays the estimated cointegrating vector (long-run relationship)
#' and the speed of adjustment back to equilibrium.
#'
#' Both sections use \code{\link{westerlund_test_reg}} to generate the statistical
#' inference (standard errors, z-stats, p-values, and 95% confidence intervals)
#' and format the results using \code{stats::printCoefmat()}.
#'
#'
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{mg_model}: A numeric matrix containing the Mean-group ECM results.
#'   \item \code{long_run}: A numeric matrix containing the Long-run relationship and adjustment results.
#' }
#'
#' @section Intended Use:
#' This is an internal display helper used in reporting flows. It is not
#' required for the calculation of the Westerlund test statistics but provides
#' the descriptive regression context for the mean-group estimator.
#'
#' @examples
#' ## Long-run relationship table inputs
#' b <- c(alpha = -0.20, beta = 1.05)
#' V <- diag(c(0.03^2, 0.10^2))
#' rownames(V) <- colnames(V) <- names(b)
#'
#' ## Mean-group ECM table inputs
#' b2 <- c(ec = -0.18)
#' V2 <- matrix(0.04^2, nrow = 1, ncol = 1)
#' rownames(V2) <- colnames(V2) <- names(b2)
#'
#' ## Print both sections (auto = TRUE prints the omission note)
#' westerlund_test_mg(b, V, b2, V2, auto = TRUE, verbose = FALSE)
#'
#' @seealso \code{\link{westerlund_test_reg}}, \code{\link{westerlund_test}}
#'
#' @export
westerlund_test_mg <- function(b, V, b2, V2, auto, verbose = FALSE) {

  # ----------------------------------------------------------------------------
  # 1. Display Second Set (Short Run / Mean Group ECM)
  # Corresponds to: eret post `b2' `V2' -> eret disp
  # ----------------------------------------------------------------------------

  if (verbose) {
    cat("\n\n")
    cat("Mean-group error-correction model\n")

    if (auto) {
      cat("Short run coefficients apart from the error-correction term are omitted as lag and lengths might differ between cross-sectional units\n")
      cat("\n")
    }
  }

  # Call the previously defined display function
  # Ensure b2 is vector and V2 is matrix with matching names
  mg_model <- westerlund_test_reg(b2, V2, verbose)

  # ----------------------------------------------------------------------------
  # 2. Display First Set (Long Run Relationship)
  # Corresponds to: eret post `b' `V' -> eret disp
  # ----------------------------------------------------------------------------

  if (verbose) {
    cat("\n")
    cat("Estimated long-run relationship and short run adjustment\n")
  }

  # Call the display function
  lr = westerlund_test_reg(b, V, verbose)

  return(list(mg_model = mg_model, long_run = lr))
}

# ============================================================================== 
# S3 Method: plot.westerlund_test
# ==============================================================================
#' Plot Bootstrap Distributions for Westerlund ECM Panel Cointegration Tests
#'
#' @param x An object of class \code{westerlund_test}.
#' @param title Main plot title.
#' @param conf_level Lower-tail bootstrap critical probability. Defaults to 0.05.
#' @param save_path Optional file path. If supplied, the plot is saved with \code{ggsave}.
#' @param dpi Resolution used when \code{save_path} is supplied.
#' @param figsize Numeric length-2 vector giving width and height in inches.
#' @param colors Named list with \code{obs}, \code{crit}, \code{fill}, and \code{density}.
#' @param lwd Named list with line widths for \code{obs}, \code{crit}, and \code{density}.
#' @param alpha Density fill transparency.
#' @param show_grid Logical; show major panel grid lines.
#' @param show_robust_p Logical; annotate the bootstrap p-value in each facet.
#' @param ... Additional arguments reserved for compatibility.
#' @return A \code{ggplot} object, invisibly saved as well when \code{save_path} is supplied.
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom scales percent
#' @importFrom stats quantile aggregate
#' @importFrom dplyr everything
#' @method plot westerlund_test
#' @export
plot.westerlund_test <- function(x,
                                 title = "Westerlund Test: Bootstrap Distributions",
                                 conf_level = 0.05,
                                 save_path = NULL,
                                 dpi = 300,
                                 figsize = c(12, 10),
                                 colors = list(obs = "#D55E00", crit = "#0072B2",
                                               fill = "grey80", density = "grey30"),
                                 lwd = list(obs = 1, crit = 0.8, density = 0.5),
                                 alpha = 0.5,
                                 show_grid = TRUE,
                                 show_robust_p = TRUE, ...) {
  if (is.null(x$bootstrap_distributions))
    stop("No bootstrap results found. Run westerlund_test() with bootstrap > 0 first.")
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be between 0 and 1.")
  if (length(figsize) != 2) stop("figsize must contain width and height.")

  boot_df <- as.data.frame(x$bootstrap_distributions)
  colnames(boot_df) <- c("Gt", "Ga", "Pt", "Pa")
  plot_data <- tidyr::pivot_longer(boot_df, cols = dplyr::everything(),
                                   names_to = "Statistic", values_to = "Value")
  plot_data <- plot_data[is.finite(plot_data$Value), ]
  obs_df <- data.frame(Statistic = c("Gt", "Ga", "Pt", "Pa"),
                       Observed = unlist(x$test_stats[c("Gt", "Ga", "Pt", "Pa")]),
                       stringsAsFactors = FALSE)
  crit_df <- stats::aggregate(Value ~ Statistic, data = plot_data,
                              FUN = function(v) stats::quantile(v, conf_level, na.rm = TRUE))
  names(crit_df)[2] <- "CriticalValue"
  ref_lines <- merge(obs_df, crit_df, by = "Statistic", all.x = TRUE)
  ref_lines$RobustP <- NA_real_
  if (!is.null(x$boot_pvals) && length(x$boot_pvals)) {
    bp <- unlist(x$boot_pvals)
    ref_lines$RobustP <- bp[match(ref_lines$Statistic, names(bp))]
  }
  ref_lines$Label <- paste0("Obs: ", sprintf("%.3f", ref_lines$Observed),
                            "\nCV: ", sprintf("%.3f", ref_lines$CriticalValue))
  if (show_robust_p)
    ref_lines$Label <- ifelse(is.finite(ref_lines$RobustP),
                              paste0(ref_lines$Label, "\nRobust p: ", sprintf("%.3f", ref_lines$RobustP)),
                              ref_lines$Label)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Value)) +
    ggplot2::geom_density(fill = colors$fill, color = colors$density,
                          alpha = alpha, linewidth = lwd$density) +
    ggplot2::geom_vline(data = ref_lines,
                        ggplot2::aes(xintercept = .data$Observed, color = "Observed"),
                        linetype = "solid", linewidth = lwd$obs) +
    ggplot2::geom_vline(data = ref_lines,
                        ggplot2::aes(xintercept = .data$CriticalValue, color = "Critical"),
                        linetype = "dashed", linewidth = lwd$crit) +
    ggplot2::geom_text(data = ref_lines,
                       ggplot2::aes(x = -Inf, y = Inf, label = .data$Label),
                       hjust = -0.1, vjust = 1.4, size = 3,
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~Statistic, scales = "free", ncol = 2) +
    ggplot2::scale_color_manual(name = NULL,
                                values = c("Observed" = colors$obs, "Critical" = colors$crit),
                                labels = c("Observed Statistic",
                                           paste0(scales::percent(conf_level), " Bootstrap CV"))) +
    ggplot2::labs(title = title,
                  subtitle = paste("H0: No Cointegration | Replications:", nrow(boot_df)),
                  x = "Statistic Value", y = "Kernel Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  if (!show_grid) p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  if (!is.null(save_path)) {
    out_dir <- dirname(save_path)
    if (!identical(out_dir, ".") && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = save_path, plot = p, width = figsize[1], height = figsize[2], dpi = dpi)
  }
  return(p)
}

#' Print Method for westerlund_test
#' @param x An object of class \code{westerlund_test}.
#' @param ... Unused.
#' @method print westerlund_test
#' @export
print.westerlund_test <- function(x, ...) {
  cat("\n--- Westerlund (2007) Panel Cointegration Test ---\n")
  tab <- data.frame(
    Statistic = c("Gt", "Ga", "Pt", "Pa"),
    Value = unlist(x$test_stats[c("Gt", "Ga", "Pt", "Pa")]),
    Z_score = unlist(x$z_scores[c("Gt", "Ga", "Pt", "Pa")]),
    P_val_asymp = unlist(x$p_values[c("Gt", "Ga", "Pt", "Pa")]),
    row.names = NULL
  )
  if (!is.null(x$boot_pvals) && length(x$boot_pvals))
    tab$P_val_boot <- unlist(x$boot_pvals[tab$Statistic])
  print(tab, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Summary Method for Westerlund Test
#' @param object An object of class \code{westerlund_test}.
#' @param ... Additional arguments.
#' @method summary westerlund_test
#' @export
summary.westerlund_test <- function(object, ...) {
  stats_table <- data.frame(
    Statistic = c("Gt", "Ga", "Pt", "Pa"),
    Value = unlist(object$test_stats[c("Gt", "Ga", "Pt", "Pa")]),
    Z_score = unlist(object$z_scores[c("Gt", "Ga", "Pt", "Pa")]),
    P_val_asymp = unlist(object$p_values[c("Gt", "Ga", "Pt", "Pa")]),
    stringsAsFactors = FALSE, row.names = NULL
  )
  if (!is.null(object$boot_pvals) && length(object$boot_pvals))
    stats_table$P_val_boot <- unlist(object$boot_pvals[stats_table$Statistic])
  out <- list(stats_table = stats_table,
              settings = object$settings,
              mean_group = object$mean_group,
              mg_results = object$mg_results,
              mg_tables = object$mg_tables,
              n_units = object$metadata$n_groups)
  class(out) <- "summary.westerlund_test"
  out
}

#' Print Summary Method
#' @param x An object of class \code{summary.westerlund_test}.
#' @param ... Unused.
#' @method print summary.westerlund_test
#' @export
print.summary.westerlund_test <- function(x, ...) {
  cat("\n======================================================================\n")
  cat("Westerlund (2007) Panel Cointegration Test Summary\n")
  cat("======================================================================\n")
  cat(sprintf("Units: %s | Average time periods: %.2f\n", x$n_units, x$settings$T))
  det <- if (isTRUE(x$settings$trend)) "trend" else if (isTRUE(x$settings$constant)) "constant" else "none"
  cat("Deterministic terms:", det, "\n")
  cat("Lag/lead selection:", x$settings$selection_method, "\n\n")
  print(x$stats_table, row.names = FALSE, digits = 4)
  cat("\nMean Group Estimates:\n")
  print(x$mean_group)
  invisible(x)
}

