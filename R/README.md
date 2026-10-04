# Westerlund: Panel Cointegration Testing in R

`Westerlund` is an R package implementing the four error-correction-based panel cointegration tests developed by **Westerlund (2007)**. The test evaluates the null hypothesis of **no cointegration** by testing whether the error-correction term in a conditional panel ECM is equal to zero. Rejection of the null provides evidence of a long-run equilibrium relationship between the variables.

The implementation is designed to reproduce the logic of the Stata `xtwest` command while providing an R-native interface, reusable result objects, plotting, and summary methods.

## Key Features

The package provides:

* **Four Test Statistics**: Computes $G_t$, $G_a$, $P_t$, and $P_a$.
* **Flexible Dynamics**: Supports fixed or unit-specific lag and lead lengths.
* **Automated Selection**: Built-in AIC/BIC selection and the Westerlund-specific information criterion.
* **Bootstrap Inference**: Residual-based bootstrap inference with synchronized time-cluster resampling to preserve contemporaneous cross-sectional dependence.
* **Reproducible Bootstrap Runs**: Optional `seed` argument.
* **Unbalanced Panels**: Supports unequal panel lengths while requiring continuous time within each unit.
* **Kernel Estimation**: Bartlett-kernel long-run variance estimation.
* **Strict Time-Series Handling**: Gap-aware lags, leads, and differences.
* **Individual ECM Output**: Optional unit-level regression output with `indiv.ecm = TRUE`.
* **Rich Result Objects**: Returns raw statistics, standardized Z-scores, asymptotic and bootstrap p-values, unit-level estimates, mean-group estimates, settings, and bootstrap distributions.
* **S3 Methods**: Includes `print()`, `summary()`, and `plot()` methods for fitted results.
* **Bootstrap Visualization**: Configurable critical level, robust-p annotations, and optional direct file export.

## Installation

You can install the development version of `Westerlund` from GitHub using the `devtools` package:

```r
# Install devtools if you haven't already
# install.packages("devtools")

# Install Westerlund
devtools::install_github("bosco-hung/Westerlund/R", build_vignettes = TRUE)
```

## Quick Start

```r
library(Westerlund)

results <- westerlund_test(
  data = my_data,
  yvar = "ln_gdp",
  xvars = c("ln_energy", "ln_capital"),
  idvar = "country",
  timevar = "year",
  constant = TRUE,
  trend = FALSE,
  lags = c(1, 3),      # AIC search between 1 and 3 lags
  leads = c(0, 2),     # AIC search between 0 and 2 leads
  bootstrap = 100,     # Run 100 bootstrap replications
  lrwindow = 3,        # Bartlett kernel window size
  seed = 123           # Reproducible bootstrap
)

# Main results
print(results)

# Extended summary
summary(results)

# Bootstrap distributions and critical values
plot(
  results,
  conf_level = 0.05,
  show_robust_p = TRUE
)
```

If no lag specification is supplied, `lags = 1` is used. If `leads = NULL`, the lead order defaults to zero.

For individual-unit ECM output, set:

```r
results <- westerlund_test(
  data = my_data,
  yvar = "ln_gdp",
  xvars = "ln_energy",
  idvar = "country",
  timevar = "year",
  constant = TRUE,
  indiv.ecm = TRUE
)

results$unit_data
results$indiv_reg
```

## References

Westerlund, J. (2007). Testing for Error Correction in Panel Data. *Oxford Bulletin of Economics and Statistics*, 69(6), 709-748.

Persyn, D., & Westerlund, J. (2008). Error-Correction-Based Cointegration Tests for Panel Data. *Stata Journal*, 8(2), 232-241.
