# Westerlund: Panel Cointegration Testing in Python

`Westerlund` is a Python implementation of the four error-correction-based panel cointegration tests developed by **Westerlund (2007)**. The test evaluates the null hypothesis of **no cointegration** by testing whether the error-correction term in a conditional panel ECM is equal to zero. Rejection of the null provides evidence of a long-run equilibrium relationship among the variables.

The implementation is designed to follow the logic of the Stata `xtwest` command closely while providing a Python-native interface, reusable result objects, summary methods, optional individual ECM output, and bootstrap visualization. The Python and R versions are maintained with aligned user-facing functionality.

## Key Features

The package provides:

* **Four Test Statistics**: Computes $G_t$, $G_a$, $P_t$, and $P_a$.
* **Flexible Dynamics**: Supports fixed or unit-specific lag and lead lengths.
* **Automated Selection**: Supports AIC/BIC lag/lead selection and the Westerlund-specific information criterion.
* **Stata-Aligned Bootstrap**: Uses the restricted null model, synchronized actual-time cluster resampling, unit-specific panel lengths, and the original `xtwest` bootstrap ordering logic.
* **Cross-Sectional Dependence**: Bootstrap time clusters are resampled jointly across units to preserve contemporaneous dependence.
* **Unbalanced Panels**: Unequal panel lengths are supported, subject to continuous time within each unit.
* **Reproducible Bootstrap Runs**: Optional `seed` argument.
* **Kernel Estimation**: Bartlett-kernel long-run variance estimation.
* **Strict Time-Series Handling**: Gap-aware lags, leads, and differences.
* **Input Validation**: Enforces the six-regressor asymptotic-table limit, `trend=True` → `constant=True`, and the restrictions associated with `westerlund=True`.
* **Individual ECM Output**: Optional unit-specific ECM regression tables through `indiv_ecm=True`.
* **Mean-Group Results**: Returns mean-group alpha and beta estimates together with their standard errors.
* **Rich Result Bundle**: Returns raw statistics, standardized Z-scores, asymptotic p-values, bootstrap p-values and distributions, unit-level results, model settings, metadata, and reporting tables.
* **Summary Interface**: Provides `summary()` and `print_summary()`.
* **Bootstrap Visualization**: Configurable critical level, robust-p annotations, line widths, colors, grid display, figure size, and direct file export.
* **R/Python Feature Parity**: The public functionality is kept aligned with the companion R implementation.

## Installation

Install the package from PyPI:

```bash
pip install Westerlund
```

## Quick Start

```python
import pandas as pd
from westerlund_test import WesterlundTest

# Panel data in long format:
# ID | Time | Y | X1 | X2 | ...
df = pd.read_csv("your_data.csv")

test = WesterlundTest(
    data=df,
    y_var="log_gdp",
    x_vars=["log_energy", "log_capital"],
    id_var="country_id",
    time_var="year",
    constant=True,
    trend=False,
    lags=(0, 2),        # Select unit-specific lags between 0 and 2
    leads=(0, 1),       # Select unit-specific leads between 0 and 1
    bootstrap=100,      # Number of bootstrap replications
    lrwindow=2,         # Bartlett-kernel window
    seed=42             # Reproducible bootstrap
)

results = test.run()
```

If `lags` is not supplied, it defaults to `1`. `leads` defaults to `0`.

## Inspecting Results

`run()` returns a dictionary containing the main test results and post-estimation information:

```python
results["test_stats"]          # Gt, Ga, Pt, Pa
results["z_scores"]            # Standardized Z statistics
results["p_values"]            # Asymptotic p-values
results["boot_pvals"]          # Bootstrap p-values, if requested
results["boot_distributions"]  # Bootstrap draws
results["unit_data"]           # Unit-level estimates and selected dynamics
results["indiv_data"]          # Detailed unit-level internal results
results["mean_group"]          # MG alpha/betas and standard errors
results["settings"]            # Model and selection settings
results["metadata"]            # Group count, average T, model type, etc.
results["mg_results"]          # Mean-group reporting results
results["mg_tables"]           # Formatted MG reporting tables
results["indiv_reg"]           # Individual ECM tables when requested
```

A compact programmatic summary is available with:

```python
summary = test.summary()
print(summary["stats_table"])
```

For a formatted console summary:

```python
test.print_summary()
```

Printing the fitted object also displays the main results:

```python
print(test)
```

## Individual ECM Results

To retain unit-specific ECM regression tables:

```python
test = WesterlundTest(
    data=df,
    y_var="log_gdp",
    x_vars=["log_energy"],
    id_var="country_id",
    time_var="year",
    constant=True,
    indiv_ecm=True
)

results = test.run()

unit_results = results["unit_data"]
individual_regressions = results["indiv_reg"]
```

Individual OLS regression tables use residual degrees of freedom and t inference. Mean-group reporting uses large-sample normal (z) inference.

## Bootstrap Inference

When `bootstrap > 0`, the implementation constructs the bootstrap distribution under the null of no error correction. In broad terms it:

1. estimates the restricted short-run model under $H_0$ for each unit;
2. selects unit-specific lag/lead orders from that restricted model when ranges are supplied;
3. centers residuals and differenced regressors within unit;
4. resamples **actual time clusters jointly across units** to preserve contemporaneous cross-sectional dependence;
5. reproduces the Stata-style duplication and common ordering procedure used by `xtwest`;
6. preserves each unit's original panel length $T_i$;
7. recursively generates bootstrap differences and integrates them to levels; and
8. reruns the Westerlund test on each simulated panel.

For compatibility with the original 2010 `xtwest` source, the bootstrap also reproduces its `currlead` behavior when automatically selected lead orders vary across units.

The package reports the finite-sample-corrected bootstrap p-value

\[
p^* = \frac{r+1}{B+1},
\]

where $r$ is the number of finite bootstrap statistics less than or equal to the observed statistic and $B$ is the number of finite bootstrap draws.

## Bootstrap Visualization

After running a bootstrap:

```python
fig = test.plot(
    conf_level=0.05,
    show_robust_p=True,
    figsize=(12, 10)
)
```

`plot()` is an alias of `plot_bootstrap()`. The critical probability is configurable:

```python
test.plot(conf_level=0.10)
```

The plot can also be saved directly:

```python
test.plot(
    conf_level=0.05,
    save_path="figures/westerlund_bootstrap.png",
    dpi=300
)
```

Additional options include `colors`, `lwd`, `alpha`, `show_grid`, `show_robust_p`, and `show`.

## Input Requirements and Validation

The data must contain the dependent variable, all regressors, the panel identifier, and the time identifier.

The implementation requires:

* a maximum of six long-run regressors;
* `constant=True` whenever `trend=True`;
* when `westerlund=True`, at least a constant and at most one regressor;
* non-negative integer lag, lead, and long-run-variance window settings; and
* a continuous time index within each unit after rows with missing model variables are excluded.

Panels may be unbalanced; different units can have different values of $T_i$.

## Reproducibility Across Python, R, and Stata

Supplying `seed` makes repeated Python bootstrap runs reproducible. The R implementation similarly provides a seed option.

However, the random-number generators and numerical linear-algebra implementations differ across Python, R, and Stata. Therefore, using the same integer seed does **not** imply replication-by-replication identical bootstrap draws across software. Cross-software validation should focus on the estimation/bootstrap algorithm and the resulting distribution rather than on matching each random draw exactly.

Non-bootstrap results should agree up to numerical precision when the same data and specification are used.

## References

Westerlund, J. (2007). Testing for Error Correction in Panel Data. *Oxford Bulletin of Economics and Statistics*, 69(6), 709-748.

Persyn, D., & Westerlund, J. (2008). Error-Correction-Based Cointegration Tests for Panel Data. *Stata Journal*, 8(2), 232-241.
