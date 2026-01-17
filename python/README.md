# Westerlund: Panel Cointegration Testing in Python

`Westerlund` is a Python package implementing a functional approximation of the four panel cointegration tests developed by **Westerlund (2007)**. The Westerlund test evaluates the null hypothesis of **no cointegration** by testing whether the error-correction term in a conditional panel ECM is equal to zero. If the null is rejected, there is evidence of a long-run equilibrium relationship between the variables.

## Key Features

The package replicates the logic of the Westerlund (2007) methodology, including:
* **Four Test Statistics**: Computes $G_t$, $G_a$, $P_t$, and $P_a$.
* **Flexible Dynamics**: Allows for unit-specific lag and lead lengths.
* **Automated Selection**: Built-in AIC/BIC selection logic for optimal lag and lead lengths.
* **Bootstrap Procedure**: Robust p-values to handle cross-sectional dependence.
* **Kernel Estimation**: Bartlett kernel long-run variance estimation.
* **Gap Handling**: Strict time-series continuity checks to ensure valid econometric results.



## Usage Examples
```python
import pandas as pd
from westerlund_test import WesterlundTest

# 1. Prepare your panel data (Long format)
# Required columns: ID, Time, Y, X1, X2...
df = pd.read_csv("your_data.csv")

# 2. Initialize the test
test = WesterlundTest(
    data=df, 
    y_var='log_gdp', 
    x_vars=['log_energy', 'log_capital'], 
    id_var='country_id', 
    time_var='year',
    lags=(0, 2),        # Auto-select lags between 0 and 2
    leads=(0, 1),       # Auto-select leads between 0 and 1
    constant=True,      # Include intercept
    trend=True,         # Include time trend
    bootstrap=100,      # Perform 100 bootstrap replications
    seed=42
)

# 3. Run the estimation
results = test.run()

# 4. Visualize the results
test.plot_bootstrap()
```

## References
Westerlund, J. (2007). Testing for Error Correction in Panel Data. Oxford Bulletin of Economics and Statistics, 69(6), 709-748.

Persyn, D., & Westerlund, J. (2008). Error-Correction-Based Cointegration Tests for Panel Data. Stata Journal, 8(2), 232-241.