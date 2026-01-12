# Westerlund: Westerlund ECM Panel Cointegration Test (Python & R)

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository provides both a Python and R implementation of a functional approximation of the Error Correction Model (ECM) based panel cointegration tests proposed by **Westerlund (2007)**. 

---

## 📖 Overview

The Westerlund test evaluates the null hypothesis of **no cointegration** by testing whether the error-correction term in a conditional panel ECM is equal to zero. If the null is rejected, there is evidence of a long-run equilibrium relationship between the variables.



### Key Features
The package replicates the logic of the Westerlund (2007) methodology, including:
* **Four Test Statistics**: Computes $G_t$, $G_a$, $P_t$, and $P_a$.
* **Flexible Dynamics**: Allows for unit-specific lag and lead lengths.
* **Automated Selection**: Built-in AIC/BIC selection logic for optimal lag and lead lengths. (BIC to be implemented)
* **Bootstrap Procedure**: Robust p-values to handle cross-sectional dependence.
* **Kernel Estimation**: Bartlett kernel long-run variance estimation.
* **Gap Handling**: Strict time-series continuity checks to ensure valid econometric results.

---

## Project Repository Structure
```text
Westerlund/
├── python/                   # Python Package Root
│   ├── LICENSE               # Python-specific license file
│   ├── pyproject.toml        # Build configuration
│   └── src/                  # Source Code
│       └── westerlund/           
│           ├── __init__.py       # API promotion (exporting WesterlundTest)
│           └── main.py           # WesterlundTest Class Implementation
├── R/                        # R Package Root
│   ├── DESCRIPTION           # Package Metadata, Dependencies, and MIT License ref
│   ├── NAMESPACE             # Exported functions and imported R packages
│   ├── LICENSE               # R-specific license file
│   ├── R/                    # R Source Code
│   │   └── Westerlund.R      # Main logic and R6/S3 class definitions
│   ├── man/                  # Documentation (Help files)
│   └── vignettes/            # Long-form Documentation
├── LICENSE                   # Root MIT License (Full Text)
└── README.md                 # Integrated Monorepo README
```

## References
Westerlund, J. (2007). Testing for Error Correction in Panel Data. Oxford Bulletin of Economics and Statistics, 69(6), 709-748.

Persyn, D., & Westerlund, J. (2008). Error-Correction-Based Cointegration Tests for Panel Data. Stata Journal, 8(2), 232-241.
