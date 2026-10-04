import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.stats import gaussian_kde
from numpy.random import Generator, MT19937

import matplotlib.pyplot as plt
import seaborn as sns

import os

class WesterlundTest:
    """
    Implements the Westerlund (2007) ECM-based panel cointegration tests.
    
    This class provides a comprehensive framework for testing the null hypothesis of 
    no cointegration in panel data by evaluating the presence of error correction. 
    It computes four main statistics: two mean-group tests (Gt, Ga) and two 
    pooled/panel tests (Pt, Pa). 
    
    The implementation is designed to replicate the logic and strict time-series 
    handling of the Stata 'xtwest' command, including support for unbalanced panels, 
    automatic lag/lead selection, and residual-based bootstrapping.
    """

    def __init__(self, data, y_var, x_vars, id_var, time_var, 
                 lags=1, leads=0, lrwindow=2, constant=False, 
                 trend=False, westerlund=False, aic=True,
                 bootstrap=-1, indiv_ecm=False, seed=None, 
                 verbose=False):
        """
        Initializes the WesterlundTest class with data and model specifications.

        Parameters
        ----------
        data : pandas.DataFrame
            The dataset containing the variables for the panel cointegration test.
        y_var : str
            The name of the dependent variable (levels).
        x_vars : str or list of str
            The name(s) of the regressors entering the long-run relationship (levels).
        id_var : str
            The column name identifying the cross-sectional units (groups).
        time_var : str
            The column name identifying the time index.
        lags : int or tuple of (int, int), optional
            Short-run lag length (p). If a tuple is provided, the model performs 
            automatic selection within the range (min, max). Defaults to 1.
        leads : int, tuple of (int, int), or None, optional
            Short-run lead length (q). If a tuple is provided, the model performs 
            automatic selection within the range (min, max). If None, defaults to 0.
        lrwindow : int, optional
            The bandwidth parameter for the Bartlett kernel used in long-run 
            variance estimation. Defaults to 2.
        constant : bool, optional
            Whether to include an intercept in the cointegrating relationship. 
            Defaults to False.
        trend : bool, optional
            Whether to include a linear time trend. If True, `constant=True` is
            required. Defaults to False.
        aic : bool, optional
            Whether to use Akaike Information Criterion for lag/lead selection. If
            False, uses Bayesian Information Criterion. Defaults to True.
        westerlund : bool, optional
            If True, enforces specific constraints from the 2007 paper, including 
            specialized information criteria and trimming logic for variance. 
            Defaults to False.
        bootstrap : int, optional
            Number of bootstrap replications. If > 0, robust p-values are 
            calculated. Defaults to -1 (no bootstrap).
        indiv_ecm : bool, optional
            Whether to store unit-specific ECM regression results.
        seed : int, optional
            Random seed for the bootstrap procedure. If None, the bootstrap RNG is
            initialized without a fixed seed; pass an integer for reproducibility.
            Defaults to None.
        verbose : bool, optional
            If True, prints detailed messages to the console. Defaults to False.

        Notes
        -----
        **Lag and Lead Logic:**
        The model uses these parameters to capture short-run dynamics. If the 
        `auto` flag is triggered (via range input), the class utilizes an 
        Information Criterion (AIC, BIC, or Westerlund-specific) to find the optimal 
        balance between model fit and parsimony for each cross-sectional unit.

        

        **Input Constraints:**
        - A maximum of 6 regressors is supported for asymptotic moment lookup.
        - If `trend=True`, the model must also include a `constant`.
        - When `westerlund=True`, the specification is restricted to at most one 
          regressor and requires at least a constant.
        """
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data must be a pandas DataFrame.")

        self.raw_data = data.copy()
        self.y_var = y_var
        self.x_vars = x_vars if isinstance(x_vars, list) else [x_vars]
        self.id_var = id_var
        self.time_var = time_var

        if len(self.x_vars) == 0:
            raise ValueError("At least one regressor must be specified.")
        if len(self.x_vars) > 6:
            raise ValueError("No more than 6 covariates can be specified.")
        required_cols = [self.y_var, self.id_var, self.time_var] + self.x_vars
        missing_cols = [c for c in required_cols if c not in self.raw_data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        if trend and not constant:
            raise ValueError("If a trend is included, a constant must be included as well.")
        if westerlund and not constant:
            raise ValueError("If westerlund=True, at least a constant must be included.")
        if westerlund and len(self.x_vars) > 1:
            raise ValueError("If westerlund=True, at most one x-variable may be included.")
        if not isinstance(lrwindow, (int, np.integer)) or lrwindow < 0:
            raise ValueError("lrwindow must be a non-negative integer.")
        if not isinstance(bootstrap, (int, np.integer)):
            raise ValueError("bootstrap must be an integer.")

        # Parse lags/leads
        if isinstance(lags, int):
            self.minlag, self.maxlag = lags, lags
        else:
            self.minlag, self.maxlag = min(lags), max(lags)
            
        if leads is None:
            self.minlead, self.maxlead = 0, 0
        elif isinstance(leads, int):
            self.minlead, self.maxlead = leads, leads
        else:
            self.minlead, self.maxlead = min(leads), max(leads)

        for nm, val in [("minlag", self.minlag), ("maxlag", self.maxlag),
                        ("minlead", self.minlead), ("maxlead", self.maxlead)]:
            if not isinstance(val, (int, np.integer)) or val < 0:
                raise ValueError(f"{nm} must be a non-negative integer.")

        self.lrwindow = lrwindow
        self.constant = constant
        self.trend = trend
        self.westerlund = westerlund
        self.aic = aic
        self.bootstrap = bootstrap
        self.indiv_ecm = indiv_ecm
        self.seed = seed
        self.verbose = verbose
        self.nox = len(self.x_vars)
        self.auto = (self.minlag != self.maxlag) or (self.minlead != self.maxlead)

        self.results_bundle = None
        self.indiv_reg = []
        self.mg_results = None
        self.mg_tables = None

    def _get_ts_map(self, series, time_vec):
        """
        Constructs a mapping between time indices and observed values to enable 
        gap-aware time-series indexing.

        This helper function creates a dictionary where the keys are the time indices 
        and the values are the corresponding observations. This map is used by 
        downstream methods (like `_get_lag_lead`) to perform lookups based on actual 
        time values rather than row positions.

        Parameters
        ----------
        series : array-like
            A numeric vector of observations for a specific variable.
        time_vec : array-like
            A vector of time indices corresponding one-to-one with the series.

        Returns
        -------
        dict
            A dictionary mapping time indices to values.

        Notes
        -----
        By indexing via a map, the code ensures that if a lag is requested for 
        time T, it retrieves the value associated with time T-1. If T-1 is missing 
        from the keys, a gap is correctly identified and an NA can be returned 
        rather than incorrectly shifting the vector.
        """
        return dict(zip(time_vec, series))

    def _get_lag_lead(self, ts_map, time_vec, k):
        """
        Retrieves lagged or led values from a time-indexed mapping while explicitly
        respecting gaps in the time series.

        This method replicates the behavior of Stata's `L.` (lag) and `F.` (lead) 
        operators. Instead of shifting the vector by position, it calculates the 
        target time (T - k) for every observation and looks it up in the provided 
        time-series map.

        Parameters
        ----------
        ts_map : dict
            A dictionary mapping time indices to observed values, typically 
            generated by `_get_ts_map`.
        time_vec : array-like
            The vector of current time indices for which lags or leads are required.
        k : int
            The offset to apply. 
            - If k > 0: Performs a lag (e.g., k=1 retrieves values for t-1).
            - If k < 0: Performs a lead (e.g., k=-1 retrieves values for t+1).

        Returns
        -------
        numpy.ndarray
            A vector of the same length as `time_vec` containing the shifted values. 
            If a target time index does not exist in `ts_map`, `np.nan` is returned 
            for that specific observation.

        Notes
        -----
        This approach is robust to "holes" in the data. For instance, if the 
        time index jumps from 5 to 7, a lag of 1 for time 7 will correctly return 
        `NaN` (since time 6 is missing), whereas a standard position-based shift 
        would incorrectly return the value from time 5.
        """
        return np.array([ts_map.get(t - k, np.nan) for t in time_vec])

    def _lrvar_kernel(self, series, maxlag):
        """
        Computes a Bartlett-kernel long-run variance estimate for a univariate 
        series.

        This method estimates the spectral density at frequency zero using a 
        weighted sum of autocovariances. It replicates the specific scaling 
        and centering behavior used in the Westerlund (2007) test implementation.

        Parameters
        ----------
        series : array-like
            The input numeric series. Missing values (NaNs) are removed prior 
            to calculation.
        maxlag : int
            The maximum lag (bandwidth) for the Bartlett kernel. If 0, the 
            function returns the sample variance (gamma_0).

        Returns
        -------
        float
            The long-run variance estimate. Returns 0.0 if the series contains 
            no non-missing observations.

        Notes
        -----
        **Mathematical Implementation:**
        1.  **Autocovariance Calculation:** For a cleaned series of length $n$, 
            the lag-$j$ autocovariance ($\gamma_j$) is calculated as:
            $$\gamma_j = \frac{1}{n} \sum_{t=j+1}^{n} x_t x_{t-j}$$
            Note the use of $1/n$ as the denominator for all lags.
        2.  **Bartlett Weights:** The weights decrease linearly as the lag 
            increases:
            $$w_j = 1 - \frac{j}{maxlag + 1}$$
        3.  **Aggregation:** The final estimate is:
            $$\hat{\Omega} = \gamma_0 + 2 \sum_{j=1}^{maxlag} w_j \gamma_j$$

        

        This helper is primarily used to calculate the unit-specific $S_i$ and 
        $S_{i1}$ terms required for the $G_a$ and $P_a$ statistics.
        """
        x = series[~np.isnan(series)]
        n = len(x)
        if n <= 0: return 0.0
        gamma_0 = np.sum(x**2) / n
        if maxlag == 0: return gamma_0
        sum_weighted = 0
        for j in range(1, maxlag + 1):
            if n > j:
                gamma_j = np.sum(x[j:] * x[:-j]) / n
                weight = 1.0 - (j / (maxlag + 1))
                sum_weighted += weight * gamma_j
        return gamma_0 + 2.0 * sum_weighted

    def _selection_loop(self, y_map, x_maps, dy_map, dx_maps, time_vec):
        """
        Performs unit-specific optimal lag and lead selection for the Error Correction Model 
        using an information criterion search.

        This method replicates the Stata selection logic by iterating through candidate 
        lag (p) and lead (q) combinations. It constructs a temporary regression matrix 
        for each pair and evaluates it using either the standard Akaike Information 
        Criterion (AIC), Bayesian Information Criterion (BIC), or the specialized Westerlund
        (2007) criterion.

        Parameters
        ----------
        y_map : dict
            Time-series mapping for the dependent variable (levels).
        x_maps : list of dict
            List of time-series mappings for the regressors (levels).
        dy_map : dict
            Time-series mapping for the first difference of the dependent variable.
        dx_maps : list of dict
            List of time-series mappings for the first differences of the regressors.
        time_vec : array-like
            The actual time indices for the specific unit being evaluated.

        Returns
        -------
        tuple
            A tuple (best_l, best_ld) representing the optimal lag order and lead 
            order, respectively, that minimized the chosen information criterion.

        Notes
        -----
        **Information Criteria:**
        - **Standard Mode:** Uses `mod.aic` or `mod.bic`.
        - **Westerlund Mode:** Implements the specific penalization structure:
          $$IC = \ln\left(\frac{RSS}{T-p-q-1}\right) + \frac{2(p+q+d+1)}{T - p_{max} - q_{max}}$$
          where $d$ represents the count of deterministic terms.

        **Search Strategy:**
        The loop iterates from the maximum requested orders down to the minimum 
        orders (`range(max, min-1, -1)`). This ensures that the most restrictive 
        data requirements are tested first.

        

        **Selection Stability:**
        For each candidate $(p, q)$, the method rebuilds the full Right-Hand Side (RHS) 
        matrix including deterministic terms, lagged levels, and the short-run dynamics 
        of all variables to ensure the IC comparison is valid and based on the same 
        model structure used in the final estimation.
        """
        best_ic = np.inf
        best_l, best_ld = self.maxlag, self.maxlead
        
        # Determine total T for this group from time_vec
        thisti = len(time_vec)
        dy = self._get_lag_lead(dy_map, time_vec, 0)
        d_count = int(self.constant) + int(self.trend)

        for l in range(self.maxlag, self.minlag - 1, -1):
            for ld in range(self.maxlead, self.minlead - 1, -1):
                rhs = []
                # 1. Deterministics
                if self.constant: rhs.append(np.ones(thisti))
                if self.trend:    rhs.append(np.arange(1, thisti + 1))
                
                # 2. Levels (Selection MUST include these for correct IC)
                rhs.append(self._get_lag_lead(y_map, time_vec, 1))
                for x_map in x_maps:
                    rhs.append(self._get_lag_lead(x_map, time_vec, 1))
                
                # 3. dy lags
                if l > 0:
                    for k in range(1, l + 1):
                        rhs.append(self._get_lag_lead(dy_map, time_vec, k))
                
                # 4. dx dynamics (Current + Leads + Lags)
                for dx_map in dx_maps:
                    for k in range(-ld, l + 1):
                        rhs.append(self._get_lag_lead(dx_map, time_vec, k))
                
                RHS = np.column_stack(rhs)
                valid = ~np.isnan(dy) & ~np.isnan(RHS).any(axis=1)
                
                if valid.sum() > RHS.shape[1]:
                    mod = sm.OLS(dy[valid], RHS[valid]).fit()
                    rss = mod.ssr
                    
                    if self.westerlund:
                        # Stata Westerlund IC
                        term1 = np.log(rss / (thisti - l - ld - 1))
                        penalty = 2 * (l + ld + d_count + 1) / (thisti - self.maxlag - self.maxlead)
                        ic = term1 + penalty
                    else:
                        ic = mod.aic if self.aic else mod.bic
                    
                    if ic < best_ic:
                        best_ic, best_l, best_ld = ic, l, ld
                        
        return best_l, best_ld

    def _selection_loop_null(self, dy_map, dx_maps, time_vec):
        """
        Select the unit-specific lag/lead orders for the bootstrap DGP under H0.

        This mirrors the restricted regression used inside Stata's
        ``WesterlundBootstrap`` routine.  Unlike the selection regression used by
        ``WesterlundPlain``, the bootstrap selection excludes lagged levels of y
        and x and estimates only deterministic terms, lagged differences of y,
        and current/lagged/led differences of x.
        """
        best_ic = np.inf
        best_l, best_ld = self.maxlag, self.maxlead

        thisti = len(time_vec)
        dy = self._get_lag_lead(dy_map, time_vec, 0)
        d_count = int(self.constant) + int(self.trend)

        for l in range(self.maxlag, self.minlag - 1, -1):
            for ld in range(self.maxlead, self.minlead - 1, -1):
                rhs = []

                # Stata bootstrap code uses a within-unit observation counter for
                # the trend (``by id: gen tren = trend*_n``).
                if self.constant:
                    rhs.append(np.ones(thisti))
                if self.trend:
                    rhs.append(np.arange(1, thisti + 1, dtype=float))

                if l > 0:
                    for k in range(1, l + 1):
                        rhs.append(self._get_lag_lead(dy_map, time_vec, k))

                for dx_map in dx_maps:
                    for k in range(-ld, l + 1):
                        rhs.append(self._get_lag_lead(dx_map, time_vec, k))

                if not rhs:
                    continue

                RHS = np.column_stack(rhs)
                valid = ~np.isnan(dy) & ~np.isnan(RHS).any(axis=1)
                if valid.sum() <= RHS.shape[1]:
                    continue

                mod = sm.OLS(dy[valid], RHS[valid]).fit()
                rss = mod.ssr

                if self.westerlund:
                    # Exact IC expression used by WesterlundBootstrap.
                    denom1 = thisti - l - ld - 1
                    denom2 = thisti - self.maxlag - self.maxlead
                    if denom1 <= 0 or denom2 <= 0 or rss <= 0:
                        continue
                    term1 = np.log(rss / denom1)
                    penalty = 2 * (l + ld + d_count + 1) / denom2
                    ic = term1 + penalty
                else:
                    # xtwest uses estat ic/AIC here.  Retain the class's BIC
                    # option for callers that deliberately request it.
                    ic = mod.aic if self.aic else mod.bic

                if ic < best_ic:
                    best_ic, best_l, best_ld = ic, l, ld

        return best_l, best_ld

    def _get_optimal_model(self, y, X, time_vec, null_model=False, return_full_model=False):
        """
        Coordinates the unit-specific model selection and final Error Correction Model (ECM) 
        estimation.

        This method follows the two-phase logic established by Westerlund (2007) and 
        implemented in Stata's 'xtwest'. It first identifies the optimal short-run 
        dynamics (lags and leads) and then estimates the full model to extract the 
        error correction coefficients.

        Parameters
        ----------
        y : array-like
            The dependent variable series for a single cross-sectional unit.
        X : numpy.ndarray
            The matrix of regressors for a single cross-sectional unit.
        time_vec : array-like
            The time index vector corresponding to the observations.
        null_model : bool, optional
            If True, the method estimates the model under the null hypothesis of no 
            cointegration by excluding the lagged level terms (y_{t-1} and x_{t-1}). 
            This is primarily used for generating bootstrap residuals. Defaults to False.
        return_full_model : bool, optional
            If True, returns an extended set of results including residuals, the 
            validity mask, and the design matrix. Defaults to False.

        Returns
        -------
        tuple
            Depending on `return_full_model`, returns:
            - (model_final, best_lag, best_lead) 
            OR
            - (model_final, resid, best_lag, best_lead, valid_F, RHS_F, dy_vec)

        Process
        -------
        **1. Map Generation:**
        Constructs time-series maps for levels and first-differences of all variables 
        to ensure gap-aware lagging and leading.

        **2. Phase 1: Selection:**
        If automatic selection is enabled, it calls `_selection_loop` to find the 
        lag (p) and lead (q) orders that minimize the information criterion. 

        

        **3. Phase 2: Final Estimation:**
        Constructs the final Right-Hand Side (RHS) design matrix including:
        - Deterministic terms (constant and/or trend).
        - Lagged levels of y and X (if not `null_model`).
        - Lagged differences of y.
        - Leads, lags, and current differences of X.
        
        The model is then fitted using Ordinary Least Squares (OLS) on the subset of 
        observations where all lagged/led terms are non-missing.

        Notes
        -----
        - The use of `null_model=True` is critical for the residual-based bootstrap 
          procedure to ensure the simulated data correctly represents the 
          non-cointegrated case.
        """
        T_len = len(y)
        
        # Generate basic TS maps
        y_map = self._get_ts_map(y, time_vec)
        x_maps = [self._get_ts_map(X[:, j], time_vec) for j in range(self.nox)]
        
        # Generate Delta variables (dy and dx)
        dy_vec = y - self._get_lag_lead(y_map, time_vec, 1)
        dy_map = self._get_ts_map(dy_vec, time_vec)
        
        dx_maps = [
            self._get_ts_map(X[:, j] - self._get_lag_lead(x_maps[j], time_vec, 1), time_vec)
            for j in range(self.nox)
        ]
        
        # ---------- Phase 1: Unit-specific lag/lead selection ----------
        if self.auto:
            if null_model:
                # WesterlundBootstrap selects p_i and q_i from the restricted
                # H0 regression, which excludes lagged level terms.
                best_lag, best_lead = self._selection_loop_null(dy_map, dx_maps, time_vec)
            else:
                best_lag, best_lead = self._selection_loop(
                    y_map, x_maps, dy_map, dx_maps, time_vec
                )
        else:
            best_lag, best_lead = self.maxlag, self.maxlead

        # ---------- Phase 2: Final Estimation (Full ECM) ----------
        rhs_fin = []
        d_count = 0
        
        # 1. Deterministic terms
        if self.constant:
            rhs_fin.append(np.ones(T_len))
            d_count += 1
        if self.trend:
            rhs_fin.append(np.arange(1, T_len + 1))
            d_count += 1

        # 2. Levels (These are excluded under H0 null_model for bootstrapping)
        if not null_model:
            # Lagged Level of Y: y_{t-1}
            rhs_fin.append(self._get_lag_lead(y_map, time_vec, 1))
            # Lagged Levels of X: x_{t-1}
            for j in range(self.nox):
                rhs_fin.append(self._get_lag_lead(x_maps[j], time_vec, 1))

        # 3. Delta Y Lags: L(1/best_lag)D.y
        if best_lag > 0:
            for k in range(1, best_lag + 1):
                rhs_fin.append(self._get_lag_lead(dy_map, time_vec, k))

        # 4. Delta X Dynamics: L(-best_lead/best_lag)D.x
        # Stata: L(-lead/lag)D.x includes the current difference (lag 0)
        for k in range(-best_lead, best_lag + 1):
            for j in range(self.nox):
                rhs_fin.append(self._get_lag_lead(dx_maps[j], time_vec, k))

        RHS_F = np.column_stack(rhs_fin)
        
        # Identify valid support for the final model
        valid_F = ~np.isnan(dy_vec) & ~np.isnan(RHS_F).any(axis=1)
        
        # Fit OLS
        model_final = sm.OLS(dy_vec[valid_F], RHS_F[valid_F]).fit()

        if return_full_model:
            resid = np.full(T_len, np.nan)
            resid[valid_F] = model_final.resid
            return model_final, resid, best_lag, best_lead, valid_F, RHS_F, dy_vec

        return model_final, best_lag, best_lead

    def _run_westerlund_plain(self, df, is_boot=None):
        """
        Internal plain (non-bootstrap) routine for computing the four Westerlund (2007)
        ECM-based panel cointegration test statistics: Gt, Ga, Pt, and Pa.

        This method implements a two-pass estimation procedure that closely replicates 
        the logic of the Stata 'xtwest' command. It handles unbalanced panels, 
        automatic lag/lead selection, and time-continuity checks implicitly through 
        time-indexed mapping helpers.

        Parameters
        ----------
        df : pandas.DataFrame
            The input panel data. Must contain the cross-sectional identifier (id_var) 
            and time index (time_var).
        is_boot : bool, optional
            A flag to indicate if the function is being called as part of a 
            bootstrap replication. If True, it suppresses diagnostic printing and 
            Mean-Group (MG) summary tables. Defaults to None.

        Returns
        -------
        tuple
            A tuple containing:
            - summary_stats (dict): A dictionary with the raw values of Gt, Ga, 
              Pt, and Pa.
            - indiv_stats (list): A list of dictionaries containing granular 
              unit-specific results (alpha_i, standard errors, optimal lags, 
              residuals, and long-run variance components).

        Process
        -------
        **Pass 1: Unit-Specific Regressions (Mean-Group Statistics)**
        For each cross-sectional unit, an Error Correction Model (ECM) is estimated:
        1.  Determines optimal lag (p) and lead (q) lengths using information 
            criteria (AIC or Westerlund-specific).
        2.  Regresses Δy on deterministic terms, y[t-1], x[t-1], and short-run 
            dynamics (Δy lags and Δx leads/lags).
        3.  Computes Gt (the mean of individual t-statistics of the error 
            correction coefficient) and Ga (the mean of scaled error 
            correction coefficients).
        4.  Computes individual Long-Run (LR) coefficients (beta) and Long-Run 
            variances using a Bartlett kernel.
        5.  Aggregates Mean-Group results for the error correction term and 
            LR relationship, displaying them if not in bootstrap mode.

        **Pass 2: Pooled Estimation (Panel Statistics)**
        1.  Calculates the average lag and lead lengths across all units.
        2.  Partials out the deterministic terms and short-run dynamics from 
            the levels of y[t-1] and Δy using the average lag/lead orders.
        3.  Aggregates the products of these filtered residuals across all units 
            to estimate a common error correction coefficient.
        4.  Computes Pt (pooled t-ratio) and Pa (pooled scaled coefficient) 
            statistics.

        Notes
        -----
        - The `westerlund` flag enables specific finite-sample adjustments and 
          trimming rules for variance estimation as proposed in the 2007 paper.
        - `is_boot` logic ensures that internal bootstrap iterations remain 
          silent and efficient.
        - Standard errors for mean-group coefficients are calculated based on 
          the variance of the individual estimates across groups.
        """
        is_boot = True if (self.bootstrap > 0 and is_boot is None) else False

        # 1. Setup
        df = df.sort_values([self.id_var, self.time_var]).reset_index(drop=True)
        ids = df[self.id_var].unique()
        n_groups = len(ids)
        
        indiv_stats = []
        opt_lags = []
        opt_leads = []

        # ============================================================================
        # PASS 1: Individual Regressions (Gt and Ga)
        # ============================================================================
        for gid in ids:            
            gdata = df[df[self.id_var] == gid].copy()
            y = gdata[self.y_var].values
            X = gdata[self.x_vars].values
            tv = gdata[self.time_var].values
            ti_orig = len(y)

            # Determine Optimal Lags/Leads for this specific group
            model_opt, blag, blead = self._get_optimal_model(y, X, tv)
            opt_lags.append(blag)
            opt_leads.append(blead)

            # Build Local Maps for Series Construction
            y_map = self._get_ts_map(y, tv)
            dy = y - self._get_lag_lead(y_map, tv, 1)
            ly = self._get_lag_lead(y_map, tv, 1)

            # Construct RHS Matrix following the exact order of lr_count
            # Order: [constant, trend, ly, lx_1...lx_n, dy_lag1...dy_lagp, dx_dynamics]
            rhs_list, col_names = self._build_rhs_ordered(gdata, tv, blag, blead)
            RHS = np.column_stack(rhs_list)
            y_map = self._get_ts_map(y, tv)
            dy = y - self._get_lag_lead(y_map, tv, 1)
            ly = self._get_lag_lead(y_map, tv, 1)

            valid = ~np.isnan(dy) & ~np.isnan(RHS).any(axis=1)
            
            # Regression
            model = sm.OLS(dy[valid], RHS[valid]).fit(method="qr")
            # Call the display function for this specific group (Unit i)
            # We use the column names from _build_rhs_ordered to label the coefficients
            
            if is_boot == False and self.indiv_ecm:
                indiv_reg = self._reg_display(
                    params=model.params,
                    bse=model.bse,
                    tvalues=model.tvalues,
                    pvalues=model.pvalues,
                    title=f"Individual Regression Results for ID: {gid} (Lags={blag}, Leads={blead})"
                )
                self.indiv_reg.append({gid: indiv_reg})
            
            # Locate alpha_i (coefficient of L.y)
            alpha_idx = int(self.constant) + int(self.trend)
            alpha_i = model.params[alpha_idx]
            se_alpha_i = model.bse[alpha_idx]

            # Extract Long-Run Coefficients (beta_i)
            # beta_i = -(coefficient of lx) / alpha_i
            beta_is = []
            for j in range(self.nox):
                gamma_idx = alpha_idx + 1 + j
                gamma_ij = model.params[gamma_idx]
                # Stata xtwest logic: beta = -gamma / alpha
                beta_is.append(-gamma_ij / alpha_i if alpha_i != 0 else 0.0)

            # --- Long Run Variance (wysq) ---
            dytmp = dy.copy()
            if self.westerlund:
                # Replicate R/Stata trimming: trim start by lags, end by leads
                if blag > 0: 
                    dytmp[:blag + 1] = np.nan
                if blead > 0: 
                    dytmp[-blead:] = np.nan
            elif self.constant and self.trend:
                # Standard detrending (centering) if not using Westerlund-specific logic
                v_mask = ~np.isnan(dytmp)
                if v_mask.sum() > 1:
                    dytmp[v_mask] -= np.nanmean(dytmp)

            wysq = self._lrvar_kernel(dytmp, self.lrwindow)

            # --- Residual u calculation (Dimension Mismatch Fix) ---
            # lr_count = count of [constant, trend, ly, lx..., dy_lags...]
            lr_count = alpha_idx + 1 + self.nox + blag
            
            # long-run-only validity: dy and only the first lr_count columns
            valid_lr = ~np.isnan(dy) & ~np.isnan(RHS[:, :lr_count]).any(axis=1)

            u_vec = np.full(ti_orig, np.nan)
            u_vec[valid_lr] = dy[valid_lr] - (RHS[valid_lr, :lr_count] @ model.params[:lr_count])

            
            if self.westerlund:
                if blag > 0: u_vec[:blag + 1] = np.nan
                if blead > 0: u_vec[-blead:] = np.nan
            
            wusq = self._lrvar_kernel(u_vec, self.lrwindow)
            aonesemi = np.sqrt(wusq / wysq) if wysq > 0 else 0

            # Stata kp: cons + trend + nox + lags + nox*(lags + leads + 1)
            kp_stata = int(self.constant) + int(self.trend) + self.nox + blag + self.nox*(blag + blead + 1)

            if self.westerlund:
                tnorm = ti_orig - blag - blead - 1
                alt = tnorm - (kp_stata + 1)
                se_alpha_i = se_alpha_i * np.sqrt(max(0, alt)) / np.sqrt(tnorm)
            else:
                # Standard Mode uses the "Small Sample" adjustment for tnorm
                # matches: gen `tnorm' = `ti' - `lags' - `leads' - 1 - (kp_stata) - 1
                tnorm = ti_orig - blag - blead - 1 - kp_stata - 1

            indiv_stats.append({
                'ai': alpha_i, 'seai': se_alpha_i, 'betai': beta_is,
                'aonesemi': aonesemi, 'wysq': wysq, 'wusq': wusq,
                'tnorm': tnorm, 'gid': gid, 'dy': dy, 'ly': ly,
                'sub': gdata, 'tv': tv, 'ti': ti_orig, 'blag': blag, 'blead': blead
            })

        Gt = np.mean([s['ai'] / s['seai'] for s in indiv_stats])
        Ga = np.mean([s['tnorm'] * s['ai'] / s['aonesemi'] for s in indiv_stats])

        # ============================================================================
        # PASS 1.5: Aggregate Mean-Group Coefficients (Stata mgdisplay logic)
        # ============================================================================
        # 1. Error Correction (Alpha) Aggregation
        if is_boot == False:
            # 1. Error Correction (Alpha) - Always included
            all_ais = [s['ai'] for s in indiv_stats]
            mg_alpha_val = np.mean(all_ais)
            se_mg_alpha = np.std(all_ais, ddof=1) / np.sqrt(n_groups)

            # 2. Long Run (Beta) - Always included
            lr_results_list = []
            for j, xname in enumerate(self.x_vars):
                all_betas_j = [s['betai'][j] for s in indiv_stats]
                beta_mean = np.mean(all_betas_j)
                beta_se = np.std(all_betas_j, ddof=1) / np.sqrt(n_groups)
                beta_z = beta_mean / beta_se
                lr_results_list.append({
                    'Variable': xname,
                    'Coef.': beta_mean,
                    'Std. Err.': beta_se,
                    'z': beta_z,
                    'P>|z|': 2 * stats.norm.sf(np.abs(beta_z))
                })

            # 3. Call MG Display with the auto flag
            # If self.auto is True, the display method will print the warning about omitted SR terms
            self.mg_results = self._mg_display(
                mg_results={
                    'Variable': 'ec (alpha)',
                    'Coef.': mg_alpha_val,
                    'Std. Err.': se_mg_alpha,
                    'z': mg_alpha_val / se_mg_alpha,
                    'P>|z|': 2 * stats.norm.sf(np.abs(mg_alpha_val / se_mg_alpha))
                },
                lr_results=lr_results_list,
                auto=self.auto # Replicates Stata's conditional display
            )
        
        # ============================================================================
        # PASS 2: Pooled Regressions
        # ============================================================================
        mean_lag = int(np.mean(opt_lags))
        mean_lead = int(np.mean(opt_leads))
        mean_big_T = df.groupby(self.id_var)[self.y_var].size().mean()

        pooled_alpha_top = 0.0
        pooled_alpha_bottom = 0.0
        sum_sisq = 0.0

        for s in indiv_stats:
            ti = s['ti']
            tv = s['tv']
            dy = s['dy']
            ly = s['ly']
            gdata = df[df[self.id_var] == s['gid']]

            # Build Pooled Z_mat (Regressors to partial out, excludes L.y)
            z_rhs, _ = self._build_rhs_ordered(gdata, tv, mean_lag, mean_lead, include_ly=False)
            Z_mat = np.column_stack(z_rhs)
            
            # Partialling out
            v_dy = ~np.isnan(dy) & ~np.isnan(Z_mat).any(axis=1)
            v_ly = ~np.isnan(ly) & ~np.isnan(Z_mat).any(axis=1)
            
            dy_resid = np.full(ti, np.nan)
            dy_resid[v_dy] = sm.OLS(dy[v_dy], Z_mat[v_dy]).fit(method="qr").resid
            
            y_resid = np.full(ti, np.nan)
            y_resid[v_ly] = sm.OLS(ly[v_ly], Z_mat[v_ly]).fit(method="qr").resid

            # Pooled u and sigma calculation
            full_X = np.column_stack([ly, Z_mat])
            v_f = ~np.isnan(dy) & ~np.isnan(full_X).any(axis=1)
            mod_full = sm.OLS(dy[v_f], full_X[v_f]).fit(method="qr")
            
            # Use mean_lag for the Long Run count in pooled pass
            lr_c_pool = 1 + int(self.constant) + int(self.trend) + self.nox + mean_lag
            valid_lr_pool = ~np.isnan(dy) & ~np.isnan(full_X[:, :lr_c_pool]).any(axis=1)

            u_pool = np.full(ti, np.nan)
            u_pool[valid_lr_pool] = dy[valid_lr_pool] - (full_X[valid_lr_pool, :lr_c_pool] @ mod_full.params[:lr_c_pool])
    
            
            if self.westerlund:
                if mean_lag > 0: u_pool[:mean_lag+1] = np.nan
                if mean_lead > 0: u_pool[-mean_lead:] = np.nan
            wusq_p = self._lrvar_kernel(u_pool, self.lrwindow)

            dytmp_p = s['dy'].copy()
            if self.westerlund:
                if mean_lag > 0: dytmp_p[:mean_lag+1] = np.nan
                if mean_lead > 0: dytmp_p[-mean_lead:] = np.nan
            elif self.constant and self.trend:
                v = ~np.isnan(dytmp_p)
                if v.sum() > 1:
                    dytmp_p[v] -= np.nanmean(dytmp_p[v])
            
            wysq_p = self._lrvar_kernel(dytmp_p, self.lrwindow)
            aonesemi_p = np.sqrt(wusq_p / wysq_p) if wysq_p > 0 else 0

            # Accumulate Alpha Components
            vp = ~np.isnan(dy_resid) & ~np.isnan(y_resid)
            pooled_alpha_top += np.sum((1.0 / aonesemi_p) * y_resid[vp] * dy_resid[vp])
            pooled_alpha_bottom += np.sum(y_resid[vp]**2)

            # Accumulate sum_sisq for SE
            sigmasqi = np.sum(mod_full.resid**2)

            if self.westerlund:
                tnorm_p = mean_big_T - mean_lag - mean_lead - 1
            else:
                kp_pool = int(self.constant) + int(self.trend) + self.nox + mean_lag + self.nox*(mean_lag + mean_lead + 1)
                tnorm_p = mean_big_T - mean_lag - mean_lead - 1 - kp_pool - 1

            si = (np.sqrt(sigmasqi / tnorm_p) / aonesemi_p)
            sum_sisq += si**2

        pooled_alpha = pooled_alpha_top / pooled_alpha_bottom
        se_pooled = (np.sqrt(sum_sisq / n_groups)) / np.sqrt(pooled_alpha_bottom)
        Pt = pooled_alpha / se_pooled
        
        # Pa calculation
        kp_pa = (int(self.constant) + int(self.trend) + self.nox + mean_lag + 
                 self.nox * (mean_lag + mean_lead + 1))
        
        if self.westerlund:
            tnorm_pa = mean_big_T - mean_lag - mean_lead - 1
        else:
            tnorm_pa = mean_big_T - mean_lag - mean_lead - 1 - kp_pa - 1
        
        Pa = tnorm_pa * pooled_alpha

        summary_stats = {'Gt': Gt, 'Ga': Ga, 'Pt': Pt, 'Pa': Pa}

        # Return both the summary stats and the granular individual data
        return summary_stats, indiv_stats
    
    def _build_rhs_ordered(self, sub, tv, l, ld, include_ly=True):
        """
        Constructs the ordered Right-Hand Side (RHS) design matrix and corresponding 
        column names for the Error Correction Model (ECM) regression.

        This internal helper ensures that the variables are stacked in a specific 
        sequence that aligns with the indexing logic used for coefficient extraction 
        (e.g., identifying the error correction coefficient $\alpha$ and the 
        long-run relationship $\beta$).

        Parameters
        ----------
        sub : pandas.DataFrame
            A subset of the data for a single cross-sectional unit.
        tv : array-like
            The time index vector for the specific unit.
        l : int
            The number of lags ($p$) to include for the short-run dynamics.
        ld : int
            The number of leads ($q$) to include for the short-run dynamics.
        include_ly : bool, optional
            Whether to include the lagged level of the dependent variable ($y_{t-1}$). 
            Setting this to False is useful for partialling-out steps in pooled 
            regressions. Defaults to True.

        Returns
        -------
        tuple
            A tuple (rhs, names) where:
            - rhs (list of numpy.ndarray): A list of vectors representing the columns 
               of the design matrix.
            - names (list of str): The descriptive labels for each column.

        Process
        -------
        The design matrix is constructed in the following strict order:
        1.  **Deterministics:** Adds a constant (vector of ones) and/or a linear 
            trend (1 to $T$) based on class settings.
        2.  **Lagged Levels (The Cointegrating Vector):** Includes $y_{t-1}$ 
            (if `include_ly` is True) followed by $x_{t-1}$ for each regressor. 
            These terms represent the long-run equilibrium.
        3.  **Lagged Differences of Y:** Adds $\Delta y_{t-k}$ for $k = 1, \dots, l$. 
            These capture the persistence in short-run dynamics of the dependent variable.
        4.  **Dynamics of X:** For each regressor $X$, it adds the lead, current, 
            and lagged differences: $\Delta x_{t-k}$ for $k = -ld, \dots, l$. 

        

        Notes
        -----
        - This function relies on `_get_ts_map` and `_get_lag_lead` to ensure that 
          lags and leads correctly handle any gaps in the time index `tv`.
        - The column names generated here are used for labeling results in 
          `_reg_display` when the `verbose` flag is active.
        - The ordering is critical because the main estimation routine assumes 
          that the error correction coefficient ($\alpha$) is located at the index 
          immediately following the deterministic terms.
        """
        rhs, names = [], []
        y_map = self._get_ts_map(sub[self.y_var], tv)
        
        # 1. Deterministics
        if self.constant: 
            rhs.append(np.ones(len(tv))); names.append('cons')
        if self.trend:    
            rhs.append(np.arange(1, len(tv) + 1)); names.append('trend')
        
        # 2. Lagged Levels (The ECM terms)
        if include_ly:    
            rhs.append(self._get_lag_lead(y_map, tv, 1)); names.append('ly')
        
        for xv in self.x_vars:
            x_map = self._get_ts_map(sub[xv], tv)
            rhs.append(self._get_lag_lead(x_map, tv, 1)); names.append(f'lx_{xv}')
            
        # 3. Lagged Differences of Y
        dy_vec = sub[self.y_var].values - self._get_lag_lead(y_map, tv, 1)
        dy_map = self._get_ts_map(dy_vec, tv)
        if l > 0:
            for k in range(1, l + 1): 
                rhs.append(self._get_lag_lead(dy_map, tv, k)); names.append(f'ldy_{k}')
            
        # 4. Dynamics of X (Leads, Current, and Lags of Delta X)
        for xv in self.x_vars:
            # Build dx_map for dynamics
            x_vals = sub[xv].values
            xm = self._get_ts_map(x_vals, tv)
            dx_vec = x_vals - self._get_lag_lead(xm, tv, 1)
            dx_map = self._get_ts_map(dx_vec, tv)
            
            # Leads (e.g., FD1.x, FD2.x)
            for k in range(-ld, l + 1): # Leads (negative), Current (0), Lags (positive)
                rhs.append(self._get_lag_lead(dx_map, tv, k))
                names.append(f'dx_{xv}_k{k}')
                
        return rhs, names
        
    def _tsset_and_clean(self):
        """
        Validates panel structure, enforces time-series continuity, and checks for 
        observation sufficiency.

        This method replicates the data-preparation guardrails found in Stata's 
        `tsset` and `marksample` routines. It ensures that the dataset is properly 
        sorted and that each cross-sectional unit (group) meets the minimum 
        requirements for estimating the Error Correction Model (ECM).

        Returns
        -------
        pandas.DataFrame
            A cleaned and sorted copy of the data containing only valid 
            observations (no NaNs in model variables) and meeting all 
            continuity requirements.

        Raises
        ------
        ValueError
            - If any unit contains "holes" (gaps) in the time index. A difference 
              of more than 1 between consecutive observations in a sorted unit 
              is considered a hole.
            - If any unit has fewer observations than the threshold required 
              to estimate the requested lag/lead orders and deterministic terms.

        Process
        -------
        1.  **Sorting and Subset Selection:** Sorts the data by the ID and Time 
            identifiers and removes any rows containing missing values in the 
            dependent variable, regressors, or panel identifiers.
        2.  **Continuity Check:** Within each group, it calculates the first 
            difference of the time variable. If any difference is greater than 1, 
            it identifies the specific IDs with gaps and raises a `ValueError`.
        3.  **Observation Sufficiency Check:** Calculates the minimum required 
            observations per unit using the formula:
            $$minobs = (p_{max} + q_{max} + 1) + (constant + trend) + 1 + nox + p_{max} + nox(p_{max} + q_{max} + 1) + 1$$
            This ensures there are enough degrees of freedom to fit the most 
            complex model allowed by the user's settings.

        Notes
        -----
        The continuity check is performed *after* dropping missing values. This 
        means that if a row is dropped because of a missing regressor, it may 
        create a "hole" in the time series, which is strictly prohibited for 
        the Westerlund test as it breaks the lag/lead structure.
        """
        df = self.raw_data.copy()
        df = df.sort_values(by=[self.id_var, self.time_var])
        cols = [self.y_var] + self.x_vars + [self.id_var, self.time_var]
        df = df.dropna(subset=cols)
        
        df['__time_diff'] = df.groupby(self.id_var)[self.time_var].diff().fillna(1.0)
        if (df['__time_diff'] > 1).any():
             ids_with_holes = df.loc[df['__time_diff'] > 1, self.id_var].unique()
             raise ValueError(f"Continuous time-series required. Holes in IDs: {ids_with_holes}")
        df = df.drop(columns=['__time_diff'])
        
        counts = df.groupby(self.id_var)[self.y_var].count()
        minobs = (self.maxlag + self.maxlead + 1 + int(self.constant) + int(self.trend) + 
                  1 + self.nox + self.maxlag + self.nox * (self.maxlag + self.maxlead + 1) + 1)
        
        if counts.min() < minobs:
            raise ValueError(f"At least {minobs} observations are required.")
            
        return df

    def _roll_nan(self, a, shift):
        """
        Performs vector-based lagging or leading with NaN padding.

        This method shifts an array by a specified number of positions and fills 
        the resulting empty slots with `np.nan`. While `_get_lag_lead` is used 
        for gap-aware lookups across a full panel, `_roll_nan` is a high-performance 
        alternative used for internal operations where the input array is already 
        known to be a continuous time-series segment (such as within bootstrap 
        recursions).

        Parameters
        ----------
        a : numpy.ndarray
            The input numeric array to be shifted.
        shift : int
            The number of positions to shift the array.
            - If shift > 0: Performs a lag (L.), shifting elements forward and 
              padding the start with NaNs.
            - If shift < 0: Performs a lead (F.), shifting elements backward and 
              padding the end with NaNs.
            - If shift == 0: Returns the original array.

        Returns
        -------
        numpy.ndarray
            A shifted copy of the input array with the same shape and type, 
            padded with `np.nan`.

        Notes
        -----
        This method is purely position-based and does not check a time index. 
        It is used extensively in the bootstrap engine to generate lagged 
        differenced terms ($\Delta y_{t-k}$) and during the integration process 
        to recover levels from differences.

        

        **Comparison to `np.roll`:**
        Unlike `numpy.roll`, which performs a circular shift (where elements 
        wrap around from the end to the beginning), `_roll_nan` treats the 
        boundaries as "lost" data and pads them with missing values, 
        maintaining the integrity of time-series causality.
        """
        out = np.full_like(a, np.nan)

        if shift == 0:
            out[:] = a
        elif shift > 0:  # lag
            out[shift:] = a[:-shift]
        else:            # lead
            k = -shift
            out[:-k] = a[k:]

        return out
    
    def _bootstrap_run(self, df):
        """
        Execute the residual-based bootstrap used by Stata's ``xtwest``.

        The bootstrap DGP is estimated under H0 (no lagged level terms).  When
        lag/lead ranges are supplied, the restricted regression performs its own
        unit-specific IC selection, as in ``WesterlundBootstrap``.  Resampling is
        then performed on *actual time clusters* so that a selected time period is
        carried across all panel units that have a valid null residual at that
        period.

        The resampling block mirrors the Stata sequence

            expandcl 2, cluster(id)
            bsample if e < ., cluster(t)
            ... construct tussent/newtt ...

        including the original xtwest ``currlead`` quirk: the lead loop used in
        the bootstrap reconstruction is controlled by the selected lead of the
        last panel unit.  Coefficients for leads not selected by a particular unit
        remain zero.
        """
        rng = Generator(MT19937(self.seed))
        idv, tv = self.id_var, self.time_var
        df = df.sort_values([idv, tv]).reset_index(drop=True)
        all_ids = sorted(df[idv].unique())

        # ------------------------------------------------------------------
        # Phase 1: estimate the restricted (H0) bootstrap DGP by panel unit.
        # ------------------------------------------------------------------
        meta_map = {}
        pool_rows = []
        stata_currlead = self.maxlead

        for gid in all_ids:
            g = df[df[idv] == gid].sort_values(tv)
            y = g[self.y_var].to_numpy(dtype=float)
            X = g[self.x_vars].to_numpy(dtype=float)
            time_vec = g[tv].to_numpy()
            Ti = len(y)

            model, resid, blag, blead, valid_F, RHS_F, dy_vec = self._get_optimal_model(
                y, X, time_vec, null_model=True, return_full_model=True
            )

            # Stata: by id: egen meane = mean(e); replace e = e - meane
            e_raw = resid.astype(float)
            valid_e = np.isfinite(e_raw)
            if not valid_e.any():
                raise ValueError(f"No valid bootstrap residuals for panel unit {gid}.")
            e_centered = e_raw - np.nanmean(e_raw)

            # Stata: d.x followed by within-id demeaning.
            dx_full = np.diff(X, axis=0, prepend=np.full((1, self.nox), np.nan))
            cdX_full = dx_full - np.nanmean(dx_full, axis=0)

            # Coefficient order in the restricted regression is
            # [constant, trend, L1..Lp D.y, Fq..F1/current/L1..Lp D.x].
            coeffs = np.asarray(model.params, dtype=float)
            c_idx = int(self.constant) + int(self.trend)

            phi = np.zeros(self.maxlag + 1, dtype=float)
            if blag > 0:
                phi[1:blag + 1] = coeffs[c_idx:c_idx + blag]
                c_idx += blag

            bFdx = np.zeros((self.maxlead + 1, self.nox), dtype=float)
            bLdx = np.zeros((self.maxlag + 1, self.nox), dtype=float)
            for k in range(-blead, blag + 1):
                block = coeffs[c_idx:c_idx + self.nox]
                c_idx += self.nox
                if k < 0:
                    bFdx[-k, :] = block
                else:
                    bLdx[k, :] = block

            meta_map[gid] = {
                "Ti": Ti,
                "phi": phi,
                "bFdx": bFdx,
                "bLdx": bLdx,
                "e": e_centered,
                "cdX": cdX_full,
                "blead": blead,
                "blag": blag,
                "time": time_vec,
            }

            # ``bsample if e < ., cluster(t)`` samples only observations whose
            # null residual is nonmissing.  Store their actual time labels now.
            for pos in np.flatnonzero(valid_e):
                rec = {
                    idv: gid,
                    tv: time_vec[pos],
                    "__e": e_centered[pos],
                    "__Ti": Ti,
                }
                for j in range(self.nox):
                    rec[f"__cdx{j}"] = cdX_full[pos, j]
                pool_rows.append(rec)

            # Literal xtwest behavior: after the foreach-id loop, currlead holds
            # the selected lead of the final panel unit and is subsequently used
            # to control the bootstrap lead loop for every unit.
            stata_currlead = blead

        if not pool_rows:
            raise ValueError("No valid observations are available for bootstrap resampling.")

        eligible_pool = pd.DataFrame(pool_rows)

        # Stata's expandcl 2, cluster(id): two copies of every eligible row are
        # present when the subsequent time-cluster bootstrap is drawn.  The new
        # cluster identifier created by expandcl is not used later by xtwest.
        expanded_pool = pd.concat(
            [eligible_pool.assign(__expanded=1), eligible_pool.assign(__expanded=2)],
            ignore_index=True,
        )

        # bsample ..., cluster(t) draws Nc time clusters by default, where Nc is
        # the number of distinct eligible time clusters.  Work with positions in
        # this list so datetime/string/numeric time labels are all handled safely.
        time_keys = list(eligible_pool[tv].drop_duplicates())
        if not time_keys:
            raise ValueError("No eligible time clusters are available for bootstrap resampling.")
        rows_by_time = {
            key: expanded_pool.loc[expanded_pool[tv] == key].copy()
            for key in time_keys
        }

        results_dist = {k: [] for k in ["Gt", "Ga", "Pt", "Pa"]}

        # ------------------------------------------------------------------
        # Phase 2: Stata-style time-cluster bootstrap and data reconstruction.
        # ------------------------------------------------------------------
        for b in range(self.bootstrap):
            # Equivalent sampling unit to: bsample if e < ., cluster(t)
            draw_pos = rng.choice(len(time_keys), size=len(time_keys), replace=True)
            sampled_parts = [rows_by_time[time_keys[i]] for i in draw_pos]
            sampled = pd.concat(sampled_parts, ignore_index=True)

            # Stata: gen newttt = ceil(uniform()*(ti-maxlag-maxlead-1))
            upper = (
                sampled["__Ti"].to_numpy(dtype=int)
                - self.maxlag
                - self.maxlead
                - 1
            )
            if np.any(upper <= 0):
                raise ValueError("Insufficient observations for xtwest bootstrap shuffle.")
            sampled["__newttt"] = np.ceil(rng.random(len(sampled)) * upper).astype(int)

            # sort id, stable; by id: gen tussent = _n
            sampled["__stable0"] = np.arange(len(sampled))
            sampled = sampled.sort_values([idv, "__stable0"], kind="mergesort").reset_index(drop=True)
            sampled["__tussent"] = sampled.groupby(idv, sort=False).cumcount() + 1

            # sort tussent, stable; by tussent: egen newtt = mean(newttt)
            sampled["__stable1"] = np.arange(len(sampled))
            sampled = sampled.sort_values(
                ["__tussent", "__stable1"], kind="mergesort"
            ).reset_index(drop=True)
            sampled["__newtt"] = sampled.groupby("__tussent", sort=False)["__newttt"].transform("mean")

            # sort id newtt, stable.  __stable2 explicitly preserves the order
            # established by the preceding stable sort when newtt ties occur.
            sampled["__stable2"] = np.arange(len(sampled))
            sampled = sampled.sort_values(
                [idv, "__newtt", "__stable2"], kind="mergesort"
            ).reset_index(drop=True)

            boot_data_list = []

            for gid in all_ids:
                meta = meta_map[gid]
                required_len = int(meta["Ti"] + self.maxlag + self.maxlead + 2)

                # Stata: by id: keep if _n <= ti + maxlag + maxlead + 2
                sg = sampled.loc[sampled[idv] == gid].head(required_len).copy()
                if sg.empty:
                    raise ValueError(
                        f"Bootstrap replication {b + 1} contains no sampled rows for panel unit {gid}."
                    )
                sg["__newt"] = np.arange(1, len(sg) + 1)

                e_b = sg["__e"].to_numpy(dtype=float)
                cdX_b = np.column_stack([
                    sg[f"__cdx{j}"].to_numpy(dtype=float)
                    for j in range(self.nox)
                ])

                # -------------------------- e -> u --------------------------
                # Stata first zeroes missing u after blanking the first maxlag
                # observations, then adds all lag terms up to maxlag.
                u = e_b.copy()
                if self.maxlag > 0:
                    u[:self.maxlag] = np.nan
                u = np.nan_to_num(u, nan=0.0)

                # Stata loops 0/maxlag globally.  Coefficients outside a unit's
                # selected lag are zero, matching the pre-bootstrap fill step.
                for k in range(0, self.maxlag + 1):
                    Lk = self._roll_nan(cdX_b, k)
                    u += (Lk * meta["bLdx"][k, :]).sum(axis=1)

                # Literal currlead quirk from xtwest: use the final unit's
                # selected lead as the global upper limit.  Unit-specific unused
                # lead coefficients are already zero in bFdx.
                if stata_currlead > 0:
                    for k in range(1, stata_currlead + 1):
                        Fk = self._roll_nan(cdX_b, -k)
                        u += (Fk * meta["bFdx"][k, :]).sum(axis=1)

                # ------------------------- u -> dy --------------------------
                # Stata turns any missing u into zero when dy is initialized.
                u_for_dy = np.nan_to_num(u, nan=0.0)
                dy = np.zeros(len(u_for_dy), dtype=float)
                phi = meta["phi"]
                for t in range(self.maxlag, len(dy)):
                    ar_term = 0.0
                    for k in range(1, self.maxlag + 1):
                        ar_term += dy[t - k] * phi[k]
                    dy[t] = u_for_dy[t] + ar_term

                # Stata subsequently sets the first maxlag dy values missing.
                # For cumulative sums, treating them as zero gives the same
                # running total before those boundary observations are masked.
                if self.maxlag > 0:
                    dy[:self.maxlag] = 0.0

                # -------------------- integrate to levels -------------------
                booty = np.cumsum(dy)
                bootx = np.cumsum(np.nan_to_num(cdX_b, nan=0.0), axis=0)

                mask = np.ones(len(booty), dtype=bool)
                mask[:self.maxlag] = False
                mask[int(meta["Ti"] + self.maxlag):] = False

                df_b = pd.DataFrame({
                    idv: gid,
                    tv: sg["__newt"].to_numpy(),
                    self.y_var: np.where(mask, booty, np.nan),
                })
                for j, xname in enumerate(self.x_vars):
                    df_b[xname] = np.where(mask, bootx[:, j], np.nan)

                # _run_westerlund_plain expects its input already cleaned.  The
                # Stata routine keeps boundary rows but marks them out because
                # booty/bootx are missing; dropping them here preserves the same
                # usable time support while retaining the generated newt labels.
                boot_data_list.append(
                    df_b.dropna(subset=[self.y_var] + self.x_vars)
                )

            boot_panel = pd.concat(boot_data_list, ignore_index=True)
            stats_b, indiv_stats_b = self._run_westerlund_plain(boot_panel)
            for k in results_dist:
                results_dist[k].append(stats_b[k])

        return results_dist

    def _asymptotic_inference(self, results):
        """Return standardized Z-scores and asymptotic left-tail p-values."""
        ridx = 0
        if self.constant and not self.trend:
            ridx = 1
        if self.constant and self.trend:
            ridx = 2
        cidx = self.nox - 1

        moments = {
            'gt': {'mean': np.array([[-0.9763,-1.3816,-1.7093,-1.9789,-2.1985,-2.4262],[-1.7776,-2.0349,-2.2332,-2.4453,-2.6462,-2.8358],[-2.3664,-2.5284,-2.7040,-2.8639,-3.0146,-3.1710]]),
                   'var':  np.array([[1.0823,1.0981,1.0489,1.0576,1.0351,1.0409],[0.8071,0.8481,0.8886,0.9119,0.9083,0.9236],[0.6603,0.7070,0.7586,0.8228,0.8477,0.8599]])},
            'ga': {'mean': np.array([[-3.8022,-5.8239,-7.8108,-9.8791,-11.7239,-13.8581],[-7.1423,-9.1249,-10.9667,-12.9561,-14.9752,-17.0673],[-12.0116,-13.6324,-15.5262,-17.3648,-19.2533,-21.2479]]),
                   'var':  np.array([[20.6868,29.9016,39.0109,50.5741,58.9595,69.5967],[29.6336,39.3428,49.4880,58.7035,67.9499,79.1093],[46.2420,53.7428,64.5591,74.7403,84.7990,94.0024]])},
            'pt': {'mean': np.array([[-0.5105,-0.9370,-1.3169,-1.6167,-1.8815,-2.1256],[-1.4476,-1.7131,-1.9206,-2.1484,-2.3730,-2.5765],[-2.1124,-2.2876,-2.4633,-2.6275,-2.7858,-2.9537]]),
                   'var':  np.array([[1.3624,1.7657,1.7177,1.6051,1.4935,1.4244],[0.9885,1.0663,1.1168,1.1735,1.1684,1.1589],[0.7649,0.8137,0.8857,0.9985,0.9918,0.9898]])},
            'pa': {'mean': np.array([[-1.0263,-2.4988,-4.2699,-6.1141,-8.0317,-10.0074],[-4.2303,-5.8650,-7.4599,-9.3057,-11.3152,-13.3180],[-8.9326,-10.4874,-12.1672,-13.8889,-15.6815,-17.6515]]),
                   'var':  np.array([[8.3827,24.0223,39.8827,53.4518,63.2406,76.6757],[19.7090,31.2637,42.9975,57.4844,69.4374,81.0384],[37.5948,45.6890,57.9985,74.1258,81.3934,91.2392]])}
        }
        required = [self.y_var] + self.x_vars + [self.id_var, self.time_var]
        valid_rows = self.raw_data[required].notna().all(axis=1)
        N = self.raw_data.loc[valid_rows, self.id_var].nunique()
        sqrt_N = np.sqrt(N)
        z_scores, p_values = {}, {}
        for name in ['Gt', 'Ga', 'Pt', 'Pa']:
            val = results[name]
            k = name.lower()
            if self.westerlund:
                if not self.trend:
                    params = {'gt':(-1.793, 0.7904), 'ga':(-7.2014, 29.3677), 'pt':(-1.4746, 1.0262), 'pa':(-4.3559, 21.0535)}
                else:
                    params = {'gt':(-2.356, 0.6450), 'ga':(-11.8978, 44.2471), 'pt':(-2.1128, 0.7371), 'pa':(-8.9536, 35.6802)}
                mu, var = params[k]
            else:
                mu = moments[k]['mean'][ridx, cidx]
                var = moments[k]['var'][ridx, cidx]
            z = ((val - sqrt_N * mu) / np.sqrt(var)) if k == 'pt' else ((sqrt_N * val - sqrt_N * mu) / np.sqrt(var))
            z_scores[name] = float(z)
            p_values[name] = float(stats.norm.cdf(z))
        return z_scores, p_values

    def run(self):
        """Run the Westerlund test and return a feature-complete results bundle."""
        df = self._tsset_and_clean()
        main_stats, indiv_data = self._run_westerlund_plain(df, is_boot=False)

        unit_results = pd.DataFrame([{
            'id': s['gid'],
            'alpha_i': s['ai'], 'ai': s['ai'],
            'se_alpha_i': s['seai'], 'seai': s['seai'],
            'beta_i': s['betai'],
            'aonesemi': s['aonesemi'], 'wysq': s['wysq'], 'wusq': s['wusq'],
            'lag': s['blag'], 'lags': s['blag'],
            'lead': s['blead'], 'leads': s['blead'],
            'tnorm': s['tnorm'], 'obs': s['ti'], 'ti': s['ti']
        } for s in indiv_data])

        boot_dist, boot_pvals = {}, {}
        if self.bootstrap > 0:
            if self.verbose:
                print(f"Bootstrapping {self.bootstrap} replications...")
            boot_dist = self._bootstrap_run(df)
            for k in ['Gt', 'Ga', 'Pt', 'Pa']:
                obs = main_stats[k]
                boots = np.asarray(boot_dist[k], dtype=float)
                boots = boots[np.isfinite(boots)]
                if boots.size == 0:
                    boot_pvals[k] = np.nan
                else:
                    r = np.sum(boots <= obs)
                    boot_pvals[k] = (r + 1) / (boots.size + 1)

        z_scores, asymp_pvals = self._asymptotic_inference(main_stats)

        all_betas = np.asarray(unit_results['beta_i'].tolist(), dtype=float)
        n_groups = len(unit_results)
        mg_alpha = float(unit_results['alpha_i'].mean())
        se_mg_alpha = float(unit_results['alpha_i'].std(ddof=1) / np.sqrt(n_groups)) if n_groups > 1 else np.nan
        mg_beta_vec = np.mean(all_betas, axis=0)
        se_mg_beta_vec = (np.std(all_betas, axis=0, ddof=1) / np.sqrt(n_groups)) if n_groups > 1 else np.full(self.nox, np.nan)

        selected_lags = unit_results['lag'].to_numpy(dtype=float)
        selected_leads = unit_results['lead'].to_numpy(dtype=float)
        realmeanlag = float(np.mean(selected_lags))
        realmeanlead = float(np.mean(selected_leads))
        meanlag = int(realmeanlag)
        meanlead = int(realmeanlead)
        selection_method = 'Westerlund IC' if self.westerlund else ('AIC' if self.aic else 'BIC')
        model_type = 'trend' if self.trend else ('constant' if self.constant else 'none')
        avg_obs = float(unit_results['obs'].mean())

        settings = {
            'constant': self.constant, 'trend': self.trend,
            'minlag': self.minlag, 'maxlag': self.maxlag,
            'minlead': self.minlead, 'maxlead': self.maxlead,
            'meanlag': meanlag, 'meanlead': meanlead,
            'realmeanlag': realmeanlag, 'realmeanlead': realmeanlead,
            'auto': self.auto, 'selection_method': selection_method,
            'lrwindow': self.lrwindow, 'bootstrap': self.bootstrap,
            'seed': self.seed, 'T': avg_obs
        }
        metadata = {
            'bandwidth': self.lrwindow, 'n_groups': n_groups,
            'avg_obs': avg_obs, 'model_type': model_type,
            'selection_method': selection_method, 'seed': self.seed
        }

        self.results_bundle = {
            'test_stats': main_stats,
            'z_scores': z_scores,
            'p_values': asymp_pvals,
            'boot_pvals': boot_pvals,
            'boot_distributions': boot_dist,
            'unit_data': unit_results,
            'indiv_data': indiv_data,
            'mean_group': {
                'mg_alpha': mg_alpha,
                'se_mg_alpha': se_mg_alpha,
                'mg_betas': dict(zip(self.x_vars, mg_beta_vec)),
                'se_mg_betas': dict(zip(self.x_vars, se_mg_beta_vec))
            },
            'settings': settings,
            'metadata': metadata,
            'mg_results': self.mg_results,
            'mg_tables': self.mg_tables,
            'indiv_reg': self.indiv_reg
        }

        if self.verbose:
            self._display_final(main_stats, boot_pvals, z_scores=z_scores, p_values=asymp_pvals)
        return self.results_bundle

    def _display_final(self, results, boot_pvals, z_scores=None, p_values=None):
        """Print the raw, standardized, asymptotic and bootstrap test results."""
        if z_scores is None or p_values is None:
            z_scores, p_values = self._asymptotic_inference(results)
        required = [self.y_var] + self.x_vars + [self.id_var, self.time_var]
        valid_rows = self.raw_data[required].notna().all(axis=1)
        N = self.raw_data.loc[valid_rows, self.id_var].nunique()
        if self.verbose:
            print("\n" + "="*75)
            print("Westerlund ECM Panel Cointegration Tests")
            print(f"Series: {self.y_var} ~ {', '.join(self.x_vars)}")
            print(f"N (Groups): {N}")
            print(f"Lags: {self.minlag}-{self.maxlag} | Leads: {self.minlead}-{self.maxlead} | Window: {self.lrwindow}")
            print("="*75)
            print(f"{'Statistic':<10} | {'Value':<10} | {'Z-score':<10} | {'P-value':<10} | {'Robust P':<10}")
            print("-"*75)
            for name in ['Gt', 'Ga', 'Pt', 'Pa']:
                rob = boot_pvals.get(name, np.nan) if boot_pvals else np.nan
                rob_txt = f"{rob:.3f}" if np.isfinite(rob) else "-"
                print(f"{name:<10} | {results[name]:<10.3f} | {z_scores[name]:<10.3f} | {p_values[name]:<10.3f} | {rob_txt:<10}")
            print("="*75)

    def summary(self, verbose=False):
        """Return a compact summary analogous to ``summary.westerlund_test`` in R."""
        if self.results_bundle is None:
            raise RuntimeError("No results available. Call run() first.")
        res = self.results_bundle
        rows = []
        for name in ['Gt', 'Ga', 'Pt', 'Pa']:
            row = {
                'Statistic': name,
                'Value': res['test_stats'][name],
                'Z_score': res['z_scores'][name],
                'P_val_asymp': res['p_values'][name]
            }
            if res['boot_pvals']:
                row['P_val_boot'] = res['boot_pvals'].get(name, np.nan)
            rows.append(row)
        stats_table = pd.DataFrame(rows)
        out = {
            'stats_table': stats_table,
            'settings': res['settings'],
            'mean_group': res['mean_group'],
            'mg_results': res['mg_results'],
            'mg_tables': res.get('mg_tables'),
            'n_units': res['metadata']['n_groups']
        }
        if verbose:
            print("\n" + "="*70)
            print("Westerlund (2007) Panel Cointegration Test Summary")
            print("="*70)
            print(f"Units: {out['n_units']} | Average time periods: {res['settings']['T']:.2f}")
            print(f"Deterministic terms: {res['metadata']['model_type']}")
            print(f"Lag/lead selection: {res['settings']['selection_method']}")
            print("\nTest Statistics:")
            print(stats_table.to_string(index=False))
            print("\nMean Group Estimates:")
            print(pd.Series(res['mean_group']).to_string())
        return out

    def print_summary(self):
        """Print and return :meth:`summary`."""
        return self.summary(verbose=True)

    def __repr__(self):
        if self.results_bundle is None:
            return (f"WesterlundTest(y_var={self.y_var!r}, x_vars={self.x_vars!r}, "
                    f"groups={self.raw_data.dropna(subset=[self.y_var] + self.x_vars + [self.id_var, self.time_var])[self.id_var].nunique()}, bootstrap={self.bootstrap})")
        tab = self.summary(verbose=False)['stats_table']
        return "WesterlundTest results\n" + tab.to_string(index=False)

    def _reg_display(self, params, bse, tvalues, pvalues, title="Regression Results"):
        """
        Formats and prints a detailed coefficient table for unit-specific regressions.

        This method is primarily used when the `indiv_ecm` flag is active to provide 
        transparency into the individual Error Correction Model (ECM) estimations 
        that form the basis of the Mean-Group statistics ($G_t$ and $G_a$).

        Parameters
        ----------
        params : array-like
            The estimated coefficients (parameters) from the OLS regression.
        bse : array-like
            The standard errors associated with the estimated parameters.
        tvalues : array-like
            The t-statistics for each coefficient.
        pvalues : array-like
            The two-sided p-values for the hypothesis test that the coefficient 
            is equal to zero.
        title : str, optional
            The header title for the regression table, typically identifying 
            the cross-sectional unit (ID) and the lag/lead specification used. 
            Defaults to "Regression Results".

        Returns
        -------
        results_df : pandas.DataFrame
            A DataFrame containing the coefficients, standard errors, t-statistics,
            and p-values.

        Process
        -------
        1.  **Data Organization:** Constructs a temporary pandas DataFrame to 
            align the statistical metrics (Coefficients, Standard Errors, 
            t-stats, and p-values) side-by-side.
        2.  **Formatting:** Uses fixed-width string conversion to ensure the 
            printed table maintains alignment, regardless of the magnitude 
            of the values.
        3.  **Console Output:** Prints a structured block to the console, 
            bracketed by horizontal separators to match standard econometric 
            software output.

        

        Notes
        -----
        - The column names for the rows are derived from the indices of the 
          input `params`, which are typically assigned by the `_build_rhs_ordered` 
          helper (e.g., 'cons', 'ly', 'lx_var1', 'ldy_1').
        - This display is essential for identifying potential issues in specific 
          units, such as near-singular design matrices or coefficients that 
          exhibit non-standard behavior before they are aggregated into 
          panel-wide statistics.
        """
        results_df = pd.DataFrame({
            'Coef.': params,
            'Std. Err.': bse,
            't': tvalues,
            'P>|t|': pvalues
        })
        if self.verbose:
            print(f"\n{title}")
            print("-" * 60)
            print(results_df.to_string())
            print("-" * 60)
        return results_df
    
    def _mg_display(self, mg_results, lr_results, auto=False):
        """
        Formats and prints the aggregated Mean-Group (MG) results.

        This method displays the average speed of adjustment (the error correction term) 
        and the average long-run equilibrium relationship ($\beta$ coefficients) across 
        all cross-sectional units. It provides a summary of the economic relationship 
        estimated by the panel model.

        Parameters
        ----------
        mg_results : dict
            A dictionary containing the aggregated statistics for the error 
            correction coefficient ($\alpha$). Expected keys include 'Variable', 'Coef.', 'Std. Err.', 'z', and 'P>|z|'.
        lr_results : list of dict
            A list where each element is a dictionary representing a regressor's 
            long-run relationship ($\beta$). This represents the cointegrating vector.
        auto : bool, optional
            A flag indicating whether automatic lag/lead selection was used. 
            If True, the method prints a warning that short-run dynamics are 
            omitted from the display because they vary by unit. Defaults to False.

        Process
        -------
        1.  **Header Generation:** Prints a clear section header for the Mean-Group 
            ECM output.
        2.  **Conditional Caveats:** If `auto` is True, it explicitly notes that 
            non-essential short-run coefficients are not shown, maintaining 
            consistency with Stata's reporting constraints for heterogeneous 
            specifications.
        3.  **Error Correction Display:** Renders a focused table for the 
            aggregated $\alpha$ coefficient, which measures the panel's average 
            speed of adjustment toward the long-run equilibrium.
        4.  **Long-Run Relationship Display:** Renders a second table showing the 
            estimated cointegrating vector (the long-run relationship between y and X).

        

        Notes
        -----
        - **Inference:** Standard errors in this table are calculated based on the 
          cross-sectional variance of the individual unit estimates, which is the 
          hallmark of Mean-Group estimation.
        - **Interpretation:** The Long-Run table represents the target equilibrium 
          level, while the Error Correction Term indicates how much of the 
          disequilibrium is corrected in each period.
        """
        
        if self.verbose:
            print("\n" + "="*60)
            print("Mean-group error-correction model")
            if auto:
                print("Short run coefficients apart from the error-correction term are omitted")
                print("as lag and lengths might differ between cross-sectional units.")
            
            # Display the Error Correction / MG results
            print("\nError Correction Term (Alpha):")
            print(pd.DataFrame([mg_results]).to_string(index=False))
            
            print("\nEstimated long-run relationship and short run adjustment")
        
        # Store both reporting tables while preserving the historical long-run
        # DataFrame return value for backwards compatibility.
        alpha_df = pd.DataFrame([mg_results])
        lr_df = pd.DataFrame(lr_results)
        self.mg_tables = {'mg_model': alpha_df, 'long_run': lr_df}

        if self.verbose:
            print(lr_df.to_string(index=False))
            print("="*60 + "\n")
        return lr_df
    
    def plot_bootstrap(self,
                       title="Westerlund Test: Bootstrap Distributions",
                       conf_level=0.05,
                       save_path=None,
                       dpi=300,
                       figsize=(12, 10),
                       colors=None,
                       lwd=None,
                       alpha=0.5,
                       show_grid=True,
                       show_robust_p=True,
                       show=True):
        """Plot bootstrap distributions with configurable critical level and export."""
        if self.results_bundle is None or not self.results_bundle.get('boot_distributions'):
            raise RuntimeError("No bootstrap results found. Call run() with bootstrap > 0 first.")
        if not (0 < conf_level < 1):
            raise ValueError("conf_level must be between 0 and 1.")
        colors = colors or {'obs': '#D55E00', 'crit': '#0072B2', 'fill': '#CCCCCC', 'density': '#4D4D4D'}
        lwd = lwd or {'obs': 2.5, 'crit': 1.5, 'density': 1.5}
        # Accept the older Python key names as aliases.
        obs_color = colors.get('obs', '#D55E00')
        crit_color = colors.get('crit', colors.get('cv', '#0072B2'))
        density_color = colors.get('density', colors.get('kde', '#4D4D4D'))
        fill_color = colors.get('fill', density_color)

        res = self.results_bundle
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
        for i, stat in enumerate(['Gt', 'Ga', 'Pt', 'Pa']):
            ax = axes[i]
            data = np.asarray(res['boot_distributions'][stat], dtype=float)
            data = data[np.isfinite(data)]
            if data.size < 2:
                ax.set_title(f"Statistic: {stat}")
                ax.text(0.5, 0.5, "Insufficient finite bootstrap draws", ha='center', va='center', transform=ax.transAxes)
                continue
            kde = gaussian_kde(data)
            sd = np.std(data)
            x_range = np.linspace(data.min() - sd, data.max() + sd, 200)
            y_kde = kde(x_range)
            ax.plot(x_range, y_kde, color=density_color, lw=lwd.get('density', 1.5), label='Bootstrap H0 Dist.')
            ax.fill_between(x_range, y_kde, color=fill_color, alpha=alpha)

            obs_val = res['test_stats'][stat]
            cv_val = float(np.quantile(data, conf_level))
            ax.axvline(obs_val, color=obs_color, linestyle='-', linewidth=lwd.get('obs', 2.5), label='Observed Statistic')
            ax.axvline(cv_val, color=crit_color, linestyle='--', linewidth=lwd.get('crit', 1.5), label=f'{conf_level:.0%} Bootstrap CV')
            ann = f"Obs: {obs_val:.3f}\nCV: {cv_val:.3f}"
            if show_robust_p:
                p = res['boot_pvals'].get(stat, np.nan) if res['boot_pvals'] else np.nan
                if np.isfinite(p):
                    ann += f"\nRobust p: {p:.3f}"
            ax.text(0.05, 0.95, ann, transform=ax.transAxes, va='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            ax.set_title(f"Statistic: {stat}", fontweight='bold')
            if show_grid:
                ax.grid(True, which='major', alpha=0.3)
            else:
                ax.grid(False)

        fig.suptitle(title, fontsize=16, fontweight='bold')
        fig.text(0.5, 0.93, "H0: No cointegration | reject when observed statistic is left of the bootstrap critical value",
                 ha='center', fontsize=10)
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True)
        fig.tight_layout(rect=[0, 0.05, 1, 0.92])
        if save_path:
            directory = os.path.dirname(save_path)
            if directory:
                os.makedirs(directory, exist_ok=True)
            fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
        if show:
            plt.show()
        return fig

    def plot(self, *args, **kwargs):
        """Plotting interface."""
        return self.plot_bootstrap(*args, **kwargs)
