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
    Implements a functional approximation of the Westerlund (2007) ECM-based panel 
    cointegration tests.
    
    This class provides a framework for testing the null hypothesis of no cointegration
    in panel data by evaluating the presence of error correction. It computes four main 
    statistics: two mean-group tests (Gt, Ga) and two pooled/panel tests (Pt, Pa). 
    """

    def __init__(self, data, y_var, x_vars, id_var, time_var, 
                 lags=1, leads=0, lrwindow=2, 
                 constant=False, trend=False, westerlund=False, 
                 bootstrap=-1, seed=None, nosily=False):
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
            Whether to include a linear time trend. If True, 'constant' is 
            automatically treated as True. Defaults to False.
        westerlund : bool, optional
            If True, enforces specific constraints from the 2007 paper, including 
            specialized information criteria and trimming logic for variance. 
            Defaults to False.
        bootstrap : int, optional
            Number of bootstrap replications. If > 0, robust p-values are 
            calculated. Defaults to -1 (no bootstrap).
        seed : int, optional
            Random seed for the bootstrap procedure to ensure reproducibility. 
            Defaults to None.
        nosily : bool, optional
            If True, prints unit-specific regression results to the console 
            during the estimation process. Defaults to False.

        Notes
        -----
        **Lag and Lead Logic:**
        The model uses these parameters to capture short-run dynamics. If the 
        `auto` flag is triggered (via range input), the class utilizes an 
        Information Criterion (AIC or Westerlund-specific) to find the optimal 
        balance between model fit and parsimony for each cross-sectional unit.

        

        **Input Constraints:**
        - A maximum of 6 regressors is supported for asymptotic moment lookup.
        - If `trend=True`, the model must also include a `constant`.
        - When `westerlund=True`, the specification is restricted to at most one 
          regressor and requires at least a constant.
        """
        self.raw_data = data.copy()
        self.y_var = y_var
        self.x_vars = x_vars if isinstance(x_vars, list) else [x_vars]
        self.id_var = id_var
        self.time_var = time_var
        
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

        self.lrwindow = lrwindow
        self.constant = constant
        self.trend = trend
        self.westerlund = westerlund
        self.bootstrap = bootstrap
        self.seed = seed
        self.nosily = nosily
        self.nox = len(self.x_vars)
        self.auto = (self.minlag != self.maxlag) or (self.minlead != self.maxlead)

        self.results_bundle = None

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
        weighted sum of autocovariances. 

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

        This method iterates through candidate lag (p) and lead (q) combinations. It constructs 
        a temporary regression matrix for each pair and evaluates it using either the standard 
        Akaike Information Criterion (AIC) or the specialized Westerlund (2007) criterion.

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
        - **Standard Mode:** Uses `mod.aic`.
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
                        ic = mod.aic
                    
                    if ic < best_ic:
                        best_ic, best_l, best_ld = ic, l, ld
                        
        return best_l, best_ld

    def _get_optimal_model(self, y, X, time_vec, null_model=False, return_full_model=False):
        """
        Coordinates the unit-specific model selection and final Error Correction Model (ECM) 
        estimation.

        This method first identifies the optimal short-run dynamics (lags and leads) and then 
        estimates the full model to extract the error correction coefficients.

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
        
        # ---------- Phase 1: Selection (Equation 9 restricted) ----------
        if self.auto and not null_model:
            # The selection regression does NOT include lagged levels
            best_lag, best_lead = self._selection_loop(y_map, x_maps, dy_map, dx_maps, time_vec)
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
            
            if is_boot == False and self.nosily:
                self._reg_display(
                    params=model.params,
                    bse=model.bse,
                    tvalues=model.tvalues,
                    pvalues=model.pvalues,
                    title=f"Individual Regression Results for ID: {gid} (Lags={blag}, Leads={blead})"
                )
            
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
                # beta = -gamma / alpha
                beta_is.append(-gamma_ij / alpha_i if alpha_i != 0 else 0.0)

            # --- Long Run Variance (wysq) ---
            dytmp = dy.copy()
            if self.westerlund:
                # trim start by lags, end by leads
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

            # --- Residual u calculation ---
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

            # kp: cons + trend + nox + lags + nox*(lags + leads + 1)
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
                'ai': alpha_i, 'seai': se_alpha_i, 'betai': beta_is, 'aonesemi': aonesemi, 'tnorm': tnorm, 'gid': gid,
                'dy': dy, 'ly': ly, 'sub': gdata, 'tv': tv, 'ti': ti_orig, 'blag': blag, 'blead': blead
            })

        Gt = np.mean([s['ai'] / s['seai'] for s in indiv_stats])
        Ga = np.mean([s['tnorm'] * s['ai'] / s['aonesemi'] for s in indiv_stats])

        # ============================================================================
        # PASS 1.5: Aggregate Mean-Group Coefficients
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
                lr_results_list.append({
                    'Variable': xname,
                    'Coef.': np.mean(all_betas_j),
                    'Std. Err.': np.std(all_betas_j, ddof=1) / np.sqrt(n_groups),
                    't': np.mean(all_betas_j) / (np.std(all_betas_j, ddof=1) / np.sqrt(n_groups)),
                    'P>|t|': 2 * (1 - stats.t.cdf(np.abs(np.mean(all_betas_j) / (np.std(all_betas_j, ddof=1) / np.sqrt(n_groups))), df=n_groups-1))
                })

            # 3. Call MG Display with the auto flag
            # If self.auto is True, the display method will print the warning about omitted SR terms
            self._mg_display(
                mg_results={
                    'Variable': 'ec (alpha)',
                    'Coef.': mg_alpha_val,
                    'Std. Err.': se_mg_alpha,
                    't': mg_alpha_val / se_mg_alpha,
                    'P>|t|': 2 * (1 - stats.t.cdf(np.abs(mg_alpha_val / se_mg_alpha), df=n_groups-1))
                },
                lr_results=lr_results_list,
                auto=self.auto
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
          `_reg_display` when the `noisily` flag is active.
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

        This method ensures that the dataset is properly sorted and that each cross-sectional 
        unit (group) meets the minimum requirements for estimating the Error Correction Model (ECM).

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
        Executes the residual-based bootstrap procedure to obtain the empirical 
        distribution of the Westerlund test statistics under the null hypothesis 
        of no cointegration.

        Parameters
        ----------
        df : pandas.DataFrame
            The cleaned and sorted panel data from `_tsset_and_clean`.

        Returns
        -------
        dict
            A dictionary containing lists of simulated statistics for "Gt", 
            "Ga", "Pt", and "Pa", representing the bootstrap distribution.

        Process
        -------
        **Phase 1: Metadata Extraction and Null-Model Estimation**
        For each cross-sectional unit:
        1.  Estimates the ECM under the null hypothesis ($H_0: \alpha_i = 0$) 
            using `_get_optimal_model(null_model=True)`.
        2.  Extracts and centers the residuals ($e$) and the first-differences 
            of the regressors ($\Delta X$).
        3.  Stores the short-run autoregressive coefficients ($\phi$) and 
            dynamic regressor coefficients ($bLdx, bFdx$).

        **Phase 2: The Bootstrap Loop**
        For each replication:
        1.  **Cluster-Bootstrap Time:** Performs a cluster bootstrap on the 
            time dimension to maintain the contemporaneous correlation between 
            units.
        2.  **Cross-Sectional Correlation (The 'Newtt' Shuffle):** Implements 
            a stable-sort permutation logic (shuffling) that ensures every unit 
            in the panel is subjected to the same temporal realignment.
        3.  **Innovation Construction:** Generates bootstrap innovations ($u$) 
            by combining resampled residuals with the previously estimated 
            short-run dynamics of $\Delta X$.
        4.  **AR Recursion and Integration:** - Recursively generates the bootstrap $\Delta y$ series using the 
              estimated AR coefficients ($\phi$).
            - Integrates $\Delta y$ and $\Delta X$ (via cumulative sums) to 
              generate bootstrap levels ($y^*$ and $X^*$).
        5.  **Re-Estimation:** Runs the full Westerlund test on the newly 
            simulated panel using `_run_westerlund_plain` and records the 
            resulting statistics.

        

        Notes
        -----
        - **Null Hypothesis Enforcement:** By estimating the model with 
          `null_model=True`, we ensure that the generated data does not 
          contain a cointegrating relationship, making it a valid reference 
          for p-value calculation.
        - **Reproducibility:** Uses a NumPy `Generator` with a Mersenne-Twister 
          (MT19937) algorithm to ensure consistent results across runs if a 
          seed is provided.
        """
        # Ensure seed consistency
        rng = Generator(MT19937(self.seed if self.seed is not None else 123))
        idv, tv = self.id_var, self.time_var
        df = df.sort_values([idv, tv]).reset_index(drop=True)
        all_ids = sorted(df[idv].unique())
        num_ids = len(all_ids)
        
        # ---------- Phase 1: Per-ID Metadata & Stata IC ----------
        meta_map = {}
        
        for gid in all_ids:
            g = df[df[idv] == gid].sort_values(tv)
            y = g[self.y_var].to_numpy()
            X = g[self.x_vars].to_numpy()
            time_vec = g[tv].to_numpy()
            Ti = len(y)

            # ic = ln(RSS / (Ti - p - q - 1)) + 2*(p + q + constant + trend + 1) / (Ti - maxlag - maxlead)
            model, resid, blag, blead, valid_F, RHS_F, dy_vec = self._get_optimal_model(
                y, X, time_vec, null_model=True, return_full_model=True
            )
            
            # Centering residuals: Stata's 'replace e = e - meane'
            e_raw = resid.astype(float)
            valid_indices = np.where(~np.isnan(e_raw))[0]
            e_mean = np.mean(e_raw[valid_indices])
            e_centered = e_raw - e_mean 
            
            # Centering dX: Stata's 'egen meandx = mean(dx)'
            dx_full = np.diff(X, axis=0, prepend=np.full((1, self.nox), np.nan))
            dx_mean = np.nanmean(dx_full, axis=0)
            cdX_full = dx_full - dx_mean
            
            coeffs = model.params
            c_idx = int(self.constant) + int(self.trend)
            
            phi = np.zeros(self.maxlag + 1)
            if blag > 0:
                phi[1:blag+1] = coeffs[c_idx : c_idx + blag]
                c_idx += blag
                
            bFdx = np.zeros((self.maxlead + 1, self.nox))
            bLdx = np.zeros((self.maxlag + 1, self.nox))
            
            for k in range(-blead, blag + 1):
                block = coeffs[c_idx : c_idx + self.nox]
                c_idx += self.nox
                if k < 0:
                    bFdx[-k, :] = block
                else:
                    bLdx[k, :] = block
                        
            meta_map[gid] = {
                "Ti": Ti, "phi": phi, "bFdx": bFdx, "bLdx": bLdx,
                "e": e_centered, "cdX": cdX_full, "blead": blead, "blag": blag
            }

        results_dist = {k: [] for k in ["Gt", "Ga", "Pt", "Pa"]}
        
        # Available Innovation Pool (where e is not missing)
        # Stata: bsample if e < ., cluster(t)
        # We take the indices where the first ID has valid residuals
        U_indices = valid_indices
        U_len = len(U_indices)

        # ---------- Phase 2: The Bootstrap Loop ----------
        for b in range(self.bootstrap):
            # 1. Cluster-Bootstrap Time: bsample ..., cluster(t)
            t_draw = rng.choice(U_indices, size=U_len, replace=True)
            
            # 2. Replicate 'expandcl 2' logic
            # Stata duplicates observations to allow for the 'newtt' shuffle
            pool_expanded = np.repeat(t_draw, 2)
            
            # 3. The "Newtt" Key Generation (The Core of Cross-Sectional Correlation)
            # Stata: by tussent: egen newtt = mean(newttt)
            newttt = rng.uniform(0, 1, size=len(pool_expanded))
            
            # Stable sort to get the permutation
            perm = np.argsort(newttt, kind="mergesort")
            shuffled_t_idx = pool_expanded[perm]
            
            boot_data_list = []
            for gid in all_ids:
                meta = meta_map[gid]
                # Stata: keep if _n <= ti + maxlag + maxlead + 2
                required_len = int(meta["Ti"] + self.maxlag + self.maxlead + 2)
                idx = shuffled_t_idx[:required_len]
                
                e_b = meta["e"][idx]
                cdX_b = meta["cdX"][idx, :]
                
                # 1. Start with residuals (contains NaNs from support/truncation)
                u = e_b.copy()
                
                # 2. Accumulate components (NaNs will propagate)
                for k in range(0, meta["blag"] + 1):
                    Lk = self._roll_nan(cdX_b, k)
                    # Note: (Lk * beta) will be NaN if Lk is NaN
                    u += (Lk * meta["bLdx"][k, :]).sum(axis=1)
                
                if meta["blead"] > 0:
                    for k in range(1, meta["blead"] + 1):
                        Fk = self._roll_nan(cdX_b, -k)
                        u += (Fk * meta["bFdx"][k, :]).sum(axis=1)
                    
                # 3. Apply Stata truncation: replace u = . if _n < maxlag + 1
                u[:self.maxlag] = np.nan
                
                # 4. Convert to zero for the recursion
                u = np.nan_to_num(u, nan=0.0)
                    
                # --- Recursion dy (Stata: command string loop) ---
                dy = np.zeros(len(u))
                phi = meta["phi"]
                for t in range(self.maxlag, len(dy)):
                    ar_term = 0.0
                    for k in range(1, self.maxlag + 1):
                        ar_term += dy[t-k] * phi[k]
                    dy[t] = u[t] + ar_term
                
                # --- Integration to Levels (Stata: sum) ---
                booty = np.cumsum(dy)
                bootx = np.cumsum(np.nan_to_num(cdX_b, nan=0.0), axis=0)
                
                # Masking based on Stata's _n constraints
                mask = np.ones(len(booty), dtype=bool)
                mask[:self.maxlag] = False
                mask[int(meta["Ti"] + self.maxlag):] = False
                
                df_b = pd.DataFrame({
                    idv: gid,
                    tv: np.arange(1, len(booty) + 1),
                    self.y_var: np.where(mask, booty, np.nan)
                })
                for j, xname in enumerate(self.x_vars):
                    df_b[xname] = np.where(mask, bootx[:, j], np.nan)
                
                boot_data_list.append(df_b.dropna(subset=[self.y_var]))

            # Calculate statistics for this replication
            boot_panel = pd.concat(boot_data_list).reset_index(drop=True)
            stats_b, indiv_stats_b = self._run_westerlund_plain(boot_panel)
            for k in results_dist:
                results_dist[k].append(stats_b[k])
                    
        return results_dist

        
    def run(self):
        """
        Executes the full Westerlund (2007) ECM-based panel cointegration testing 
        workflow, including data cleaning, primary estimation, optional bootstrapping, 
        and results aggregation.

        This is the primary user-facing method. It coordinates the transition from 
        raw data to standardized test statistics and robust p-values by orchestrating 
        the internal estimation and simulation engines.

        Returns
        -------
        dict
            A comprehensive results bundle containing:
            - **test_stats**: A dictionary of the four raw observed statistics (Gt, Ga, Pt, Pa).
            - **boot_pvals**: Empirical p-values derived from the bootstrap distribution (if enabled).
            - **boot_distributions**: The full set of simulated statistics for each replication.
            - **unit_data**: A pandas DataFrame containing unit-specific estimates (alpha, beta, lags, leads).
            - **mean_group**: Aggregated coefficients for the long-run relationship and speed of adjustment.
            - **metadata**: Information on group counts, average observations, and model specifications.

        Process
        -------
        **1. Data Preparation:**
        Calls `_tsset_and_clean` to enforce time-series continuity and verify that 
        every cross-sectional unit has sufficient observations to support the 
        requested lag/lead structure.

        **2. Primary Estimation:**
        Executes `_run_westerlund_plain` to calculate the observed statistics. 
        This step estimates unit-specific Error Correction Models (ECMs) to derive 
        Mean-Group statistics and performs partialling-out regressions for the 
        pooled Panel statistics.

        **3. Bootstrap Inference (Optional):**
        If `bootstrap > 0`, the method triggers `_bootstrap_run`. It then 
        calculates "robust" p-values by comparing the observed statistics 
        against the simulated null distribution using the finite-sample correction:
        $$p = \frac{r + 1}{B + 1}$$
        where $r$ is the number of bootstrap replicates less than or equal to the 
        observed statistic, and $B$ is the number of valid replications.

        

        **4. Results Bundling:**
        Extracts and organizes granular unit-level data into a structured 
        format. This includes calculating the Mean-Group $\beta$ coefficients 
        which represent the average long-run equilibrium relationship across the panel.

        **5. Reporting:**
        Calls `_display_final` to print the formatted test results, standardized 
        Z-scores, and p-values to the console, following the standard output 
        conventions of econometric software.

        Notes
        -----
        - The method automatically handles the calculation of the "effective" 
          sample size and bandwidth parameters used in the standardization process.
        - By returning the `results_bundle`, it allows for further programmatic 
          analysis, such as plotting the bootstrap distributions or conducting 
          post-estimation hypothesis tests on the $\beta$ coefficients.
        """
        # 1. Primary Estimation
        df = self._tsset_and_clean()
        # Ensure _run_westerlund_plain returns both stats AND indiv_stats
        main_stats, indiv_data = self._run_westerlund_plain(df, is_boot=False)
        
        # 2. Extract Unit-Level Data into a DataFrame
        # This captures alpha_i, beta_i, specific lags, leads, and tnorm for every group
        unit_results = pd.DataFrame([{
            'id': s['gid'],
            'alpha_i': s['ai'],
            'se_alpha_i': s['seai'],
            'beta_i': s['betai'], # List of betas for x_vars
            'lag': s['blag'],
            'lead': s['blead'],
            'tnorm': s['tnorm'],
            'obs': s['ti']
        } for s in indiv_data])

        # 3. Process Bootstrap if requested
        boot_dist = {}
        boot_pvals = {}
        if self.bootstrap > 0:
            print(f"Bootstrapping {self.bootstrap} replications...")
            boot_dist = self._bootstrap_run(df) # Returns {Gt: [...], Ga: [...], ...}
            
            for k in ["Gt", "Ga", "Pt", "Pa"]:
                obs = main_stats[k]
                boots = np.array(boot_dist[k])
                boots = boots[np.isfinite(boots)]
                B_valid = boots.size
                if B_valid == 0:
                    boot_pvals[k] = np.nan
                else:
                    r = np.sum(boots <= obs)
                    boot_pvals[k] = (r + 1) / (B_valid + 1)

        # 4. Bundle Everything for further processing
        # We calculate MG coefficients directly from the unit_results DataFrame
        all_betas = np.array(unit_results['beta_i'].tolist())
        
        self.results_bundle = {
            'test_stats': main_stats,
            'boot_pvals': boot_pvals,
            'boot_distributions': boot_dist,
            'unit_data': unit_results,
            'mean_group': {
                'mg_alpha': unit_results['alpha_i'].mean(),
                'mg_betas': dict(zip(self.x_vars, all_betas.mean(axis=0)))
            },
            'metadata': {
                'bandwidth': self.lrwindow,
                'n_groups': len(unit_results),
                'avg_obs': unit_results['obs'].mean(),
                'model_type': 'constant' if self.constant and not self.trend else 'trend' if self.trend else 'none'
            }
        }

        # 5. Output to Console
        self._display_final(main_stats, boot_pvals)
        
        return self.results_bundle

    def _display_final(self, results, boot_pvals):
        """
        Calculates standardized Z-scores and formats the final results table for the 
        Westerlund ECM panel cointegration tests.

        This method performs the final statistical inference by comparing raw 
        test statistics against asymptotic moments. It generates a comprehensive 
        console output that includes raw values, Z-scores, asymptotic p-values, 
        and robust bootstrap p-values (if available).

        Parameters
        ----------
        results : dict
            A dictionary containing the raw values of the four test statistics: 
            'Gt', 'Ga', 'Pt', and 'Pa'.
        boot_pvals : dict or None
            A dictionary containing the empirical p-values calculated from the 
            bootstrap distribution. If None or empty, the robust p-value column 
            will display as "-".

        Process
        -------
        **1. Deterministic and Covariate Mapping:**
        Determines the appropriate index for the asymptotic moments table based on:
        - `ridx`: The deterministic specification (0: None, 1: Constant, 2: Trend).
        - `cidx`: The number of regressors (capped at 6 as per standard 
          Westerlund tables).

        **2. Standardization:**
        Converts raw statistics to standardized Z-scores using the formula:
        $$Z = \frac{\sqrt{N}(S - \mu)}{\sqrt{\sigma^2}}$$
        The mean ($\mu$) and variance ($\sigma^2$) are retrieved from either:
        - **Lookup Tables:** Hard-coded asymptotic moments indexed by 
          deterministic case and covariate count.
        - **Westerlund Constants:** Specific constants used when the 
          `westerlund` flag is active, varying by trend inclusion.

        **3. P-value Calculation:**
        - **Asymptotic:** Calculated as the left-tail probability from a 
          standard normal distribution ($\Phi(Z)$).
        - **Robust:** Retrieved from the `boot_pvals` dictionary.

        

        **4. Formatted Output:**
        Prints a Stata-style summary header including:
        - Series names and group counts ($N$).
        - Lag, lead, and bandwidth specifications.
        - A structured table comparing all four statistics side-by-side.

        Notes
        -----
        - The lookup moments provided in the code are derived from the 
          simulations in Westerlund (2007).
        - Rejection of the null hypothesis (no cointegration) occurs for 
          large negative Z-scores or small p-values (typically $< 0.05$).
        """
        ridx = 0
        if self.constant and not self.trend: ridx = 1
        if self.constant and self.trend: ridx = 2
        cidx = min(self.nox - 1, 5) 
        
        moments = {
            'gt': {'mean': np.array([[-0.976,-1.382,-1.709,-1.979,-2.199,-2.426],[-1.778,-2.035,-2.233,-2.445,-2.646,-2.836],[-2.366,-2.528,-2.704,-2.864,-3.015,-3.171]]),
                   'var':  np.array([[1.082,1.098,1.049,1.058,1.035,1.041],[0.807,0.848,0.889,0.912,0.908,0.924],[0.660,0.707,0.759,0.823,0.848,0.860]])},
            'ga': {'mean': np.array([[-3.802,-5.824,-7.811,-9.879,-11.724,-13.858],[-7.142,-9.125,-10.967,-12.956,-14.975,-17.067],[-12.012,-13.632,-15.526,-17.365,-19.253,-21.248]]),
                   'var':  np.array([[20.687,29.902,39.011,50.574,58.960,69.597],[29.634,39.343,49.488,58.704,67.950,79.109],[46.242,53.743,64.559,74.740,84.799,94.002]])},
            'pt': {'mean': np.array([[-0.511,-0.937,-1.317,-1.617,-1.882,-2.126],[-1.448,-1.713,-1.921,-2.148,-2.373,-2.577],[-2.112,-2.288,-2.463,-2.628,-2.786,-2.954]]),
                   'var':  np.array([[1.362,1.766,1.718,1.605,1.494,1.424],[0.989,1.066,1.117,1.174,1.168,1.159],[0.765,0.814,0.886,0.999,0.992,0.990]])},
            'pa': {'mean': np.array([[-1.026,-2.499,-4.270,-6.114,-8.032,-10.007],[-4.230,-5.865,-7.460,-9.306,-11.315,-13.318],[-8.933,-10.487,-12.167,-13.889,-15.682,-17.652]]),
                   'var':  np.array([[8.383,24.022,39.883,53.452,63.241,76.676],[19.709,31.264,42.998,57.484,69.437,81.038],[37.595,45.689,57.999,74.126,81.393,91.239]])}
        }

        N = len(self.raw_data[self.id_var].unique())
        sqrt_N = np.sqrt(N)
        
        print("\n" + "="*75)
        print(f"Westerlund ECM Panel Cointegration Tests")
        print(f"Series: {self.y_var} ~ {', '.join(self.x_vars)}")
        print(f"N (Groups): {N}")
        print(f"Lags: {self.minlag}-{self.maxlag} | Leads: {self.minlead}-{self.maxlead} | Window: {self.lrwindow}")
        print("="*75)
        print(f"{'Statistic':<10} | {'Value':<10} | {'Z-score':<10} | {'P-value':<10} | {'Robust P':<10}")
        print("-" * 75)

        for name in ['Gt', 'Ga', 'Pt', 'Pa']:
            val = results[name]
            k = name.lower()
            if self.westerlund:
                if not self.trend:
                    params = {'gt':(-1.793, 0.7904), 'ga':(-7.2014, 29.3677), 'pt':(-1.4746, 1.0262), 'pa':(-4.3559, 21.0535)}
                else:
                    params = {'gt':(-2.356, 0.6450), 'ga':(-11.8978, 44.2471), 'pt':(-2.1128, 0.7371), 'pa':(-8.9536, 35.6802)}
                mu, var = params[k]
                z = (val - sqrt_N * mu) / np.sqrt(var) if k=='pt' else (sqrt_N * val - sqrt_N * mu) / np.sqrt(var)
            else:
                mu = moments[k]['mean'][ridx, cidx]
                var = moments[k]['var'][ridx, cidx]
                z = (val - sqrt_N * mu) / np.sqrt(var) if k=='pt' else (sqrt_N * val - sqrt_N * mu) / np.sqrt(var)
            
            pval = stats.norm.cdf(z)
            rob_p = f"{boot_pvals[name]:.3f}" if boot_pvals else "-"
            print(f"{name:<10} | {val:<10.3f} | {z:<10.3f} | {pval:<10.3f} | {rob_p:<10}")
        print("="*75)
    
    def _reg_display(self, params, bse, tvalues, pvalues, title="Regression Results"):
        """
        Formats and prints a detailed coefficient table for unit-specific regressions.

        This method is primarily used when the `noisily` flag is active to provide 
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
        print(f"\n{title}")
        print("-" * 60)
        print(results_df.to_string())
        print("-" * 60)
    
    def _mg_display(self, mg_results, lr_results, auto=False):
        """
        Formats and prints the aggregated Mean-Group (MG) results, replicating the 
        reporting structure of Stata's 'xtwest' MG output.

        This method displays the average speed of adjustment (the error correction term) 
        and the average long-run equilibrium relationship ($\beta$ coefficients) across 
        all cross-sectional units. It provides a summary of the economic relationship 
        estimated by the panel model.

        Parameters
        ----------
        mg_results : dict
            A dictionary containing the aggregated statistics for the error 
            correction coefficient ($\alpha$). Expected keys include 'Variable', 
            'Coef.', 'Std. Err.', 't', and 'P>|t|'.
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
        print("\n" + "="*60)
        print("Mean-group error-correction model")
        if auto:
            print("Short run coefficients apart from the error-correction term are omitted")
            print("as lag and lengths might differ between cross-sectional units.")
        
        # Display the Error Correction / MG results
        print("\nError Correction Term (Alpha):")
        print(pd.DataFrame([mg_results]).to_string(index=False))
        
        print("\nEstimated long-run relationship and short run adjustment")
        # This matches the second 'eret disp' in Stata
        lr_df = pd.DataFrame(lr_results)
        print(lr_df.to_string(index=False))
        print("="*60 + "\n")
    
    def plot_bootstrap(self, 
                       title="Westerlund Panel Cointegration Test (Bootstrap)",
                       save_path=None, 
                       dpi=300,
                       figsize=(12, 10),
                       colors={'obs': '#D55E00', 'cv': '#0072B2', 'kde': 'grey'},
                       alpha=0.5):
        """
        Visualizes the bootstrap distributions of the Westerlund test statistics using 
        Kernel Density Estimation (KDE).

        This method generates a 2x2 grid of subplots for Gt, Ga, Pt, and Pa. Each plot 
        compares the simulated distribution under the null hypothesis (no cointegration) 
        against the observed test statistic. It serves as a diagnostic tool to assess 
        the significance and stability of the results.

        Parameters
        ----------
        title : str, optional
            The main title of the figure. Defaults to "Westerlund Panel Cointegration 
            Test (Bootstrap)".
        save_path : str, optional
            The filesystem path (e.g., 'path/to/plot.png') where the figure should 
            be saved. If None, the plot is only displayed in the console/notebook.
        dpi : int, optional
            The resolution of the saved image in dots per inch. Defaults to 300.
        figsize : tuple, optional
            The width and height of the figure in inches. Defaults to (12, 10).
        colors : dict, optional
            A dictionary defining the hex colors for different plot elements:
            - 'obs': The vertical line for the observed test statistic.
            - 'cv': The dashed vertical line for the 5% critical value.
            - 'kde': The color of the density curve and fill area.
        alpha : float, optional
            The transparency level (0 to 1) for the KDE shaded fill area. 
            Defaults to 0.5.

        Process
        -------
        1.  **Data Retrieval:** Extracts the bootstrap distributions and observed 
            statistics from the `results_bundle`. It requires that `.run()` has been 
            called with a positive bootstrap value.
        2.  **KDE Calculation:** Uses `scipy.stats.gaussian_kde` to compute the 
            probability density function for each of the four statistics.
        3.  **Critical Value Identification:** Calculates the 5th percentile of each 
            distribution to determine the lower-tail critical value.
        4.  **Statistical Overlay:** - Draws a solid vertical line for the **Observed Statistic**.
            - Draws a dashed vertical line for the **5% Critical Value**.
            - Displays a text box containing the observed value, critical value, 
              and the calculated robust p-value.
        5.  **Layout & Export:** Standardizes the formatting with a unified legend 
            and saves the output to the specified `save_path` if provided.

        

        Interpretation
        --------------
        The Westerlund test is a lower-tail test. Cointegration is suggested if the 
        **Observed Statistic** (solid line) is located to the **LEFT** of the 
        **Critical Value** (dashed line). In such cases, the robust p-value will 
        typically be less than 0.05.

        Notes
        -----
        - Non-finite (NaN/Inf) bootstrap draws are automatically filtered out 
          before plotting.
        - The x-axis range is automatically expanded by one standard deviation 
          to ensure the tails of the distribution are visible.
        """
        
        # 1. Check if results exist
        if not hasattr(self, 'results_bundle') or not self.results_bundle.get('boot_distributions'):
            print("Error: No bootstrap results found. Please call .run() with bootstrap > 0 first.")
            return

        # 2. Extract Data from Bundle
        res = self.results_bundle
        boot_dist = res['boot_distributions']
        obs_stats = res['test_stats']
        p_vals = res['boot_pvals']
        stats_to_plot = ['Gt', 'Ga', 'Pt', 'Pa']

        # 3. Setup Plot
        sns.set_theme(style="whitegrid")
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()

        for i, stat in enumerate(stats_to_plot):
            ax = axes[i]
            
            # Clean data (Pure NumPy)
            data = np.array(boot_dist[stat])
            data = data[np.isfinite(data)]
            
            if len(data) > 1:
                # --- Manual KDE Logic ---
                kde = gaussian_kde(data)
                # Expand range slightly for better visual tail representation
                x_range = np.linspace(data.min() - np.std(data), data.max() + np.std(data), 200)
                y_kde = kde(x_range)
                
                ax.plot(x_range, y_kde, color=colors.get('kde', 'grey'), lw=1.5, label='Bootstrap H0 Dist.')
                ax.fill_between(x_range, y_kde, color=colors.get('kde', 'grey'), alpha=alpha)

                # Observed Statistic (Red Solid)
                obs_val = obs_stats[stat]
                ax.axvline(obs_val, color=colors.get('obs', '#D55E00'), linestyle='-', linewidth=2.5, label='Observed Stat')

                # 5% Critical Value (Blue Dashed)
                # Westerlund tests are typically lower-tail tests (reject if stat < CV)
                cv_val = np.percentile(data, 5)
                ax.axvline(cv_val, color=colors.get('cv', '#0072B2'), linestyle='--', linewidth=1.5, label='5% Crit Value')

                # Display Stats
                p_val = p_vals.get(stat, np.nan)
                ax.set_title(f"Statistic: {stat}", fontweight='bold')
                ax.text(0.05, 0.95, f"Obs: {obs_val:.3f}\nCV: {cv_val:.3f}\nRobust p: {p_val:.3f}", 
                        transform=ax.transAxes, verticalalignment='top', fontsize=9,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # 4. Global Formatting
        plt.suptitle(title, fontsize=16, fontweight='bold')
        plt.figtext(0.5, 0.93, "Null Rejected if Observed (solid) is to the LEFT of Critical Value (dashed)", 
                    ha="center", fontsize=10, color="dimgrey")
        
        # Unified Legend
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True)
        
        plt.tight_layout(rect=[0, 0.05, 1, 0.92])

        # 5. Save Logic
        if save_path:
            try:
                # Create directory if it doesn't exist
                directory = os.path.dirname(save_path)
                if directory and not os.path.exists(directory):
                    os.makedirs(directory)
                
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
                print(f"Figure successfully saved to: {save_path}")
            except Exception as e:
                print(f"Warning: Could not save figure. Error: {e}")

        plt.show()
