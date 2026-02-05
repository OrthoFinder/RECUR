import os
import math
from typing import Dict, List, Optional, Tuple, Union, Sequence

import numpy as np 
from scipy.stats import beta
from statsmodels.stats.multitest import multipletests



"""
statsmodels.stats.multitest.multipletests provides the following methods:

bonferroni : one-step correction. 
sidak : one-step correction. 
holm-sidak : step down method using Sidak adjustments
holm : step-down method using Bonferroni adjustments
simes-hochberg : step-up method (independent)
hommel : closed method based on Simes tests (non-negative)
fdr_bh : Benjamini/Hochberg (non-negative)
fdr_by : Benjamini/Yekutieli (negative)
fdr_tsbh : two stage fdr correction (non-negative)
fdr_tsbky : two stage fdr correction (non-negative)

Since
sidak / holm-sidak: Only beats Bonferroni/Holm if tests are strictly independent; the power gain is tiny, so not worth clutter.
simes-hochberg / hommel: Optimal for ≤ 20 tests; computationally and conceptually overkill for thousands.
fdr_by: Two-stage BY is so conservative that it's almost never used in practice.

Therefore, in RECUR we use only provide for the following methods.

bonferroni (FWER): bullet-proof, simplest, most conservative. 
                   Use only if you truly need “no false positives at all.”
holm (FWER): same guarantee as Bonferroni but strictly more powerful. 
             If you need FWER, prefer Holm.
fdr_bh (FDR): baseline for exploratory/many tests under independence/PRDS. 
              Simple and widely accepted.
fdr_tsbh (FDR, two-stage BH): more power when there are many non-nulls (π₀<1). 
                              Good when dependence is mild.
fdr_tsbky (FDR, two-stage BKY): like tsbh but a hair more conservative (pilot at q/(1+q) 
                                and a tiny safety factor). Safer default if dependence is moderate/uncertain.
"""


METHODS_AVAILABLE = {
    "bonferroni": "Bonferroni (FWER, one-step)",
    "holm": "Holm-Bonferroni (FWER, step-down)",
    "fdr_bh": "Benjamini-Hochberg (FDR)", # Most recommended, should be used by default
    "fdr_tsbh": "two-stage BH (adaptive FDR)",
    "fdr_tsbky": "Benjamini–Krieger–Yekutieli (FDR, correlated safe)"
}

METHODS_FWER = {"bonferroni", "holm"}
METHODS_FDR  = {"fdr_bh", "fdr_tsbky", "fdr_tsbh"}   # alias for Storey


def cp_ci_vec(R: List[int], B: int, alpha=0.05) -> Tuple[float, float]:
    """
    Exact Clopper-Pearson confidence interval (CI) for a Monte-Carlo p-value.
    
    p_hat = (R + 1) / (B + 1)
    R = #simulated statistics >= observed
    R ~ Binomial(B, p_hat) follows a binomial distribution, so we can use
    exact Clopper-Pearson or any binomial CI

    lower bound: Prob(R <= r_obs) >= alpha / 2
    upper bound: Prob(R >= r_obs) >= alpha / 2
    Solving those inequalities yields the closed-form "Beta-inverse" expressions.
     
    After you adjust p-values (Bonferroni, BH, …), 
    compare both ends of the interval to the family-wise or FDR threshold.
    If the whole interval is below the threshold, the conclusion is robust to MC error.

    """
    R = np.asarray(R)
    lower = np.where(
        R == 0,
        0.0,                                     # CI lower bound is exactly 0
        beta.ppf(alpha / 2,  R, B - R + 1),
    )

    upper = np.where(
        R == B,
        1.0,                                     # CI upper bound is exactly 1
        beta.ppf(1 - alpha / 2, R + 1, B - R),
    )

    return lower, upper

def sitewise_decision(R, B, alpha=0.05, q=0.05, method="fdr_bh"):
    """
    Perform site-wise hypothesis testing with Monte Carlo p-values,
    multiple testing correction, and confidence interval analysis.

    Parameters
    ----------
    R : array-like of int
        Vector of exceedance counts per site. Each R[i] is the number of
        Monte Carlo statistics ≥ observed at site i.

    B : int
        Number of Monte Carlo simulations. Each p-value is estimated as (R + 1) / (B + 1).

    alpha : float, default=0.05
        Family-wise error rate (FWER) threshold used for Bonferroni and Holm correction.

    q : float, default=0.05
        False discovery rate (FDR) threshold used for FDR-based methods like BH, BY, etc.

    method : str, default="fdr_bh"
        Multiple testing correction method. One of:
        {"bonferroni", "holm", "fdr_bh", "fdr_by", "fdr_tsbh"}.

    Returns
    -------
    p_hat : ndarray of float
        Monte Carlo p-value estimates: (R + 1) / (B + 1).

    p_adj : ndarray of float
        Multiple-testing adjusted p-values.

    ci_lo_adj : ndarray of float
        Adjusted lower bounds of the Clopper–Pearson confidence intervals for p_hat.

    ci_hi_adj : ndarray of float
        Adjusted upper bounds of the Clopper–Pearson confidence intervals for p_hat.

    decision_sig : ndarray of bool
        Boolean mask indicating which hypotheses are significant after adjustment.

    robust : ndarray of bool
        Boolean mask indicating which decisions are robust to Monte Carlo error
        (i.e., the entire CI is either above or below the threshold).

    The Phipson & Smyth method is used for adjustment
    R = #simulated statistics >= observed

    p_hat = (R + 1) / (B + 1) 

    The +1 in numerator and denominator is a bias correction 
    (ensures p_hat >0) and is standard (see Davison & Hinkley 1997).

    """

    method = method.lower()
    if method not in {"bonferroni", "holm", "fdr_bh", "fdr_by", "fdr_tsbh"}:
        raise ValueError("method must be bonferroni, holm, fdr_bh, fdr_by, or fdr_tsbh")

    R_vec = np.asarray(R)
    M  = len(R_vec)

    p_hat = (R_vec + 1) / (B + 1) 
    ci_lo, ci_hi = cp_ci_vec(R_vec, B, alpha)

    if method == "bonferroni":
        p_adj = np.minimum(1.0, M * p_hat)
        ci_lo_adj = np.minimum(1.0, M * ci_lo)
        ci_hi_adj = np.minimum(1.0, M * ci_hi)
        thresh = alpha
    else:
        thresh = q 
        pvals = np.clip(p_hat, 0.0, 1.0)
        lo_clipped = np.clip(ci_lo, 0.0, 1.0)
        hi_clipped = np.clip(ci_hi, 0.0, 1.0)

        _, p_adj,  *_ = multipletests(pvals, alpha=thresh, method=method)
        _, ci_lo_adj, *_ = multipletests(lo_clipped, alpha=thresh, method=method)
        _, ci_hi_adj, *_ = multipletests(hi_clipped, alpha=thresh, method=method)

    decision_sig = p_adj < thresh  # point-estimate verdict

    robust_sig = ci_hi_adj < thresh  # always significant
    robust_nonsig= ci_lo_adj > thresh    # always non-significant
    robust = robust_sig | robust_nonsig
    # unstable = ~(robust_sig | robust_nonsig)   # CI crosses threshold

    return p_hat, p_adj, ci_lo_adj, ci_hi_adj, decision_sig, robust


def method_selection(M: int, suspect_dependence: bool = False) -> str:
    """
    Heuristic chooser for multiple-testing methods (keeps Bonferroni).

    Goals:
    - Keep Bonferroni for very small M (ultra-conservative, simple).
    - Allow Holm for small/moderate M when you still want FWER but Bonferroni is too harsh.
    - For larger M, switch to FDR methods to avoid Monte Carlo simulation blow-up.
    - If dependence is suspected, prefer BKY (two-stage FDR) rather than reverting to FWER.

    Rules:
    - 1–50: Bonferroni (FWER, avoid any false positive)
    - 51–200: Holm (FWER, less harsh than Bonferroni)
    - 201–2000:
        * if suspect_dependence: BKY two-stage FDR
        * else: BH FDR
    - >2000:
        * if suspect_dependence: BKY two-stage FDR
        * else: two-stage BH (extra power)

    Note: BY is not auto-selected; use it only if worst-case arbitrary dependence
    guarantees are explicitly required (very conservative).
    """
    if not isinstance(M, int) or M < 1:
        raise ValueError("M must be a positive integer.")

    if M <= 50:
        return "bonferroni"          # strong FWER

    if M <= 200:
        return "holm"                # strong FWER, less harsh than Bonferroni

    if M <= 2000:
        return "fdr_tsbky" if suspect_dependence else "fdr_bh"   # FDR

    return "fdr_tsbky" if suspect_dependence else "fdr_tsbh"     # FDR


def min_mcs(
        M,
        method="fdr_bh",
        *,
        alpha=0.05,
        q=0.05,
        pi0=1.0, 
        p_expected=None,
        eps=None,
        rel_tol=0.1,
        extra_grid_cushion=False,
        suspect_dependence=False,
        mc_error_control=False,
    ):
    
    """
    Estimate the minimum number of Monte Carlo simulations (B) required to control
    Monte Carlo standard error (MC-SE) or to ensure grid resolution is fine enough
    for multiple-testing correction under various methods.

    Parameters
    ----------
    M : int
        Total number of hypothesis tests.

    method : str, default="fdr_bh"
        Multiple testing correction method. One of:
        {"bonferroni", "holm", "fdr_bh", "fdr_tsbky", "fdr_tsbh"}

    alpha : float, default=0.05
        Family-wise error rate (FWER) threshold. Used for "bonferroni" and "holm".

    q : float, default=0.05
        False discovery rate (FDR) threshold. Used for FDR methods.

    pi0 : float, default=1.0
        Estimate of the proportion of true null hypotheses (π₀). Used in Storey-type methods.
        Set to 1.0 for a conservative (worst-case) estimate.

    p_expected : float or None, default=None
        The p-value level where MC error should be tightly controlled.
        Required if `eps` is provided.

    eps : float or None, default=None
        Desired maximum standard error of the Monte Carlo p-value estimate at `p_expected`.
        If provided, the function will ensure MC-SE(p̂) ≤ eps.
        You must also supply `p_expected`.

    rel_tol : float, default=0.1
        Relative tolerance when determining the minimum B for grid resolution or
        MC error control. A value of 0.1 allows 10% slack over the strict minimum.

    extra_grid_cushion : bool, default=False
        If True, adds a small cushion to the estimated B to ensure stable
        multiple-testing correction thresholds (useful in threshold-sensitive procedures).

    suspect_dependence : bool, default=False
        If True, switches from "fdr_bh" to "fdr_by", or from "storey" to a conservative variant,
        to account for potentially dependent test statistics.

    mc_error_control : bool, default=False
        If True, enforces Monte Carlo standard error (MC-SE) control using `p_expected` and `eps`.

    Returns
    -------
    B_required : int
        Minimum number of Monte Carlo simulations needed to meet the grid and/or
        MC error criteria under the specified correction method.

    """
    if method is not None:
        if method not in (METHODS_FWER | METHODS_FDR):
            raise ValueError("method must be one of: "
                             "bonferroni, holm, hochberg, hommel, "
                             "sidak, holm-sidak, fdr_bh, fdr_tsbh, "
                             "fdr_tsbky, or fdr_by")
        selected = method
    else:
        selected = method_selection(M, suspect_dependence)

    if selected in {"bonferroni", "holm", "hochberg", "hommel"}:
        threshold = alpha / M    # smallest cutoff among these FWER methods

    elif selected == "sidak":
        threshold = 1.0 - (1.0 - alpha)**(1.0 / M)

    elif selected == "holm-sidak":
        threshold = 1.0 - (1.0 - alpha)**(1.0 / M)

    elif selected == "fdr_bh":
        threshold = q / M

    elif selected in ["fdr_tsbh", "fdr_tsbky"]:
        threshold = q / (pi0 * M) # pi0 from pilot or worst-case 1.0

    elif selected == "fdr_by":
        c_m = np.sum(1.0 / np.arange(1, M + 1)) # harmonic sum
        threshold = q / (M * c_m)

    if extra_grid_cushion:
        threshold /= 2.0

    # grid-resolution
    B_grid = math.ceil(1 / threshold)
    B_required = B_grid
    
    if mc_error_control:
        # Monte-Carlo standard-error requirement (optional)
        if p_expected is None:
            p_expected = threshold 
        
        if eps is None:                  
            eps = rel_tol * p_expected 

        B_error = math.ceil(p_expected * (1 - p_expected) / eps**2)
        B_error = max(B_error, 0)      # just in case

        B_required = max(B_grid, B_error)

    return selected, B_required
