from __future__ import annotations

import math
from typing import Literal, Optional, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import beta
from statsmodels.stats.multitest import multipletests


MethodFWER = Literal["bonferroni", "holm"]
MethodFDR = Literal["fdr_bh", "fdr_tsbh", "fdr_tsbky"]
Method = Union[MethodFWER, MethodFDR]


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


def cp_ci_vec(
    R: ArrayLike,
    B: int,
    alpha: float = 0.05,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
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
    Rv = np.asarray(R, dtype=np.int64)

    lower = np.where(
        Rv == 0,
        0.0,
        beta.ppf(alpha / 2.0, Rv, B - Rv + 1),
    ).astype(np.float64)

    upper = np.where(
        Rv == B,
        1.0,
        beta.ppf(1.0 - alpha / 2.0, Rv + 1, B - Rv),
    ).astype(np.float64)

    return lower, upper


def sitewise_decision(
    R: ArrayLike,
    B: int,
    alpha: float = 0.05,
    q: float = 0.05,
    method: Method = "fdr_bh",
) -> Tuple[
    NDArray[np.float64],  # p_hat
    NDArray[np.float64],  # p_adj
    NDArray[np.float64],  # ci_lo_adj
    NDArray[np.float64],  # ci_hi_adj
    NDArray[np.bool_],    # decision_sig
    NDArray[np.bool_],    # robust
]:
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

    method_lc = method.lower()
    if method_lc not in METHODS_FWER | METHODS_FDR:
        raise ValueError(f"Unsupported method: {method}")

    R_vec = np.asarray(R, dtype=np.int64)
    M = int(R_vec.size)

    p_hat = (R_vec + 1) / (B + 1)
    ci_lo, ci_hi = cp_ci_vec(R_vec, B, alpha)

    if method_lc == "bonferroni":
        p_adj = np.minimum(1.0, M * p_hat)
        ci_lo_adj = np.minimum(1.0, M * ci_lo)
        ci_hi_adj = np.minimum(1.0, M * ci_hi)
        thresh = alpha
    elif method_lc == "holm":
        # Holm is FWER with familywise alpha, not q
        thresh = alpha
        pvals = np.clip(p_hat, 0.0, 1.0)
        lo = np.clip(ci_lo, 0.0, 1.0)
        hi = np.clip(ci_hi, 0.0, 1.0)
        _, p_adj, *_ = multipletests(pvals, alpha=thresh, method="holm")
        _, ci_lo_adj, *_ = multipletests(lo, alpha=thresh, method="holm")
        _, ci_hi_adj, *_ = multipletests(hi, alpha=thresh, method="holm")
    else:
        # FDR methods use q
        thresh = q
        pvals = np.clip(p_hat, 0.0, 1.0)
        lo = np.clip(ci_lo, 0.0, 1.0)
        hi = np.clip(ci_hi, 0.0, 1.0)
        _, p_adj, *_ = multipletests(pvals, alpha=thresh, method=method_lc)
        _, ci_lo_adj, *_ = multipletests(lo, alpha=thresh, method=method_lc)
        _, ci_hi_adj, *_ = multipletests(hi, alpha=thresh, method=method_lc)

    decision_sig = (p_adj < thresh)
    robust_sig = (ci_hi_adj < thresh)
    robust_nonsig = (ci_lo_adj > thresh)
    robust = robust_sig | robust_nonsig

    # Casts for stable dtypes
    return (
        np.asarray(p_hat, dtype=np.float64),
        np.asarray(p_adj, dtype=np.float64),
        np.asarray(ci_lo_adj, dtype=np.float64),
        np.asarray(ci_hi_adj, dtype=np.float64),
        np.asarray(decision_sig, dtype=np.bool_),
        np.asarray(robust, dtype=np.bool_),
    )



def method_selection(M: int, suspect_dependence: bool = False) -> Method:
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
    M: int,
    method: Optional[Method] = "fdr_bh",
    *,
    alpha: float = 0.05,
    q: float = 0.05,
    pi0: float = 1.0,
    p_expected: Optional[float] = None,
    eps: Optional[float] = None,
    rel_tol: float = 0.1,
    extra_grid_cushion: bool = False,
    suspect_dependence: bool = False,
    mc_error_control: bool = False,
) -> Tuple[Method, int]:
    
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
    selected: str
        Selected p-value adjustment method
    B_required : int
        Minimum number of Monte Carlo simulations needed to meet the grid and/or
        MC error criteria under the specified correction method.
    

    """
    if method is None:
        selected: Method = method_selection(M, suspect_dependence)
    else:
        method_lc = method.lower()
        if method_lc not in METHODS_FWER | METHODS_FDR:
            raise ValueError(
                "method must be one of: bonferroni, holm, fdr_bh, fdr_tsbh, fdr_tsbky"
            )
        selected = method_lc  # type: ignore[assignment]

    # smallest relevant cutoff (grid resolution target)
    if selected in {"bonferroni", "holm"}:
        threshold = alpha / M
    elif selected == "fdr_bh":
        threshold = q / M
    else:  # fdr_tsbh or fdr_tsbky
        threshold = q / (pi0 * M)

    if extra_grid_cushion:
        threshold /= 2.0

    B_grid = math.ceil(1.0 / threshold)
    B_required = B_grid

    if mc_error_control:
        pe = threshold if p_expected is None else float(p_expected)
        e = (rel_tol * pe) if eps is None else float(eps)
        B_error = math.ceil(pe * (1.0 - pe) / (e * e))
        B_required = max(B_grid, B_error)

    return selected, int(B_required)