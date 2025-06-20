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
fdr_tsbky: Two-stage BY is so conservative that it's almost never used in practice.

Therefore, in RECUR we use only provide for the following methods.
"""


METHODS_AVAILABLE = {
    "bonferroni": "Bonferroni (FWER, one-step)",
    "holm": "Holm-Bonferroni (FWER, step-down)",
    "fdr_bh": "Benjamini-Hochberg (FDR)", # Most recommended, should be used by default
    "fdr_tsbh": "Storey 2-stage BH (adaptive FDR)",
    "fdr_by": "Benjamini-Yekutieli (FDR, correlated safe)"
}

METHODS_FWER = {"bonferroni", "holm"}
METHODS_FDR  = {"fdr_bh", "fdr_by", "fdr_tsbh"}   # alias for Storey

def mc_pvalues(R: int, B: int) -> float:
    
    """
    The Phipson & Smyth method is used for adjustment
    R = #simulated statistics >= observed

    The +1 in numerator and denominator is a bias correction 
    (ensures p_hat >0) and is standard (see Davison & Hinkley 1997).

    """
    p_hat = (R + 1) / (B + 1)

    return p_hat


def adjusted_pvalue_test(pvals, alpha, method: str = "bonferroni"):

    """
    alpha: family-wise significance level (0.05, 0.1, 0.001)
    M = len(pvals): test size
    H0i: the i-th null hypothesis 
    q: the target False Discovery Rate (FDR) level you wish to control.
       It bounds the expected proportion of mistakes rather than the chance of any mistake
       It's the the FDR-analogue of the familiar family-wise significance level alpha.
       (0.05, 0.1, 0.001)

    M0: true null
    M1: true alternative, M1 = M- M0
    pi0 = M0 / M: 
    pi1 = M1 / M: 

    1. Family-Wise Error Rate (FWER):
    Probability of making even one false positive in the entire family of tests.

        A. Bonferroni (classic)
        Reject the null hypothesis if p-value <= alpha / M
        Very conservative; power drops quickly asMgrows.
        You only have a handful of tests or a false positive is very costly (e.g., clinical trial safety).
        
        B. Holm-Bonferroni (step-down Bonferroni)
        sorted_idx = np.argsort(pvals)
        p_sorted = pvals[sorted_idx]

        reject = np.zeros(M, dtype=bool)
        for i in range(M):
            threshold = alpha / (M - i)
            if p_sorted[i] > threshold:
                break
            reject[i] = True

        reject_final = np.zeros(M, dtype=bool)
        reject_final[sorted_idx[:np.sum(reject)]] = True

        Stop at the first non-significant, declare none beyond it. 
        Uniformly more powerful than plain Bonferroni; still works under any dependence.
        Still conservative for large m; sequential procedure slightly more work.
        Use it when you need FWER control but want a power boost over Bonferroni.

        Early (small) p-values are compared to a stricter threshold, 
        but as you move down the list the denominator shrinks, 
        so later tests get looser cut-offs.

    2. 	False Discovery Rate (FDR):
    Expected proportion of false positives among the discoveries you call significant.

        A. Benjamini-Hochberg (BH, "step-up FDR")
        sorted_idx = np.argsort(pvals)
        p_sorted = pvals[sorted_idx]
        thresholds = np.arange(1, M + 1) / M * q

        below = p_sorted <= thresholds
        if not np.any(below):
            return np.zeros(M, dtype=bool) 

        k_max = np.max(np.where(below)[0])
        cutoff_p = p_sorted[k_max]

        Much higher power for large-scale studies; simple to implement. 
        Formal FDR guarantee assumes independent or positively dependent tests (still robust in practice). 
        Use it anytimeMis big and you can tolerate a few false calls.
        
        B. Storey's q-value method (adaptive BH)
        thresholds = np.arange(1, M + 1) / M* (q / pi0) 
        For each test, the minimum FDR at which it becomes significant. 
        Gains extra power when many tests are truly null (common in genomics). 
        Provides a per-test q-value you can treat like an FDR-based p-value. 
        Slightly more complex.
        pi0 can be unstable in small samples.
        Same dependence caveats as BH.

    FWER control is stricter. it tries to keep any false positive from slipping through.
    FDR control relaxes that to "some mistakes are fine, as long as the proportion is small," which gives you more power (i.e., you keep more true signals).
    
    - Need almost zero risk of any false positive?
      Bonferroni (small M) or Holm (medium-large M).

    - Large-scale study where a few false hits are acceptable if you keep lots of true ones?
      BH 

      if you believe the vast majority of tests are null and want maximum power,
      use Storey's q-values.

    M <= 50: Bonferroni
    You need to avoid any false positive at all costs.

    50 < M <= 500: Holm
    You still want strong FWER control, but Bonferroni is too harsh.

    500 < M <= 5000: BH
    You're running many tests (e.g. microarrays), and FDR is acceptable.

    M > 5000 and you believe most hypothesis are null: Storey's q-values
    You want maximum power and don't mind a small number of false positives.
    """

    reject, padj, *_ = multipletests(pvals, alpha, method=method)

    return reject, padj

# def p_value_confidence_interval():

# def p_value_relative_error():


def grid_resolution(B: int):
    return 1 / (B + 1)

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
    lower = beta.ppf(alpha/2, R, B - R + 1)
    upper = beta.ppf(1-alpha/2, R+1, B - R)

    return  lower, upper

def sitewise_decision(R, B, alpha=0.05, q=0.05, method="fdr_bh"):
    """
    Parameters
    ----------
    R_vec  : array-like  counts of exceedances (length M)
    B      : int         number of Monte-Carlo draws
    alpha  : float       FWER threshold for Bonferroni / Holm
    q      : float       FDR threshold for BH / Storey
    method : {"bonferroni","holm","fdr_bh","fdr_tsbh"}

    Returns
    -------
    dict  with keys p_hat, p_adj, ci_lo_adj, ci_hi_adj,
          decision, robust, unstable

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
        p_adj = np.minimum(1, M * p_hat)
        ci_lo_adj = np.minimum(1, M * ci_lo)
        ci_hi_adj = np.minimum(1, M * ci_hi)
        thresh = alpha
    else:
        thresh = q 
        _, p_adj, *_ = multipletests(p_hat, thresh, method=method)
        _, ci_lo_adj, *_ = multipletests(ci_lo, thresh, method=method)
        _, ci_hi_adj, *_ = multipletests(ci_hi, thresh, method=method)

    decision_sig = p_adj <= thresh  # point-estimate verdict

    robust_sig = ci_hi_adj < thresh  # always significant
    robust_nonsig= ci_lo_adj > thresh    # always non-significant
    robust = robust_sig | robust_nonsig
    # unstable = ~(robust_sig | robust_nonsig)   # CI crosses threshold

    return p_hat, p_adj, ci_lo_adj, ci_hi_adj, decision_sig, robust


def method_selection(M: int, suspect_dependence: bool = False) -> str:

    if M <= 50:
        selected = "bonferroni"
    elif M <= 500:
        selected = "holm"
    elif M <= 5000:
        selected = "fdr_by" if suspect_dependence else "fdr_bh"
    else:
        # > 5 k tests, adaptive FDR gives best power
        selected = "fdr_tsbh" if not suspect_dependence else "fdr_by"

    return selected

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
        # return_details=False,
    ):
    
    """
    Parameters
    ----------
    M         : int      total hypotheses
    method     : str      {"bonferroni","holm","bh","storey"}
    alpha      : float    family-wise α (Bonferroni/Holm)
    q          : float    target FDR q  (BH/Storey)
    pi0        : float    π₀ estimate   (Storey; default 1 = worst-case)
    p_expected : float    p-value level you want the MC-SE control for
    eps        : float    maximum allowed Monte-Carlo SE at p_expected
    return_details : bool if True, also return the two individual B's
    extra_grid_cushion : bool if True, 
    
    Returns
    -------
    B_required : int                minimum simulations
    details    : dict  (optional)   {"grid_B":…, "error_B":…, "threshold":…}
    
    ===========================================================================================
    Importance concepts:
    - M:
        It is used in statistical correction (how strict the thresholds are)
        It sets how small a p-value you might ever need to resolve.
        For example, for Bonferroni / Holm: smallest critical value is α / M
    - B:
      It controls the numerical accuracy of every Monte-Carlo p-value (how finely you can measure p).
      If B is too small, you cannot reach below α / M (grid too coarse), 
      it will not make SE small enough that borderline decisions are stable

    ----------------------------------------------------------------------------------------
    Minimum number of Monte-Carlo simulations (B) needed for
    1. Resolution:  1/(B+1) ≤ smallest adjusted-p threshold

       For Bonferroni:
       smallest adjusted-p threshold = α / M

       1/(B+1) represents the smallest possible nonzero p-value 
       you can obtain from a Monte Carlo hypothesis test based on B simulated null samples.
       This is also called the grid resolution of the Monte Carlo p-values.
       
       Define R = #{T* >= T_obs}, where T* is simulated recurrence, T_obs is the observed recurrence.
       It says that number of simulated test statistics ≥ observed statistic
       In Monte-Carlo test, p-values are calculated as 
       p_hat = (R+1) / (B+1)

       So your p-values can only take values on the grid
       {1 / (B+1), 1 / (B+2), ..., 1}
       This limits the precision of your significance testing,
       particularly critical in multiple-testing scenarios like Bonferroni or BH where adjusted thresholds
       (e.g. α/M) may be very small

       B >= M/α - 1 

    2. MC error:  sqrt(p(1-p)/(B+1) ≤ eps (optional)
       Since we are trying to estimate a true (ideal) p-value by sampling,
       and the count R is random, which in fact, follows a Binomial distribution

       R ~ Binomial(B, p)
       so the estimate p_hat = (R+1) / (B+1) is is a function of a Binomial random variable, 
       and we can approximate its standard error using the delta method or basic variance propagation.

       Since Var(R) = B * p * (1-p)
       Var(p_hat) = Var((R+1) / (B+1))
                  = Var(R / (B+1) + 1 / (B+1))
                  = Var(R) / (B+1)**2
                  = B * p * (1-p) / (B+1)**2
                  = (1 / (B+1)) * (1 - 1 / (B+1)) * p * (1-p)
                  ≈ p * (1-p) / (B+1)  where 1 / (B+1) ≪ 1

       SE(p_hat) = sqrt(Var(p_hat)) 
       SE (standard error) tells you how much the Monte-Carlo p-value fluctuates from run to run.
       If that noise is a big chunk of the p-value itself, 
       a site can slide above or below the threshold just by RNG luck.

       We want Monte-Carlo noise to be small compared with the p-value itself,
       so that numerical jitter doesn't flip "significant" to "not significant". 

       A common criterion
       SE(p_hat) / p_hat <= 0.1, i.e. sqrt(p * (1-p) / (B+1)) <= 0.1 * p, eps = 0.1 * p

       Solve for B, we can get:
       B >= (1-p) / (0.01*p) - 1

       This tells you how large B needs to be to keep Monte Carlo noise within 10% of the estimate.
       A 10 % relative SE means the 95 % CI is roughly ±20 % of p. 
       It's tight enough that the decision is unlikely to flip, yet still affordable in CPU time.

       Here, the 10% is just a rule-of-thumb, not a law. 
       You can also set SE / p ≤ 5% or 20% depending on your tolerance.

       p_hat ~ approximately normal(p, p*(1-p)/(B + 1))
       Var(p_hat) = p_hat * (1 - p_hat) / B
       For a normal distribution, a 95% confidence interval is:
       p_hat ± z_0.975 * SE = p_hat ± sqrt(p_hat * (1 - p_hat) / (B + 1))
       when SE = 0.1 * p_hat, the 95% confidence interval becomes
       CI width = z_0.975 * SE = 1.96 * 0.1 * p_hat = 0.196 * p_hat ≈ 0.2 * p_hat
       Similarly, when the relative tolerence is 0.05, 
       CI width = 1.96 * 0.05 * p_hat ≈ 0.1 * p_hat 

    -------------------------------------------------------------------------------------    
    alpha: family-wise significance level (0.05, 0.1, 0.001)
    M = len(pvals): test size
    H0i: the i-th null hypothesis 
    q: the target False Discovery Rate (FDR) level you wish to control.
       It bounds the expected proportion of mistakes rather than the chance of any mistake
       It's the the FDR-analogue of the familiar family-wise significance level alpha.
       (0.05, 0.1, 0.001)

    M0: true null
    M1: true alternative, M1 = M- M0
    pi0 = π₀ = M0 / M: the fraction of hypotheses whose null model is true.
    pi1 = M1 / M = 1 - pi0
    π₀ = 1: every test is null (no signal). 
            very conservative, meaning you have no idea about the test beforehand, 
            generate large B.
    π₀ ≈ 0: almost everything is a true effect (rare in practice).

    1. Family-Wise Error Rate (FWER):
        Probability of making even one false positive in the entire family of tests.

        A. Bonferroni (classic)
        Reject the null hypothesis if p-value <= alpha / M
        Very conservative; power drops quickly asMgrows.
        You only have a handful of tests or a false positive is very costly (e.g., clinical trial safety).
        
        B. Holm-Bonferroni (step-down Bonferroni)
        sorted_idx = np.argsort(pvals)
        p_sorted = pvals[sorted_idx]

        reject = np.zeros(M, dtype=bool)
        for i in range(M):
            threshold = alpha / (M - i)
            if p_sorted[i] > threshold:
                break
            reject[i] = True

        reject_final = np.zeros(M, dtype=bool)
        reject_final[sorted_idx[:np.sum(reject)]] = True

        Stop at the first non-significant, declare none beyond it. 
        Uniformly more powerful than plain Bonferroni; still works under any dependence.
        Still conservative for large m; sequential procedure slightly more work.
        Use it when you need FWER control but want a power boost over Bonferroni.

        Early (small) p-values are compared to a stricter threshold, 
        but as you move down the list the denominator shrinks, 
        so later tests get looser cut-offs.

    2. 	False Discovery Rate (FDR):
        Expected proportion of false positives among the discoveries you call significant.

        A. Benjamini-Hochberg (BH, "step-up FDR")
        sorted_idx = np.argsort(pvals)
        p_sorted = pvals[sorted_idx]
        thresholds = np.arange(1, M + 1) / M * q

        below = p_sorted <= thresholds
        if not np.any(below):
            return np.zeros(M, dtype=bool) 

        k_max = np.max(np.where(below)[0])
        cutoff_p = p_sorted[k_max]

        Much higher power for large-scale studies; simple to implement. 
        Formal FDR guarantee assumes independent or positively dependent tests (still robust in practice). 
        Use it anytimeMis big and you can tolerate a few false calls.

        Rule of thumb for choosing q:
            If a single false positive is costly but not catastrophic, stick with 0.05;
            if discovery is cheap and you'd rather miss nothing, relax to 0.10;
            if each follow-up is expensive or safety-critical, tighten to 0.01.
        By default, we set q = α = 0.05 in RECUR.
        
        B. Storey's q-value method (adaptive BH)
        thresholds = np.arange(1, M + 1) / M* (q / pi0) 
        For each test, the minimum FDR at which it becomes significant. 
        Gains extra power when many tests are truly null (common in genomics). 
        Provides a per-test q-value you can treat like an FDR-based p-value. 
        Slightly more complex.
        pi0 can be unstable in small samples.
        Same dependence caveats as BH.

        C. Benjamini-Yekutieli (BY)
        order = np.argsort(pvals)
        p_sorted = pvals[order]

        c_m = np.sum(1.0 / np.arange(1, M + 1))
        threshold = (np.arange(1, M + 1) / M) * (q / c_m)

        below = p_sorted <= threshold
        k_max = np.max(np.where(below)[0]) if np.any(below) else -1

        reject_sorted = np.zeros(m, dtype=bool)
        if k_max >= 0:
            reject_sorted[:k_max + 1] = True

        reject = np.zeros(m, dtype=bool)
        reject[order] = reject_sorted

        Benjamini-Hochberg (BH) controls the false-discovery rate (FDR) only when the test statistics are independent or positively dependent.
        When you have unknown/complex correlation (e.g. linkage blocks in a genome, long-range co-evolution among sites), 
        BH can slightly under-shoot its target q.
        BY fixes that by multiplying the BH cut-off by a harmonic-series factor
        c_m = sum([1, 1/2, ..., 1/M]) 
            ≈log(M) + γ   where γ≈0.577

        The price: the threshold line is lower by roughly log(M), so power is reduced.
        Therefore, you choose it, If 
        You suspect strong, unknown correlation across sites even after accounting for the tree (e.g. pervasive co-evolution).
        You can afford the power loss — often BY yields far fewer significant sites than BH/Storey.

    FWER control is stricter. it tries to keep any false positive from slipping through.
    FDR control relaxes that to "some mistakes are fine, as long as the proportion is small," which gives you more power (i.e., you keep more true signals).
    
    - Need almost zero risk of any false positive?
      Bonferroni (small M) or Holm (medium-large M).

    - Large-scale study where a few false hits are acceptable if you keep lots of true ones?
      BH 

      if you believe the vast majority of tests are null and want maximum power,
      use Storey's q-values.

    M <= 50: Bonferroni
    You need to avoid any false positive at all costs.

    50 < M <= 500: Holm
    You still want strong FWER control, but Bonferroni is too harsh.

    500 < M <= 5000: BH
    You're running many tests (e.g. microarrays), and FDR is acceptable.

    M > 5000 and you believe most hypothesis are null: Storey's q-values
    You want maximum power and don't mind a small number of false positives.

    """
     
    if method is not None:
        method = method.lower()
        if method not in METHODS_FWER | METHODS_FDR:
            raise ValueError("method must be bonferroni, holm, fdr_bh, "
                             "fdr_by, or fdr_tsbh")
        selected = method
    else:
        selected = method_selection(M, suspect_dependence)

    if selected in {"bonferroni", "holm"}:
        threshold = alpha / M   # one-sided
        if extra_grid_cushion:
            threshold /= 2
    elif selected == "fdr_bh":
        threshold = q / M
    elif selected == "fdr_tsbh":

        threshold = q / (pi0 * M) # pi0 from pilot or worst-case 1.0
    else:   # "fdr_by"
        c_m = np.sum(1.0 / np.arange(1, M + 1)) # harmonic sum
        threshold = q / (M * c_m)

    # grid-resolution
    B_grid = math.ceil(1 / threshold) - 1

    # Monte-Carlo standard-error requirement (optional)
    if p_expected is None:
        p_expected = threshold 
    
    if eps is None:                  
        eps = rel_tol * p_expected 

    B_error = math.ceil(p_expected * (1 - p_expected) / eps**2) - 1
    B_error = max(B_error, 0)      # just in case

    B_required = max(B_grid, B_error)
    return selected, B_required
