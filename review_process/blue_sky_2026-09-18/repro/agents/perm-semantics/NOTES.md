# perm-semantics probe notes (agent: perm-semantics)

TARGET: every permute() p-value and interval is the object it is documented to be.
METHOD: semantic conformance -- recompute every reported quantity from the saved null
distribution, using ONLY the help file's definitions, and compare.

## Machinery built
- `dgp.do`      self-contained staggered-DiD generator (knobs dN dT dons dnever dtau dsd dseed)
- `recheck.do`  Mata `recheck(tag, lvl, absp, exact)` -- independent recomputation of all 18
                kernel rows (n, s_obs, p, qlo, qhi, sd, lb, ub, se_shift, offcentre, bandlo,
                bandhi, live, p_lo, p_hi, se_p, cp_lo, cp_hi) from e(permdist_b) +
                e(perm_comparator_b) + e(b), compared to e(perm_kernel) exactly
                (p-values bit-exact, quantile-derived numbers rtol 1e-12).
- `harvest.do`  copies the surfaces out of e() and cross-checks
                e(perm_p)/e(perm_ci)/e(p)/e(ci_l)/e(ci_u)/e(wdidmainresults) against the kernel.

## Documented definitions being tested (sthlp 589-650, 430-545)
- coefficient rows, default: p = min(1, 2 min(P_lo,P_hi)); P_lo = (1+#{s_r<=s_obs})/(n+1)
- permabsp:  p = (1+#{|s_r|>=|s_obs|})/(n+1)          floor 1/(n+1)
- permrawp:  p = #{|s_r|>=|s_obs|}/n                  can be 0 (ruled design)
- permside(left|right): P_lo or P_hi alone            floor 1/(n+1)
- permexact: plain shares over M                      floor 1/M or 2/M, no MC error
- joint/equality/adversarial: (1+#{W_r>=W_obs})/(n+1)
- intervals: shift = b - equal-tailed quantiles (x se_analytic under permtstat);
  under permabsp/permrawp the symmetric dual b +/- q_level(|s_r|)
- permtci: b +/- q_level(|t_r|) * se_obs (NOT test inversion)
- permstepdown: max-|studentized| step-down, monotone

## HAZARD FOUND IN MY OWN HARNESS (not a wooldid defect)
`matrix A = e(nonexistent_matrix)` in Stata silently yields a 1x1 matrix of `.` rather than
erroring -- guard with `cap confirm matrix e(...)` first.
Also: Mata refuses `{ stmt; stmt }` inside a function body with "illegal arglist".

## RESULTS SO FAR
(see FINDINGS.md; running log below)

## RESULTS (running)

| probe | what | verdict |
|---|---|---|
| probe01 | 5 conventions (default / permabsp / permrawp / permside(left) / permside(right)) on one design and seed; full 18-row kernel recomputed | CLEAN, bit-exact. p=0 only under permrawp (ruled). LEFT p == p_lo, RIGHT p == p_hi, DEF p == min(1,2min). p_lo/p_hi identical across all five runs; draws identical across conventions. |
| probe02 | permtstat, permtci, level(90) | CLEAN. permtci == b +/- q95(&#124;t_r&#124;)*se_obs to 1 ulp; coincides with the permabsp table interval (reldif 3e-16). Footer says "on the t-ratio scale" explicitly. |
| probe03 | permexact (M=28) x {plain, permabsp, permanalyticse, permtstat} | CLEAN. Plain shares, floors 2/28 and 1/28 exactly, e(perm_p_se)/e(perm_p_ci) missing, no add-one fallback, comparator bitwise == e(b). |
| probe04/04b | permstepdown, post AND pre families | CLEAN. Recomputed t_obs/p_raw/p_singlestep/p_stepdown exactly; monotone; p_raw <= p_stepdown <= p_singlestep. |
| probe05 | failures: permute(probabilities) 19/200 failed; long pre-window | CLEAN. n = completed; every kernel row conformant on the reduced n. |
| probe06 | ARE THE DRAWS WHAT THEY CLAIM? rebuilt ttre from the dump's onset_ columns and refit | CLEAN. 5 of 5 replications (incl. rep 0) BIT-IDENTICAL to an independent refit. |
| probe07 | joint / equality / adversarial p from (W_r, W_obs); pairwise contrasts | CLEAN. All 8 test p-values reconstruct exactly; 2 contrasts identified by comparator difference and matched exactly. |
| probe08 | failure accounting off the dump | CLEAN. 19 zero-treated draws == 19 all-missing rows == e(perm_failed); 201 rows kept; _dta chars present. |
| probe09 | permonsetcohorts (§7, shipped 2026-09-16, unverified) | CLEAN. e(perm_irebuilt)=1; kernel conformant; 5 of 5 draws bit-identical to a hand-rebuilt-i refit. |
| probe10 | permstrata, sampled and permexact | CLEAN. npossible = 90^4 = 65,610,000 and 36 exactly; zero composition violations in 100 and 36 replications; 36 distinct assignments, no duplicates. |
| probe11 | element-level holes (subgroup level on 3 of 24 clusters) | CLEAN. n=193 vs 199 per column, all kernel rows conformant; all six step-down matrices posted. |
| probe13 | permexact MEMBERSHIP, bitwise | CLEAN. In all four regimes exactly one of the 28 enumerated replications carries the observed assignment and its coefficients are bitwise equal to rep 0, so p=0 is unreachable. |
| probe15/15b | permexact with elements some assignments cannot estimate | **FINDING**: exact shares over the reduced n are right (2/27, 2/13) and p_se is withheld -- but e(perm_sd) switches denominator mid-table. |
| probe14 | row alignment of the draw matrices when replications fail | **FINDING**: e(permdist_b)/e(permdist) are compacted (111) while e(perm_jointW)/e(perm_advW) keep full length (120). |
| probe12/12b | e(gteffects) blanking; contreat cfx cell tables | CLEAN. gteffects cols se/p/lb/ub all missing, e(perm_gt_blank)=1. contreat: ame rows conformant; lattice(all) permuted, lattice(baseTreated) blanked with e(perm_cfx_blank)=1. |
| probe16 | is the sampled shuffle UNIFORM over the 28 assignments? | CLEAN. chi2 = 21.24, df = 27, p = 0.77; all 28 drawn. |
| probe17 | permsave reconstructability under permtstat | CLEAN. e(perm_p) rebuilt from the FILE ALONE, 0 mismatches; e(perm_ase_obs) == e(perm_comparator_se) bitwise. |

## FINDINGS (both low)
F1 e(permdist_b)/e(permdist) compacted (111) vs e(perm_jointW)/e(perm_advW) full length (120)
   when replications fail -> repro_permdist_rowalign.do
F2 e(perm_sd) switches denominator (1/n vs 1/(n-1)) inside one permexact table
   -> repro_permsd_denominator.do

## UNTESTED
permtci when se_obs is missing (could not construct); FWER SIZE of permstepdown (algorithm
verified, not its coverage); the >10% failure caution and the max(20,half) refusal
thresholds (9.5% and 14% observed, neither boundary crossed); permute(probabilities) +
permonsetcohorts together; a ritest cross-check of permrawp.
