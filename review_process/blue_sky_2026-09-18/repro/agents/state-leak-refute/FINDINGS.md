# SKEPTIC verdicts on the "state-leak" findings (wooldid blue-sky review)

Build: /private/tmp/claude-501/.../scratchpad/bluesky/build/wooldid.ado (sha c1d8979), loaded by
`adopath ++`. Every leg below was rebuilt from the CLAIM on my own generated data, in my own
do-files. I never executed the original prober's scripts.

---

## state-leak-01-missing-control-outcome -> **PARTIAL**, severity revised **medium -> low**

My repro: `rA_missctrl.do` (8 legs) and `rB_verbose_and_remedy.do` part 1 (exact hex + the
second remedy). DGP: 40 states x 12 years, cohorts 2009/2011/2013 of 10 each, 10 never-treated.

What holds:
* y missing on every never-treated row -> r(459) with exactly the quoted message. TRUE.
* The message does not name missing data. TRUE. This is the only surviving half.

What does NOT hold:

1. **"safety(off) ... proceeds silently ... and reports it with NO caution of any kind."**
   False as stated. The safety(off) run's own header reads
   `N clusters = 30; N = 240 (240 obs dropped from initial sample)` on a 480-row dataset.
   That IS the dropped-observation disclosure the help promises (sthlp:898/927/1411); half the
   sample being gone is stated on the first line of output.
   Separately, the *suppression of cautions* is the documented contract of the option the user
   typed: sthlp:1159 - "When safety(off) is not explicitly specified, wooldid will caution the
   user or halt estimation ... Specifying safety(off) removes them ... suppressed along with
   everything else." The ado gate is literal (`if "`safety'" != "off"` at ado:838).

2. **"substitutes a different sample and estimand" (implying a wrong number).**
   The number is right. e(wdidmainresults) from
   (a) y missing on all never-treated + safety(off), and
   (b) the never-treated rows PHYSICALLY DELETED + safety(off)
   are **bit-identical**: max|a-b| = +0.0000000000000X-3ff (exactly 0), BITIDENTICAL=1.
   safety(off) returns precisely the answer the honest hand-cleaned dataset gives.

3. **"both remedies the message offers are wrong for this data."**
   False. The second remedy - "modify the sample window manually" - works on exactly this data:
   `wooldid y st yr onset if yr <= 2012, cluster(st) att` -> rc 0, N = 240, N_dropped = 80,
   att = .21855514 (the 2013 cohort becomes never-treated and supplies the control group).
   A hand-restricted `keep if yr<=2012` dataset gives the identical table. This is also the
   author's standing guidance for the no-never-treated design (brief, KNOWN BY DESIGN).

4. **"the diagnosis is not [correct]".** Deleting the never-treated rows outright produces the
   IDENTICAL r(459) and the IDENTICAL message (leg A4 == leg A2). So the message is not a
   mis-diagnosis of a missing-data problem; it is the same accurate statement about the screened
   sample that a genuinely control-less design gets. It is also preceded by
   "Caution: estimability checks found time periods exist in sample without a control unit,
   causing them to be dropped."

5. **"the documented disclosure never happens [on this path]".** The run errored; there is no
   results header on any failed run, of any command. That is not a doc-vs-code divergence.

Not outcome-specific: zero weight on the never-treated rows (leg A7) and missing `i` on them
(leg A8) give the same r(459) - it is the general markout->screen ordering, not a y bug.
Partial missingness is handled (leg A6: 30% of control rows missing -> rc 0, N=436, N_dropped=44).

6. **Why the two routes coincide (measured, not theorised).** The safety(off) run's realized
   sample is e(N) = 240 of 480 - exactly 30 treated states x the 8 calendar years 2005-2012 -
   and the `if yr <= 2012` run's realized sample is also e(N) = 240, with an identical table.
   The estimability screen that drops control-less periods (ado:830-836,
   `qui replace ifvar = 0 if controlcheck2 == 0`) is NOT gated on safety(); only the caution at
   ado:838 is. So safety(off) does not "proceed on a different sample" of its own invention -
   it lands on precisely the window the message's second remedy names. The prober's 48% move on
   their richer DGP is consistent with a dynamic treatment effect averaged over a shorter
   post-window (inference, not something I measured on their data).

CORRECTION: the live issue is one sentence of message text - the r(459) message could add
"(note: rows are also removed for missing y/i/t/weight before this check)". That is a low-severity
UX nit, not a medium-severity estimand substitution.

---

## state-leak-02-verbose-leaves-variable -> **PARTIAL**, severity **low -> low (doc nit only)**

My repro: `rB_verbose_and_remedy.do` part 2.

Observation CONFIRMED: plain c(k) 8 -> 8 and no `_est_` variable; verbose c(k) 8 -> 9 with
`_est____wooldidfullmodel` present and `estimates dir` = `___wooldidfullmodel`.

But the framing "leaves it there" and the proposed remedy are both wrong:

1. **It is not a leak; it is the storage the documented feature requires.** sthlp:1132 documents
   that verbose "stores the underlying regression model ... to ___wooldidfullmodel; this model
   can be accessed using estimates restore". `_est_<name>` is Stata's own esample marker for a
   stored estimate. wooldid did not invent it and cannot keep the estimate without it.

2. **The prober's remedy destroys the documented feature, silently.** I ran it:
   verbose run (marker marks 480 rows) -> `drop _est____wooldidfullmodel` ->
   `est restore ___wooldidfullmodel` returns **rc 0** and `estimates describe` returns **rc 0**,
   but `count if e(sample)` returns **N = 0** (was 480). Dropping the marker silently empties
   e(sample) of the very model the help tells users to restore. "the variable can be re-created"
   is false.

3. **It self-heals through the documented object and does not accumulate.**
   `est drop ___wooldidfullmodel` removes the variable too (c(k) back to 8).
   A subsequent PLAIN wooldid run removes it (ado:4539): c(k) back to 8.
   Three back-to-back verbose runs: c(k) = 9, 9, 9 - no growth.

CORRECTION: severity stands at low and the only actionable item is one clause in the verbose help
entry ("...and, with it, Stata's `_est____wooldidfullmodel` sample marker variable"), exactly as
sthlp:712 already does for `_ppmlhdfe_d`. Not a state leak.

---

## state-leak-03-mata-namespace-leftovers -> **REFUTED** (attribution), severity **low -> none**

My repro: three SEPARATE Stata processes, each starting from `mata: mata clear` and planting
`mata: k = 12345`: `rD_mata_wooldid.do`, `rE_mata_reghdfe.do`, `rF_mata_gtools.do`.

| arm | k after | Mata globals added | Stata scalars added |
|---|---|---|---|
| wooldid alone | **3** (clobbered) | `__00000L`, `__00000M`, `k` | `__gtools_gfile_ghdfeabs`, `__gtools_gfile_gregvcov` |
| reghdfe alone | **12345** (intact) | none | none |
| **gtools alone** (`gegen`+`gcollapse`, no wooldid, no reghdfe) | **2** (clobbered) | `__000000`, `__000001`, `k` | `__gtools_gfile_ghdfeabs`, `__gtools_gfile_gregvcov` |

**gtools owns it.** A bare `gegen ... , by(st)` + `gcollapse ..., by(st yr)` in a session that has
never loaded wooldid clobbers Mata `k` and creates both `__gtools_*` scalars. The value left in
`k` tracks the by-variable count (2 under `by(st yr)`, 3 under wooldid's internal collapse), which
is gtools' own by-key machinery. reghdfe is innocent. wooldid's only involvement is calling a
documented hard dependency.

This is exactly the outcome the prober's own code reading predicted and the condition they set
for striking it ("If it is a dependency's, STRIKE this finding"). Refuted as a wooldid defect;
it belongs upstream in gtools if anywhere.

---

## Negative results from my own catch-all (`rC_state_catchall.do`)

A plain wooldid run on a session deliberately holding NON-DEFAULT settings
(`set emptycells keep`, `set level 80`, `set varabbrev off`, `set cformat %9.4f`) changes
**nothing** except the two gtools scalars above. Specifically unchanged across the run:
c(level), c(emptycells), c(coeftabresults), c(type), c(varabbrev), c(matsize), c(dp), c(cformat),
c(pformat), c(sformat), c(showbaselevels), c(showemptycells), c(showomitted), c(fvlabel),
c(linesize), c(k), c(N), `: sortedby`, `c(rngstate)`, `c(rng)`, `c(frame)`, the datasignature,
`estimates dir`, `matrix dir`, `frames dir`, the global namespace, and physical row order
(`assert myrow == _n` passes).

Attacked and clean:
* **`set emptycells` program-scoping** (the ado RELIES on it; the record calls B18 "verified
  illusory"). Caller at `keep`, wooldid sets `drop` internally, caller is still at `keep`
  afterwards - on the SUCCESS path and on the r(459) ERROR path. The assumption holds.
* **`set level` restore (F054/WOOLDID_LEVEL0).** Caller at 80, `level(90)` requested:
  c(level) = 80 afterwards on the success path, on the r(459) estimability-error path, and on a
  parse-stage r(198) path. No `WOOLDID_LEVEL0` global survives any of them.
* **RNG.** A plain run (and a `gteffects makeplots` run) leaves `c(rngstate)` byte-identical - an
  estimator that silently advanced the user's stream would break their simulations. It does not.
  `bootstrap()` does advance it, which is what every replication command does.
* **Error-path residue.** The only object a failed run leaves is the documented
  `WOOLDID_LASTRC` global (NFRPL-01/02, by design). Frame, level, emptycells, data and estimates
  are all as they were.

---

## Additional negative results from my own probes (`rG`, `rH`, `rI`, `rJ`)

* **Microdata / fast-collapse path (12,000 rows, 40 individuals per state-year).**
  `rG_microdata_state.do`: full session snapshot before and after -> **NO DIFFERENCES AT ALL**
  (level, emptycells, c(k), c(N), sortedby, frame, datasignature, estimates dir, global
  namespace, frames dir, varlist), row order preserved. Fast collapse vs `nofastcollapse`:
  max|diff| over `e(wdidmainresults)` exactly 0 - bit-identical.
* **`xtset` survives** a wooldid call: `r(panelvar)`/`r(timevar)` and `_dta[iis]`/`_dta[tis]`
  unchanged (`rH`, H1).
* **A caller's own `preserve` straddling a wooldid call** works and restores byte-exactly
  (`rH`, H2: datasignature identical).
* **The verbose marker is on the right scale** even when the collapse fires: on 12,000-row
  microdata it marks 12000 rows = c(N) = e(N), and `est restore` + `count if e(sample)` = 12000
  (`rH`, H3). `est sto` fires before the collapse's `preserve`, so there is no scale mismatch.
* **`e(sample)` is sound on every path I tried** (`rI`): unit-level 300/300, microdata with the
  collapse 12000/12000, `nofastcollapse` 12000/12000, with 40 rows made missing 11960/11960
  (e(N_dropped)=40), with an `if()` 290/290 and 11600/11600. This was the natural place for a
  collapsed-scale esample to leak into the caller's data. It does not.
* **`permute(shuffle)` and `bootstrap(standard)`** leave the caller's datasignature unchanged,
  row order pinned, sortedby preserved, `frames dir` = default only, and no wooldid globals
  (`rI`, I8/I9). `rseed()` reproduces bit-exactly across a deliberately moved RNG stream
  (`rH` H6, max diff 0), and `set seed` alone reproduces bit-exactly without `rseed()`
  (`rH` H7, max diff 0).
* **Row-order invariance** (`rJ_roworder.do`, same rows in shuffled and reverse order):
  unweighted att/att_i/att_it/att_itime 3.6e-15 (shuffled) and 8.2e-15 (reversed);
  microdata fast-collapse 9.9e-13 and `nofastcollapse` 9.9e-13; bootstrap with `rseed()`
  5.3e-15; `permute(shuffle)` 7.5e-15; and `set sortseed` makes no difference at all
  (bit-identical). All at or below the project's 1e-10 indifference threshold.
  **One arm is above it** - `[aw=wt]` with an event-study window and `jointtests` moved
  3.73e-08 relative. Isolation in `rK_isolate_aw.do`; see the note there.

---

## NEW (mine, not seeded): `skept-roworder-01` -- the same rows in a different physical order
## change an event-study SE by 3.7e-08 relative, but only on the HDFE path with continuous weights

Repro: `rL_attribute.do` (self-contained; `rJ_roworder.do` found it, `rK_isolate_aw.do` bisected it).

CLAIM: with a weight whose values are not exactly representable, wooldid's reghdfe path returns a
different event-study standard error depending on the caller's physical row order. The relative
move is 3.7e-08, which is 370x the project's own 1e-10 indifference threshold and about 0.3 of a
single-precision ulp. Estimates do not move; one horizon's `se`, and its `t`, `p`, `lb`, `ub`, do.

Two files, identical rows, one sorted `st yr` and one in random order. 30 states x 10 years,
cohorts 2008/2011/2013, `att esprelength(3) espostlength(3)`, max relative difference over
`e(wdidmainresults)`:

| arm | max rel diff |
|---|---|
| L1 `[aw=wt]`, default `regtol(9)` | **3.732e-08** |
| L2 `[aw=wt]`, `regtol(12)` | **3.732e-08** (unchanged) |
| L3 `[aw=wt]`, `regtol(14)` | **3.732e-08** (unchanged) |
| L4 `[aw=wt]`, `regtol(6)` | 4.793e-09 |
| L7 `[pw=wt]` | **3.732e-08** |
| L11 `[aw=wt] robust` | **2.496e-08** |
| L12 `[aw=wt] fe(st)` | **1.662e-08** |
| L5 `[aw=wt] unconditionalse` (regress, no HDFE) | 1.407e-13 |
| L6 `[aw=wt] customvce(cluster st)` (regress, no HDFE) | 1.790e-14 |
| L8 `[aw=iwt]` (integer, exactly representable) | 1.019e-13 |
| L9 `[aw=one]` (constant) | 4.870e-15 |
| L10 unweighted | 1.448e-14 |

Diagnosis, as far as I took it:
* It is **not** the absorption convergence tolerance. `regtol(9)`, `regtol(12)` and `regtol(14)`
  give the identical 3.732e-08; the documented mechanism at sthlp:1117 does not explain it.
* It is **HDFE-side**: both paths that bypass reghdfe (`unconditionalse`, `customvce()`, which run
  `regress` on explicit dummies) are clean at 1e-13/1e-14 on the same wooldid weight construction.
* It needs **non-representable weight values**: exactly-representable integer aweights are clean at
  1e-13, constant weights at 5e-15, unweighted at 1e-14.
* **It is NOT a surviving B21-class float store in wooldid** (checked, `rM_floatcheck.do`).
  ado:1619-1621 build the cell weight denominators with three UNTYPED `gegen ... = total()` calls,
  which under `c(type)` = float would be the exact B21 signature. They are not: gtools' `gegen`
  defaults to **double** here (M0: untyped type = double, bit-identical to an explicit
  `gegen double`; plain `egen` would have given float). And running the firing arm under
  `set type double` changes nothing (M1 6.270e-13 == M2 6.270e-13; M3 == M4).
* **It IS reghdfe's weighted clustered VCE, isolated directly** (`rM`, M7-M9). Bare
  `reghdfe y post i.st#c.post [aw=wt], absorb(st yr) cluster(st) tol(1e-9)` on the two row
  orderings, with no wooldid anywhere: max relative difference in `e(V)` = **8.079e-10** with
  continuous aweights, **exactly 0** with integer aweights, **exactly 0** unweighted. `e(b)` is
  ~1e-13 in all three. So reghdfe's VCE is order-dependent at ~1e-9 precisely when weights are
  continuous, and exactly order-invariant otherwise.
* **The 3.7e-08 is wooldid amplifying that ~1e-09 input, and is DGP-specific.** The same option
  set on a different draw of the same DGP shape (`rM` M1) moves only 6.270e-13. So the size of the
  move depends on how well conditioned the particular horizon is, not on the option combination
  alone; a weakly identified pre-horizon amplifies the reghdfe perturbation by ~10^2.
* Which cells move (from `rK`): exactly the `pre3` row's `se`, `t`, `p`, `lb_95ci`, `ub_95ci`. The
  `att` point estimate and the other horizons stay below 1e-10.

Why it is worth the author's time despite being the 8th significant digit: the project blesses
baselines and adjudicates battery diffs on bit-identity, so a weighted event-study config that is
only reproducible to ~1e-8 across row orderings is a false-positive generator there, and the
workaround is knowable (`unconditionalse`/`customvce()`, or integer weights).
Severity **low**; realistic **yes** (weighted admin/survey panels with an event-study window are an
ordinary wooldid use, and re-sorting a dataset between runs is ordinary); confidence **probable** -
the measurement is certain and the root cause is isolated to reghdfe, but how large the amplified
move gets is dataset-dependent (3.7e-08 on one draw, 6.3e-13 on another).
