# SKEPTIC pass on FW-01 / FW-02 / FW-03 (weights surface)

Build asserted in every probe: `which wooldid` ->
`/private/tmp/.../a575a3bcab42/scratchpad/bluesky/build/wooldid.ado`, version 1.4 2026-08-19.
All data generated fresh in each do-file. The project's test suite was not touched.

Every repro below is self-contained: run it with
`"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do <file>` from this directory.

---
## FW-01 -- VERDICT: CONFIRMED (severity medium; magnitudes LARGER than reported)

**My own repro:** `p2_mechanism.do` (silent-wrong branch), `p1_exprgrid.do` (error branch).

I did not run the author's script. I rebuilt the case from the claim with an *asymmetric*
prediction that only the text-paste mechanism can satisfy:

| call | att | att_it | att_i | att_itime |
|---|---|---|---|---|
| `[aw=wf+zero]`, zero==0 | 1.464771975439270 | **1.464771975439270** | **1.464771975439270** | **1.464771975439270** |
| `[aw=wf]` (same data) | 1.464771975439270 | 1.458344350902463 | 1.415271428230021 | 1.417598312200727 |

All four bitwise equal (`+1.76fb4bd29be62X+000`) under `[aw=wf+zero]`, SEs too. Mathematically
`wf+zero` IS `wf`. Only verbatim pasting explains it: `gen double w_att_it = wf+zero/w_per_it`
-> `wf + 0` -> the cell division is annihilated and every cell-weighted estimand silently
becomes `att`. The mirror image also holds: `[aw=zero+wm]` is bitwise identical to `[aw=wm]`
on all four estimands and all four SEs, because a leading zero term is harmless.

Magnitudes on the finding's own setup (`wtot == wf+wm` asserted row by row,
`customvce(robust)`), all r(0), no error, no caution:

| estimand | `[aw=wtot]` | `[aw=wf+wm]` | reldiff |
|---|---|---|---|
| att | 1.470710677532787 | 1.470710677532787 | 0 (exact) |
| att_it | 1.461823093540426 | 1.468148949125649 | 4.33e-03 |
| att_i | 1.425732957679875 | 1.468059193713792 | **2.97e-02** |
| att_itime | 1.420357846860070 | 1.468041600513423 | **3.36e-02** |

That is 3.4%, vs the 5.1e-03 max the finding reported. On the `unconditionalse robust
safety(off)` leg the SE of att_it moves 0.126915925 -> 0.143077097, a **12.7%** SE error.
`customweightmultiplier(cwm)` breaks `att` too: 1.455847603 -> 1.464679018 (6.07e-03), as claimed.

**Corrections / additions to the finding:**
1. The error branch is INHERITED FROM reghdfe VERBATIM, not produced by wooldid, and the
   fix site is not the cell collapse. Grid of 9 expressions x 5 engines (`p1_exprgrid.do`):
   wooldid-default, wooldid-`nofastcollapse` and bare `reghdfe` return the **identical** rc
   and message for every expression (`2*k`->111 "variable 2*k not found"; `k/2`, `k*2`,
   `wf+wm`, `(wf+wm)`, `1/pscore`->198 "+ invalid name"/"1 invalid name"; `int(k)`,
   `cond(k>2,k,1)`->101). Bare `regress` accepts all nine. `nofastcollapse` changes nothing,
   so the collapse is not implicated.
2. The *realistic* half of this finding is the error branch, and the finding under-sells it:
   **`[pw=1/pscore]` -- inverse-probability weighting, the single most common weight
   expression in applied work -- dies r(198) "1 invalid name" on the DEFAULT path**, with no
   mention of weights. The silent-wrong branch needs `customvce()`/`unconditionalse` AND a
   top-level `+`/`-`, which is narrow.
3. New symptom worth naming: when the second term is small or zero the three cell-weighted
   estimands silently *converge on att*. In the `wf+wm` run all three land within 1e-4 of each
   other (1.46815/1.46806/1.46804) instead of spanning 1.4204-1.4618. "All my estimands agree"
   reads as reassuring, not as a bug.
4. I did not separately test `safety(off)` on the default path; it cannot help, since bare
   reghdfe with the same weight produces the same rc.

---
## FW-02 -- VERDICT: PARTIAL (claim true and bitwise reproduced; stated consequence is not)

**My own repro:** `p3_pw_scale_collapse.do` (section 1).

Reproduced exactly, on my own 30-state x 12-year x 8-county panel:
- `[pw=pop]`, no vce option: se 0.10884885, `e(vce)=[ols]`, `e(vcetype)=[.]`, header
  "Standard Errors: " blank.
- `[aw=pop], robust`: `mreldif(e(V)_pw, e(V)_awrobust)` = **exactly 0** (`+0.0X-3ff`).
- `[aw=pop]` iid: `mreldif` = 4.5e-03; se 0.08545442 (the pw se is 1.2738x this).
- `reghdfe y post [pw=pop], absorb(i t)` -> `e(vce)=[robust]`, `e(vcetype)=[Robust]`.
  `regress ... [pw=pop]` -> same. So this is wooldid's own posting, as the finding says.
- Distinctness from PS-33/PS-1 holds: `rulings_2026-09-07.md:68` rules "the iid default and
  the unlabelled header stand", and the project's own PS-1 audit note says e(vce)=="ols" is
  posted "for estout-class consumers" on a *bare unclustered* run. Neither covers a run whose
  e(V) is not iid. `grep` shows wooldid never READS e(vce)/e(vcetype) itself, so the harm is
  entirely downstream.

**Correction (why PARTIAL, and why I lower severity to low):** the finding's stated harm --
"estout/esttab read e(vcetype) to decide whether to print 'Robust standard errors' beneath a
table, so a pweighted wooldid result is published labelled OLS" -- is not what esttab does.
Measured in `p7_esttab.do` with SINGLE-model tables (so a mixed-vcetype table cannot confound
the answer): `esttab rclu, se nostar` on a `cluster(state)` model with `e(vcetype)=="Robust"`
prints "**Standard errors in parentheses**"; so does the `robust` model; so does the iid model.
esttab's default footer does not read e(vcetype) at all. What DOES surface it is an explicit
request: `estout rrob, stats(vcetype N)` prints a `vcetype  Robust` row, and `stats(vce)` would
print the token. So the false label is real and contradicts ado:4283's own stated contract, and
a user who asks for the vce row, reads `ereturn list`, or writes code keying on it gets "ols"
for a robust covariance -- but the automatic "published labelled OLS" harm the finding names
does not happen.

---
## FW-03 -- VERDICT: PARTIAL (core fact true; the "expected" and one component are wrong)

**My own repro:** `p4_fw_collapsed.do`. 8717 microdata rows -> 360 (state,year) cells.

Plain-`regress` oracle legs, SE relative to the microdata SE:

| leg | b reldiff vs micro | se / micro se |
|---|---|---|
| `[fw=n]` cells, `cluster(state)` | 2e-14 | **1.000000000** |
| `[aw=n]` cells, `cluster(state)` | 2e-14 | 1.060009592 |
| `[fw=n]` cells, iid | 2e-14 | **0.192690664** |
| `[aw=n]` cells, iid | 2e-14 | **1.006426015** |

wooldid legs: `[fw=n]` -> r(198) "Error: wooldid does not accept fweights. Specify aweights or
pweights."; `[aw=n]` == `[pw=n]` bitwise; att b matches the microdata wooldid run to ~1e-14;
att se ratio cells/micro under `cluster(state)` = 1.258867.

**What survives:** yes, the prescribed substitute reproduces the point estimate exactly and not
the clustered SE, and neither the message nor the help says so (`sthlp:57` documents only
`[aw/pw=weights]`). That is a low-severity message/help gap.

**Three corrections:**
1. **The iid leg is attributed to the wrong thing.** The finding says the substitute's SE is
   "9% too small iid" as though fweights would have fixed it. They would not: `[fw=n]` on cell
   means is **0.1927** of the microdata SE (a factor of 5 too small -- fw treats each cell mean
   as n identical observations, destroying the within-cell SSR), while `[aw=n]` is **1.0064**,
   i.e. nearly exact. Under iid the prescribed remedy is the RIGHT one and fweights are the
   disastrous one. Whatever the finding measured under iid is a cells-vs-microdata artefact of
   its DGP, not an aw-vs-fw gap.
2. **The clustered gap is entirely Stata's small-sample multiplier, not a difference in the
   estimated variance.** Meat and bread are identical (b agrees to 2e-14 and the per-cluster
   score sums coincide because x is cell-constant). The observed 1.060009592 equals
   `sqrt( [(360-1)/(360-42)] / [(8717-1)/(8717-42)] ) = 1.060009...` to 7 digits -- pure
   `(N-1)/(N-k)`. This is generic Stata aweight behaviour in every command, not anything
   wooldid does. The wooldid ratio is larger (1.2589) only because reghdfe's K is far larger
   relative to N=360 (an implied K of 137 reproduces 1.25887 to 5 digits; I solved for K, so
   treat that as a consistency check rather than a proof).
3. **`e(N)=360` is not a defect.** Every Stata command reports the number of observations
   under aweights. Reporting 8941/8717 there would be the anomaly.

---
## NEW FINDING (mine) -- SKEP-01: on the `unconditionalse` path e(vce)/e(vcetype) name the STAGE-1 regression's vce, not e(V)

**Repro:** `r_skep01_uncond_label.do`. Severity low, confidence certain, realistic.

ado:4291's ladder tests `robust` and `cluster()` BEFORE `unconditionalse`, so the
`unconditionalse` branch at ado:4300 is only reachable by the one call that supplies neither --
a bare `[pw=]` run (ado:4676 exempts pweights from the "requires robust or clustered" refusal).
**The load-bearing evidence is (i) vs (iii) below: ONE covariance, bitwise identical, carrying
TWO different names depending only on whether the user typed the (redundant) `robust`.**
Measured on one 24-state panel, `att` only:

| call | se | e(vce) | e(vcetype) |
|---|---|---|---|
| (i) `[pw=w], unconditionalse robust` | 0.2855386999 | robust | Robust |
| (ii) `[pw=w], robust` (no unconditionalse) | 0.2437086534 | robust | Robust |
| (iii) `[pw=w], unconditionalse` | 0.2855386999 | **unconditional** | . |
| (iv) `[aw=w], unconditionalse robust` | 0.2865183827 | robust | Robust |
| (v) `[aw=w], robust` | 0.2437086534 | robust | Robust |

(i) and (iii) are the SAME matrix bitwise (`+1.246441c7f884cX-002`) under two different
names. That is the claim. ado:4283 states the contract being broken: "e(vce) and e(vcetype)
name what e(V) IS", and on this path they name the stage-1 regression instead. It reaches
EVERY `unconditionalse` user, since every one of them must supply `robust` or `cluster()`
unless they use a bare pweight.

A weaker secondary observation, offered but not leaned on: (i) and (ii) are different
covariances 17.2% apart in SE (0.2855 vs 0.2437) that carry identical labels. I do NOT rate
this as false the way "ols" is false under a pweight -- `margins, vce(unconditional)` is itself
a sandwich, so "Robust" is defensible for it, and wooldid's own "unconditional" token is not a
standard Stata `e(vce)` value either. Severity stays low; this is not a second FW-02.

---
## NEGATIVE RESULTS (attacked, clean)

1. **Weight-scale invariance** (`p3`, section 2). ado:1655 asserts "weight scale here actually
   doesn't matter". `[aw=pop]` vs `[aw=pop1000]` (=1000*pop, a pre-made variable so no
   expression is involved): all four estimands and all four SEs agree to 1-2 ulp
   (`+1.473577ab01cc8` vs `+1.473577ab01cc6`). Same under `[pw]`. Holds.
2. **Default cell collapse vs `nofastcollapse` with a weight that varies WITHIN the
   cohort-by-period cell** (`p3`, section 3). i=state, 8 counties per state-year, lognormal
   county population. Under iid, `robust` AND `cluster(i)`, all four estimands and all four
   SEs agree to <=1 ulp. If the collapse were averaging the cell before forming the robust
   meat this would have shown up as a visible SE gap. Clean.
3. **`poisson` weight-type gates match sthlp:710** (`p5`, section B). `[pw=]`+poisson rc 0;
   `[aw=]`+poisson r(198) "Error: poisson only accepts pweights, not aweights."; `safety(off)`
   does not lift it.
4. **Zero survey weights are counted, not swallowed** (`p5`, section C). ado:683 drops rows
   with `w <= 0`. With 1760 of 8717 rows at weight 0, wooldid prints "N = 6957 (1760 obs
   dropped from initial sample)", posts `e(N_dropped)=1760`, and the run is BITWISE identical
   to the hand-restricted `if w0>0` run. `regress [aw=w0]` also reports N=6957.
5. **sthlp:268's pweight exemption is really implemented** (`p6`). `[pw=w], unconditionalse`
   with no robust/cluster runs; `[aw=w]` in the same call is refused r(198). (Killed a
   candidate finding of mine before it cost a Stata run -- the gate is at ado:4676.)
6. **pweights survive the replication engines** (`p6`). `[pw=w] + bootstrap(standard) reps(25)
   cluster(state)` and `[pw=w] + permute(shuffle) permreps(25)` both run.
