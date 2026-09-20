# LENS-1 (DOCUMENTED-BEHAVIOUR) tiebreak on FP-CTZERO-ENGINE — CONCLUDED

Probe: `lens_uz_documented.do` / `lens_uz_documented.log` (this dir). Build asserted
`*! version 1.4 2026-08-19 wooldid` from the pinned bluesky build (log:83).

## Verdict: PARTIAL on the prober's claim; the underlying defect is REAL, UNDOCUMENTED,
## UNRULED, and HIGH severity. It is not a correct-but-surprising design.

### What I measured (seed 20260917, my own DGP, conforming untreatedZero: `d2 = dose*post`)

PART A — homogeneous true dose effect 3.0, att truth 4.3206:
| arm | att | ame |
|---|---|---|
| hand oracle (cell slopes, NO pooled dose control, no wooldid) | 4.2729432 | **2.9632049** |
| wooldid DEFAULT (reghdfe) | 4.3822937 | **2.9632049** |
| wooldid DEFAULT + nofastcollapse | 4.3822937 | 2.9632049 |
| wooldid customvce(robust) | 0.1935925 | 0.0548165 |
| wooldid unconditionalse robust | 0.1935925 | 0.0548165 |

- DEFAULT ame is **bit-identical to my independent hand oracle**. The default arm is RIGHT.
- `nofastcollapse` changes nothing ⇒ the cell collapse is NOT the axis; the engine is.
- customvce and unconditionalse agree to the last digit ⇒ it is the SPECIFICATION, not the VCE.
- Hand `reg y i.i i.t c.d2 i.cell i.cell#c.d2` puts **b(d2) = 2.9083884** in the pooled term:
  the dose effect migrates into the pooled control that margins never switches off.

PART B — **heterogeneous** slopes (cohorts 1-3 β=3.0, cohorts 4-6 β=1.0; truth ame 2.2500,
att 3.1259). This is my discriminating test and it OVERTURNS the skeptic's invariant:
| arm | att | ame |
|---|---|---|
| DEFAULT | 3.0727451 | 2.1925650 |
| customvce(robust) | 1.7079442 | **1.2142114** |

Predicted a priori: reported ame = mean(β) − β_k for the one dropped reference cell k.
β_k = 1 ⇒ 2.25 − 1 = **1.25**; observed 1.214. The explicit arm is not "collapsed to ~0";
it is **arbitrarily re-referenced to one dropped cell's effect**. "~0" is only the
homogeneous special case. A plausible nonzero wrong number is WORSE than a null.

PART C — **poisson hits the DEFAULT engine** (truth xb-scale ame 0.300):
`ppmlhdfe y ... , absorb(i.i i.t **i.__00000A##c.(__000009)**)` (log:781-783) — the ringer is
IN the absorb; ppmlhdfe then fails ("absorbed categorical fixed effect with only 1 level"),
falls back to the explicit-regressor path, drops a cell slope, and reports
**att −0.0343963 (p 0.162), ame −0.0198475 (p 0.208)** against a truth of +0.30.
So the trigger set is untreatedZero × {customvce() | unconditionalse | poisson} and the
default engine is included via poisson.

Verbose commands (the mechanism, printed by wooldid itself):
- DEFAULT: `reghdfe y i._i#i._t#i.TF#c.d2 i._i#i._t#i.TF, absorb(i.i i.t)` — no ringer (log:491).
- customvce: `reg y i.i i.t **i.__00000A##c.(__000009)** i._i#i._t#i.TF#c.d2 i._i#i._t#i.TF, vce(robust)` (log:376).

### Documented-behaviour findings (my lens)

1. **Not documented.** sthlp:1033 states untreatedZero ⇒ "g() = 0 ... By not including a g()
   function at all". No engine caveat anywhere. sthlp:1034's "makes other technical changes to
   the estimation process" is about untreatedZero changing things, not engines undoing it.
   Grep of the whole sthlp for `untreatedZero` (98, 472, 703, 903, 986, 993, 1033, 1086, 1235,
   1418, 1551) finds no customvce/unconditionalse/poisson interaction anywhere; sthlp:1235's
   compat table covers permute()/bootstrap() only.
2. **CONTRA-documented, twice.**
   - sthlp:1086 tells the *user* that dose-function controls should go with
     "contreatcontrols(basic), or any option other than untreatedZero". The code does
     internally exactly what the help warns the user not to do under untreatedZero.
   - sthlp:1096-1097 recommends, for contreat designs where ame inference fails under
     cluster containment, "robust (or a robust **customvce()**)" and notes **unconditionalse**
     "reports ame standard errors on such a design without this caution". The help routes an
     untreatedZero user straight into the broken arms. This is why `realistic = true`.
3. **The governing ruling's premise is false here, and the record says the combination was
   never examined.** The DF-1 comment (author ruling 2026-09-13, ado:1263-1298) justifies the
   ringer as "a numerical twin of columns the margins-visible block already carries through its
   ## chain (c.contreat, ...). It adds no rank", and says the explicit paths "already report the
   correct residual count". Both are true under `basic` and FALSE under untreatedZero, where
   contreatinj (ado:1313-1314) deliberately removes the base dose level, so the ringer adds REAL
   rank — and, because the dose is exactly 0 outside treated cells, it is exactly collinear with
   the set of cell dose slopes, forcing a drop. `review_process/cell_collapse_contreat_2026-09-02.md:118-120`
   correctly says untreatedZero *alone* puts the ringer nowhere and never looked at the combination.
   Per the coordinator note, a ruling whose premise fails on measurement is a first-class finding.
4. **No gate, no targeted caution.** ado grep: the only untreatedZero compat screen is permute()
   (5841); nothing for customvce/unconditionalse/poisson. The SI-A4/R6 dose-slope screen
   (ado:11949) fires only when EVERY cell slope drops and sthlp:903 explicitly exempts
   untreatedZero, so its silence is by design, not a divergence. The one caution that does fire
   names the dropped cell slope (`6o.__000006#3o.__000007#1o.__000008#co.d2`) and never the
   pooled control that caused the drop, so it points the user at the wrong object.
5. **Not "poorly identified" in the documented sense.** sthlp:993-998's disclaimer covers data
   that violates the always-zero premise. Both the prober's and the skeptic's DGPs, and mine,
   set `d2 = dose*post` — conforming. The design is fully identified; the rank deficiency is
   manufactured by wooldid's own added column, and the default engine proves identification.

### Where the prober was wrong (why PARTIAL, not CONFIRMED)
- "opposite sign": the sign is residual noise. Under heterogeneity the wrong number tracks
  −β_k, so it can be any sign or magnitude. The invariant is re-referencing, not sign.
- "could not determine which side is right": determined. DEFAULT correct (bit-identical to an
  independent oracle on two DGPs); explicit arms wrong.
- "the screen should fire on both paths": moot — the default path has nothing to warn about.
  The right remedy is not to insert the control under untreatedZero at all (or to refuse
  the combination), which is a different fix from the one the prober asked for.
- The skeptic's "collapsed to ~0" is right only for homogeneous effects (Part B refutes it as
  the general statement).

### Two reading notes on the Part A table
- Only **ame** is bit-checked against the hand oracle (2.9632049 both). The oracle's att
  (4.2729432) sums `slope*d2` alone, while wooldid's att also carries the treated-cell
  intercept terms `i._i#i._t#i.TF` (zero in the DGP, estimated with noise), so att is checked
  against the DGP truth (4.3822937 vs 4.3206415), not against the oracle. Not a discrepancy.
- `A_RC_DEFAULT_NOFC=0` doubles as the v1.3 exclusion: the stale v1.3 in the PLUS adopath
  refuses `nofastcollapse` with r(198), so rc=0 there proves the pinned build ran.
