# oob-untreatedZero — independent verification of the "basic control leaks into untreatedZero" claim

STATUS: designing my own exact-fit oracle (do NOT reuse the prober's design).

## Code facts established by reading (pinned build c1d8979)
- ado:1234 gate: `if strpos("contreatcontrols","basic") != 0 | unconditionalse | customvce != "" | poisson`
  -> builds `ccringerhead`/`ccringerheadabs` = `i.null##c.(ringer)` (a pooled dose slope; `null` is a
  constant tempvar, so this collapses to a single pooled `c.ringer` column).
- ado:1300-1302: `contreatcontrolinserterFE = ccringerhead ccringertail`;
  `contreatcontrolinserterABS = ccringerheadabs ccringertail` (headabs blanked iff a byX was asked for).
- ado:1311-1315 untreatedZero branch: rewrites `contreatinj` ONLY. It never clears
  `contreatcontrolinserterFE` / `...ABS`. Comment there literally says
  "if always zero when not treated, the only thing to do is not control for anything!".
- ado:1485 (reg/poisson path used by unconditionalse & customvce) order:
  `regtype y i.i i.t fe trend inserter inserterFE ccc i._i#i._t#i.targetflag#c.contreat i._i#i._t#i.targetflag`
  => the pooled ringer column PRECEDES the cell#dose columns, so Stata's sweep keeps the pooled
  one and drops the LAST cell#dose column.
- ado:1476 (reghdfe default) puts the same set in absorb() instead.
- untreatedZero contreatinj (ado:1314) = `#c.contreat i._i#i._t#i.targetflag`, i.e. cell#dose with NO
  pooled base dose term + cell intercepts. Dose is exactly 0 off treated cells, so
  pooled_dose == sum(cell#dose) EXACTLY -> perfect collinearity.
- Arithmetic prediction: if cell j's dose column is dropped, surviving margins-visible slopes are
  beta_c - beta_j, and the ame (any weights summing to 1) = mean_w(beta) - beta_j.
- `_t` on the DEFAULT path (wooldid_reltime) = t - ttre + 1, so relyr 1 == onset period and
  targetflag==1 <=> t >= ttre. The prober's oracle sample (t>=ttre) is therefore CORRECT.
- sthlp:986 verbatim: "conttreatcontroltype: untreatedZero, which cannot be paired with any other
  type; or one or more of the following: basic, byTreated, byTtre, byCohort, and byTime"
- sthlp (untreatedZero prose): "If the continuous treatment variable is not meaningfully defined
  outside of treated cohort-periods, the control type untreatedZero can be used to specify g() = 0.
  By not including a g() function at all, identification of marginal effects of Contreat_jt will
  occur using only within treated cohort-period variation."

## My design (DIFFERENT from the prober's)
- 6 groups i=1..6, 7 units each (42 units), t=1..9 (378 obs).  [prober: 4 groups, 10 units, 8 periods]
- onsets ttre: i=1 -> 3, i=2 -> 5, i=3 -> 7; i=4,5,6 never treated.  [prober: 2 treated cohorts]
- 15 treated cells (7+5+3).  [prober: 8]
- dose = exp(N(0.2,0.55)) on treated cells, EXACTLY 0 elsewhere.  [prober: 1+2*U(0,1)]
- beta_it = 1.45 - 0.31*i + 0.22*r - 0.06*r^2, r = t-ttre (quadratic in r; prober used linear)
- y = a_i + g_t + kappa_it + beta_it*dose, exact (no noise).
- PREDICTED dropped cell = last in factor order = (_i=3, _t=3) i.e. i=3, t=9.

## RUN 1 RESULT (oob_uz.do / oob_uz.log) -- CLAIM CONFIRMED, AND WIDER THAN CLAIMED
TRUE ame = 0.907333333333 (hex +1.d08dfea27983cX-001). beta(i=3,t=9) = 0.72.
Predicted broken value = 0.187333333333.

My own replication of ado:1485 (plain reg): R2 = 1.0000, Residual SS = 0 (exact fit).
  note: 3._ii#3._tt#1.tf#c.dose omitted because of collinearity   <-- the LAST treated cell
  ringer coef = 0.7200000000                                     <-- == beta of the dropped cell
  every surviving cell dose coef == beta_c - 0.72 EXACTLY (15/15 cells checked)

wooldid ame, contreatcontrols(untreatedZero), att:
  att                         0.907333333333   relerr vs truth 5.8e-15   OK
  att robust                  0.907333333333   OK
  att cluster(unit)           0.907333333333   OK
  att ccnoabsorb              0.907333333333   OK  <- reghdfe drops the RINGER, not a cell
  att ccnoabsorb robust       0.907333333333   OK
  att customvce(robust)       0.187333333333   relerr vs truth 7.9e-01 ; vs prediction 3.9e-13
  att customvce(cluster unit) 0.187333333333   same
  att unconditionalse robust  0.187333333333   same
  att unconditionalse         RC=198 (needs a vce string)
  implied beta_dropped on all three broken paths = 0.7200000000 == beta(i=3,t=9) to 1.0e-13

ATT row moves too (prober only checked ame): att = 2.1175542 (good) vs 1.1486734 (broken)

contreatcontrols(basic) on the SAME data returns 0.187333333333 on BOTH engines
  => on customvce/unconditionalse, untreatedZero silently delivers the basic answer -- the
     very pairing wooldid_validate (ado:4797) REFUSES when the user asks for it.

POISSON -- WORSE: every poisson path is broken, INCLUDING the default one.
  ppmlhdfe's first attempt always fails ("absorbed categorical fixed effect with only 1 level"
  = the i.null term the gate inserted), so wooldid takes the rc-3200 fallback at ado:1468,
  which carries contreatcontrolinserterFE as an EXPLICIT regressor -> same drop.
  homogeneous twin, true semi-elasticity 0.25 everywhere:
    att poisson / robust / ccnoabsorb    ame = -0.000000000000   (truth 0.25!)
    att poisson customvce(robust)        ame =  0.000000000003
    att poisson unconditionalse robust   ame =  0.000000000003
  heterogeneous twin, truth 0.226833333333, gam(3,9)=0.18, prediction 0.046833333333:
    all four poisson paths               ame =  0.0468333...  EXACT match to prediction

USER-VISIBLE TEXT
  default / ccnoabsorb: NO caution of any kind.
  customvce()/unconditionalse/poisson: only the GENERIC multicollinearity caution
    "Caution: a treatment period coefficient interacted with something (likely, a term
     implementing a continuous treatment variable or an interactive control variable) was
     dropped for multicollinearity."
    "This could imply your results are not identified, depending on exactly what was dropped."
    "Re-estimate with option verbose for details on what was dropped."
    "Note: the multicollinearity caution above was triggered 1 times total."
  Nothing names untreatedZero, nothing says a pooled baseline dose control was inserted.

## REGRESSION CHECK
`git log -L1234,1236:wooldidinstall/wooldid.ado` touches exactly two commits:
  fc17051 2023-08-15 "Add files via upload"  -- ORIGINAL
  1df8613 2026-09-05 "SI-20: reference-cell margin - basic always absorbs the baseline dose control"
In the ORIGINAL 2023 upload the line already read (its line 935):
  if (strpos("contreatcontrols","basic") != 0 & "latticepart2hash" == "" ) | "unconditionalse" ==
     "unconditionalse" | "customvce" != "" | "poisson" == "poisson" {
=> the `| unconditionalse | customvce | poisson` disjunction is ORIGINAL. SI-20 only deleted the
   `& latticepart2hash == ""` qualifier on the `basic` arm. So this is NOT a defect introduced by
   a previous repair; it has been there since the first upload.

The valuable record is instead review_process/cell_collapse_contreat_2026-09-02.md, which read
this exact block and wrote down the OPPOSITE conclusion (s1, and again in s7 item 2):
  "`contreatcontrols(untreatedZero)` **alone** puts the ringer nowhere - it only rewrites
   `contreatinj` (1152-1154) - so the prototype's own comment, which names `untreatedZero` as a
   ringer path, is wrong on that point (section 7)."
The sentence immediately before it correctly notes the block fires "unconditionally under
unconditionalse, customvce() or poisson" -- the assessment saw the disjunction and still recorded
that untreatedZero puts the ringer nowhere, and the prototype comment was corrected in the WRONG
direction as a result. Nobody then crossed untreatedZero with those three engines.

Near-miss: review_process/amese/verifier1_adversarial.md P4 probed
  "dose zeroed outside treated cells + contreatcontrols(untreatedZero) + safety(off) + cluster(unitid)"
i.e. the DEFAULT engine (which is correct) -- never customvce/unconditionalse/poisson.

Existing machinery that ALMOST catches it: SI-A4/d3 (wooldid_standing_issues.md ~300-315,
wooldid_amescreen mode(detect), ado:1542) scans e(b) colnames for omitted treated-cell dose terms
(`co.` marker) and voids the ame -- but only "on all-dropped + none-live". Here exactly 1 of 15 is
dropped, so it stays silent. The dropped term IS in the scan's scope (targetflag==1).

The shipped comment at ado:1230-1233 justifies the i.null##ringer spelling as ensuring "the base
level of contreat gets partialled out... if not specified... one of the contreat control
coefficients gets dropped for colinearity instead". Under untreatedZero contreatinj uses `#`
(ado:1314), not `##`, so there IS no base level to partial out: the manipulation has nothing to do
and the outcome it was written to PREVENT is exactly what happens.

## RUN 2 RESULT (oob_uz2.do / oob_uz2.log)
TRUE att = 2.117554176197 -- exactly what the default path prints (2.1175542).
The broken paths print att = 1.1486734, i.e. the ATT is understated by 0.9689 (46%).
THE PROBER ONLY CHECKED THE ame ROW; the att row is corrupted too.

VERBOSE proof of the code half (the echoed regression under customvce(robust)):
  reg y i.i i.t   i.__000009##c.(__000008)   i.__000005#i.__000006#i.__000007#c.dose
      i.__000005#i.__000006#i.__000007   if __000001 == 1 ,  vce(robust)
  __000009 = null, __000008 = ringer -> the pooled baseline dose control IS in the model under
  contreatcontrols(untreatedZero). verbose then names the casualty:
  "The specific problematic coefficient is: 3o.__000005#3o.__000006#1o.__000007#co.dose"
  "Given we are looking at something with an _i _t (other items here) structure,"
  "the problem values of i and relative-t are 3, 2."
  == exactly the cell my arithmetic predicted (i=3, t=9).

reghdfe control experiment (B1 vs B2 in the log) -- WHY the default path is safe:
  B2 absorb(i.i i.t)                       -> df_a 14, cell slopes = 1.14 .99 ... .72 = TRUTH
  B1 absorb(i.i i.t i.nullv##c.ringer)     -> df_a 15, cell slopes = 0.42 .27 ... 0.00
  So the default path is correct ONLY because the gate at ado:1234 is false there and the ringer
  set is never built. Had it been built, reghdfe would break identically. Same reason ccnoabsorb
  is safe: with the gate false, contreatcontrolinserterFE is empty, so ccnoabsorb has nothing
  to move out of absorb().

EVENT STUDY: every horizon is off by the same 0.72 on the broken path.
  ES_ame post0/1/2 : 0.83 0.99 1.03  (default)  ->  0.11 0.27 0.31  (customvce(robust))
  ES_att post0/1/2 : 1.52224 2.09894 2.29665    ->  0.64856 1.04817 1.31716

LATTICE (the shape of the help's own example 8, lattice(baseTreatedPost 3)):
  att                        att = 2.1175541762   ame = 0.9899399972
  att customvce(robust)      att = 0.6846554498   ame = 0.1625714...
  att unconditionalse robust att = 0.6846554498   ame = 0.1625714...
  With a lattice the head carries i.lattice##ringerpoly too, so more columns are destroyed.
  (The lattice arm is not a clean oracle - a lattice legitimately redefines the ame - but the
   default-vs-customvce gap of 3.2x in att and 6x in ame is the same defect, larger.)

REALISTIC noisy panel (240 units, 9 periods, unit effects sd .6, noise sd .4, cluster(unit)):
  truth 0.90733333
  att cluster(unit)                 ame = 0.80421677  se 0.03053  ci [0.744, 0.864]
  att customvce(cluster unit)       ame = 0.11917189  se 0.14136  ci [-0.159, 0.398]
  att unconditionalse cluster(unit) ame = 0.11917189  se 0.14185  ci [-0.160, 0.399]
  => the broken paths return a statistically INSIGNIFICANT near-zero where the truth is 0.91.
  (realistic poisson skipped there: y was continuous, rc 411 -- rerun as a count in run 3)

unconditionalse WITHOUT a vce is an informative refusal, NOT a defect:
  "Error: the option unconditionalse requires specifying either robust or clustered standard errors."
  "Specify robust or cluster(), or provide a robust/clustered vce via customvce()."

## RUN 3 RESULT (oob_uz3.do / oob_uz3.log) -- the arithmetic holds EXACTLY on noisy data
Same realistic noisy panel as run 2 part D. Population ame 0.90733333.
  ORACLE A (correct spec, reghdfe, NO ringer): mean of 15 cell slopes = 0.80421677
      == wooldid "att cluster(unit)" to 8 d.p.  -> the default path IS the correct OLS answer
         on this draw; its 0.10 gap from the population value is ordinary sampling error.
  estimated slope of the cell wooldid drops (i=3,relyr=3) = 0.68504489
  predicted broken value = 0.80421677 - 0.68504489 = 0.11917189
  ORACLE B (same spec WITH the pooled ringer) : mean = 0.11917189, ringer coef = 0.68504489
  wooldid att customvce(cluster unit)          : 0.11917189  se .14136  ci [-0.159, 0.398] p .400
  wooldid att unconditionalse cluster(unit)    : 0.11917189  se .14185  ci [-0.160, 0.399] p .402
  => reported == correct - slope_of_dropped_cell to 8 decimals on realistic noisy data.
  => the user is handed an insignificant near-zero (p = .40) where the truth is 0.91.

REALISTIC COUNT OUTCOME + poisson (mean 7.68, median 7, max 32, 11 zeros of 2160):
true dose semi-elasticity 0.25 EVERYWHERE (homogeneous), so the mechanism predicts ~0.
  att poisson cluster(unit)            ame = 0.00505342  se .05910  ci [-0.1108, 0.1209]
  att poisson robust                   ame = 0.00505342  se .05954  ci [-0.1116, 0.1217]
  att poisson customvce(cluster unit)  ame = 0.00505342  se .05910  ci [-0.1108, 0.1209]
  => 98% attenuation, ON THE DEFAULT POISSON ENGINE, on an ordinary count panel. The user sees
     only the generic multicollinearity caution plus the ppmlhdfe-fallback caution.

## VERDICT: CONFIRMED (and broader than the prober claimed)

## INNOCENT-EXPLANATION TEST -- closed by the help's OWN contract
sthlp:266 unconditionalse: "Compute standard errors that account for sampling variation in
  covariates... This option also forces wooldid to abandon high dimensional fixed effect
  estimation packages for the underlying two-way fixed effect regression in favor of the
  commands regress and poisson. As such, this option should be slower to use on average."
sthlp:270 customvce(): "Forces wooldid to switch to using either regress or poisson for
  underlying regression estimation, and then passes the user specified text without
  modification into the vce() option for these commands. This allows specifying a broader set
  of standard error calculation methods..."
sthlp:1146: "standard errors are the regression's coefficient-level (conditional) VCE --
  unconditionalse and impvalues do not alter them"
ado DF-1 comment (1270s): "Point estimates, ll, rss, r2 and N do not move."
=> both options are documented purely as VARIANCE/engine options. Neither is documented to
   redefine the estimand, and nothing in the contreatcontrols() entry conditions untreatedZero
   on the engine. There is no reading on which 0.187333 is the requested quantity.
sthlp:986 (not :995 in the pinned build) is where "untreatedZero, which cannot be paired with
  any other type" lives; wooldid_validate ado:4797-4798 enforces it with
  "Error: untreatedZero cannot be paired with other control types".
The help's worked example 8 at sthlp:1551 is
  wooldid y i t ttre, contreat(intensity) contreatcontrols(untreatedZero) makecfxplot lattice(baseTreatedPost 5)
  -- an intensity that is 0 outside the treatment period, NOT a count outcome and NOT poisson.
  (The brief's characterisation of that line is wrong on both counts.) It is still the lattice
  arm, which is the worse one.
