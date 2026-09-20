*==========================================================================
* gen_mutate.do -- in-memory transformations of a gen_core dataset that
*                  have a KNOWN effect on the answer.  Each one is half of
*                  an invariance test: run wooldid, mutate, run again,
*                  compare.
*
* USAGE:
*     local MUT dupall
*     include "`GENDIR'/gen_mutate.do"
*
* CATALOG (what must / must not move; see data_shapes.md section "invariances")
*   dupall      every row duplicated. POINT ESTIMATES must be identical;
*               SEs must NOT be (N doubles).  Cluster/unit counts double
*               only if you also relabel.
*   dupcohort   cohort 1's rows duplicated.  att MUST move (cell weights),
*               att_i MUST NOT (every unit still counts once).
*   halvew      w = w/2.  Combined with dupall this is the aweight-vs-
*               physical-duplication identity for point estimates.
*   aw2         collapses each duplicated pair back to one row with w*2.
*   relabel     unit / hh / cohort ids replaced by a random bijection.
*               EVERYTHING must be identical (ids carry no information).
*   reorder     row order shuffled, ids untouched.  EVERYTHING must be
*               identical.  Any movement is an order dependence.
*   shift_t     t and ttre shifted by +500.  EVERYTHING must be identical,
*               event-study rows included.
*   scale_t     t and ttre multiplied by 7 (after centring).  The att
*               family must be identical; horizons are now on a 7-grid,
*               so compare with espre(0) espost(0).
*   scale_w     w multiplied by 1e6.  Point estimates AND standard errors
*               must be identical under aw and pw.
*   tofloat     y, w, dose, x1 recast to float.  Differences here are the
*               fcol_sg precision class; expect small, not zero.
*   extmiss     extended missings (.a .b) sprinkled into y, x1, w, sggap.
*   ifflag      adds `flag' in {0,1,2} for the exact-one if-rule.
*   cellmean    collapses repeated cross-sections to one row per unit-
*               period cell (y = weighted mean, w = summed).  With NO
*               controls every regressor is cell-constant, so the att
*               family must be identical to the uncollapsed run.
*   negw        flips the sign of 5% of the weights.
*   zerow       zeroes 10% of the weights.
*   dropmid     deletes one interior calendar period for half the units.
*==========================================================================

if "`MUT'" == "dupall" {
  qui expand 2
}
else if "`MUT'" == "dupcohort" {
  qui expand 2 if cohort == 1
}
else if "`MUT'" == "halvew" {
  qui replace w = w / 2
}
else if "`MUT'" == "aw2" {
  qui duplicates drop
  qui replace w = w * 2
}
else if "`MUT'" == "relabel" {
  preserve
  qui keep unit
  qui duplicates drop
  set seed 424242
  set sortseed 424242
  qui gen double _r = runiform()
  sort _r unit
  qui gen long newunit = _n * 3 + 7
  qui keep unit newunit
  tempfile _map
  qui save `_map'
  restore
  qui merge m:1 unit using `_map', nogen
  qui gen long _newhh = cond(hh == unit, newunit, hh - unit * 100 + newunit * 100)
  drop hh
  rename _newhh hh
  drop unit
  rename newunit unit
  qui replace cohort = cohort * 5 + 11
  order unit cohort hh
}
else if "`MUT'" == "reorder" {
  set seed 515151
  set sortseed 515151
  qui gen double _shuf = runiform()
  sort _shuf
  drop _shuf
}
else if "`MUT'" == "shift_t" {
  qui replace t    = t + 500
  qui replace ttre = ttre + 500
}
else if "`MUT'" == "scale_t" {
  qui summarize t, meanonly
  local _t0 = r(min)
  qui replace t    = (t - `_t0') * 7
  qui replace ttre = (ttre - `_t0') * 7
  qui replace h    = t - ttre
}
else if "`MUT'" == "scale_w" {
  qui replace w = w * 1e6
}
else if "`MUT'" == "tofloat" {
  qui recast float y
  qui recast float ync
  qui recast float w
  qui recast float dose
  qui recast float x1
}
else if "`MUT'" == "extmiss" {
  qui replace y     = .a if mod(_n, 23) == 0
  qui replace ync   = .a if mod(_n, 23) == 0
  qui replace x1    = .b if mod(_n, 29) == 0
  qui replace w     = .c if mod(_n, 31) == 0
  qui replace sggap = .a if mod(_n, 37) == 0
}
else if "`MUT'" == "ifflag" {
  cap drop flag
  qui gen byte flag = 1
  qui replace flag = 2 if mod(unit, 5) == 0
  qui replace flag = 0 if mod(unit, 7) == 0
}
else if "`MUT'" == "cellmean" {
  qui gen double _wy   = w * y
  qui gen double _wync = w * ync
  qui gen double _wd   = w * dose
  qui gen double _wx1  = w * x1
  qui gen double _wx2  = w * x2
  collapse (sum) _wy _wync _wd _wx1 _wx2 w (first) cohort region ttre h post truetau xcell ///
           sgall sg2 sggap sgtiny sgbig sgmiss sgvary, by(unit t)
  qui gen double y    = _wy   / w
  qui gen double ync  = _wync / w
  qui gen double dose = _wd   / w
  qui gen double x1   = _wx1  / w
  qui gen double x2   = _wx2  / w
  drop _wy _wync _wd _wx1 _wx2
  qui gen long hh = unit
  order unit cohort hh region t ttre h post y ync truetau dose x1 x2
}
else if "`MUT'" == "negw" {
  qui replace w = -w if mod(_n, 20) == 0
}
else if "`MUT'" == "zerow" {
  qui replace w = 0 if mod(_n, 10) == 0
}
else if "`MUT'" == "dropmid" {
  qui summarize t, meanonly
  local _tm = int((r(min) + r(max)) / 2)
  qui drop if t == `_tm' & mod(unit, 2) == 0
}
else {
  di as error "gen_mutate: unknown MUT `MUT'"
  error 198
}

di as text "MUT " as result "`MUT'" as text " applied; N = " as result _N
