*=====================================================================
* p4b_thinvsthick.do -- SELF-FOUND.  Exact column counts by counting the
* NONEMPTY cells directly (fvexpand builds the full grid before dropping,
* so it cannot be used at unit-level i).  No estimation; seconds.
*
* Calibration: for i=state, fvexpand gave 2579 columns for
* `cell##i.lat5##c.exposure'.  The identity below reproduces it exactly:
*     cols = 2*(n_cell + n_cell_lat) + 2*L + 1
* SELF-CONTAINED.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"

clear
set seed 20260916
set sortseed 20260916
set obs 25
gen long state = _n
gen int ttre = cond(_n <= 18, 2005 + mod(_n-1,8), .)
expand 20
sort state
by state: gen int cnum = _n
gen long county = state*1000 + cnum
gen double exposure = exp(rnormal()*0.7)
expand 20
sort county
by county: gen int t = 1999 + _n
gen byte post = (t >= ttre) & ttre < .
gen int h = t - ttre
qui gquantiles lat5 = exposure , xtile nquantiles(5)

foreach IV in state county {
  preserve
  gen double _ii = cond(post, `IV', 0)
  egen long _hh = group(h) if post
  qui replace _hh = 0 if missing(_hh)
  gen byte _tf = post
  tempvar t1 t2 t3
  qui gegen byte `t1' = tag(_ii _hh _tf)
  qui count if `t1' == 1
  local ncell = r(N)
  qui gegen byte `t2' = tag(_ii _hh _tf lat5)
  qui count if `t2' == 1
  local ncelllat = r(N)
  qui gegen byte `t3' = tag(`IV' h) if post
  qui count if `t3' == 1
  local cells = r(N)
  qui levelsof `IV', local(L)
  local ilev : word count `L'
  local bound = `cells'*2 + `ilev' + 20 + 200
  local a0 = 2*`ncell' + 1
  local a5 = 2*(`ncell' + `ncelllat') + 2*5 + 1
  dis "P4B i=`IV'  si8cells=" `cells' "  si8demand=" `bound' ///
      "  ncell=" `ncell' " ncell_x_lat=" `ncelllat' ///
      "  cols_nolat=" `a0' "  cols_lat5=" `a5' ///
      "  RATIO_actual_over_bound=" %6.3f `a5'/`bound'
  restore
}
dis "P4B_DONE"
