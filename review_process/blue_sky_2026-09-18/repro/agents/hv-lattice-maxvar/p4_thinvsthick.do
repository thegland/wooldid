*=====================================================================
* p4_thinvsthick.do -- SELF-FOUND. Pure counting, no estimation, seconds.
* Where is the preflight's lattice blind spot WORST?
* Answer: on COARSE i with THICK cells, i.e. exactly the shape ado:1371
* and sthlp:1163 tell the user to move to (fe() + a coarser i).
* With i = unit (one row per cohort-period cell) `emptycells drop'
* collapses the lattice interaction and the bound is fine.
* SELF-CONTAINED.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set maxvar 40000
set emptycells drop

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
  tempvar tg
  egen `tg' = tag(`IV' h) if post
  qui count if `tg' == 1
  local cells = r(N)
  qui levelsof `IV', local(L)
  local ilev : word count `L'
  qui levelsof t, local(L2)
  local tlev : word count `L2'
  local bound = `cells'*2 + `ilev' + `tlev' + 200
  qui fvexpand i._ii#i._hh#i._tf##c.exposure
  local a0 = wordcount("`r(varlist)'")
  qui fvexpand i._ii#i._hh#i._tf##i.lat5##c.exposure
  local a5 = wordcount("`r(varlist)'")
  dis "P4 i=`IV'  cells=" `cells' "  si8demand=" `bound' ///
      "  actual_nolat=" `a0' "  actual_lat5=" `a5' ///
      "  ratio_actual_over_bound=" %6.3f `a5'/`bound'
  restore
}
dis "P4_DONE"
