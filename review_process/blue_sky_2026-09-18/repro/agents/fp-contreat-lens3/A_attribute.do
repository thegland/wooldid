* ==========================================================================
* LENS 3 / probe A -- ATTRIBUTION + WARNING AUDIT, single draw
* Rebuilds the prober's exact fixture (seed 121212) so the claimed numbers
* (ame 18.091, customvce se .43558, unconditionalse se 10.04) are on the table,
* then asks: is wooldid's customvce(robust) ame SE the SAME object that a hand
* regress + margins delta method produces on the identical model?  And does the
* default path print ANY caution when its ame SE is missing?
* ==========================================================================
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 121212
set sortseed 121212

clear
set obs 12
gen int st = _n
gen int nc = 5 + mod(_n*7,8)
gen double ttre = .
qui replace ttre =  8 if inlist(st,1,2,3)
qui replace ttre = 12 if inlist(st,4,5,6)
qui replace ttre = 15 if inlist(st,7,8)
expand nc
gen long unit = _n
expand 20
bysort unit: gen int year = _n
bysort unit (year): gen double ubase = exp(rnormal(0,0.6)) if _n==1
bysort unit (year): replace ubase = ubase[1]
gen double dose_it = ubase*(1 + 0.10*rnormal())
bysort unit (year): gen double ufe = rnormal(0,1) if _n==1
bysort unit (year): replace ufe = ufe[1]
gen double y_dose = 3 + ufe + 0.2*year + 0.5*dose_it + rnormal(0,1)
qui replace y_dose = y_dose + 0.4*dose_it + 0.3 if ttre<. & year>=ttre
tempfile F
qui save `F'

*--- 1. default path, NOISY: does anything caution? -----------------------
dis "=====BEGIN_DEFAULT_NOISY====="
use `F', clear
cap noi wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) robust
dis "DEFRC=" _rc
mat M = e(wdidmainresults)
local rn : rownames M
dis "ROWNAMES: `rn'"
mat V = e(V)
mata: st_local("mv", strofreal(max(abs(st_matrix("V")))))
dis "DEFAULT max|e(V)|=`mv'"
dis "=====END_DEFAULT_NOISY====="

*--- 2. customvce(robust) and unconditionalse, values by rowname ----------
foreach sp in "customvce(robust)" "unconditionalse robust" {
  use `F', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) `sp'
  dis "SPEC=`sp' RC=" _rc
  if _rc==0 {
    mat M = e(wdidmainresults)
    local rn : rownames M
    local k = 0
    foreach r of local rn {
      local ++k
      dis "   ROW `r'  b=" %14.8f M[`k',1] "  se=" %14.8g M[`k',2]
    }
  }
}

*--- 3. HAND model on the identical data: is 0.4356 the delta-method HC1? -
use `F', clear
gen byte post = (ttre<. & year>=ttre)
qui egen cell = group(st year) if post==1
qui replace cell = 0 if cell>=.
qui xtile lat = dose_it , nq(3)
gen double ringer = dose_it
gen byte nullv = 1

foreach v in robust ols hc3 {
  local vceopt "vce(`v')"
  if "`v'"=="ols" local vceopt ""
  cap qui reg y_dose i.st i.year i.nullv##c.ringer i.lat##c.ringer ib0.cell##i.lat##c.dose_it , `vceopt'
  dis "HANDREG `v' RC=" _rc "  N=" e(N) " rank=" e(rank) " df_r=" e(df_r) " rss=" %14.4f e(rss)
  if _rc==0 {
    cap qui margins , dydx(dose_it) subpop(if post==1) noestimcheck
    if _rc==0 {
      mat b = r(b)
      mat Vm = r(V)
      dis "   HANDMARG `v' (noestimcheck):  dydx=" %14.8f b[1,1] "  se=" %14.8g sqrt(Vm[1,1])
    }
    else dis "   HANDMARG `v' (noestimcheck) RC=" _rc
    cap qui margins , dydx(dose_it) subpop(if post==1)
    if _rc==0 {
      mat b = r(b)
      mat Vm = r(V)
      dis "   HANDMARG `v' (ESTIMCHECK ON):  dydx=" %14.8f b[1,1] "  se=" %14.8g sqrt(Vm[1,1])
    }
    else dis "   HANDMARG `v' (ESTIMCHECK ON) RC=" _rc
  }
}

*--- 4. how bad is the conditioning?  within-(cell,bin) dose spread --------
use `F', clear
gen byte post = (ttre<. & year>=ttre)
qui xtile lat = dose_it , nq(3)
qui keep if post==1
qui egen grp = group(st year lat)
qui by grp, sort: gen int ng = _N
qui by grp: gen byte first = (_n==1)
qui egen double sdd = sd(dose_it), by(grp)
qui egen double mnd = mean(dose_it), by(grp)
qui sum ng if first, detail
dis "CELLBIN groups=" r(N) "  n per group: p10=" r(p10) " p50=" r(p50) " p90=" r(p90) " min=" r(min)
qui count if first & ng<=2
dis "CELLBIN groups with n<=2: " r(N)
qui gen double cv = sdd/mnd
qui sum cv if first, detail
dis "CELLBIN within-group dose CV: p10=" %8.5f r(p10) " p50=" %8.5f r(p50) " p90=" %8.5f r(p90)
dis "DONE_A"
