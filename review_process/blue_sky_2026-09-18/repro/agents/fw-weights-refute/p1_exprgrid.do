clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

* ---- BUILD ASSERTION ----
which wooldid
capture noisily wooldid
dis "BUILDCHECK_RC=" _rc

* ---- DATA: 40 units x 10 periods, staggered onset, plausible survey weights ----
set obs 40
gen long i = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 4 if grp==1
replace ttre = 6 if grp==2
replace ttre = 8 if grp==3
gen double ai = rnormal()
expand 10
bys i: gen int t = _n
gen double y = ai + 0.3*t + rnormal() + 1.5*(t>=ttre & ttre<.)
* two plausible sub-weights (male/female headcount) that sum to the total
gen double wf = 20 + 10*runiform()
gen double wm = 15 + 12*runiform()
gen double wtot = wf + wm
gen double zero = 0
gen double pscore = 0.2 + 0.6*runiform()
gen double k = 1 + 4*runiform()
compress
save "p1dat.dta", replace

* ---- GRID: each expression, three engines ----
local exprs `" "wtot" "2*k" "k/2" "k*2" "int(k)" "wf+wm" "(wf+wm)" "1/pscore" "cond(k>2,k,1)" "'
local n 0
foreach e of local exprs {
  local ++n
  dis _n "================ EXPR`n' = [`e'] ================"
  capture noisily wooldid y i t ttre [aw=`e'], att att_it
  dis "WD_DEFAULT_`n'_RC=" _rc
  capture noisily wooldid y i t ttre [aw=`e'], att att_it nofastcollapse
  dis "WD_NOFC_`n'_RC=" _rc
  capture noisily wooldid y i t ttre [aw=`e'], att customvce(robust)
  dis "WD_CVCE_`n'_RC=" _rc
  qui gen byte post`n' = (t>=ttre & ttre<.)
  capture noisily reghdfe y post`n' [aw=`e'], absorb(i t)
  dis "RHDFE_`n'_RC=" _rc
  capture noisily regress y post`n' i.i i.t [aw=`e']
  dis "REG_`n'_RC=" _rc
}
dis "P1_DONE"
