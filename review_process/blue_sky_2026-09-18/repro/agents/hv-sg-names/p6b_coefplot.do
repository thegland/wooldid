*==========================================================================
* p6b_coefplot.do -- exactly how many coefficients does coefplot plot when
*   handed the help's non-matching drop(*_s1*) pattern?  SELF-CONTAINED.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 120
gen long unit = _n
gen byte grp = 1 + mod(_n, 12)
gen byte sgA = 1 + mod(grp, 2)
gen int ttre = .
replace ttre = 2003 if inrange(grp,1,3)
replace ttre = 2006 if inrange(grp,4,6)
gen double g_fe = rnormal()
bysort grp (unit): replace g_fe = g_fe[1]
expand 10
bysort unit: gen int year = 1999 + _n
gen double post = (year >= ttre) & ttre < .
gen double y = g_fe + 0.2*(year-2000) + post*1.0 + 0.25*rnormal()

wooldid y grp year ttre, att subgroups(sgA) esprelength(1) espostlength(1) robust
local nb = colsof(e(b))
di "NCOEF=" `nb'
local cn : colnames e(b)
di "COLNAMES=`cn'"
est store W

* --- coefplot with the HELP's pattern ---
preserve
cap noi coefplot W, drop(*_s1*) nodraw generate(cpA)
di "DROP_HELPNAME_RC=" _rc
qui ds cpA*
di "DROP_HELPNAME_VARS=`r(varlist)'"
cap qui count if !mi(cpAb)
di "DROP_HELPNAME_NPLOTTED=" r(N)
restore

* --- coefplot with the TRUE pattern ---
preserve
cap noi coefplot W, drop(*SG1*) nodraw generate(cpB)
di "DROP_TRUENAME_RC=" _rc
qui ds cpB*
di "DROP_TRUENAME_VARS=`r(varlist)'"
cap qui count if !mi(cpBb)
di "DROP_TRUENAME_NPLOTTED=" r(N)
restore

* --- coefplot with nothing dropped, as the benchmark ---
preserve
cap noi coefplot W, nodraw generate(cpC)
di "NOOPT_RC=" _rc
qui ds cpC*
di "NOOPT_VARS=`r(varlist)'"
cap qui count if !mi(cpCb)
di "NOOPT_NPLOTTED=" r(N)
restore

di "PROBE6B_DONE"
