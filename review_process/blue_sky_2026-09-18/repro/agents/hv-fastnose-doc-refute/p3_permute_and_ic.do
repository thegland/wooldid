* =========================================================================
* P3 -- the two follow-ups the FN2 claim turns on:
*  (A) permute() fast-mode default on the SAME d3 design: does the
*      structural-zero ame reach a REPORTABLE inference table, and is it named?
*  (B) the interactivecontrols() twin of the same screen (BS-3, ado:12005),
*      which carries the identical `fastnose' gate: analytic vs fastnose.
* Self-contained.
* =========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 12800
gen long row = _n
gen int  st  = ceil(row/320)
gen int  yr  = 2005 + mod(ceil(row/20)-1, 16)
gen int adopt = .
replace adopt = 2010 if inrange(st, 1, 10)
replace adopt = 2013 if inrange(st, 11, 20)
replace adopt = 2016 if inrange(st, 21, 28)
bysort st: gen double u_st = runiform() if _n == 1
bysort st (u_st): replace u_st = u_st[1]
sort row
gen double dose = round(0.5 + 2.5*u_st, 0.01)
gen byte D = (yr >= adopt) & !missing(adopt)
bysort st: gen double afe = rnormal() if _n == 1
bysort st (afe): replace afe = afe[1]
sort row
gen double y = 2 + afe + 0.03*(yr-2005) + 0.40*dose*D + rnormal(0,1)

* a control that is CONSTANT within each cohort-period cell (the BS-3 geometry):
* state-level unemployment, identical for every individual in a state-year
bysort st yr: gen double urate = 5 + rnormal(0,1) if _n == 1
bysort st yr (urate): replace urate = urate[1]
sort row

dis "=== (A1) permute() fast-mode default, d3 design ==="
cap noi wooldid y st yr adopt, att cluster(st) contreat(dose) contreatcontrols(basic) ///
    esprelength(3) espostlength(4) esfixedbaseperiod permute(20)
dis "PA_RC=" _rc
cap confirm matrix e(b_es_ame)
dis "PA_BESAME_RC=" _rc
dis "PA_NUMPRE_AME=[" e(numpre_ame) "]"
cap noi mat li e(wdidmainresults), format(%12.6f)

dis "=== (B1) interactivecontrols, analytic ==="
cap noi wooldid y st yr adopt, att robust interactivecontrols(urate) esprelength(3) espostlength(4) esfixedbaseperiod
dis "B1_RC=" _rc

dis "=== (B2) interactivecontrols, fastnose ==="
cap noi wooldid y st yr adopt, att robust interactivecontrols(urate) esprelength(3) espostlength(4) esfixedbaseperiod fastnose
dis "B2_RC=" _rc

dis "P3_DONE"
