*==========================================================================
* e1_exact.do -- SKEPTIC's own repro of PI-1, built from the CLAIM, with the
* Monte-Carlo confound REMOVED.
*
* PI-1's evidence was 30 random replications on two datasets whose cluster
* COUNT differs (10 vs 9), so their shuffle streams differ: a p-value gap at
* 30 draws is not evidence that the reference DISTRIBUTION differs.
* permexact enumerates every distinct assignment, so there is no RNG at all
* and any difference is structural.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do e1_exact.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*--- my own panel: 8 states x 3 counties, 8 annual periods, state-level onset
*    states 1..7 carry data; state 8 is the PHANTOM.
clear
set seed 31337
set sortseed 31337
qui set obs 24
qui gen long unit  = _n
qui gen long state = ceil(unit/3)
qui gen double ai  = rnormal()
qui gen double ttre = .
qui replace ttre = 2004 if state == 1
qui replace ttre = 2006 if state == 2
qui expand 8
qui bysort unit: gen int t = 1999 + _n
qui gen double x1   = rnormal()
qui gen double eps  = rnormal()*0.4
qui gen byte   post = (t >= ttre) & !missing(ttre)
qui gen double y    = ai + 0.15*(t-2000) + 0.5*x1 + 0.9*post + eps
qui drop eps
tempfile base
save `base'

local OPT att cluster(state) permute(shuffle) permexact dots(0)

prog def rpt
  args tag
  dis "RES `tag' rc=" _rc
  if _rc == 0 {
    dis "RES `tag' N="        e(N)      " Nclust=" e(N_clust)
    dis "RES `tag' npossible=" e(perm_npossible) " completed=" e(perm_completed) " failed=" e(perm_failed)
    dis "RES `tag' b="  %21x e(wdidmainresults)[1,1]
    dis "RES `tag' p="  %21x e(perm_p_att)
    dis "RES `tag' lo=" %21x e(perm_lb_att) " hi=" %21x e(perm_ub_att)
    dis "RES `tag' pdec=" e(perm_p_att) " lodec=" e(perm_lb_att) " hidec=" e(perm_ub_att)
  }
end

*========== ARM 1: phantom is NEVER-TREATED, y missing on every row ==========
* A : phantom present, y missing everywhere on it
use `base', clear
qui replace y = . if state == 8
cap noi wooldid y unit t ttre , `OPT'
rpt A_ymiss

* B : phantom's ROWS DROPPED  (identical estimation sample to A)
use `base', clear
qui drop if state == 8
cap noi wooldid y unit t ttre , `OPT'
rpt B_rowdrop

* C : phantom's CLUSTER value set missing (rows kept, y intact)
use `base', clear
qui replace state = . if state == 8
cap noi wooldid y unit t ttre , `OPT'
rpt C_clusmiss

* D : phantom removed by the user's if (rows kept, y intact)
use `base', clear
cap noi wooldid y unit t ttre if state != 8 , `OPT'
rpt D_ifdrop

*========== ARM 2: phantom killed by a MISSING CONTROL, not a missing y =====
* realistic: a covariate that one state does not collect.
* E : phantom present, x1 missing on every one of its rows, controls(x1)
use `base', clear
qui replace x1 = . if state == 8
cap noi wooldid y unit t ttre , `OPT' controls(x1)
rpt E_xmiss

* Edrop : same estimation sample, phantom rows dropped
use `base', clear
qui drop if state == 8
cap noi wooldid y unit t ttre , `OPT' controls(x1)
rpt E_xdrop

*========== ARM 3: phantom is the SOLE carrier of an onset level ===========
* J : state 8 is the ONLY 2008 adopter and its y is missing everywhere, so the
*     estimation sample contains no 2008 cohort at all.
use `base', clear
qui replace ttre = 2008 if state == 8
qui replace y = . if state == 8
cap noi wooldid y unit t ttre , `OPT'
rpt J_soleonset

* Jdrop : same estimation sample, phantom rows dropped
use `base', clear
qui drop if state == 8
cap noi wooldid y unit t ttre , `OPT'
rpt J_drop

dis "DONE_E1"
