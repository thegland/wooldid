* probe09: permtstat and weights under the rebuild -- identity against hand recompute
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* ---------- A. permtstat ----------
use "dataB.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permtstat permsave(T_A) replace
dis "T_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " M=" e(perm_npossible) " p=" e(perm_p_att)
cap mat CB = e(perm_comparator_b)
cap mat list CB, format(%21x)
cap mat OB = e(b)
cap mat list OB, format(%21x)

use "T_A.dta", clear
desc
* hand-check three replications
use "T_A.dta", clear
forval r = 1/3 {
  local rr = `r' + 1
  forval s = 1/6 {
    local a`r'_`s' = onset_`s'[`rr']
  }
  local sb`r' : di %21x main_att[`rr']
  cap local ss`r' : di %21x se_main_att[`rr']
}
forval r = 1/3 {
  use "dataB.dta", clear
  gen double repons = .
  forval s = 1/6 {
    qui replace repons = `a`r'_`s'' if state == `s'
  }
  gen double repi = cond(repons==4,1,cond(repons==6,2,3))
  cap qui wooldid y repi t repons, att cluster(state) esrelativeto(-1)
  dis "TCHK rep `r' rc=" _rc
  dis "   saved b=`sb`r''  ext b=" %21x _b[main_att]
  dis "   saved se=`ss`r''  ext se=" %21x _se[main_att]
}

* ---------- B. weights ----------
use "dataB.dta", clear
gen double popw = 1 + 0.5*state
cap noi wooldid y coh t onset [aw=popw], att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(W_A) replace
dis "W_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att)
use "W_A.dta", clear
forval r = 1/3 {
  local rr = `r' + 1
  forval s = 1/6 {
    local w`r'_`s' = onset_`s'[`rr']
  }
  local wb`r' : di %21x main_att[`rr']
}
forval r = 1/3 {
  use "dataB.dta", clear
  gen double popw = 1 + 0.5*state
  gen double repons = .
  forval s = 1/6 {
    qui replace repons = `w`r'_`s'' if state == `s'
  }
  gen double repi = cond(repons==4,1,cond(repons==6,2,3))
  cap qui wooldid y repi t repons [aw=popw], att fastnose cluster(state) esrelativeto(-1)
  dis "WCHK rep `r' rc=" _rc "  saved b=`wb`r''  ext b=" %21x _b[main_att]
}

* ---------- C. customweightmultiplier ----------
use "dataB.dta", clear
gen double cwm = 1 + 0.25*state
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts customweightmultiplier(cwm) permsave(CW_A) replace
dis "CW_RC=" _rc " done=" e(perm_completed) " p=" e(perm_p_att)
use "CW_A.dta", clear
forval r = 1/2 {
  local rr = `r' + 1
  forval s = 1/6 {
    local c`r'_`s' = onset_`s'[`rr']
  }
  local cb`r' : di %21x main_att[`rr']
}
forval r = 1/2 {
  use "dataB.dta", clear
  gen double cwm = 1 + 0.25*state
  gen double repons = .
  forval s = 1/6 {
    qui replace repons = `c`r'_`s'' if state == `s'
  }
  gen double repi = cond(repons==4,1,cond(repons==6,2,3))
  cap qui wooldid y repi t repons, att fastnose cluster(state) esrelativeto(-1) customweightmultiplier(cwm)
  dis "CWCHK rep `r' rc=" _rc "  saved b=`cb`r''  ext b=" %21x _b[main_att]
}
