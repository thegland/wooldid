* probe19: permstepdown under the rebuild; e(perm_n) completeness; timetrends coherence
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

use "dataB.dta", clear
cap noi qui wooldid y coh t onset, att esprelength(2) espostlength(2) cluster(state) permute(shuffle) permexact permonsetcohorts permstepdown permtstat
dis "SD_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed)
cap mat PN = e(perm_n)
cap noi mat list PN
cap mat SD = e(perm_p_stepdown)
cap noi mat list SD
cap mat PP = e(perm_p)
cap noi mat list PP

* timetrends coherence: is the internal replication the same as the hand-built one?
use "dataB.dta", clear
cap noi wooldid y coh t onset, att cluster(state) timetrends permute(shuffle) permexact permonsetcohorts permsave(TT_A) replace
dis "TT_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att)
use "TT_A.dta", clear
forval r = 1/3 {
  local rr = `r' + 1
  forval s = 1/6 {
    local q`r'_`s' = onset_`s'[`rr']
  }
  local tb`r' : di %21x main_att[`rr']
}
forval r = 1/3 {
  use "dataB.dta", clear
  gen double repons = .
  forval s = 1/6 {
    qui replace repons = `q`r'_`s'' if state == `s'
  }
  gen double repi = cond(repons==4,1,cond(repons==6,2,3))
  cap qui wooldid y repi t repons, att fastnose cluster(state) timetrends esrelativeto(-1)
  dis "TTCHK rep `r' rc=" _rc " saved=`tb`r''  ext=" %21x _b[main_att]
}
