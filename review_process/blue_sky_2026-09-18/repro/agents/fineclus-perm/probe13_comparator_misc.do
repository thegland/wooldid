* probe13: (a) e(perm_comparator_b) vs e(b) under the rebuild
*          (b) rows with a MISSING cluster value (the ado calls the sentinel unreachable)
*          (c) odd i codes (negative / non-contiguous)
*          (d) jointtests + event study under the rebuild
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

prog def cmpb
  mat CB = e(perm_comparator_b)
  mat OB = e(b)
  local k = colsof(OB)
  local worst = 0
  forval j = 1/`k' {
    local d = abs(CB[1,`j'] - OB[1,`j'])
    if `d' > `worst' local worst = `d'
  }
  dis "   COMPARATOR_MAXABSDIFF=" %12.4e `worst' "   K=" `k'
  local n1 : colfullnames OB
  local n2 : colfullnames CB
  dis "   COMPARATOR_NAMES_MATCH=" cond("`n1'"=="`n2'","YES","NO")
end

* ---------- (a) comparator identity, fast mode and permanalyticse ----------
use "dataB.dta", clear
cap noi qui wooldid y coh t onset, att att_i esprelength(2) espostlength(2) cluster(state) permute(shuffle) permexact permonsetcohorts
dis "A1 rc=" _rc
cmpb
use "dataB.dta", clear
cap noi qui wooldid y coh t onset, att att_i esprelength(2) espostlength(2) cluster(state) permute(shuffle) permexact permonsetcohorts permanalyticse
dis "A2 (permanalyticse) rc=" _rc
cmpb
use "dataB.dta", clear
cap noi qui wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permtstat
dis "A3 (permtstat) rc=" _rc
cmpb

* ---------- (b) missing cluster values ----------
use "dataB.dta", clear
replace state = . if state == 1
cap noi wooldid y coh t onset, att cluster(state) fastnose
dis "B_OBS rc=" _rc " N=" e(N) " b=" %21x _b[main_att]
use "dataB.dta", clear
replace state = . if state == 1
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
dis "B_PERM rc=" _rc " N=" e(N) " b=" %21x _b[main_att] " p=" e(perm_p_att) " irebuilt=" e(perm_irebuilt)
cap cmpb
* the hand-built rep-0 the ado constructs for such rows (repi = ncode for the unlinked rows)
use "dataB.dta", clear
replace state = . if state == 1
gen double repons = onset
replace repons = . if missing(state)
gen double repi = coh
replace repi = 3 if missing(state)
cap noi wooldid y repi t repons, att cluster(state) fastnose esrelativeto(-1)
dis "B_HAND rc=" _rc " N=" e(N) " b=" %21x _b[main_att]

* ---------- (c) odd i codes ----------
use "dataB.dta", clear
replace coh = cond(coh==1,-7,cond(coh==2,0,55))
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
dis "C_ODDI rc=" _rc " p=" e(perm_p_att) " b=" %21x _b[main_att] " irebuilt=" e(perm_irebuilt)
* never-treated code NOT the minimum (55) - wdp_ncode is min(i) among the never-treated
use "dataB.dta", clear
replace coh = cond(coh==1,20,cond(coh==2,30,10))
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
dis "C_ODDI2 rc=" _rc " p=" e(perm_p_att) " b=" %21x _b[main_att]

* ---------- (d) jointtests + event study ----------
use "dataB.dta", clear
cap noi wooldid y coh t onset, att esprelength(2) espostlength(2) jointtests cluster(state) permute(shuffle) permexact permonsetcohorts
dis "D_JT rc=" _rc " done=" e(perm_completed) " fail=" e(perm_failed)
