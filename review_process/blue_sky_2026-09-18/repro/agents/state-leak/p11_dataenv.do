* p11: the user's DATA ENVIRONMENT -- xtset/tsset/svyset state, notes, labels,
* nested preserve, repeated makeplots
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: xtset survives a run (default, fastcollapse ON)
wdgen
qui xtset unit year
dis "A_TSPANEL_BEFORE=" "`: char _dta[_TSpanel]'"
dis "A_TSTVAR_BEFORE=" "`: char _dta[_TStvar]'"
capture noisily wooldid y unit year ttre, cluster(state) att
dis "A_RC=" _rc
dis "A_TSPANEL_AFTER=" "`: char _dta[_TSpanel]'"
dis "A_TSTVAR_AFTER=" "`: char _dta[_TStvar]'"
capture noisily xtset
dis "A_XTSET_RC=" _rc
capture noisily xtreg y post, fe
dis "A_XTREG_RC=" _rc

* ============ ARM B: same under nofastcollapse (no preserve/restore)
wdgen
qui xtset unit year
capture noisily wooldid y unit year ttre, cluster(state) att nofastcollapse
dis "B_RC=" _rc
dis "B_TSPANEL_AFTER=" "`: char _dta[_TSpanel]'"
capture noisily xtset
dis "B_XTSET_RC=" _rc

* ============ ARM C: notes on _dta and on a variable
wdgen
note _dta: "this panel came from the 2019 extract"
note y: "winsorised at the 99th percentile"
dis "C_DTANOTES_BEFORE=" "`: char _dta[note1]'"
dis "C_YNOTE_BEFORE=" "`: char y[note1]'"
capture noisily wooldid y unit year ttre, cluster(state) att
dis "C_RC=" _rc
dis "C_DTANOTES_AFTER=" "`: char _dta[note1]'"
dis "C_YNOTE_AFTER=" "`: char y[note1]'"

* ============ ARM D: value labels, variable labels, formats survive
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att subgroups(sg)
dis "D_RC=" _rc
dis "D_SGVALLBL=" "`: value label sg'"
dis "D_YLBL=" "`: variable label y'"
dis "D_YFMT=" "`: format y'"
dis "D_TTREFMT=" "`: format ttre'"
capture noisily label dir
capture noisily label list sglbl

* ============ ARM E: a pre-existing _ppmlhdfe_d on a NON-poisson run
wdgen
qui gen double _ppmlhdfe_d = 12345
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E_RC=" _rc
capture confirm variable _ppmlhdfe_d
dis "E_PPML_STILL_THERE=" (_rc == 0)
if _rc == 0 {
  qui sum _ppmlhdfe_d, meanonly
  dis "E_PPML_MEAN=" r(mean)
}

* ============ ARM F: makeplots twice back to back (graph-name collisions)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) makeplots
dis "F1_RC=" _rc
capture noisily graph dir
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) makeplots
dis "F2_RC=" _rc
capture noisily graph dir
* and with saveplots to the SAME name twice, no replace
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) makeplots saveplots(p11plots)
dis "F3_RC=" _rc
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(2) espostlength(2) makeplots saveplots(p11plots)
dis "F4_RC=" _rc
capture noisily dir "p11*"

* ============ ARM G: wooldid called from inside a user program that has
* already issued `preserve' (nested preserve across program boundaries)
cap prog drop userwrap
prog def userwrap
  preserve
  qui keep if year >= 2003
  capture noisily wooldid y unit year ttre, cluster(state) att
  dis "G_INNER_RC=" _rc
  restore
end
wdgen
capture noisily userwrap
dis "G_RC=" _rc
wdline G_post
qui count
dis "G_ROWS_BACK=" r(N)

* ============ ARM H: string variables present in the data
wdgen
qui gen str12 countyname = "county" + string(unit)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "H_RC=" _rc
capture confirm variable countyname
dis "H_STRVAR_INTACT=" (_rc == 0)
wdline H_post

dis "P11_DONE"
