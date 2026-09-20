* p04: DIFFERENT specifications back to back -- does the first change the second?
* Every arm regenerates byte-identical data before the reference spec S, so any
* difference between arms is session state, not data.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* the reference spec S, dumped under a tag
cap prog drop wdS
prog def wdS
  args tag
  wdgen
  capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3)
  dis "S_RC_`tag'=" _rc
  wddumpe ds_`tag'
  dis "S|`tag'|vcvstate=" "`e(wdid_vcvstate)'"
  dis "S|`tag'|vblocks=" "`e(wdid_vblocks)'"
  dis "S|`tag'|fcol=" "`e(wdid_fcol)'" "|fcncell=" e(wdid_fcncell)
  wdline S_`tag'
end

* ============ baseline: S first, in a clean session
wdS 00_first

* ---- X1: verbose
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att verbose
dis "X1_RC=" _rc
capture noisily estimates dir
wdS 01_after_verbose
capture noisily estimates dir

* ---- X2: subgroups()
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att subgroups(sg)
dis "X2_RC=" _rc
wdS 02_after_subgroups

* ---- X3: poisson (leaves _ppmlhdfe_d scaffolding behind on failure paths)
wdgen
qui replace y = exp(y/3)
capture noisily wooldid y unit year ttre, cluster(state) att poisson
dis "X3_RC=" _rc
capture confirm variable _ppmlhdfe_d
dis "X3_PPML_LEFT=" (_rc == 0)
wdS 03_after_poisson

* ---- X4: nofastcollapse (the non-preserve path)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att nofastcollapse
dis "X4_RC=" _rc
wdline X4_post
wdS 04_after_nofastcollapse

* ---- X5: impvalues
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att impvalues
dis "X5_RC=" _rc
wdS 05_after_impvalues

* ---- X6: fastnose
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att fastnose
dis "X6_RC=" _rc
* known_bydesign s7 item 11: does fastnose post e(b_es_att) / e(numpre_att)?
capture confirm matrix e(b_es_att)
dis "X6_b_es_att_present=" (_rc == 0)
capture confirm matrix e(V_es_att)
dis "X6_V_es_att_present=" (_rc == 0)
dis "X6_numpre_att=" e(numpre_att)
wdS 06_after_fastnose

* ---- X7: level(80)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att level(80)
dis "X7_RC=" _rc
dis "X7_LEVEL_AFTER=" c(level)
dis "X7_e_level=" e(level)
wdS 07_after_level80

* ---- X8: a customvce run
wdgen
capture noisily wooldid y unit year ttre, att customvce(cluster state)
dis "X8_RC=" _rc
wdS 08_after_customvce

* ---- X9: unconditionalse
wdgen
capture noisily wooldid y unit year ttre, att unconditionalse
dis "X9_RC=" _rc
wdS 09_after_unconditionalse

* ---- X10: contreat
wdgen
qui gen double dose = 0
qui replace dose = 0.5 + 0.5*mod(unit,3) if !missing(ttre) & year >= ttre
capture noisily wooldid y unit year ttre, cluster(state) att contreat(dose)
dis "X10_RC=" _rc
wdS 10_after_contreat

dis "P04_DONE"
