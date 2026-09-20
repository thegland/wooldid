* p17: known_bydesign s7 item 11 -- sthlp:1325 says that when NO horizon has an
* identified variance (naming fastnose, timetrends, and an all-zero covariance)
* "the covariance is not returned but the coefficient vector still is, under a
* Caution". The 2026-09-13 ruling was to REMOVE fastnose (and impvalues) from that
* list; sthlp:1325 still opens it with {cmd:fastnose}. Settle it by measurement,
* WITH an event-study window so the absence cannot be blamed on "no event study".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

cap prog drop wdes
prog def wdes
  args tag
  foreach m in b_es_att V_es_att b_es_att_it V_es_att_it {
    capture confirm matrix e(`m')
    dis "ES|`tag'|`m'_present=" (_rc == 0)
  }
  dis "ES|`tag'|numpre_att=" e(numpre_att)
  dis "ES|`tag'|numpost_att=" e(numpost_att)
  dis "ES|`tag'|es_ndegen_att=" e(es_ndegen_att)
  dis "ES|`tag'|vcvstate=" "`e(wdid_vcvstate)'"
  dis "ES|`tag'|ES_cols_in_b=" strpos("`: colfullnames e(b)'", "ES_")
end

* ---- 1 reference: full analytic run with an ES window
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3)
dis "R1_RC=" _rc
wdes 1_reference

* ---- 2 fastnose, SAME window
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) fastnose
dis "R2_RC=" _rc
wdes 2_fastnose

* ---- 3 impvalues, SAME window
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) impvalues
dis "R3_RC=" _rc
wdes 3_impvalues

* ---- 4 timetrends, SAME window (the other member of the sthlp:1325 list)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) timetrends
dis "R4_RC=" _rc
wdes 4_timetrends

* ---- 5 fastnose with att_it as well
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att att_it esprelength(3) espostlength(3) fastnose
dis "R5_RC=" _rc
wdes 5_fastnose_attit

* ---- 6 fastnose under bootstrap is refused (documented); permute + fastnose?
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) permute(shuffle) permreps(20) dots(0) rseed(7)
dis "R6_RC=" _rc
wdes 6_permute

dis "P17_DONE"
