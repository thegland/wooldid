* p09: simultaneousbands RNG neutrality, e(sample) sanity, level() interplay
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* ============ ARM A: simultaneousbands is documented (ado:10888) as saving and
* restoring the caller's RNG state.  Measure it.
wdgen
global WDH_R1 "`c(rngstate)'"
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) simultaneousbands
dis "A_RC=" _rc
global WDH_R2 "`c(rngstate)'"
dis "A_RNG_PRESERVED=" ("$WDH_R1" == "$WDH_R2")
dis "A_esband_engine=" "`e(esband_engine)'"
dis "A_esband_families=" "`e(esband_families)'"
dis "A_esband_ndraw=" e(esband_ndraw)

* the band must not depend on how many draws the caller had consumed beforehand
wdgen
set seed 987654
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) simultaneousbands
dis "A2_RC=" _rc
capture confirm matrix e(esband_att)
if _rc == 0 {
  tempname EB
  matrix `EB' = e(esband_att)
  forvalues a = 1/`=rowsof(`EB')' {
    forvalues b = 1/`=colsof(`EB')' {
      dis "A2_esband[" `a' "," `b' "]=" %21x `EB'[`a',`b']
    }
  }
}
dis "A2_crit_att=" %21x e(esband_crit_att)

wdgen
set seed 13579
capture noisily wooldid y unit year ttre, cluster(state) att esprelength(3) espostlength(3) simultaneousbands
dis "A3_RC=" _rc
dis "A3_crit_att=" %21x e(esband_crit_att)

* ============ ARM B: e(sample) sanity after the fastcollapse preserve/restore
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att
dis "B_RC=" _rc
dis "B_e_N=" e(N)
qui count if e(sample)
dis "B_count_esample=" r(N)
qui count
dis "B_rows_in_memory=" r(N)
* the same under nofastcollapse
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att nofastcollapse
dis "B2_RC=" _rc
dis "B2_e_N=" e(N)
qui count if e(sample)
dis "B2_count_esample=" r(N)

* ============ ARM C: level() interplay with a non-default session level
wdgen
set level 80
dis "C_LEVEL_SESSION=" c(level)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "C_RC=" _rc
dis "C_e_level=" e(level)
dis "C_LEVEL_AFTER=" c(level)
capture noisily wooldid y unit year ttre, cluster(state) att level(99)
dis "C2_RC=" _rc
dis "C2_e_level=" e(level)
dis "C2_LEVEL_AFTER=" c(level)
set level 95

* ============ ARM D: estimand changed between back-to-back runs
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att att_i att_it att_itime
dis "D_RC=" _rc
wddumpe misc_all4
capture noisily wooldid y unit year ttre, cluster(state) att
dis "D2_RC=" _rc
wddumpe misc_att_after_all4
capture noisily ereturn list

* ============ ARM E: the outcome is missing on every never-treated row.
* A very common administrative-data shape (outcome only measured for participants).
* What does wooldid tell the user?
wdgen
qui replace y = . if missing(ttre)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E_RC=" _rc

* milder: missing on 30% of the control rows only
wdgen
set seed 24680
qui replace y = . if missing(ttre) & runiform() < 0.3
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E2_RC=" _rc
dis "E2_N=" e(N)
dis "E2_N_dropped=" e(N_dropped)

* and: outcome missing on every row of ONE calendar year
wdgen
qui replace y = . if year == 2006
capture noisily wooldid y unit year ttre, cluster(state) att
dis "E3_RC=" _rc
dis "E3_N=" e(N)
dis "E3_N_dropped=" e(N_dropped)

dis "P09_DONE"
