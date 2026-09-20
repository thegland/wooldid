* probe08: permstrata + permonsetcohorts composition; permprobs + rebuild; guard checks
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* ---------- A. permstrata() x permonsetcohorts ----------
* design B: s1,s2 onset 4 ; s3,s4 onset 6 ; s5,s6 never
* strata: region 1 = {s1,s3,s5}, region 2 = {s2,s4,s6}
* each stratum holds one of each onset -> 3! x 3! = 36 assignments
use "dataB.dta", clear
gen byte region = cond(inlist(state,1,3,5),1,2)
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permstrata(region) permonsetcohorts permsave(ST_A) replace
dis "ST_A_RC=" _rc " M=" e(perm_npossible) " (hand count 36) done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att)
* verify the strata invariant survives into the drawn assignments
use "ST_A.dta", clear
gen byte ok = 1
forval r = 1/`=_N' {
  local o1 = onset_1[`r']
  local o2 = onset_2[`r']
  local o3 = onset_3[`r']
  local o4 = onset_4[`r']
  local o5 = onset_5[`r']
  local o6 = onset_6[`r']
  * region 1 = clusters 1,3,5 ; region 2 = 2,4,6 ; each must hold exactly one of {4,6,.}
  local n4a = (`o1'==4) + (`o3'==4) + (`o5'==4)
  local n6a = (`o1'==6) + (`o3'==6) + (`o5'==6)
  local nva = missing(`o1') + missing(`o3') + missing(`o5')
  if `n4a' != 1 | `n6a' != 1 | `nva' != 1 qui replace ok = 0 in `r'
}
count if ok == 0
dis "ST_STRATA_VIOLATIONS=" r(N)

* ---------- B. permprobs() x permonsetcohorts ----------
use "dataB.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(probabilities) permprobs(0.34 0.33 0.33) permreps(200) rseed(99) permonsetcohorts permsave(PR_A) replace
dis "PR_RC=" _rc " done=" e(perm_completed) " fail=" e(perm_failed) " p=" e(perm_p_att)
cap mat list e(perm_n)
cap use "PR_A.dta", clear
cap count if missing(main_att)
dis "PR_NMISSING=" r(N)

* ---------- C. guard contradiction checks ----------
* C1: i splits the never-treated across two codes -> must be refused
use "dataB.dta", clear
replace coh = 4 if state == 6
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
dis "C1_RC=" _rc

* C2: onset varies within i -> must be refused
use "dataB.dta", clear
replace onset = 5 if state == 2
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
dis "C2_RC=" _rc

* C3: cluster finer than i WITHOUT permonsetcohorts -> must be refused (459)
use "dataB.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact
dis "C3_RC=" _rc

* C4: subgroups() with permonsetcohorts -> must be refused (198)
use "dataB.dta", clear
gen byte sg = mod(state,2)
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts subgroups(sg)
dis "C4_RC=" _rc

* C5: permonsetcohorts when cluster() IS i -> rebuild must be a no-op, irebuilt=0
use "dataB.dta", clear
cap noi wooldid y coh t onset, att cluster(coh) permute(shuffle) permexact permonsetcohorts
dis "C5_RC=" _rc " irebuilt=" e(perm_irebuilt) " p=" e(perm_p_att) " b=" %21x _b[main_att]
use "dataB.dta", clear
cap noi wooldid y coh t onset, att cluster(coh) permute(shuffle) permexact
dis "C5b_RC=" _rc " p=" e(perm_p_att) " b=" %21x _b[main_att]
