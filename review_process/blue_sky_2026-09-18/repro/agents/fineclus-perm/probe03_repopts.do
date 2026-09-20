* probe03: does repopts carry ORIGINAL-i-referencing options into a rebuilt-i replication?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

use "dataA.dta", clear
gen int cohcopy = coh
* a cohort-level control: constant within the OBSERVED cohort, varies across cohorts
gen double cohctrl = 10*coh
* a state-level control (finer than i) - the realistic case
gen double statectrl = state*1.0
save "dataA2.dta", replace

* ---------- R1 baseline ----------
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(p_base) replace
dis "R1_RC=" _rc
dis "R1_B=" %21x _b[main_att] "  R1_P=" e(perm_p_att)

* ---------- R2 fe(state): finer than i, the realistic spec ----------
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fe(state) permute(shuffle) permexact permonsetcohorts permsave(p_fest) replace
dis "R2_RC=" _rc
dis "R2_B=" %21x _b[main_att] "  R2_P=" e(perm_p_att)

* ---------- R3 fe(cohcopy): a copy of i, collinear in the observed run ----------
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fe(cohcopy) permute(shuffle) permexact permonsetcohorts permsave(p_fecoh) replace
dis "R3_RC=" _rc
dis "R3_B=" %21x _b[main_att] "  R3_P=" e(perm_p_att)

* ---------- R3b the same fe(cohcopy) WITHOUT permute: is the point estimate a no-op? ----------
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fastnose
dis "R3b_NOOPT_B=" %21x _b[main_att]
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) fe(cohcopy) fastnose
dis "R3b_FECOH_B=" %21x _b[main_att]

* ---------- R4 controls(cohctrl): constant within observed cohort ----------
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) controls(cohctrl) permute(shuffle) permexact permonsetcohorts permsave(p_cctrl) replace
dis "R4_RC=" _rc
dis "R4_B=" %21x _b[main_att] "  R4_P=" e(perm_p_att)
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) controls(cohctrl) fastnose
dis "R4b_CCTRL_B=" %21x _b[main_att]

* ---------- R5 controls(statectrl): the realistic finer-than-i control ----------
use "dataA2.dta", clear
cap noi wooldid y coh t onset, att cluster(state) controls(statectrl) permute(shuffle) permexact permonsetcohorts permsave(p_sctrl) replace
dis "R5_RC=" _rc
dis "R5_B=" %21x _b[main_att] "  R5_P=" e(perm_p_att)

* ---------- compare the saved reference distributions ----------
foreach f in p_base p_fest p_fecoh p_cctrl p_sctrl {
  use "`f'.dta", clear
  qui sum main_att if _rep > 0
  dis "DIST `f': n=" r(N) " mean=" %12.8f r(mean) " sd=" %12.8f r(sd) " min=" %12.8f r(min) " max=" %12.8f r(max)
  dis "   obs row = " %21x main_att[1]
}
use "p_base.dta", clear
rename main_att b_base
merge 1:1 _rep using "p_fecoh.dta", nogen
rename main_att b_fecoh
merge 1:1 _rep using "p_fest.dta", nogen
rename main_att b_fest
merge 1:1 _rep using "p_cctrl.dta", nogen
rename main_att b_cctrl
gen double d_fecoh = b_fecoh - b_base
gen double d_fest  = b_fest  - b_base
gen double d_cctrl = b_cctrl - b_base
list _rep b_base b_fecoh d_fecoh d_fest d_cctrl, sep(0)
