* probe07: all four estimands + event study under the rebuild; column-name claim;
*          attrition design where a drawn onset leaves a cluster with no post period
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* ---------- A. four estimands + event study on the ragged design B ----------
use "dataB.dta", clear
cap noi wooldid y coh t onset, att att_i att_it att_itime esprelength(2) espostlength(2) cluster(state) permute(shuffle) permexact permonsetcohorts
dis "A_RC=" _rc
dis "A done=" e(perm_completed) " fail=" e(perm_failed) " M=" e(perm_npossible)
cap mat B = e(b)
cap dis "A colnames: " "`: colfullnames B'"
cap mat list e(perm_n)
cap mat list e(perm_p)

* ---------- B. the "relabelling i leaves colfullnames identical" claim ----------
use "dataB.dta", clear
cap qui wooldid y coh t onset, att att_i esprelength(2) espostlength(2) cluster(state) fastnose
mat B1 = e(b)
local n1 : colfullnames B1
use "dataB.dta", clear
replace coh = cond(coh==1,3,cond(coh==3,1,2))
cap qui wooldid y coh t onset, att att_i esprelength(2) espostlength(2) cluster(state) fastnose
mat B2 = e(b)
local n2 : colfullnames B2
dis "COLNAMES_MATCH=" cond("`n1'" == "`n2'", "YES", "NO")
dis "n1: `n1'"
dis "n2: `n2'"

* ---------- C. attrition design: a cluster that draws a late onset has no post rows ----
clear
set seed 31415
set sortseed 31415
set obs 6
gen int state = _n
gen double sfx = rnormal()*1.5
gen int thi = cond(inlist(state,5,6),7,8)
expand 8
bysort state: gen int t = _n
drop if t > thi
gen double eps = rnormal()
gen int coh = cond(inlist(state,1,2),1,cond(state==3,2,3))
gen double onset = cond(coh==1,4,cond(coh==2,8,.))
gen byte D = (onset<. & t>=onset)
gen double y = 1 + sfx + 0.3*t + 1.0*D + eps
keep state t coh onset y
save "dataC.dta", replace
tab state coh
use "dataC.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(C_att) replace
dis "C_RC=" _rc
dis "C done=" e(perm_completed) " fail=" e(perm_failed) " M=" e(perm_npossible) " p=" e(perm_p_att)
cap mat list e(perm_n)
use "C_att.dta", clear
count if missing(main_att)
dis "C_NMISSING_ROWS=" r(N)
list, sep(0)
