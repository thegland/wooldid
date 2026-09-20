clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 515151
set sortseed 515151
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid

* ---- my own panel: 22 states, 3 units, 9 years, 3 cohorts, 6 never treated
clear
set obs 22
gen int state = _n
gen int onset = .
qui replace onset = 2002 if state <= 5
qui replace onset = 2003 if state >= 6  & state <= 11
qui replace onset = 2004 if state >= 12 & state <= 16
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 100*state
gen double u_i = rnormal()*0.5
expand 9
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t-onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
save bs22, replace

*==================================================================
dis "### J: boottstat + jointtests + bootsave"
use bs22, clear
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(199) boottstat jointtests bootsave(jw1) replace dots(0)
dis "J_RC=" _rc
dis "J_joint_horiz=" %14.10f e(joint_horiz)
dis "J_reps_done=" e(boot_reps_done)
mat BES = e(b_es_att)
local esn : colnames BES
dis "J_b_es_att colnames: `esn'"
forval k = 1/`=colsof(BES)' {
  local nm : word `k' of `esn'
  dis "J_BES `nm' = " %14.10f BES[1,`k']
}
mat JW = e(boot_jointW)
local jwn : colnames JW
dis "J_jointW colnames: `jwn'"
dis "J_jointW rows=" rowsof(JW) " cols=" colsof(JW)

preserve
use jw1, clear
qui ds jointW_*
dis "J_DUMPCOLS: `r(varlist)'"
foreach v of varlist jointW_* {
  qui sum `v' if _rep==0
  local obs0 = r(mean)
  qui sum `v' if _rep>0
  dis "J_FILE `v'  rep0=" %16.10f `obs0' "  repmin=" %16.10f r(min) "  repmax=" %16.6f r(max) "  n=" r(N)
  qui count if _rep>0 & `v' >= `obs0' & !missing(`v')
  local ge = r(N)
  qui count if _rep>0 & !missing(`v')
  local nn = r(N)
  dis "J_FILE `v'  p_from_file_addone=" %10.7f min(1,(`ge'+1)/(`nn'+1))
}
restore
dis "J_NEGATIVE_W_IN_FILE_CHECK done"

*==================================================================
dis "### K: same but with adversarial (advW leg)"
use bs22, clear
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(199) boottstat jointtests adversarial bootsave(jw2) replace dots(0)
dis "K_RC=" _rc
dis "K_adv_p=" %14.10f e(adversarial_p)
qui ereturn list
mat AW = e(boot_advW)
dis "K_advW rows=" rowsof(AW) " cols=" colsof(AW)
mat BES = e(b_es_att)
local esn : colnames BES
forval k = 1/`=colsof(BES)' {
  local nm : word `k' of `esn'
  dis "K_BES `nm' = " %14.10f BES[1,`k']
}
preserve
use jw2, clear
qui ds advW*
dis "K_DUMPCOLS: `r(varlist)'"
foreach v of varlist advW* {
  qui sum `v' if _rep==0
  local obs0 = r(mean)
  qui sum `v' if _rep>0
  dis "K_FILE `v'  rep0=" %16.10f `obs0' "  repmin=" %16.10f r(min) "  repmax=" %16.6f r(max) "  n=" r(N)
  qui count if _rep>0 & `v' >= `obs0' & !missing(`v')
  local ge = r(N)
  qui count if _rep>0 & !missing(`v')
  local nn = r(N)
  dis "K_FILE `v'  p_from_file_addone=" %10.7f min(1,(`ge'+1)/(`nn'+1))
}
restore
dis "K_e_scalars:"
foreach s in adversarial_p adv_p joint_horiz {
  cap dis "`s' = " %14.10f e(`s')
}
dis "R05_DONE"
