*=====================================================================
* BT-04 VERDICT REPRO (skeptic's own; adds the advW leg the original only read)
* Claim: in a bootsave() file the _rep == 0 row of jointW_<family> -- and of advW --
*   holds an observed event-study COEFFICIENT, not the observed joint/adversarial
*   statistic, because wooldid_bootsavemk (ado:14915-14931) indexes the observed
*   matrix by the ACCUMULATOR's column position.
* Look for: BT04_JOINT_FAIL=1 and BT04_ADV_FAIL=1
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 515151
set sortseed 515151
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
tempfile base
save `base', replace

dis "### JOINT leg"
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(199) boottstat jointtests bootsave(bt04a) replace dots(0)
dis "J_RC=" _rc
local jp = e(joint_horiz)
mat BES = e(b_es_att)
local pre3 = BES[1,1]
dis "J_reported_joint_p=" %14.10f `jp'
dis "J_b_es_att_pre3   =" %14.10f `pre3'
preserve
use bt04a, clear
qui sum jointW_att if _rep==0
local o0 = r(mean)
qui count if _rep>0 & jointW_att < .
local nn = r(N)
qui count if _rep>0 & jointW_att >= `o0' & jointW_att < .
local ge = r(N)
qui sum jointW_att if _rep>0
dis "J_file_rep0=" %16.10f `o0' "  rep_min=" %16.10f r(min) "  rep_max=" %16.4f r(max) "  n=" `nn'
dis "J_p_recomputed_from_file=" %10.7f min(1,(`ge'+1)/(`nn'+1))
local jfail = (abs(`o0' - `pre3') < 1e-9) & (`o0' < 0)
restore
dis "BT04_JOINT_FAIL=" `jfail'

dis "### ADVERSARIAL leg"
use `base', clear
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(199) boottstat jointtests adversarial bootsave(bt04b) replace dots(0)
dis "K_RC=" _rc
mat BES = e(b_es_att)
local pre3 = BES[1,1]
preserve
use bt04b, clear
qui sum advW if _rep==0
local o0 = r(mean)
qui count if _rep>0 & advW < .
local nn = r(N)
qui count if _rep>0 & advW >= `o0' & advW < .
local ge = r(N)
qui sum advW if _rep>0
dis "K_file_rep0=" %16.10f `o0' "  rep_min=" %16.10f r(min) "  rep_max=" %16.4f r(max) "  n=" `nn'
dis "K_p_recomputed_from_file=" %10.7f min(1,(`ge'+1)/(`nn'+1))
local afail = (abs(`o0' - `pre3') < 1e-9) & (`o0' < 0)
restore
dis "K_b_es_att_pre3=" %14.10f `pre3'
dis "BT04_ADV_FAIL=" `afail'
dis "BT04_DONE"
