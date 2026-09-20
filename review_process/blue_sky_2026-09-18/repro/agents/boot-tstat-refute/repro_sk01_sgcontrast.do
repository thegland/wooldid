*=====================================================================
* SK-01  SELF-CONTAINED REPRO
* Claim: on an ordinary ragged staggered panel with cluster() at the level of
*   i and subgroups(), wooldid prints a subgroup-contrast standard error that is
*   an order of magnitude too small, a t of 12 and a p of 0.000, and an
*   across-subgroup equality p of 0.0000, with NO caution -- because the ruled
*   rank-deficiency guard wooldid_sgdeg (ado:11092, cut 1 - 1e-2) sees an implied
*   correlation of 0.9786 and lets the row through.
* Oracle: the remedy wooldid's OWN caution names for this design --
*   bootstrap(standard) WITHOUT boottstat, same cluster() -- plus a robust control.
* Run:  cd <this dir> ; stata-mp -e do repro_sk01_sgcontrast.do
* Look for:  SK01_FAIL=1
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 123757
set sortseed 123757

* ---- 20 states, 2-5 units each, 10 years, three staggered cohorts, 30% never
*      treated, 8 percent of rows missing at random (ragged administrative panel)
clear
set obs 20
gen int state = _n
gen int onset = .
qui replace onset = 2002 if state <= 5
qui replace onset = 2003 if state >  5 & state <= 10
qui replace onset = 2004 if state > 10 & state <= 14
gen double a_i = rnormal()
gen int nun = 2 + mod(state,4)
expand nun
bys state: gen int unit = _n + 100*state
gen double u_i = rnormal()*0.6
expand 10
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t-onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
gen double _rr = runiform()
qui drop if _rr < 0.08
drop _rr
* subgroup: a fixed, cohort-uneven split of the states
gen byte sg = mod(state*7 + 3, 5) <= 1
qui count
dis "SK01_NOBS=" r(N)
qui tab sg

cap prog drop pw1
prog def pw1
  syntax , tag(string)
  tempname P A B
  mat `P' = e(pwsgs)
  mat `A' = e(wdidsgresults_0)
  mat `B' = e(wdidsgresults_1)
  local afn : rowfullnames `A'
  local k 0
  local sea .
  local seb .
  foreach wn of local afn {
    local ++k
    if "`wn'" == "main:att" {
      local sea = `A'[`k',2]
      local seb = `B'[`k',2]
    }
  }
  local sed = `P'[1,2]
  local rho = .
  if `sea'<. & `seb'<. & `sed'<. local rho = (`sea'^2+`seb'^2-`sed'^2)/(2*`sea'*`seb')
  c_local `tag'se = `sed'
  c_local `tag'p  = `P'[1,4]
  dis "`tag': att contrast b=" %12.8f `P'[1,1] "  se=" %14.9g `sed' "  t=" %8.3f `P'[1,3] ///
      "  p=" %14.9f `P'[1,4] "  rho=" %12.9f `rho'
  dis "`tag': sgdegen_flag=" cond("`e(wdid_sgdegenact)'"=="", "NONE", "`e(wdid_sgdegenact)'") ///
      "  n_flagged=" e(wdid_sgdegen) "  allsgssame_att=" %12.9f e(joint_att_allsgssame)
end

dis "### 1. ANALYTIC, cluster(state) == i   (what wooldid prints)"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "ANA_RC=" _rc
pw1 , tag(ANA)

dis "### 2. ORACLE A: bootstrap(standard) reps(399), SAME cluster(state), no boottstat"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg) ///
   bootstrap(standard) reps(399) dots(0)
dis "BOOT_RC=" _rc
pw1 , tag(BOOT)

dis "### 3. ORACLE B: robust"
cap noi wooldid y state t ttre , robust att esprelength(2) espostlength(2) subgroups(sg)
dis "ROB_RC=" _rc
pw1 , tag(ROB)

local ratio = `BOOTse'/`ANAse'
dis "SK01_SE_RATIO_boot_over_analytic=" %10.4f `ratio'
dis "SK01_ANALYTIC_P=" %14.9f `ANAp' "   SK01_BOOTSTRAP_P=" %14.9f `BOOTp' "   SK01_ROBUST_P=" %14.9f `ROBp'
local fail = (`ratio' > 4) & (`ANAp' < 0.01) & (`BOOTp' > 0.10)
dis "SK01_FAIL=" `fail'
dis "SK01_DONE"
