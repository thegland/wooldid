clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 313131
set sortseed 313131
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid

* report: subgroup vs pooled SEs, pairwise contrast SEs, and the implied rho
cap prog drop sgrep
prog def sgrep
  syntax , tag(string)
  tempname M A B P
  mat `M' = e(wdidmainresults)
  cap mat `A' = e(wdidsgresults_0)
  if _rc != 0 exit
  cap mat `B' = e(wdidsgresults_1)
  if _rc != 0 exit
  local rfn : rowfullnames `M'
  local nr = rowsof(`M')
  forval r = 1/`nr' {
    local nm : word `r' of `rfn'
    dis "`tag'_LEV `nm'  main_se=" %21x `M'[`r',2] "  sg0_se=" %21x `A'[`r',2] "  sg1_se=" %21x `B'[`r',2] ///
        "  r0=" %18.15f `A'[`r',2]/`M'[`r',2] "  r1=" %18.15f `B'[`r',2]/`M'[`r',2]
  }
  cap mat `P' = e(pwsgs)
  if _rc != 0 exit
  local prn : rowfullnames `P'
  local pr = rowsof(`P')
  local afn : rowfullnames `A'
  forval r = 1/`pr' {
    local rn : word `r' of `prn'
    * eq of the pairwise row -> matching level row
    gettoken eqn rest : rn, parse(":")
    local ispre = strpos("`eqn'","pre") > 0
    local item = cond(`ispre', subinstr(subinstr("`eqn'","pre","",1),"_","",1), "`eqn'")
    local eqwant = cond(`ispre', "pretreatment", "main")
    local sea .
    local seb .
    local k 0
    foreach wn of local afn {
      local ++k
      if "`wn'" == "`eqwant':`item'" {
        local sea = `A'[`k',2]
        local seb = `B'[`k',2]
      }
    }
    local sed = `P'[`r',2]
    local rho = .
    if `sea' < . & `seb' < . & `sed' < . local rho = (`sea'^2 + `seb'^2 - `sed'^2)/(2*`sea'*`seb')
    dis "`tag'_PW `rn'  b=" %14.10f `P'[`r',1] "  se=" %14.9g `sed' "  p=" %10.6f `P'[`r',4] ///
        "  sea=" %14.9g `sea' "  seb=" %14.9g `seb' "  rho=" %12.9f `rho' "  flagged=" (`rho' > 1-1e-2)
  }
  dis "`tag'_sgdegen=" e(wdid_sgdegen) "  act=`e(wdid_sgdegenact)'  maxrho=`e(wdid_sgdegenr)'"
  cap dis "`tag'_allsgssame_att=" %10.6f e(joint_att_allsgssame)
  cap dis "`tag'_pre_allsgssame_att=" %10.6f e(joint_pre_att_allsgssame)
end

* ================= DGP 1: PERFECTLY BALANCED panel =================
capture drop _all
clear
set obs 24
gen int state = _n
gen int onset = .
qui replace onset = 2002 if state <= 6
qui replace onset = 2003 if state >= 7  & state <= 12
qui replace onset = 2004 if state >= 13 & state <= 18
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 100*state
gen double u_i = rnormal()*0.5
expand 8
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t-onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
gen byte sg = mod(state,2)
save bal24, replace
dis "### P1: PERFECTLY BALANCED, balanced sg, cluster(state)==i"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "P1_RC=" _rc
sgrep , tag(P1)

dis "### P2: PERFECTLY BALANCED panel, UNBALANCED sg (5 of 24 states), cluster(state)"
use bal24, clear
qui replace sg = inlist(state,1,2,7,13,19)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "P2_RC=" _rc
sgrep , tag(P2)

dis "### P3: balanced panel except units per state VARY (2..5); balanced sg"
use bal24, clear
qui drop if unit - 100*state > 2 + mod(state,2)
qui count
dis "P3_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "P3_RC=" _rc
sgrep , tag(P3)

* ================= DGP 2: RAGGED panel =================
use bal24, clear
set seed 989898
gen double _rr = runiform()
qui drop if _rr < 0.08
drop _rr
qui replace sg = inlist(state,1,2,3,7,8,13,19,20)
save rag24, replace
qui count
dis "P4_NOBS=" r(N)
dis "### P4: RAGGED panel, unbalanced sg, cluster(state)"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "P4_RC=" _rc
sgrep , tag(P4)

dis "### P5: SAME ragged data, robust (oracle-ish control)"
use rag24, clear
cap noi wooldid y state t ttre , robust att esprelength(2) espostlength(2) subgroups(sg)
dis "P5_RC=" _rc
sgrep , tag(P5)

dis "### P6: SAME ragged data, cluster(state) + bootstrap(standard) reps(399) (the remedy the caution names)"
use rag24, clear
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg) ///
   bootstrap(standard) reps(399) dots(0)
dis "P6_RC=" _rc
sgrep , tag(P6)
dis "R06_DONE"
