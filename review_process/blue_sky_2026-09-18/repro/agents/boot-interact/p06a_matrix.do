* p06a — INTERACTION MATRIX: bootstrap() x the option surface, cells 1-26.
* For each cell: rc, e(wdid_vcvstate), boot_completed, and whether the reported att SE
* actually MOVED off the analytic one (se_boot vs se_analytic from e(V) / e(V_analytic)).
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 40
gen int state = _n
gen int onset = .
replace onset = 5 if state <= 10
replace onset = 7 if state > 10 & state <= 20
replace onset = 9 if state > 20 & state <= 30
gen double sfe = rnormal()
gen byte grp  = mod(state,3)
gen byte reg  = mod(state,4)
expand 12
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double x1 = rnormal()
gen double wt = 1 + runiform()
gen int fwt = 1 + int(3*runiform())
gen double dose = cond(D, 0.5 + runiform(), 0)
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + sfe + 0.05*t + te + 0.2*x1 + rnormal()*0.8
gen double ypos = exp(y/4)
save "p06_data.dta", replace

local BASE cluster(state) esfixedbaseperiod esprelength(3) espostlength(3)
local BOOT bootstrap(standard) reps(20) rseed(7) dots(0)

local NC 26
local c1  ""
local c2  "att_i"
local c3  "att_it"
local c4  "att_itime"
local c5  "att att_i att_it att_itime"
local c6  "subgroups(grp)"
local c7  "jointtests"
local c8  "adversarial"
local c9  "simultaneousbands"
local c10 "gteffects"
local c11 "impvalues"
local c12 "fastnose"
local c13 "unconditionalse"
local c14 "customvce(cluster state)"
local c15 "permute(onset)"
local c16 "subgroups(grp) suppresssgprimaryeffects"
local c17 "oldsyntax"
local c18 "poisson"
local c19 "timetrends"
local c20 "controls(x1)"
local c21 "interactivecontrols(x1)"
local c22 "anticipation(1)"
local c23 "esrelativeto(-2)"
local c24 "boottstat"
local c25 "bootdisplay(bc)"
local c26 "bootdisplay(bca)"

forval k = 1/`NC' {
  dis "##### CELL `k' : [`c`k''] #####"
  cap noi wooldid y state t onset, `BASE' `BOOT' `c`k''
  local rc = _rc
  dis "CELL`k'_RC=" `rc'
  if `rc' == 0 {
    dis "CELL`k'_VCV=`e(wdid_vcvstate)' METHOD=`e(boot_method)' COMPLETED=" e(boot_completed) " FAILED=" e(boot_failed)
    dis "CELL`k'_VPOSTED=" e(boot_v_posted) " NV=" e(boot_N_V) " DFR=" e(df_r) " DFRAN=" e(df_r_analytic)
    dis "CELL`k'_VESFAILED=[`e(boot_ves_failed)'] BCABSENT=[`e(boot_bc_absent)'] JTABSENT=[`e(boot_jt_absent)']"
    cap confirm matrix e(V_analytic)
    local hva = (_rc==0)
    dis "CELL`k'_HASVANALYTIC=" `hva'
    if `hva' {
      matrix VB = e(V)
      matrix VA = e(V_analytic)
      dis "CELL`k'_SEBOOT=" %12.6g sqrt(VB[1,1]) " SEANA=" %12.6g sqrt(VA[1,1]) " RATIO=" %8.4f sqrt(VB[1,1]/VA[1,1])
    }
  }
  dis "##### END CELL `k' #####"
}
dis "P06A_DONE"
