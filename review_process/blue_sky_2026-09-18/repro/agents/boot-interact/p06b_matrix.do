* p06b — INTERACTION MATRIX: bootstrap() x the option surface, cells 27-52.
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
local c1  "bootdf(normal)"
local c2  "bootdf(7)"
local c3  "bootalpha(equal)"
local c4  "level(90)"
local c5  "ties"
local c6  "histogramcohorteffects"
local c7  "summarystats"
local c8  "nofastcollapse"
local c9  "permonsetcohorts"
local c10 "emptycellsoverride"
local c11 "cleanmatrices"
local c12 "pretrendpower(80)"
local c13 "safety(off)"
local c14 "dots(5)"
local c15 "contreat(dose) contreatcontrols(x1)"
local c16 "oostreatedcontrols"
local c17 "coarsecohortcontrols(reg)"
local c18 "ccc_absorb(reg)"
local c19 "fe(reg)"
local c20 "regtol(1e-6)"
local c21 "subgroups(grp) jointtests adversarial"
local c22 "subgroups(grp) makecfxplotsbysg"
local c23 "esprelength(0)"
local c24 "makeplots"
local c25 "boottstat jointtests adversarial"
local c26 "boottstat simultaneousbands pretrendpower(80)"
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
dis "P06B_DONE"
