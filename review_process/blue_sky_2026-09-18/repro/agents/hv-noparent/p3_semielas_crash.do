*=====================================================================
* p3_semielas_crash.do -- isolate r(3598) "variable __00001O not found"
* seen on:  wooldid ync unit t ttre [aw=w], att att_i robust semielasticity
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set maxvar 20000
which wooldid

cap prog drop tryit
prog def tryit
  syntax , nm(string) opts(string asis) [yv(string) wt(string asis)]
  if "`yv'" == "" local yv ync
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  dis "SENT `nm' RC=" _rc
end

local gN      80
local gT      10
local gt0     2000
local gonsets 2004 2006 2008
local gysd    0.30
local gdose   lattice
local gwtype  unit
local gwspread 2
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
qui save "pdat3.dta", replace
qui sum ync
dis "YNC min=" r(min) " max=" r(max)
qui gen double ypos = ync - r(min) + 1
qui save "pdat3.dta", replace

* ---- 1: exact repro from p2
use "pdat3.dta", clear
tryit , nm(s1_attatti_aw) opts(att att_i robust semielasticity) wt([aw=w])
* ---- 2: no att_i
use "pdat3.dta", clear
tryit , nm(s2_att_aw) opts(att robust semielasticity) wt([aw=w])
* ---- 3: no weights
use "pdat3.dta", clear
tryit , nm(s3_att_noweight) opts(att robust semielasticity)
* ---- 4: no weights, att_i
use "pdat3.dta", clear
tryit , nm(s4_attatti_noweight) opts(att att_i robust semielasticity)
* ---- 5: nofastcollapse
use "pdat3.dta", clear
tryit , nm(s5_att_nofastcollapse) opts(att robust semielasticity nofastcollapse)
* ---- 6: strictly positive outcome (elasticity domain is sane)
use "pdat3.dta", clear
tryit , nm(s6_ypos) opts(att robust semielasticity) yv(ypos)
* ---- 7: control, no semielasticity
use "pdat3.dta", clear
tryit , nm(s7_control_nosemi) opts(att robust)
* ---- 8: cluster instead of robust
use "pdat3.dta", clear
tryit , nm(s8_cluster) opts(att robust semielasticity cluster(region))
* ---- 9: with event study
use "pdat3.dta", clear
tryit , nm(s9_es) opts(att robust semielasticity esprelength(2) espostlength(2))
* ---- 10: default vce (iid)
use "pdat3.dta", clear
tryit , nm(s10_iid) opts(att semielasticity)
* ---- 11: subgroups
use "pdat3.dta", clear
tryit , nm(s11_sg) opts(att robust semielasticity subgroups(sg2))
* ---- 12: tiny clean panel built by hand, no generator at all
clear
set obs 40
gen long unit = _n
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2003 if unit <= 10
qui replace ttre = 2005 if unit > 10 & unit <= 20
set seed 424242
gen double ync = 3 + 0.1*unit + 0.05*t + 0.6*(ttre < . & t >= ttre) + rnormal()*0.2
tryit , nm(s12_handbuilt) opts(att robust semielasticity)
* ---- 13: same hand-built, control
tryit , nm(s13_handbuilt_control) opts(att robust)
* ---- 14: same hand-built, nofastcollapse
tryit , nm(s14_handbuilt_nofc) opts(att robust semielasticity nofastcollapse)
* ---- 15: same hand-built + poisson (the documented remap path)
qui gen int ycnt = max(0, round(exp(ync/3)))
tryit , nm(s15_handbuilt_pois) opts(att robust semielasticity poisson) yv(ycnt)
* ---- 16: same hand-built, semielasticity + poisexpresults
tryit , nm(s16_handbuilt_poisexp) opts(att robust semielasticity poisson poisexpresults) yv(ycnt)

dis "PROBE3_COMPLETE"
