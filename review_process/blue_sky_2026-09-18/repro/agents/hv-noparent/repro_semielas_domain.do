*=====================================================================
* repro_semielas_domain.do -- SELF-CONTAINED.
* Claim: the semielasticity / contreatelasticitytype(eyex|eydx) family
* has no screen on the DOMAIN of the reporting scale.  d ln E[y]/dx is
* undefined where the fitted outcome is not positive.  wooldid does
* three different things there and none of them names the problem:
*   (i)  returns a number, with a standard error, t, p and interval,
*        and no caution at all
*   (ii) dies with Stata's generic margins message, r(459)
*   (iii) dies inside Mata with r(3598) naming an internal tempvar,
*        "variable __00001O not found" -- the weighted leg
* Which of the three you get depends on the data and on whether a
* weight is present.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

cap prog drop T
prog def T
  syntax , nm(string) opts(string asis) yv(string) [wt(string asis)]
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  if `rc' == 0 {
    matrix Z = e(wdidmainresults)
    dis "===VERDICT `nm': rc=0  att=" string(Z[1,1],"%12.6f") "  se=" string(Z[1,2],"%12.6f") "  p=" string(Z[1,4],"%12.6f")
  }
  else dis "===VERDICT `nm': rc=`rc'"
end

*---- panel A: hand-built, same model, three outcomes that differ only
*     by an additive constant.  A constant shift cannot change the LEVEL
*     ATT at all, so any change in the reported semi-elasticity is pure
*     reporting-scale behaviour.
clear
set obs 40
gen long unit = _n
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2003 if unit <= 12
qui replace ttre = 2005 if unit > 12 & unit <= 24
set seed 987654
gen double aw = 0.5 + runiform()
gen double base = 0.1*unit + 0.05*(t-1999) + 0.8*(ttre < . & t >= ttre) + rnormal()*0.2
gen double ypos  = base + 5      // fitted values strictly positive
gen double yneg  = base - 2      // fitted values cross zero
qui sum ypos
dis "min(ypos)=" r(min)
qui sum yneg
dis "min(yneg)=" r(min)
qui save "rpse.dta", replace

dis "=== level ATT is identical for the two outcomes (shift invariance) ==="
use "rpse.dta", clear
T , nm(level_ypos)  opts(att robust) yv(ypos)
use "rpse.dta", clear
T , nm(level_yneg)  opts(att robust) yv(yneg)

dis "=== semielasticity, strictly positive outcome: fine ==="
use "rpse.dta", clear
T , nm(semi_ypos_unweighted) opts(att robust semielasticity) yv(ypos)
use "rpse.dta", clear
T , nm(semi_ypos_aweighted)  opts(att robust semielasticity) yv(ypos) wt([aw=aw])

dis "=== semielasticity, outcome crossing zero: a NUMBER, with no caution ==="
use "rpse.dta", clear
T , nm(semi_yneg_unweighted) opts(att robust semielasticity) yv(yneg)
use "rpse.dta", clear
T , nm(semi_yneg_aweighted)  opts(att robust semielasticity) yv(yneg) wt([aw=aw])
use "rpse.dta", clear
T , nm(semi_yneg_pweighted)  opts(att robust semielasticity) yv(yneg) wt([pw=aw])

*---- panel B: a different panel of the same kind, on which the SAME
*     option dies instead -- unweighted with Stata's generic margins
*     message r(459), weighted with an opaque Mata r(3598).
clear
set obs 80
gen long unit = _n
expand 10
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2004 if unit <= 16
qui replace ttre = 2006 if unit > 16 & unit <= 32
qui replace ttre = 2008 if unit > 32 & unit <= 48
set seed 20260916
gen double aunit = .
bysort unit: replace aunit = rnormal() if _n == 1
bysort unit: replace aunit = aunit[1]
gen double wu = .
bysort unit: replace wu = exp(rnormal()*1.15) if _n == 1
bysort unit: replace wu = wu[1]
gen double gt = .
bysort t: replace gt = rnormal()*0.5 if _n == 1
bysort t: replace gt = gt[1]
gen double yx = aunit + gt + 0.8*(ttre < . & t >= ttre) + rnormal()*0.3
qui sum yx
dis "min(yx)=" r(min) "  max=" r(max)
qui save "rpse2.dta", replace

dis "=== the SAME option on a panel of the same kind ==="
use "rpse2.dta", clear
T , nm(semi_yx_unweighted) opts(att robust semielasticity) yv(yx)
use "rpse2.dta", clear
T , nm(semi_yx_aweighted)  opts(att robust semielasticity) yv(yx) wt([aw=wu])
dis "    (the aweighted leg is the r(3598) 'variable __00001O not found' traceback"
dis "     printed just above this line)"
use "rpse2.dta", clear
T , nm(semi_yx_control_nosemi) opts(att robust) yv(yx) wt([aw=wu])
*---- panel C: the opaque Mata leg.  Built with the shared blue-sky
*     generator so any reviewer on this machine can rerun it verbatim.
*     UNWEIGHTED it dies r(459) with Stata's generic margins message;
*     with an aweight it dies r(3598) inside Mata, "variable __00001O
*     not found" -- an internal tempvar name, from st_data() under
*     __marg_cr_eval().
clear all
set seed 20260916
set sortseed 20260916
local gN      80
local gT      10
local gt0     2000
local gonsets 2004 2006 2008
local gysd    0.30
local gdose   lattice
local gwtype  unit
local gwspread 2
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
qui sum ync
dis "min(ync)=" r(min) "  max=" r(max)
qui save "rpse3.dta", replace
cap prog drop T
prog def T
  syntax , nm(string) opts(string asis) yv(string) [wt(string asis)]
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  if `rc' == 0 {
    matrix Z = e(wdidmainresults)
    dis "===VERDICT `nm': rc=0  att=" string(Z[1,1],"%12.6f")
  }
  else dis "===VERDICT `nm': rc=`rc'"
end
use "rpse3.dta", clear
T , nm(semi_gen_unweighted) opts(att robust semielasticity) yv(ync)
use "rpse3.dta", clear
T , nm(semi_gen_AWEIGHTED_expect3598) opts(att robust semielasticity) yv(ync) wt([aw=w])
use "rpse3.dta", clear
T , nm(semi_gen_control_nosemi) opts(att robust) yv(ync) wt([aw=w])
dis "REPRO_SEMIELAS_COMPLETE"
