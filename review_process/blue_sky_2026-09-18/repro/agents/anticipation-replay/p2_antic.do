*==============================================================================
* p2_antic.do -- EXACT EQUIVALENCE test of anticipation(k) against the
*                documented hand-edit (sthlp:821: "exactly equivalent to
*                re-running wooldid after replacing ttre with ttre - #").
*
* For each config: run A = anticipation(k) on the raw ttre
*                  run B = same options, on data where ttre has been
*                          overwritten IN PLACE with ttre - k
* Both runs use the SAME variable name, so e() macros are comparable.
* Every visible e() scalar/macro/matrix is dumped in %21x and diffed by diffe.py.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

* ------------------------------------------------------- datasets -----------
* D2: workhorse annual panel, 2000-2013, onsets 2005/2008/2010, 30% never treated
local dN 100
local dT 14
local donsets "2005 2008 2010"
local dantic 2
include "`W'/mkdata.do"
gen double ypos = exp(y/3)
save "`W'/d2.dta", replace

* D2B: biennial panel (t spacing 2) -- an ODD k lands off-grid
clear all
local dN 120
local dT 10
local dstep 2
local dt0 2000
local donsets "2006 2010 2012"
local dnevsh .30
local dysd .30
local dantic 2
local dseed 771
include "`W'/mkdata.do"
save "`W'/d2b.dta", replace

* D2C: NO never-treated units at all; latest cohort is beyond the window.
*      Only usable with oostreatedcontrols(); anticipation() moves the boundary.
clear all
local dN 120
local dT 10
local dstep 1
local dt0 2000
local donsets "2003 2005 2011"
local dnevsh 0
local dysd .30
local dantic 2
local dseed 772
include "`W'/mkdata.do"
save "`W'/d2c.dta", replace

* D2D: unbalanced -- late entry / early exit / interior gaps (administrative shape)
clear all
local dN 120
local dT 14
local dstep 1
local dt0 2000
local donsets "2004 2007 2009"
local dnevsh .30
local dysd .30
local dantic 2
local dseed 773
include "`W'/mkdata.do"
set seed 9931
gen double uu = runiform()
gen double entry = 2000 + floor(3*runiform())
bysort i (t): replace entry = entry[1]
drop if t < entry
drop if uu < .04
bysort i (t): gen byte nobs = _N
drop if nobs < 6
drop uu entry nobs
save "`W'/d2d.dta", replace

include "`W'/dumpe.do"

cap prog drop runpair
prog def runpair
  syntax , tag(string) k(integer) opts(string) w(string) [dta(string)]
  if "`dta'" == "" local dta d2
  dis as result "==== CONFIG `tag'  data=`dta'  k=`k'  opts: `opts' ===="
  qui use "`w'/`dta'.dta", clear
  cap noi wooldid `opts' anticipation(`k')
  local rcA = _rc
  dis "A_RC_`tag'=`rcA'"
  if `rcA' == 0 dumpe , f("`w'/out_`tag'_A.txt")
  qui use "`w'/`dta'.dta", clear
  qui replace ttre = ttre - `k'
  cap noi wooldid `opts'
  local rcB = _rc
  dis "B_RC_`tag'=`rcB'"
  if `rcB' == 0 dumpe , f("`w'/out_`tag'_B.txt")
  dis "PAIR_`tag'_RCS=`rcA',`rcB'"
end

* ------------------------------------------------------ core sweep ----------
runpair , tag(c01) k(2) w("`W'") opts(y i t ttre, att robust)
runpair , tag(c02) k(2) w("`W'") opts(y i t ttre, att att_it att_i att_itime cluster(i) espre(3) espost(4) jointtests)
runpair , tag(c03) k(2) w("`W'") opts(y i t ttre, att espre(3) espost(4) esfixedbaseperiod esrelativeto(-2) robust)
runpair , tag(c04) k(2) w("`W'") opts(y i t ttre, att robust oostreatedcontrols(3) espre(2) espost(2))
runpair , tag(c05) k(2) w("`W'") opts(y i t ttre, att subgroups(sg) robust espre(2) espost(2))
runpair , tag(c06) k(2) w("`W'") opts(y i t ttre, att robust contreat(dose) espre(2) espost(2))
runpair , tag(c07) k(2) w("`W'") opts(ypos i t ttre, att poisson robust espre(2) espost(2))
runpair , tag(c08) k(2) w("`W'") opts(y i t ttre, att timetrends robust espre(2) espost(2))
runpair , tag(c09) k(2) w("`W'") opts(y i t ttre, att robust gteffects summarystats espre(2) espost(2))
runpair , tag(c10) k(2) w("`W'") opts(y i t ttre, att robust adversarial espre(3) espost(3))
runpair , tag(c11) k(2) w("`W'") opts(y i t ttre, att robust controls(x1 x2) fe(i.sg) espre(2) espost(2))
runpair , tag(c12) k(1) w("`W'") opts(y i t ttre [aw=w], att robust espre(2) espost(2))
runpair , tag(c13) k(5) w("`W'") opts(y i t ttre, att robust espre(2) espost(2))
runpair , tag(c14) k(2) w("`W'") opts(y i t ttre, att robust espre(3) espost(3) pretrendpower(80))
runpair , tag(c15) k(2) w("`W'") opts(y i t ttre if sg==1, att robust espre(2) espost(2))
runpair , tag(c16) k(2) w("`W'") opts(y i t ttre, att_it robust espre(2) espost(2) impvalues cluster(i))

* -------------------------------------------- replication engines -----------
runpair , tag(c17) k(2) w("`W'") opts(y i t ttre, att cluster(i) espre(2) espost(2) permute(shuffle) permreps(25) rseed(4242) dots(0))
runpair , tag(c18) k(2) w("`W'") opts(y i t ttre, att cluster(i) espre(2) espost(2) bootstrap(standard) reps(25) rseed(4242) dots(0))
runpair , tag(c19) k(2) w("`W'") opts(y i t ttre, att_i cluster(i) espre(2) espost(2) bootstrap(standard) reps(25) rseed(4242) dots(0))

* ------------------------------------------------- edge geometries ----------
* biennial + ODD k  (off-grid effective onsets)
runpair , tag(c20) k(1) w("`W'") dta(d2b) opts(y i t ttre, att robust espre(2) espost(2))
* biennial + EVEN k (on grid) -- control for c20
runpair , tag(c21) k(2) w("`W'") dta(d2b) opts(y i t ttre, att robust espre(2) espost(2))
* no never-treated: control group exists only via oostreatedcontrols(); the shift moves it
runpair , tag(c22) k(2) w("`W'") dta(d2c) opts(y i t ttre, att robust espre(2) espost(2) oostreatedcontrols(0))
runpair , tag(c23) k(3) w("`W'") dta(d2c) opts(y i t ttre, att robust espre(2) espost(2) oostreatedcontrols(4))
* unbalanced administrative shape
runpair , tag(c24) k(2) w("`W'") dta(d2d) opts(y i t ttre, att att_it cluster(i) espre(3) espost(3) jointtests)
runpair , tag(c25) k(2) w("`W'") dta(d2d) opts(y i t ttre, att espre(3) espost(3) esfixedbaseperiod robust)
* shift large enough to push the earliest cohort's effective onset to/below period 1
runpair , tag(c26) k(6) w("`W'") opts(y i t ttre, att robust espre(2) espost(2))
* subgroups() named on the ONSET variable itself
runpair , tag(c27) k(2) w("`W'") opts(y i t ttre, att subgroups(cohort) robust espre(2) espost(2))

* --------------------------------------------- k = 0 must be a no-op --------
qui use "`W'/d2.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2) anticipation(0)
dis "K0_RC=" _rc
if _rc == 0 dumpe , f("`W'/out_k0_A.txt")
qui use "`W'/d2.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
dis "K0B_RC=" _rc
if _rc == 0 dumpe , f("`W'/out_k0_B.txt")

* --------------------------------------------- k < 0 must be refused --------
qui use "`W'/d2.dta", clear
cap noi wooldid y i t ttre, att robust anticipation(-1)
dis "KNEG_RC=" _rc

dis "P2_DONE"
