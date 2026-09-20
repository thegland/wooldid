*=====================================================================
* BT-05 VERDICT REPRO (skeptic's own)
* BT-05 claimed: under cluster() at/above i, EVERY subgroups() row's analytic SE
*   equals the pooled one to 15 significant digits.
* This file shows (a) the identity needs a PERFECTLY BALANCED panel; on a ragged
*   panel it is gone; and (b) where it does hold, wooldid ALREADY detects it,
*   blanks the pairwise contrast and prints the named caution -- i.e. it is the
*   level-side face of the documented contrast degeneracy (sthlp boottstat entry,
*   wooldid_sgdeg ado:11092), not an undiscovered property.
* Look for: BAL_MAXRELDIF (tiny) / RAG_MAXRELDIF (not tiny) and BAL_CONTRAST_BLANKED=1
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 313131
set sortseed 313131

cap prog drop cmpse
prog def cmpse
  syntax , tag(string)
  tempname M A B
  mat `M' = e(wdidmainresults)
  mat `A' = e(wdidsgresults_0)
  mat `B' = e(wdidsgresults_1)
  local nr = rowsof(`M')
  local rfn : rowfullnames `M'
  local mx 0
  forval r = 1/`nr' {
    local nm : word `r' of `rfn'
    local d0 = reldif(`A'[`r',2], `M'[`r',2])
    local d1 = reldif(`B'[`r',2], `M'[`r',2])
    if `d0' > `mx' local mx = `d0'
    if `d1' > `mx' local mx = `d1'
    dis "`tag'_LEV `nm'  main=" %21x `M'[`r',2] "  sg0=" %21x `A'[`r',2] "  sg1=" %21x `B'[`r',2] ///
        "  reldif0=" %10.3e `d0' "  reldif1=" %10.3e `d1'
  }
  dis "`tag'_MAXRELDIF=" %12.4e `mx'
  tempname P
  mat `P' = e(pwsgs)
  local blanked = (`P'[1,2] >= .)
  dis "`tag'_CONTRAST_SE=" %14.9g `P'[1,2] "   `tag'_CONTRAST_BLANKED=" `blanked'
  dis "`tag'_sgdegen=" e(wdid_sgdegen) "  act=`e(wdid_sgdegenact)'  maxrho=`e(wdid_sgdegenr)'"
end

* ---- PERFECTLY BALANCED panel: 24 states x 3 units x 8 years, no gaps
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
tempfile bal
save `bal', replace

dis "### BAL: perfectly balanced panel, cluster(state) == i"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "BAL_RC=" _rc
cmpse , tag(BAL)

dis "### RAG: same panel with 8 percent of rows deleted at random"
use `bal', clear
set seed 989898
gen double _rr = runiform()
qui drop if _rr < 0.08
drop _rr
qui replace sg = inlist(state,1,2,3,7,8,13,19,20)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
dis "RAG_RC=" _rc
cmpse , tag(RAG)
dis "BT05_DONE"
