* p14_messy.do -- messy-but-plausible administrative panel shapes.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime
local ESW esprelength(2) espostlength(2)

cap prog drop nrun
prog def nrun
  syntax anything(name=tg) , CMD(string asis)
  ereturn clear
  cap noi `cmd'
  dis "`tg'_RC=" _rc
  if _rc != 0 exit
  gtcheck , tag(`tg') pre es
  gtdefcheck , tag(`tg')
  gtwincheck , tag(`tg')
end

*---- N1: an `if' whose expression is missing on some rows -------------------
mkpanel , n(90) t(10) onsets("2004 2006") ysd(0.4) seed(2211)
gen double scr = rnormal()
qui replace scr = . in 1/200
tempfile b1
qui save `b1'
nrun N1_ifmissing , cmd(wooldid y cohp year ttre if scr > -3, `EST' gteffects `ESW')

*---- N2: fweights ----------------------------------------------------------
use `b1', clear
bys unit (year): gen int fw = 1 + mod(unit, 4) if _n == 1
bys unit (year): replace fw = fw[1]
nrun N2_fw , cmd(wooldid y cohp year ttre [fw=fw], `EST' gteffects `ESW')

*---- N3: a singleton treated cohort (one unit) ------------------------------
use `b1', clear
qui replace cohp = 4 if unit == 5
qui replace ttre = 2005 if unit == 5
nrun N3_singleton , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW')

*---- N4: a cohort treated in the LAST observed period -----------------------
use `b1', clear
qui replace ttre = 2009 if cohp == 3
nrun N4_lastperiod , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW')

*---- N5: cohorts treated after the sample end, as controls -----------------
use `b1', clear
qui replace ttre = 2012 if cohp == 3
nrun N5_postsample , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW')
use `b1', clear
qui replace ttre = 2012 if cohp == 3
nrun N5b_oos , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' oostreatedcontrols(4))

*---- N6: each unit observed in a random subset of periods ------------------
mkpanel , n(150) t(10) onsets("2004 2006") ysd(0.4) gapprob(0.45) seed(3311)
nrun N6_sparse , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW')

*---- N7: unit-level i with cluster at a coarser level, aweights -------------
mkpanel , n(90) t(10) onsets("2004 2006 2008") ysd(0.4) entryprob(.3) exitprob(.3) seed(4411)
gen byte st = 1 + mod(unit, 10)
bys unit (year): gen double wt = 0.3 + 2*runiform() if _n == 1
bys unit (year): replace wt = wt[1]
nrun N7_unit_clus , cmd(wooldid y unit year ttre [aw=wt], `EST' gteffects `ESW' cluster(st))


*---- N8: help claims about the se column -----------------------------------
* "under timetrends the coefficient-level variances are typically not identified
*  and the se column is missing (degenerate zero variances are always reported as
*  missing, never zero)"
cap prog drop secheck
prog def secheck
  syntax , TAG(string)
  cap confirm matrix e(gteffects)
  if _rc != 0 exit
  preserve
  clear
  qui svmat double e(gteffects), names(col)
  qui count
  local n = r(N)
  qui count if se >= .
  local nm = r(N)
  qui count if se == 0
  local nz = r(N)
  qui count if p >= .
  local np = r(N)
  dis as result "SE|`tag'|rows=`n'|se_missing=`nm'|se_exactly_zero=`nz'|p_missing=`np'"
  restore
end

mkpanel , n(90) t(10) onsets("2004 2006") ysd(0.4) seed(2211)
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW' timetrends
dis "N8_RC=" _rc
secheck , tag(N8_timetrends)
gtcheck , tag(N8_timetrends) pre es

ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW'
dis "N9_RC=" _rc
secheck , tag(N9_plain)

ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects `ESW' fastnose
dis "N10_RC=" _rc
secheck , tag(N10_fastnose)

dis "DONE_P14"
