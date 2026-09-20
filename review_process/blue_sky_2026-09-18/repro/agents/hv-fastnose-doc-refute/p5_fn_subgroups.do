* =========================================================================
* P5 -- my own probe, not seeded.  6fe814d (2026-09-14) made `fastnose' post
* the event-study b/numpre.  The help's fastnose paragraph (sthlp:286) states
* flatly "e(V_es_) is not returned, so honestdid and Roth-style diagnostics
* have no covariance to work from", while ado:1981 keeps a VARIANCE on the
* subgroup-by-horizon margins call under fastnose.  Does a per-subgroup
* e(V_es_att_s*) therefore survive into e() under fastnose?
* Self-contained.
* =========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 12800
gen long row = _n
gen int  st  = ceil(row/320)
gen int  yr  = 2005 + mod(ceil(row/20)-1, 16)
gen int adopt = .
replace adopt = 2010 if inrange(st, 1, 10)
replace adopt = 2013 if inrange(st, 11, 20)
replace adopt = 2016 if inrange(st, 21, 28)
gen byte D = (yr >= adopt) & !missing(adopt)
* subgroup constant within the (i,t) cell: region, a state-level attribute
gen byte region = 1 + mod(st, 2)
bysort st: gen double afe = rnormal() if _n == 1
bysort st (afe): replace afe = afe[1]
sort row
gen double y = 2 + afe + 0.03*(yr-2005) + (0.4 + 0.3*(region==2))*D + rnormal(0,1)

prog def sgrep
  args tag
  foreach nm in b_es_att V_es_att b_es_att_s1 V_es_att_s1 b_es_att_s2 V_es_att_s2 {
    cap confirm matrix e(`nm')
    dis "`tag'_`nm'_RC=" _rc
  }
  dis "`tag'_numpre_att=[" e(numpre_att) "]"
  dis "`tag'_numpre_att_s1=[" e(numpre_att_s1) "]"
  dis "`tag'_numpre_att_s2=[" e(numpre_att_s2) "]"
  dis "`tag'_VCVSTATE=[" e(wdid_vcvstate) "]"
  cap noi mat li e(V_es_att_s1), format(%12.8f)
  dis "`tag'_LISTV_S1_RC=" _rc
end

dis "=== SG ARM 1: analytic ==="
cap noi wooldid y st yr adopt, att cluster(st) subgroups(region) esprelength(3) espostlength(4) esfixedbaseperiod
dis "S1_RC=" _rc
sgrep S1

dis "=== SG ARM 2: fastnose ==="
cap noi wooldid y st yr adopt, att cluster(st) subgroups(region) esprelength(3) espostlength(4) esfixedbaseperiod fastnose
dis "S2_RC=" _rc
sgrep S2
dis "--- full ereturn matrix inventory, fastnose arm ---"
cap noi ereturn list
dis "P5_DONE"
