* =========================================================================
* P1 -- tests HV-FN1: does `fastnose' post e(b_es_att) / e(numpre_att)?
* Self-contained: regenerates its own data inline (no external .dta needed).
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
bysort st: gen double afe = rnormal() if _n == 1
bysort st (afe): replace afe = afe[1]
sort row
gen double y = 2 + afe + 0.03*(yr-2005) + 0.5*D + rnormal(0,1)

prog def rep
  args tag
  cap confirm matrix e(b_es_att)
  dis "`tag'_BESATT_RC=" _rc
  cap confirm matrix e(V_es_att)
  dis "`tag'_VESATT_RC=" _rc
  cap confirm matrix e(V)
  dis "`tag'_V_RC=" _rc
  local np = e(numpre_att)
  dis "`tag'_NUMPRE=[" "`np'" "]"
  dis "`tag'_VCVSTATE=[" e(wdid_vcvstate) "]"
  cap local cn : colnames e(b_es_att)
  dis "`tag'_COLNAMES=[" "`cn'" "]"
  cap mat li e(b_es_att), format(%21x)
end

dis "=== ARM 1: analytic ==="
cap noi wooldid y st yr adopt, att cluster(st) esprelength(3) espostlength(4) esfixedbaseperiod
dis "A1_RC=" _rc
rep A1
cap mat BA = e(b_es_att)

dis "=== ARM 2: fastnose ==="
cap noi wooldid y st yr adopt, att cluster(st) esprelength(3) espostlength(4) esfixedbaseperiod fastnose
dis "A2_RC=" _rc
rep A2
cap mat BF = e(b_es_att)

dis "=== compare the two coefficient vectors ==="
cap noi {
  local mx = 0
  forval j = 1/`=colsof(BA)' {
    local d = abs(BA[1,`j'] - BF[1,`j'])
    if `d' > `mx' local mx = `d'
  }
  dis "MAXABSDIFF=" %21x `mx'
  dis "MAXABSDIFF_dec=" `mx'
}
dis "DIFF_RC=" _rc

dis "=== ARM 3: fastnose, does e(b) still carry ES_att cells? ==="
cap noi mat li e(b)
dis "P1_DONE"
