clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
mata:
real scalar sk_add(real scalar a, real scalar b)
{
  return(a+b)
}
end
dis "PRE_WOOLDID"
mata: mata describe
cap noisily mata: x = sk_add(2,3); st_numscalar("XX", x)
dis "PRE_RC=" _rc
cap dis "XX=" XX
use "r1_panel.dta", clear
qui wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero)
dis "POST_WOOLDID rc=" _rc
mata: mata describe
cap noisily mata: y2 = sk_add(4,5); st_numscalar("YY", y2)
dis "POST_RC=" _rc
cap dis "YY=" YY
dis "MATATEST_DONE=1"
