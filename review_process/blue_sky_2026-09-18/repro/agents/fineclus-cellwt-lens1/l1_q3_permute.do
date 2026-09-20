*! l1 Q3 (re-run): the SAME phantom row under permute(). Loud, or silent?
*! The first attempt used permreps(19), below wooldid's floor of 20, and died r(198)
*! on my error in both arms. Fixed here.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cap prog drop mkpanel
prog def mkpanel
  syntax , dseed(int)
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs 15
  gen long cty = _n
  gen double a_g = rnormal(0,.4)
  gen onset = .
  replace onset = 4 if cty <= 5
  replace onset = 6 if cty > 5 & cty <= 10
  expand 4
  sort cty
  gen long unit = _n
  gen double b_i = rnormal(0,.4)
  gen double uf = rnormal(0,1)
  expand 7
  bysort unit (cty): gen int t = _n
  gen double tau = 1 + a_g + b_i
  gen byte D = (onset < . & t >= onset)
  gen double y = 2 + uf + .1*t + tau*D + rnormal(0,1)
  replace y = . if unit == 1 & t == 3
  keep y unit cty t onset
end

dis "--- PERMUTE on A data (no phantom straddler) ---"
mkpanel , dseed(7201)
cap noi wooldid y unit t onset, att att_i cluster(cty) permute(shuffle) permreps(20) rseed(41) dots(0)
dis "PERM_A_RC=" _rc

dis "--- PERMUTE on D data (phantom row recorded under cty 2) ---"
mkpanel , dseed(7201)
qui replace cty = 2 if unit == 1 & t == 3
cap noi wooldid y unit t onset, att att_i cluster(cty) permute(shuffle) permreps(20) rseed(41) dots(0)
dis "PERM_D_RC=" _rc
dis "DONE_Q3"
