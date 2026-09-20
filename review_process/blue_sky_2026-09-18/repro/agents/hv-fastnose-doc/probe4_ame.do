* probe4_ame.do -- SELF-CONTAINED: ame-family behaviour under fastnose
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
clear
set obs 48
gen long i = _n
gen byte coh = mod(_n,4)
gen int ttre = .
replace ttre = 6  if coh==1
replace ttre = 9  if coh==2
replace ttre = 12 if coh==3
gen double ui = rnormal()
gen double dose_u = 0.5 + runiform()
expand 16
bysort i: gen int t = _n
gen double et = rnormal()*0.5
gen byte D = (ttre<. & t>=ttre)
gen double dose = cond(D==1, dose_u, 0)
gen double yc = 1 + 0.30*ui + 0.05*t + 0.9*dose + et
compress
tempfile p4
save "`p4'", replace

foreach arm in fn an {
  dis "=== G1 arm `arm' (contreatcontrols(untreatedZero)) ==="
  use "`p4'", clear
  ereturn clear
  if "`arm'"=="fn" cap noi wooldid yc i t ttre, att cluster(i) contreat(dose) contreatcontrols(untreatedZero) esprelength(3) espostlength(3) esfixedbaseperiod fastnose
  if "`arm'"=="an" cap noi wooldid yc i t ttre, att cluster(i) contreat(dose) contreatcontrols(untreatedZero) esprelength(3) espostlength(3) esfixedbaseperiod
  dis "G1_`arm'_RC=" _rc
  foreach f in ame ame_it att att_it {
    cap confirm matrix e(b_es_`f')
    dis "G1_`arm'_BES_`f'_RC=" _rc
    cap confirm matrix e(V_es_`f')
    dis "G1_`arm'_VES_`f'_RC=" _rc
    dis "G1_`arm'_NUMPRE_`f'=" e(numpre_`f')
  }
  cap confirm matrix e(V)
  dis "G1_`arm'_V_RC=" _rc
  cap local bn : colnames e(b)
  dis "G1_`arm'_BNAMES=" "`bn'"
}
dis "=== PROBE4 DONE ==="
