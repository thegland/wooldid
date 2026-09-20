* ===================================================================
* probe3_surface2.do -- SELF-CONTAINED (uses panel2.dta built by probe2,
* rebuilt here so the file stands alone)
*  F1 fastnose + contreat(): ame family
*  F2 permute(shuffle) fast mode: b_es must NOT post (carve-out)
*  F3 is e(b_es_att) byte-identical to the ES_att_* cells of e(b)?
*  F4 honestdid handoff after a fastnose run
*  F5 user-supplied fastnose TOGETHER with permute()
* ===================================================================
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
gen byte sg = mod(_n,2)
gen double ui = rnormal()
gen double dose_u = 0.5 + runiform()
expand 16
bysort i: gen int t = _n
gen double et = rnormal()*0.5
gen byte D = (ttre<. & t>=ttre)
gen double rel = t - ttre
gen double dose = cond(D==1, dose_u, 0)
gen double y  = 1 + 0.30*ui + 0.05*t + cond(D==1, 0.8 + 0.15*rel, 0) + 0.4*sg + et
gen double yc = 1 + 0.30*ui + 0.05*t + 0.9*dose + et
compress
save "panel3.dta", replace

* ---------------- F1: contreat, fastnose vs analytic ---------------
foreach arm in fn an {
  dis "=== F1 arm `arm' ==="
  use "panel3.dta", clear
  ereturn clear
  if "`arm'"=="fn" cap noi wooldid yc i t ttre, att cluster(i) contreat(dose) contreatcontrols(byCohort) esprelength(3) espostlength(3) esfixedbaseperiod fastnose
  if "`arm'"=="an" cap noi wooldid yc i t ttre, att cluster(i) contreat(dose) contreatcontrols(byCohort) esprelength(3) espostlength(3) esfixedbaseperiod
  dis "F1_`arm'_RC=" _rc
  foreach f in ame att {
    cap confirm matrix e(b_es_`f')
    dis "F1_`arm'_BES_`f'_RC=" _rc
    cap confirm matrix e(V_es_`f')
    dis "F1_`arm'_VES_`f'_RC=" _rc
    dis "F1_`arm'_NUMPRE_`f'=" e(numpre_`f')
  }
  cap confirm matrix e(V)
  dis "F1_`arm'_V_RC=" _rc
}

* ---------------- F2: permute fast mode ----------------------------
dis "=== F2 permute(shuffle), no permtstat ==="
use "panel3.dta", clear
ereturn clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod permute(shuffle) permreps(20) rseed(7)
dis "F2_RC=" _rc
cap confirm matrix e(b_es_att)
dis "F2_BES_RC=" _rc
cap confirm matrix e(V_es_att)
dis "F2_VES_RC=" _rc
dis "F2_NUMPRE=" e(numpre_att)

* ---------------- F5: user fastnose + permute ----------------------
dis "=== F5 permute + explicit fastnose ==="
use "panel3.dta", clear
ereturn clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod permute(shuffle) permreps(20) rseed(7) fastnose
dis "F5_RC=" _rc
cap confirm matrix e(b_es_att)
dis "F5_BES_RC=" _rc
dis "F5_NUMPRE=" e(numpre_att)

* ---------------- F3: b_es vs e(b) cells ---------------------------
foreach arm in fn an {
  dis "=== F3 arm `arm' ==="
  use "panel3.dta", clear
  ereturn clear
  if "`arm'"=="fn" qui wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod fastnose
  if "`arm'"=="an" qui wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod
  matrix F3b   = e(b)
  matrix F3bes = e(b_es_att)
  local bn : colnames F3b
  local en : colnames F3bes
  local ndiff 0
  local jj 0
  foreach c of local en {
    local ++jj
    local kk 0
    local hit 0
    foreach cb of local bn {
      local ++kk
      if "`cb'" == "ES_att_`c'" {
        local hit 1
        if F3b[1,`kk'] != F3bes[1,`jj'] local ++ndiff
        dis "F3_`arm'_`c' eb=" %21x F3b[1,`kk'] " bes=" %21x F3bes[1,`jj'] " same=" (F3b[1,`kk']==F3bes[1,`jj'])
      }
    }
    if `hit'==0 dis "F3_`arm'_`c'_NOMATCH"
  }
  dis "F3_`arm'_NDIFF=" `ndiff'
}

* ---------------- F4: honestdid after fastnose ---------------------
dis "=== F4 honestdid after fastnose ==="
use "panel3.dta", clear
ereturn clear
qui wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod fastnose
dis "F4_NUMPRE=" e(numpre_att)
cap noi honestdid, pre(1/5) post(6/11) mvec(0.5) parallel(0)
dis "F4_HD_RC=" _rc

dis "=== PROBE3 DONE ==="
