* ============================================================================
* P1 -- reproduce recon F5: contreat() ships b=0 with e(V) diag = -1e100,
*       while e(V_es_ame) carries EXACTLY 0 on the same horizons.
* Self-contained: builds its own panel.
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

dis "=== WHICH ==="
which wooldid

* ---------------------------------------------------------------- diagnostics
cap prog drop dumpV
prog def dumpV
  syntax , tag(string)
  dis "##### DUMP_`tag' #####"
  dis "DUMP_`tag'_RC=" _rc
  cap confirm matrix e(b)
  if _rc != 0 {
    dis "DUMP_`tag'_NOB=1"
    exit
  }
  tempname B V
  matrix `B' = e(b)
  local bn : colfullnames `B'
  dis "DUMP_`tag'_K=" colsof(`B')
  cap confirm matrix e(V)
  local hasV = (_rc == 0)
  dis "DUMP_`tag'_HASV=" `hasV'
  if `hasV' matrix `V' = e(V)
  local nsent 0
  local nzero 0
  local j 0
  foreach nm of local bn {
    local ++j
    local bj = `B'[1,`j']
    if `hasV' local vj = `V'[`j',`j']
    else local vj = .
    local mark ""
    if `hasV' & `vj' < 0 & !missing(`vj') {
      local mark "SENTINEL"
      local ++nsent
    }
    if `hasV' & `vj' == 0 {
      local mark "ZERODIAG"
      local ++nzero
    }
    dis "ROW_`tag' | " %-28s "`nm'" " | b=" %12.0g `bj' " | Vdiag=" %12.0g `vj' " | `mark'"
  }
  dis "DUMP_`tag'_NSENT=" `nsent'
  dis "DUMP_`tag'_NZERO=" `nzero'
  * event-study copies
  foreach f in att att_it ame ame_it {
    cap confirm matrix e(V_es_`f')
    if _rc == 0 {
      tempname E EB
      matrix `E' = e(V_es_`f')
      local en : colfullnames `E'
      local nz 0
      local nneg 0
      local jj 0
      foreach en1 of local en {
        local ++jj
        local ejj = `E'[`jj',`jj']
        local mk ""
        if `ejj' == 0 local mk "ZERODIAG"
        if `ejj' < 0 & !missing(`ejj') local mk "NEGDIAG"
        if `ejj' == 0 local ++nz
        if `ejj' < 0 & !missing(`ejj') local ++nneg
        dis "VES_`tag'_`f' | " %-12s "`en1'" " | diag=" %12.0g `ejj' " | `mk'"
      }
      dis "VESSUM_`tag'_`f'_NZERO=" `nz' " NNEG=" `nneg'
    }
    else dis "VES_`tag'_`f'_ABSENT=1"
    cap confirm matrix e(b_es_`f')
    if _rc == 0 dis "BES_`tag'_`f'_PRESENT=1"
    else dis "BES_`tag'_`f'_ABSENT=1"
    cap dis "NDEGEN_`tag'_`f'=" e(es_ndegen_`f')
  }
  dis "##### ENDDUMP_`tag' #####"
end

* ---------------------------------------------------------------- DGP
* 5 cohorts x 20 units x 12 periods.  cohort 5 = never treated.
* dose is a UNIT-level intensity that switches on at onset (the canonical
* continuous-DiD "exposure x post" design).
clear
set obs 100
gen long unit = _n
gen byte coh = ceil(unit/20)
gen int ttre = .
replace ttre = 4  if coh==1
replace ttre = 6  if coh==2
replace ttre = 8  if coh==3
replace ttre = 10 if coh==4
gen double dose0 = 0.5 + 1.5*runiform()
gen double ufe   = rnormal()
expand 12
bys unit: gen int t = _n
gen double tfe = 0.2*t
gen double dose = dose0*(ttre<. & t>=ttre)
gen double eps = rnormal()
gen double y = ufe + tfe + 0.8*dose + eps
save "p1_panel.dta", replace

dis "=== CONFIG A: untreatedZero, default vce (iid) ==="
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero)
dis "A_RC=" _rc
dumpV , tag(A)

dis "=== CONFIG B: basic, default vce ==="
use "p1_panel.dta", clear
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(basic)
dis "B_RC=" _rc
dumpV , tag(B)

dis "=== CONFIG C: untreatedZero, robust ==="
use "p1_panel.dta", clear
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero) customvce(robust)
dis "C_RC=" _rc
dumpV , tag(C)

dis "=== CONFIG D: untreatedZero, cluster(unit) [finer than i=coh] ==="
use "p1_panel.dta", clear
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero) cluster(unit)
dis "D_RC=" _rc
dumpV , tag(D)

dis "=== CONFIG E: basic, cluster(unit), with fe(unit) [help-recommended geometry] ==="
use "p1_panel.dta", clear
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(basic) cluster(unit) fe(unit)
dis "E_RC=" _rc
dumpV , tag(E)

dis "=== ALLDONE ==="
