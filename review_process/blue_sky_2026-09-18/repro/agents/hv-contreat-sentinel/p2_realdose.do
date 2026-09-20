* ============================================================================
* P2 -- does the -1e100 sentinel / zero-diagonal fire when the continuous
*       treatment VARIES IN EVERY PERIOD (a genuine intensity measure that
*       exists before onset too)?  This is the discriminating realism probe.
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771331
set sortseed 771331
dis "=== WHICH ==="
which wooldid

cap prog drop dumpV
prog def dumpV
  syntax , tag(string)
  dis "##### DUMP_`tag' #####"
  cap confirm matrix e(b)
  if _rc != 0 {
    dis "DUMP_`tag'_NOB=1"
    dis "##### ENDDUMP_`tag' #####"
    exit
  }
  tempname B V
  matrix `B' = e(b)
  local bn : colfullnames `B'
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
    if `hasV' & `vj' == 0 local ++nzero
    dis "ROW_`tag' |" %-26s "`nm'" "| b=" %11.0g `bj' "| V=" %11.0g `vj' "|`mark'"
  }
  dis "DUMP_`tag'_NSENT=" `nsent'
  dis "DUMP_`tag'_NZERODIAG_IN_V=" `nzero'
  foreach f in ame ame_it {
    cap confirm matrix e(V_es_`f')
    if _rc == 0 {
      tempname E
      matrix `E' = e(V_es_`f')
      local en : colfullnames `E'
      local nz 0
      local jj 0
      local uniqd ""
      foreach en1 of local en {
        local ++jj
        local ejj = `E'[`jj',`jj']
        if `ejj' == 0 local ++nz
        dis "VES_`tag'_`f' |" %-10s "`en1'" "| diag=" %14.0g `ejj'
      }
      dis "VESSUM_`tag'_`f'_NZERO=" `nz'
    }
    else dis "VES_`tag'_`f'_ABSENT=1"
    local nd = e(es_ndegen_`f')
    dis "NDEGEN_`tag'_`f'=" `nd'
  }
  dis "##### ENDDUMP_`tag' #####"
end

* ---------------------------------------------------------------- DGP
* A genuine intensity measure: exposure index observed in EVERY period for
* EVERY unit (treated and never-treated), moving over time.
clear
set obs 100
gen long unit = _n
gen byte coh = ceil(unit/20)
gen int ttre = .
replace ttre = 4  if coh==1
replace ttre = 6  if coh==2
replace ttre = 8  if coh==3
replace ttre = 10 if coh==4
gen double base_i = 1 + 2*runiform()
gen double ufe   = rnormal()
expand 12
bys unit: gen int t = _n
gen double dose = base_i + 0.15*t + 0.4*rnormal()
replace dose = 0.05 if dose < 0.05
gen byte post = (ttre < . & t >= ttre)
gen double eps = rnormal()
gen double y = ufe + 0.2*t + 0.8*dose*post + eps
save "p2_panel.dta", replace
sum dose

foreach cc in basic byTime byTreated {
  foreach vv in NONE robust clusterunit {
    use "p2_panel.dta", clear
    local vopt ""
    if "`vv'" == "robust"      local vopt customvce(robust)
    if "`vv'" == "clusterunit" local vopt cluster(unit)
    dis "=== CONFIG `cc'_`vv' ==="
    cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(`cc') `vopt'
    dis "RC_`cc'_`vv'=" _rc
    if _rc == 0 dumpV , tag(`cc'_`vv')
  }
}

* and the untreatedZero comparator on the SAME varying dose
use "p2_panel.dta", clear
dis "=== CONFIG untreatedZero_NONE ==="
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero)
dis "RC_uz_NONE=" _rc
if _rc == 0 dumpV , tag(uz_NONE)

dis "=== ALLDONE ==="
