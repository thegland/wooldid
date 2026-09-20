* ===================================================================
* probe6_amevoid2.do -- SELF-CONTAINED
*
* The R6/SI-A4 dose-slope screen (ado:11958) that VOIDS the ame block
* -- withholding e(b_es_ame) because its entries are structural zeros
* rather than estimates -- is gated on `fastnose' == "".  Before
* 2026-09-14 that was harmless: fastnose posted no e(b_es_*) at all.
* Since 6fe814d, fastnose POSTS e(b_es_*).  So on a design where the
* analytic run withholds e(b_es_ame) and names the problem, does the
* fastnose twin now hand back the fabricated zeros?
*
* GEOMETRY (the help's own example of the void: "dose constant within
* each cohort").  30 cohorts x 6 units x 12 periods, i = cohort,
* robust VCE (so the cell-containment R1 predicate is NOT what fires).
* ===================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 30
gen long i = _n
gen int ttre = cond(mod(_n,2)==0, 6, .)
gen double dcoh = exp(0.2*mod(_n,5))        // CONSTANT within cohort
expand 6
bysort i: gen int unit = _n
gen double dunit = exp(0.8*runiform() + 0.2*mod(i,5))   // varies within cohort
expand 12
bysort i unit: gen int t = _n
gen double y0 = 0.1*i + 0.05*t + rnormal()
gen byte D = (ttre<. & t>=ttre)
gen double yV = y0 + cond(D==1, 0.5*dcoh , 0)
gen double yL = y0 + cond(D==1, 0.5*dunit, 0)
sort i unit t
compress
save "panel6.dta", replace

foreach arm in VOID LIVE {
  if "`arm'"=="VOID" {
    local yv yV
    local dv dcoh
  }
  if "`arm'"=="LIVE" {
    local yv yL
    local dv dunit
  }
  foreach mode in analytic fastnose {
    dis "=========== ARM `arm' MODE `mode' ==========="
    use "panel6.dta", clear
    ereturn clear
    if "`mode'"=="analytic" cap noi wooldid `yv' i t ttre, att robust contreat(`dv') contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod
    if "`mode'"=="fastnose" cap noi wooldid `yv' i t ttre, att robust contreat(`dv') contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod fastnose
    dis "Q_`arm'_`mode'_RC=" _rc
    foreach f in ame att {
      cap confirm matrix e(b_es_`f')
      dis "Q_`arm'_`mode'_BES_`f'_RC=" _rc
      cap confirm matrix e(V_es_`f')
      dis "Q_`arm'_`mode'_VES_`f'_RC=" _rc
      dis "Q_`arm'_`mode'_NUMPRE_`f'=" e(numpre_`f')
    }
    cap confirm matrix e(b_es_ame)
    if _rc == 0 {
      tempname BM
      matrix `BM' = e(b_es_ame)
      local nz 0
      forval j = 1/`=colsof(`BM')' {
        if `BM'[1,`j'] != 0 local ++nz
      }
      dis "Q_`arm'_`mode'_AME_NONZERO=" `nz' " of " colsof(`BM')
      mat list `BM', format(%21x)
    }
  }
}
dis "=== PROBE6 DONE ==="
