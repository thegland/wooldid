* ===================================================================
* probe5_amevoid.do -- SELF-CONTAINED
* The dose-slope (R6/SI-A4) screen that voids the ame block is SKIPPED
* under fastnose (ado:11958).  Since 2026-09-14 fastnose POSTS e(b_es_*).
* Question: on a design where the analytic run withholds e(b_es_ame)
* because its ame points are structural zeros, does the fastnose twin
* post those zeros -- and under what message?
*
* DGP: 12 cohort-groups x 8 units x 12 periods.  i = GROUP (so each
* cohort-period cell holds 8 units -- the help's own recommended shape).
* Arm VOID : dose is CONSTANT WITHIN COHORT when treated, 0 otherwise
*            (the help names exactly this case).
* Arm LIVE : dose varies across units within cohort (identified).
* ===================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 96
gen long unit = _n
gen long grp  = ceil(_n/8)               // 12 groups, 8 units each
gen int ttre = .
replace ttre = 5  if inrange(grp,1,3)
replace ttre = 8  if inrange(grp,4,6)
replace ttre = 10 if inrange(grp,10,12)   // grp 7-9 never treated
gen double au = rnormal()
gen double dose_unit = 0.5 + runiform()   // unit-level dose
gen double dose_grp  = 0.5 + 0.1*grp      // group-level dose
expand 12
bysort unit: gen int t = _n
gen double e = rnormal()*0.4
gen byte D = (ttre<. & t>=ttre)
gen double doseV = cond(D==1, dose_grp , 0)   // constant within cohort
gen double doseL = cond(D==1, dose_unit, 0)   // varies within cohort
gen double yV = 1 + 0.3*au + 0.04*t + 0.9*doseV + e
gen double yL = 1 + 0.3*au + 0.04*t + 0.9*doseL + e
assert !missing(yV) & !missing(yL)
compress
tempfile p5
save "`p5'", replace
save "panel5.dta", replace

foreach arm in VOID LIVE {
  if "`arm'"=="VOID" {
    local yv yV
    local dv doseV
  }
  if "`arm'"=="LIVE" {
    local yv yL
    local dv doseL
  }
  foreach mode in analytic fastnose {
    dis "=========== ARM `arm' MODE `mode' ==========="
    use "`p5'", clear
    ereturn clear
    if "`mode'"=="analytic" cap noi wooldid `yv' grp t ttre, att cluster(grp) fe(unit) contreat(`dv') contreatcontrols(untreatedZero) esprelength(3) espostlength(3) esfixedbaseperiod
    if "`mode'"=="fastnose" cap noi wooldid `yv' grp t ttre, att cluster(grp) fe(unit) contreat(`dv') contreatcontrols(untreatedZero) esprelength(3) espostlength(3) esfixedbaseperiod fastnose
    dis "R_`arm'_`mode'_RC=" _rc
    foreach f in ame att {
      cap confirm matrix e(b_es_`f')
      dis "R_`arm'_`mode'_BES_`f'_RC=" _rc
      cap confirm matrix e(V_es_`f')
      dis "R_`arm'_`mode'_VES_`f'_RC=" _rc
      dis "R_`arm'_`mode'_NUMPRE_`f'=" e(numpre_`f')
    }
    cap confirm matrix e(b_es_ame)
    if _rc == 0 {
      cap mat list e(b_es_ame), format(%21x)
      tempname BM
      matrix `BM' = e(b_es_ame)
      local nz 0
      forval j = 1/`=colsof(`BM')' {
        if `BM'[1,`j'] != 0 local ++nz
      }
      dis "R_`arm'_`mode'_AME_NONZERO=" `nz' " of " colsof(`BM')
    }
  }
}
dis "=== PROBE5 DONE ==="
