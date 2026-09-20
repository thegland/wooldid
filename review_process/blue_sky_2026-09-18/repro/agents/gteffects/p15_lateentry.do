* p15_lateentry.do -- ado:12056-12059 drops, WITHOUT ANY MESSAGE, every `i' whose
* effective onset is at or before ITS OWN first observed period:
*     gegen byte nopre = max(ttre <= mint) if ifstatement, by(i)
*     replace ifvar = 0 if nopre == 1
* `mint' is gegen min(t) BY(i), i.e. unit-specific when i is the unit.  The help
* (sthlp:73) describes this as "cohorts treated before the sample begins are dropped",
* which does not describe a unit that ENTERS the panel late and is already treated.
* The anticipation() version of the very same drop, twelve lines above, prints a caution
* naming the unit and observation counts.  This one prints nothing.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime

* an administrative register: 40 never-treated firms present throughout; 20 firms
* treated in 2004 present throughout; 10 firms treated in 2004 that only enter the
* register in 2004 (they exist, they are treated, they have post-treatment data).
cap prog drop mk
prog def mk
  syntax , ENTry(int)
  clear
  set obs 70
  gen long unit = _n
  gen byte late = (_n > 60)
  gen double ttre = cond(_n > 40, 2004, .)
  gen double ai = rnormal()
  qui expand 11
  bys unit: gen int year = 1999 + _n
  qui drop if late == 1 & year < `entry'
  tempvar ts
  bys year: gen double `ts' = rnormal() if _n==1
  bys year: replace `ts' = `ts'[1]
  gen double te = 0
  qui replace te = 1 + 0.25*(year - ttre) if year >= ttre & ttre < .
  * the late entrants have a LARGER effect, so dropping them moves the number
  qui replace te = te*2.5 if late == 1 & year >= ttre & ttre < .
  gen double y = 2 + ai + 0.5*`ts' + te + rnormal()*0.3
  gen byte coh = cond(ttre < ., 2, 1)
end

foreach e in 2000 2003 2004 2005 {
  dis as result _n "================ late entrants first observed in `e' ================"
  mk , entry(`e')
  qui count
  local nobs = r(N)
  qui gegen byte tgu = tag(unit)
  qui count if tgu == 1 & ttre < .
  local ntreatunits = r(N)
  drop tgu
  dis "P15_`e'_INPUT: obs=`nobs' treated_units=`ntreatunits'"
  ereturn clear
  cap noi wooldid y unit year ttre, `EST' gteffects esprelength(2) espostlength(2)
  dis "P15_`e'_RC=" _rc
  if _rc == 0 {
    dis "P15_`e'_N=" e(N) " N_dropped=" e(N_dropped)
    * how many distinct treated i appear in the table?
    preserve
    clear
    qui svmat double e(gteffects), names(col)
    qui bys i: keep if _n == 1
    qui count
    dis "P15_`e'_UNITS_IN_TABLE=" r(N)
    restore
    dis "P15_`e'_att=" %21x e(b)[1,colnumb(e(b),"main_att")]
    gtcheck , tag(P15_`e') pre es
  }
}

* the same data with i = COHORT rather than i = unit: `mint' is then the cohort-wide
* minimum, so the same late entrants are NOT dropped.
dis as result _n "================ i = cohort, entrants first observed in 2004 ================"
mk , entry(2004)
ereturn clear
cap noi wooldid y coh year ttre, `EST' gteffects esprelength(2) espostlength(2)
dis "P15_COH_RC=" _rc
dis "P15_COH_N=" e(N) " N_dropped=" e(N_dropped)
dis "P15_COH_att=" %21x e(b)[1,colnumb(e(b),"main_att")]
gtcheck , tag(P15_coh) pre es

dis "DONE_P15"
