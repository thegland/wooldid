* =============================================================================
* R4 -- SKEPTIC: the neighbours nobody probed.
*
* (A) bootstrap(): under bootstrap() the e(V_es_*) blocks are REBUILT from the
*     replications.  With an exposure x post dose the pre-period ame is 0 in
*     EVERY replication, so the rebuilt pre block is the sd of a constant =
*     exactly 0.  Does the degenerate-row Caution and e(es_ndegen_ame) -- the
*     only thing standing between a user and the honestdid handoff -- still fire
*     when the bootstrap path supplies the matrix?
* (B) contreatcontrols(basic) on the same dose: same degeneracy?
* (C) control: a dose that is defined (nonzero) in EVERY period.
* (D) subgroups(): do the per-subgroup e(V_es_ame_s#) twins get the caution?
* (E) is the ame POINT estimate right?  low-noise design, true AME = 0.8.
*
* Data: r1_panel.dta from r1_own.do (self-contained).
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r4_boot.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242
which wooldid

cap prog drop vesdump
prog def vesdump
  syntax , tag(string)
  foreach f in att ame {
    cap confirm matrix e(V_es_`f')
    if !_rc {
      tempname E
      matrix `E' = e(V_es_`f')
      local k = colsof(`E')
      local nz 0
      forval j = 1/`k' {
        local d = `E'[`j',`j']
        if `d' == 0 local ++nz
        dis "VES_`tag'_`f'_d`j'=" %14.0g `d'
      }
      dis "VESZERO_`tag'_`f'=" `nz' " of " `k'
    }
    else dis "VES_`tag'_`f'_ABSENT=1"
    dis "NDEGEN_`tag'_`f'=" e(es_ndegen_`f')
    cap confirm matrix e(b_es_`f')
    dis "BES_`tag'_`f'=" (_rc==0)
  }
  * subgroup twins
  forval s = 1/4 {
    cap confirm matrix e(V_es_ame_s`s')
    if !_rc {
      tempname ES
      matrix `ES' = e(V_es_ame_s`s')
      local k = colsof(`ES')
      local nz 0
      forval j = 1/`k' {
        if `ES'[`j',`j'] == 0 local ++nz
      }
      dis "VESZERO_`tag'_ame_s`s'=" `nz' " of " `k'
      dis "NDEGEN_`tag'_ame_s`s'=" e(es_ndegen_ame_s`s')
    }
  }
  dis "BOOTVESFAILED_`tag'=" e(boot_ves_failed)
  dis "VCVSTATE_`tag'=" e(wdid_vcvstate)
  tempname FV
  cap confirm matrix e(V)
  if !_rc {
    matrix `FV' = e(V)
    local ns 0
    forval j = 1/`=colsof(`FV')' {
      if `FV'[`j',`j'] < 0 & !missing(`FV'[`j',`j']) local ++ns
    }
    dis "NSENT_`tag'=" `ns'
  }
  else dis "NOV_`tag'=1"
end

* ============================= (A) bootstrap ===============================
use "r1_panel.dta", clear
dis "############### A: bootstrap(standard) reps(100), untreatedZero ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) cluster(hid) bootstrap(standard) reps(100) rseed(9911) dots(0)
dis "RC_A=" _rc
if _rc == 0 {
  vesdump , tag(BOOT)
  cap confirm matrix e(V_es_ame)
  if !_rc {
    matrix BB = e(b_es_ame)
    matrix VB = e(V_es_ame)
    matrix list BB
    matrix list VB, format(%12.0g)
    local NPB = e(numpre_ame)
    dis "--- honestdid on the BOOTSTRAP ame block ---"
    cap noisily honestdid, b(BB) vcov(VB) numpre(`NPB') mvec(0(0.5)1) delta(rm) parallel(0) matasave(BR)
    dis "HON_BOOT_RC=" _rc
    cap noisily mata: BR.CI
    cap noisily mata: printf("{txt}HON_BOOT_INVARIANT=%g\n", (BR.CI[2,2]==BR.CI[3,2]) & (BR.CI[3,2]==BR.CI[4,2]) & (BR.CI[2,3]==BR.CI[3,3]) & (BR.CI[3,3]==BR.CI[4,3]))
  }
}

* =========================== (B) basic controls ============================
use "r1_panel.dta", clear
dis "############### B: contreatcontrols(basic), robust ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(basic) adversarial jointtests customvce(robust)
dis "RC_B=" _rc
if _rc == 0 {
  vesdump , tag(BASIC)
  dis "ADVP_BASIC_ame=" e(adversarial_pval_ame)
  dis "ADVP_BASIC_att=" e(adversarial_pval_att)
}

* ================== (C) control: dose nonzero in all periods ===============
use "r1_panel.dta", clear
replace dose = expo + 0.3*(expyr<. & year>=expyr)
dis "############### C: dose defined in EVERY period ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(basic) adversarial jointtests customvce(robust)
dis "RC_C=" _rc
if _rc == 0 {
  vesdump , tag(ALLPER)
  dis "ADVP_ALLPER_ame=" e(adversarial_pval_ame)
}

* ============================ (D) subgroups ================================
use "r1_panel.dta", clear
gen byte sg = mod(hid, 2)
dis "############### D: subgroups(sg) + untreatedZero ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) subgroups(sg) customvce(robust)
dis "RC_D=" _rc
if _rc == 0 {
  vesdump , tag(SG)
  dis "=== e() matrix names (SG) ==="
  ereturn list
}

* ================ (E) is the ame POINT estimate right? =====================
* same shape, tiny noise, big within-cell dose spread.  True AME = 0.8 exactly.
clear
set obs 12
gen byte state = _n
gen int expyr = .
replace expyr = 2012 if inlist(state,1,2,3,4)
replace expyr = 2014 if inlist(state,5,6,7)
replace expyr = 2016 if inlist(state,8,9)
expand 60
bys state: gen int hosp = _n
gen long hid = state*1000 + hosp
gen double expo = 0.2 + 1.6*runiform()
gen double hfe = rnormal(0,1)
expand 12
bys hid: gen int year = 2007 + _n
gen double dose = expo*(expyr<. & year>=expyr)
gen double y = 4 + hfe + 0.03*(year-2007) + 0.8*dose + rnormal(0,0.05)
dis "############### E: low-noise, true AME = 0.8 ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) customvce(robust)
dis "RC_E=" _rc
if _rc == 0 {
  matrix BE = e(b)
  local bn : colfullnames BE
  local j 0
  foreach nm of local bn {
    local ++j
    dis "E_ROW|" %-22s "`nm'" "|b=" %14.10f BE[1,`j']
  }
}

dis "R4_DONE=1"
