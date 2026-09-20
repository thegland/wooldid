* ======================================================================
* LENS-0 INDEPENDENT REPRODUCTION of fineclus-perm-FP1
* Built from the CLAIM alone. No prior agent's script was read or run.
* CLAIM: under permonsetcohorts, a covariate keyed to the OBSERVED cohort
*        is inert for the point estimate but rewrites the permutation
*        reference distribution.
* ======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* ---------------------------------------------------------------- panel
* 7 administrative states, ragged reporting windows, state level + trend
* heterogeneity, common year shocks.  i = coh IS the onset grouping.
* onsets: states 1,2 -> 2006 ; states 3,4 -> 2010 ; states 5,6,7 -> never
* M = 7!/(2!2!3!) = 210 enumerated assignments (exact p grid 1/210)
clear
set obs 7
gen int state = _n
gen double onset = .
replace onset = 2006 if inlist(state,1,2)
replace onset = 2010 if inlist(state,3,4)
gen byte coh = cond(onset==2006,1,cond(onset==2010,2,3))

* ragged windows: each state reports over its own subrange
gen int y0 = 1998 + mod(state,4)
gen int y1 = 2014 + mod(state,3)

* state heterogeneity
gen double sfe    = rnormal()*2
gen double strend = rnormal()*0.10
* a STATE-level covariate (varies inside cohort) -- the contrast arm
gen double statechar = rnormal()

expand 20
bysort state: gen int year = 1996 + _n
drop if year < y0 | year > y1

* common year shocks
tempfile pnl
sort year
by year: gen byte yrtag = (_n==1)
gen double yshk = .
replace yshk = rnormal()*0.5 if yrtag==1
bysort year (yrtag): replace yshk = yshk[_N]
drop yrtag

gen byte post = !missing(onset) & year >= onset
gen double y = sfe + strend*(year-2006) + yshk + 1.10*post + rnormal()*0.90

* messy but plausible: a handful of unreported cell-years
replace y = . if state==5 & year==2003
replace y = . if state==2 & year==2011

* --- the COHORT-KEYED covariate: a cohort-level characteristic.
*     constant within coh by construction -> collinear with the i FE.
gen double cohbase = cond(coh==1, 3.2, cond(coh==2, 5.7, 1.4))
* --- an 'adoption wave' grouping, COARSER than coh
gen byte wavegrp = (coh==1)

count
tab coh, m

* ================================================================ arms
local A1 ""
local A2 "controls(cohbase)"
local A3 "fe(wavegrp)"
local A4 "coarsecohortcontrols(c.cohbase)"
local A5 "coarsecohortcontrols(c.statechar)"
local A6 "fe(coh)"

local L1 "BASE"
local L2 "CTRL_cohbase"
local L3 "FE_wavegrp"
local L4 "CCC_cohbase"
local L5 "CCC_statechar_CONTRAST"
local L6 "FE_coh_literal"

forvalues k = 1/6 {
  dis "===== ARM `k' = `L`k'' ====="
  timer clear 1
  timer on 1
  cap noi qui wooldid y coh year onset, att permute(shuffle) permexact ///
      permonsetcohorts cluster(state) `A`k''
  local rc = _rc
  timer off 1
  dis "ARM`k'_RC=" `rc'
  if `rc' == 0 {
    mat B`k' = e(b)
    mat PD`k' = e(permdist_b)
    scalar p`k'  = e(perm_p_att)
    scalar sd`k' = e(perm_sd_att)
    scalar nf`k' = e(perm_failed)
    scalar reb`k'= e(perm_irebuilt)
    scalar npos`k'= e(perm_npossible)
    local cn`k' : colnames B`k'
    dis "ARM`k'_LBL=`L`k''"
    dis "ARM`k'_BHEX=" %21x B`k'[1,1]
    dis "ARM`k'_BDEC=" %20.15f B`k'[1,1]
    dis "ARM`k'_P="    %12.9f p`k'
    dis "ARM`k'_SD="   %14.10f sd`k'
    dis "ARM`k'_FAIL=" nf`k'
    dis "ARM`k'_REB="  reb`k'
    dis "ARM`k'_NPOS=" npos`k'
    dis "ARM`k'_PDR=" rowsof(PD`k') " PDC=" colsof(PD`k')
    local pc : colnames PD`k'
    dis "ARM`k'_PDCOLS=`pc'"
    local bc : colnames B`k'
    dis "ARM`k'_BCOLS=`bc'"
  }
  timer list 1
}

* ------------------------------------------- paired per-assignment check
* find the att column of permdist_b (exact name match, fallback col 1)
local pc : colnames PD1
local j = 0
local jatt = 0
foreach c of local pc {
  local ++j
  if "`c'" == "att" & `jatt' == 0 local jatt = `j'
}
if `jatt' == 0 local jatt = 1
dis "JATT=`jatt'"

forvalues k = 2/6 {
  cap confirm matrix PD`k'
  if _rc == 0 {
    local nr = rowsof(PD1)
    local nr2 = rowsof(PD`k')
    dis "PAIR`k'_ROWS=`nr' vs `nr2'"
    if `nr' == `nr2' {
      local nd = 0
      local maxd = 0
      forvalues r = 1/`nr' {
        local d = abs(PD1[`r',`jatt'] - PD`k'[`r',`jatt'])
        if `d' > 1e-12 local ++nd
        if `d' > `maxd' local maxd = `d'
      }
      dis "PAIR`k'_NDIFF=`nd' of `nr'"
      dis "PAIR`k'_MAXD=" %12.4e `maxd'
      dis "PAIR`k'_ROW1D=" %21x (PD1[1,`jatt'] - PD`k'[1,`jatt'])
      dis "PAIR`k'_ROW2D=" %21x (PD1[2,`jatt'] - PD`k'[2,`jatt'])
    }
  }
}

* =============================== CONTROL ARM: cluster(coh) == i, no rebuild
dis "===== CONTROL ARM cluster(coh)==i (no rebuild expected) ====="
forvalues k = 1/4 {
  cap noi qui wooldid y coh year onset, att permute(shuffle) permexact ///
      permonsetcohorts cluster(coh) `A`k''
  local rc = _rc
  dis "CTLARM`k'_RC=" `rc'
  if `rc' == 0 {
    mat CB`k' = e(b)
    dis "CTLARM`k'_LBL=`L`k''"
    dis "CTLARM`k'_BHEX=" %21x CB`k'[1,1]
    dis "CTLARM`k'_P="  %12.9f e(perm_p_att)
    dis "CTLARM`k'_SD=" %14.10f e(perm_sd_att)
    dis "CTLARM`k'_REB=" e(perm_irebuilt)
    dis "CTLARM`k'_NPOS=" e(perm_npossible)
  }
}

* ================================ OBSERVED FIT ALONE (no permute) inertness
dis "===== OBSERVED-ONLY (no permute) ====="
forvalues k = 1/6 {
  cap noi qui wooldid y coh year onset, att `A`k''
  local rc = _rc
  dis "OBS`k'_RC=" `rc'
  if `rc' == 0 {
    mat OB`k' = e(b)
    dis "OBS`k'_LBL=`L`k''"
    dis "OBS`k'_BHEX=" %21x OB`k'[1,1]
  }
}

* ================================ WHAT IS PRINTED (non-quiet, arms 2 and 4)
dis "===== NONQUIET ARM 2 (controls(cohbase)) ====="
cap noi wooldid y coh year onset, att permute(shuffle) permexact ///
    permonsetcohorts cluster(state) controls(cohbase)
dis "NQ2_RC=" _rc

dis "===== NONQUIET ARM 4 (coarsecohortcontrols(c.cohbase)) ====="
cap noi wooldid y coh year onset, att permute(shuffle) permexact ///
    permonsetcohorts cluster(state) coarsecohortcontrols(c.cohbase)
dis "NQ4_RC=" _rc

dis "PROBEA_DONE"
