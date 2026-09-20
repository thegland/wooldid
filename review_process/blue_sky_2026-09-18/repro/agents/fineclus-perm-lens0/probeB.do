* ======================================================================
* LENS-0 sweep: is the reference-distribution move SYSTEMATIC, and in
* which direction, per covariate slot?  Sharp null, assignment redrawn.
* 6 states (2/2/2) -> M = 6!/(2!2!2!) = 90 enumerated assignments.
* ======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771031
set sortseed 771031

which wooldid

local NDS = 12

forvalues d = 1/`NDS' {
  clear
  set obs 6
  gen int state = _n
  * redraw WHICH states get which onset in every dataset
  gen double u = runiform()
  sort u
  gen int slot = _n
  gen double onset = .
  replace onset = 2006 if slot <= 2
  replace onset = 2010 if slot > 2 & slot <= 4
  gen byte coh = cond(onset==2006,1,cond(onset==2010,2,3))
  sort state

  gen int y0 = 1999 + mod(state,3)
  gen int y1 = 2013 + mod(state,2)
  gen double sfe    = rnormal()*2
  gen double strend = rnormal()*0.12
  gen double statechar = rnormal()

  expand 18
  bysort state: gen int year = 1997 + _n
  drop if year < y0 | year > y1

  sort year
  by year: gen byte yrtag = (_n==1)
  gen double yshk = .
  replace yshk = rnormal()*0.5 if yrtag==1
  bysort year (yrtag): replace yshk = yshk[_N]
  drop yrtag

  * SHARP NULL: onset enters the outcome nowhere
  gen double y = sfe + strend*(year-2006) + yshk + rnormal()*0.90

  gen double cohbase = cond(coh==1, 3.2, cond(coh==2, 5.7, 1.4))
  gen byte wavegrp = (coh==1)

  local A1 ""
  local A2 "controls(cohbase)"
  local A3 "fe(wavegrp)"
  local A4 "coarsecohortcontrols(c.cohbase)"

  forvalues k = 1/4 {
    cap noi qui wooldid y coh year onset, att permute(shuffle) permexact ///
        permonsetcohorts cluster(state) `A`k''
    local rc = _rc
    if `rc' == 0 {
      dis "D`d'_A`k'_OK b=" %21x e(b)[1,1] " p=" %11.8f e(perm_p_att) ///
          " sd=" %13.9f e(perm_sd_att) " fail=" e(perm_failed) " reb=" e(perm_irebuilt)
    }
    else dis "D`d'_A`k'_RC=`rc'"
  }
}

dis "PROBEB_DONE"
