* p4 -- does fastnose disagree with the analytic run about WHICH ROWS are estimable?
*  (a) binary subgroups x ES: the SI-7 estimability drop.  BS-26 carved the binary
*      sg-by-horizon margins call out of `nose' (ado:1981) so the marker survives.
*  (b) the CONTINUOUS-treatment twin of that call (ado:3290) still threads `marginsdof',
*      which IS nose under fastnose -- so the SI-7 filter should be blind there.
*  (c) wooldid_fcdeg reads r(V); under nose K==0 and the screen returns at once, so the
*      cell-collapse fallback (e(wdid_fcguard)) cannot fire in fast mode.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fastnose"
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gN 90
local gT 11
local gonsets "2003 2006 2009"
local gysd 0.30
local gdose unitconst
include "`GEN'/gen_core.do"

* a REALISTIC ragged subgroup: a policy variant that only the 2009 cohort ever adopted.
* every horizon outside that cohort's observed window is non-estimable for level 2.
gen byte sgvar = 1
qui replace sgvar = 2 if ttre == 2009
* and a second, rarer one: a handful of treated units in a single cohort
gen byte sgrare = 1
qui replace sgrare = 2 if ttre == 2003 & mod(unit,17) == 0
tab sgvar, m
tab sgrare, m
qui save "`WD'/p4panel.dta", replace

local sA "y unit t ttre , subgroups(sgvar) esprelength(4) espostlength(5)"
local sB "y unit t ttre , subgroups(sgrare) esprelength(4) espostlength(5)"
local sC "y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) subgroups(sgvar) esprelength(4) espostlength(5)"
local sD "y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) subgroups(sgrare) esprelength(4) espostlength(5)"

foreach arm in A B C D {
  foreach mode in AN FN {
    local ex
    if "`mode'" == "FN" local ex fastnose
    dis "########## ARM `arm' `mode'"
    use "`WD'/p4panel.dta", clear
    ereturn clear
    cap noi wooldid `s`arm'' `ex'
    dis "ARM`arm'_`mode'_RC=" _rc
    if _rc == 0 {
      dis "ARM`arm'_`mode'_fcguard=" e(wdid_fcguard)
      dis "ARM`arm'_`mode'_sglist=" e(wdid_sglist)
      foreach sg in 1 2 {
        cap confirm matrix e(wdidsgresults_`sg')
        if _rc == 0 {
          local rn : rowfullnames e(wdidsgresults_`sg')
          dis "ARM`arm'_`mode'_SGROWS_`sg': `rn'"
        }
        foreach fam in att ame att_it ame_it {
          cap confirm matrix e(b_es_`fam'_s`sg')
          if _rc == 0 {
            local cn : colnames e(b_es_`fam'_s`sg')
            dis "ARM`arm'_`mode'_BES_`fam'_s`sg': `cn'"
            mat list e(b_es_`fam'_s`sg'), format(%12.0g)
          }
        }
      }
      cap confirm matrix e(escomp_s2)
      if _rc == 0 {
        dis "ARM`arm'_`mode'_ESCOMP_S2:"
        mat list e(escomp_s2), format(%12.0g)
      }
    }
  }
}
dis "DONE_P4"
