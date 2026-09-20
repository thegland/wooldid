* ENGINE-INVARIANCE INVARIANT.  permute() and bootstrap() replace the INFERENCE, not the
* estimates.  So the OBSERVED point estimates (column 1 of e(wdidmainresults)) must be
* bit-identical to the analytic run's, for every engine and every engine option --
* including permonsetcohorts, which rebuilds i inside replications only.
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_col1
program define wdfz_col1, rclass
  syntax , TAG(string)
  local ra = rowsof(MA)
  local rb = rowsof(MB)
  local na : rownames MA
  local nb : rownames MB
  local ea : roweq MA
  local eb : roweq MB
  local bad ""
  local worst 0
  local nov 0
  forvalues i = 1/`ra' {
    local ni : word `i' of `na'
    local ei : word `i' of `ea'
    forvalues j = 1/`rb' {
      local nj : word `j' of `nb'
      local ej : word `j' of `eb'
      if "`ni'" != "`nj'" | "`ei'" != "`ej'" continue
      local nov = `nov' + 1
      local a = MA[`i',1]
      local b = MB[`j',1]
      if missing(`a') & missing(`b') continue
      if missing(`a') | missing(`b') {
        local bad "`bad' MISS(`ei':`ni')"
        continue
      }
      local d = abs(`a'-`b')/max(1,abs(`a'))
      if `d' > `worst' local worst = `d'
      if `d' > 1e-10 local bad "`bad' D(`ei':`ni':`a'/`b'=`d')"
    }
  }
  dis "ENG|`tag'|overlap=`nov'|worst=`worst'|`bad'"
end

local tagn 0
foreach ds in canon unbal rcs fewG10 clustonset hetero_unit {
  foreach base in "att att_i cluster(cohort) esprelength(2) espostlength(2)" "att cluster(cohort)" {
    use "`D'/ds_`ds'.dta", clear
    ereturn clear
    capture qui wooldid y unit t ttre, `base'
    local rcA = _rc
    if `rcA' == 0 matrix MA = e(wdidmainresults)
    foreach eng in "permute(shuffle) permreps(20)" ///
                   "permute(shuffle) permreps(20) permtstat" ///
                   "permute(shuffle) permreps(20) permanalyticse" ///
                   "permute(shuffle) permreps(20) permtci permabsp" ///
                   "bootstrap(standard) reps(20) rseed(5)" ///
                   "bootstrap(standard) reps(20) rseed(5) boottstat" ///
                   "bootstrap(bayes) reps(20) rseed(5) safety(off)" {
      local tagn = `tagn' + 1
      use "`D'/ds_`ds'.dta", clear
      ereturn clear
      capture qui wooldid y unit t ttre, `base' `eng'
      local rcB = _rc
      if `rcA'==0 & `rcB'==0 {
        matrix MB = e(wdidmainresults)
        wdfz_col1 , tag(`tagn')
      }
      else dis "ENG|`tagn'|RC(`rcA'/`rcB')"
      dis "ENGLEGEND|`tagn'|`ds'|`base'|`eng'"
    }
  }
}

* --- permonsetcohorts: cluster strictly finer than i (the one-day-old T1/T2/T3 surface) ---
foreach ds in rcs {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y cohort t ttre, att att_i cluster(hh) esprelength(2) espostlength(2)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  foreach eng in "permute(shuffle) permreps(20) permonsetcohorts" ///
                 "permute(shuffle) permreps(20) permonsetcohorts permtstat" ///
                 "bootstrap(standard) reps(20) rseed(5)" {
    local tagn = `tagn' + 1
    use "`D'/ds_`ds'.dta", clear
    ereturn clear
    capture qui wooldid y cohort t ttre, att att_i cluster(hh) esprelength(2) espostlength(2) `eng'
    local rcB = _rc
    if `rcA'==0 & `rcB'==0 {
      matrix MB = e(wdidmainresults)
      wdfz_col1 , tag(`tagn')
    }
    else dis "ENG|`tagn'|RC(`rcA'/`rcB')"
    dis "ENGLEGEND|`tagn'|FINE-`ds'|cluster(hh)|`eng'"
  }
}
dis "ENGINE_DONE"
