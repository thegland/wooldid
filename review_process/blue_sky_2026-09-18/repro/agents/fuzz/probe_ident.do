* Internal self-consistency identities that must hold inside / across wooldid runs.
* Compares only the rows that appear in BOTH tables, keyed by equation:name.
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_ovl
program define wdfz_ovl, rclass
  * compare the overlapping (eq:name) rows of MA and MB, columns given in cols()
  syntax , TAG(string) [ COLS(string) ]
  if "`cols'" == "" local cols 1 2 4 5 6
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
      if "`ni'" != "`nj'" continue
      if "`ei'" != "`ej'" continue
      local nov = `nov' + 1
      foreach c of local cols {
        local a = MA[`i',`c']
        local b = MB[`j',`c']
        if missing(`a') & missing(`b') continue
        if missing(`a') | missing(`b') {
          local bad "`bad' MISS(`ei':`ni'c`c')"
          continue
        }
        local d = abs(`a'-`b')/max(1,abs(`a'))
        if `d' > `worst' local worst = `d'
        if `d' > 1e-9 local bad "`bad' D(`ei':`ni'c`c':`a'/`b'=`d')"
      }
    }
  }
  dis "IDENT|`tag'|overlap=`nov'|worst=`worst'|`bad'"
end

* ---------------------------------------------------------------- 1. estimand additivity
foreach ds in canon unbal rcs longT clustonset {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre, att cluster(region) esprelength(2) espostlength(2)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i att_it att_itime cluster(region) esprelength(2) espostlength(2)
  local rcB = _rc
  if `rcB'==0 matrix MB = e(wdidmainresults)
  if `rcA'==0 & `rcB'==0 wdfz_ovl , tag(addestimand-`ds')
  else dis "IDENT|addestimand-`ds'|RCFAIL(`rcA'/`rcB')"
}

* ---------------------------------------------------------------- 2. sgall partition == main
foreach ds in canon unbal miss_y pw_survey {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2) subgroups(sgall)
  local rcS = _rc
  capture matrix MA = e(wdidmainresults)
  capture matrix MB = e(wdidsgresults_1)
  if _rc == 0 & `rcS'==0 wdfz_ovl , tag(sgall-`ds')
  else dis "IDENT|sgall-`ds'|NO_SGMAT(`=_rc') mats=`e(matrices)'"
}

* ---------------------------------------------------------------- 3. wider ES window leaves inner horizons alone
foreach ds in canon longT unbal {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre, att cluster(region) esprelength(1) espostlength(1)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  ereturn clear
  capture qui wooldid y unit t ttre, att cluster(region) esprelength(2) espostlength(3)
  local rcB = _rc
  if `rcB'==0 matrix MB = e(wdidmainresults)
  if `rcA'==0 & `rcB'==0 wdfz_ovl , tag(eswindow-`ds')
  else dis "IDENT|eswindow-`ds'|RCFAIL(`rcA'/`rcB')"
}

* ---------------------------------------------------------------- 4. subgroups() leaves the primary effects alone
foreach ds in canon unbal sg_gaps {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2) subgroups(sg2)
  local rcB = _rc
  if `rcB'==0 matrix MB = e(wdidmainresults)
  if `rcA'==0 & `rcB'==0 wdfz_ovl , tag(sgprimary-`ds')
  else dis "IDENT|sgprimary-`ds'|RCFAIL(`rcA'/`rcB')"
}

* ---------------------------------------------------------------- 5. level() moves only the interval
foreach ds in canon unbal {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2) level(90)
  local rcB = _rc
  if `rcB'==0 matrix MB = e(wdidmainresults)
  if `rcA'==0 & `rcB'==0 wdfz_ovl , tag(level-`ds') cols(1 2 4)
  else dis "IDENT|level-`ds'|RCFAIL(`rcA'/`rcB')"
}

* ---------------------------------------------------------------- 6. diagnostics do not move estimates
foreach ds in canon unbal {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  ereturn clear
  capture qui wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2) jointtests adversarial gteffects summarystats histogramcohorteffects
  local rcB = _rc
  if `rcB'==0 matrix MB = e(wdidmainresults)
  if `rcA'==0 & `rcB'==0 wdfz_ovl , tag(diagflags-`ds')
  else dis "IDENT|diagflags-`ds'|RCFAIL(`rcA'/`rcB')"
}

* ---------------------------------------------------------------- 7. nofastcollapse must not move a number
foreach ds in canon unbal rcs pw_survey miss_y clustonset {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture qui wooldid y unit t ttre [aw=w], att att_i att_it att_itime cluster(region) esprelength(2) espostlength(2)
  local rcA = _rc
  if `rcA'==0 matrix MA = e(wdidmainresults)
  ereturn clear
  capture qui wooldid y unit t ttre [aw=w], att att_i att_it att_itime cluster(region) esprelength(2) espostlength(2) nofastcollapse
  local rcB = _rc
  if `rcB'==0 matrix MB = e(wdidmainresults)
  if `rcA'==0 & `rcB'==0 wdfz_ovl , tag(nofastcollapse-`ds')
  else dis "IDENT|nofastcollapse-`ds'|RCFAIL(`rcA'/`rcB')"
}
dis "IDENT_DONE"
