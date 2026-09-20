* Permutation and bootstrap p-values must sit on their own discrete grid.
*   default permute two-sided add-one : p = 2*(m+1)/(reps+1)  -> p*(reps+1)/2 integral
*   permrawp (plain share)            : p = c/reps           -> p*reps integral
*   bootstrap percentile coefficient  : plain share          -> p*reps integral
* A p off the grid means the count, the denominator or the doubling is not what is claimed.
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

capture program drop wdfz_grid
program define wdfz_grid, rclass
  syntax , TAG(string) DEN(real) [ HALF(int 0) ]
  matrix M = e(wdidmainresults)
  local R = rowsof(M)
  local rn : rownames M
  local re : roweq M
  local off 0
  local n 0
  local det ""
  forvalues r = 1/`R' {
    local p = M[`r',4]
    if missing(`p') continue
    local n = `n' + 1
    local u = `p' * `den'
    if `half' local u = `u' / 2
    local nearest = round(`u')
    if abs(`u' - `nearest') > 1e-6 {
      local off = `off' + 1
      local nn : word `r' of `rn'
      local ee : word `r' of `re'
      if `off' <= 5 local det "`det' (`ee':`nn' p=`p' units=`u')"
    }
  }
  dis "PGRID|`tag'|den=`den'|half=`half'|rows=`n'|offgrid=`off'|`det'"
end

local RP 40
foreach ds in canon unbal fewG10 {
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture noisily qui wooldid y unit t ttre, att att_i cluster(cohort) permute(shuffle) permreps(`RP') esprelength(2) espostlength(2)
  if _rc==0 wdfz_grid , tag(perm-default-`ds') den(`=`RP'+1') half(1)
  else dis "PGRID|perm-default-`ds'|RC=`=_rc'"

  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture noisily qui wooldid y unit t ttre, att att_i cluster(cohort) permute(shuffle) permreps(`RP') permrawp esprelength(2) espostlength(2)
  if _rc==0 wdfz_grid , tag(perm-rawp-`ds') den(`RP')
  else dis "PGRID|perm-rawp-`ds'|RC=`=_rc'"

  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture noisily qui wooldid y unit t ttre, att att_i cluster(cohort) permute(shuffle) permreps(`RP') permabsp esprelength(2) espostlength(2)
  if _rc==0 wdfz_grid , tag(perm-absp-`ds') den(`=`RP'+1')
  else dis "PGRID|perm-absp-`ds'|RC=`=_rc'"

  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture noisily qui wooldid y unit t ttre, att att_i cluster(cohort) permute(shuffle) permreps(`RP') permside(right) esprelength(2) espostlength(2)
  if _rc==0 wdfz_grid , tag(perm-right-`ds') den(`=`RP'+1')
  else dis "PGRID|perm-right-`ds'|RC=`=_rc'"

  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture noisily qui wooldid y unit t ttre, att att_i cluster(cohort) bootstrap(standard) reps(`RP') rseed(11) esprelength(2) espostlength(2)
  if _rc==0 wdfz_grid , tag(boot-perc-`ds') den(`RP')
  else dis "PGRID|boot-perc-`ds'|RC=`=_rc'"

  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  capture noisily qui wooldid y unit t ttre, att att_i cluster(cohort) bootstrap(standard) reps(`RP') rseed(11) boottstat esprelength(2) espostlength(2)
  if _rc==0 wdfz_grid , tag(boot-tstat-`ds') den(`=`RP'+1') half(1)
  else dis "PGRID|boot-tstat-`ds'|RC=`=_rc'"
}
dis "PGRID_DONE"
