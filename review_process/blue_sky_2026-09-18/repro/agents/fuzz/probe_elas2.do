* Correct comparison: under contreat() the elasticity type governs the ame rows, not att.
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 60
gen long unit   = _n
gen long cohort = 1 + mod(_n, 6)
gen byte region = 1 + mod(_n, 5)
gen double a_i  = 0.3*rnormal()
gen ttre = .
replace ttre = 2003 if inlist(cohort,1,2)
replace ttre = 2005 if inlist(cohort,3,4)
gen double dose_i = runiform()*3
replace dose_i = 0 if runiform() < .4
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t  = 0.10*(t-2000)
gen byte post   = !missing(ttre) & t >= ttre
gen double dose = cond(post, dose_i, 0)
gen double ypos = 50 + a_i + g_t + 1.5*post + 0.8*dose + 0.2*rnormal()

capture program drop showame
program define showame
  syntax , TAG(string)
  matrix M = e(wdidmainresults)
  local rn : rownames M
  local re : roweq M
  local R = rowsof(M)
  forvalues r = 1/`R' {
    local nn : word `r' of `rn'
    local ee : word `r' of `re'
    if "`ee'" == "main" dis "AME|`tag'|`ee':`nn'|" M[`r',1] "|" M[`r',2]
  }
end
foreach ty in DEFAULT dydx dyex eydx eyex {
  local o ""
  if "`ty'" != "DEFAULT" local o "contreatelasticitytype(`ty')"
  ereturn clear
  capture qui wooldid ypos cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) `o'
  if _rc==0 showame , tag(`ty')
  else dis "AME|`ty'|RC=`=_rc'"
}
dis "############ binary-mode headers ############"
foreach ty in dydx dyex eydx eyex {
  dis "@@HDR|`ty'"
  capture noisily wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(`ty')
}
dis "ELAS2_DONE"
