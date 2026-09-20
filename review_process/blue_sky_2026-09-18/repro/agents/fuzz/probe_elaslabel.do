* Does contreatelasticitytype() relabel the table as elasticities while leaving the
* numbers as plain ATTs when contreat() is absent?
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

dis "############ BINARY MODE (no contreat) ############"
qui wooldid ypos cohort t ttre, att cluster(cohort)
matrix BASE = e(wdidmainresults)
dis "LAB|baseline|" BASE[1,1] "|" BASE[1,2]
foreach ty in dydx dyex eydx eyex {
  ereturn clear
  capture noisily wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(`ty')
  local rc = _rc
  if `rc'==0 {
    matrix M = e(wdidmainresults)
    dis "LAB|`ty'|" M[1,1] "|" M[1,2] "|identical_to_baseline=" (abs(M[1,1]-BASE[1,1])<1e-12)
    dis "LAB_ETYPE|`ty'|e(contreatelasticitytype)=`e(contreatelasticitytype)'"
  }
  else dis "LAB|`ty'|RC=`rc'"
}
dis "############ semielasticity in binary mode ############"
ereturn clear
qui wooldid ypos cohort t ttre, att cluster(cohort) semielasticity
matrix S = e(wdidmainresults)
dis "LAB|semielasticity|" S[1,1] "|" S[1,2] "|identical_to_baseline=" (abs(S[1,1]-BASE[1,1])<1e-12)

dis "############ CONTINUOUS MODE (contreat present) -- the option's real home ############"
ereturn clear
qui wooldid ypos cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic)
matrix C0 = e(wdidmainresults)
dis "LABC|dydx_default|" C0[1,1]
foreach ty in dyex eydx eyex {
  ereturn clear
  capture noisily qui wooldid ypos cohort t ttre, att cluster(cohort) contreat(dose) contreatcontrols(basic) contreatelasticitytype(`ty')
  if _rc==0 {
    matrix C = e(wdidmainresults)
    dis "LABC|`ty'|" C[1,1] "|identical_to_dydx=" (abs(C[1,1]-C0[1,1])<1e-12)
  }
  else dis "LABC|`ty'|RC=`=_rc'"
}
dis "############ the printed header for each binary-mode call ############"
wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(eyex)
wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(eydx)
wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(dyex)
dis "ELASLABEL_DONE"
