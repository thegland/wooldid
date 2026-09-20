clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
dis "=== plain Stata set level ==="
foreach L in 99.9 99.99 97.5 90.25 99.5 {
  capture noisily set level `L'
  dis "SETLEVEL|`L'|rc=" _rc "|c(level)=" c(level)
}
set level 95
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
gen double ypos = 50 + a_i + g_t + 1.5*post + 0.2*rnormal()
dis "=== binary-mode baseline ==="
qui wooldid ypos cohort t ttre, att cluster(cohort)
matrix A = e(wdidmainresults)
dis "BASE_att = " A[1,1] "  se=" A[1,2]
dis "=== binary + semielasticity (NO contreat) ==="
capture noisily qui wooldid ypos cohort t ttre, att cluster(cohort) semielasticity
dis "SEMI_RC=" _rc
if _rc==0 {
  matrix Bm = e(wdidmainresults)
  dis "SEMI_att = " Bm[1,1] "  se=" Bm[1,2]
  dis "RATIO = " Bm[1,1]/A[1,1]
}
dis "=== binary + contreatelasticitytype(eyex) (NO contreat) ==="
capture noisily qui wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(eyex)
dis "EYEX_RC=" _rc
if _rc==0 {
  matrix C = e(wdidmainresults)
  dis "EYEX_att = " C[1,1] "  se=" C[1,2]
}
dis "=== level(99.9) through wooldid ==="
capture noisily qui wooldid ypos cohort t ttre, att cluster(cohort) level(99.9)
dis "L999_RC=" _rc
dis "=== level(99.5) through wooldid ==="
capture noisily qui wooldid ypos cohort t ttre, att cluster(cohort) level(99.5)
dis "L995_RC=" _rc
if _rc==0 dis "colname5=" abbrev("`: word 5 of `: colnames e(wdidmainresults)''",32)
dis "FINAL_DONE"
