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
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t  = 0.10*(t-2000)
gen byte post   = !missing(ttre) & t >= ttre
gen double ypos = 50 + a_i + g_t + 1.5*post + 0.2*rnormal()
dis "#################### BASELINE ####################"
wooldid ypos cohort t ttre, att cluster(cohort)
dis "#################### + semielasticity ####################"
wooldid ypos cohort t ttre, att cluster(cohort) semielasticity
dis "#################### + contreatelasticitytype(eyex) ####################"
wooldid ypos cohort t ttre, att cluster(cohort) contreatelasticitytype(eyex)
dis "SEMIHDR_DONE"
