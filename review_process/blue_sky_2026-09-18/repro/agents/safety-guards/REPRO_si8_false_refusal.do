*=====================================================================
* REPRO_si8_false_refusal.do  --  FULLY SELF-CONTAINED (no includes,
* no repo data, no test suite).  Demonstrates that wooldid's maxvar
* preflight (wooldid.ado:1355-1377) refuses a run that completes.
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_si8_false_refusal.do
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260916
set sortseed 20260916

* ---- an ordinary wide panel: 4800 units x 8 years, 20 ever treated ----
set obs 4800
gen long unit = _n
gen double ttre = .
replace ttre = 2004 if unit <= 10
replace ttre = 2005 if unit > 10 & unit <= 20
gen double a_i = rnormal()
expand 8
bysort unit: gen int t = 1999 + _n
gen double g_t = 0.1*(t - 2000)
gen byte post = (ttre < . & t >= ttre)
gen double y = a_i + g_t + 1.5*post + 0.3*rnormal()

dis "MAXVAR = " c(maxvar)          // Stata's own default; nothing is set here
qui count
dis "ROWS = " r(N)

dis "---- (1) DEFAULT safety: wooldid refuses ----"
capture noisily wooldid y unit t ttre, att cluster(unit)
dis "RC_DEFAULT = " _rc

dis "---- (2) the SAME call with safety(off): it runs ----"
capture noisily wooldid y unit t ttre, att cluster(unit) safety(off)
dis "RC_SAFETYOFF = " _rc
dis "N = " e(N)
mat list e(wdidmainresults), format(%10.6f)

dis "---- (3) the preflight's level terms ARE real on the reg path ----"
dis "---- (reg/poisson engine, taken under customvce()/unconditionalse) ----"
capture noisily wooldid y unit t ttre, att customvce(cluster unit) safety(off)
dis "RC_REGPATH = " _rc
dis "REPRO_DONE"
