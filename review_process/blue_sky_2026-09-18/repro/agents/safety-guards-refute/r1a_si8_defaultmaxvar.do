* R1a -- SAFETY-01 refutation attempt, DEFAULT maxvar, own data.
* Claim under test: the si8 maxvar preflight refuses (r907) a run that,
* with safety(off), completes correctly on the SAME maxvar.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

* ---- HARNESS ASSERTION ----
which wooldid
dis "EDITION=" c(edition)
dis "MAXVAR_DEFAULT=" c(maxvar)

* ---- own data: 4800 units (county-scale), 8 years, 20 ever-treated ----
clear
set obs 4800
gen long unit = _n
gen byte tg = (unit <= 20)
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 + mod(unit,3) if tg == 1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
qui count
dis "NROWS=" r(N)
qui levelsof unit if tg==1, local(tu)
dis "NTREATED=" wordcount("`tu'")
save "r1a.dta", replace

* ---- LEG A: default safety ----
dis "=== LEG A: default safety, cluster(unit) ==="
cap noi wooldid y unit t ttre, att cluster(unit)
dis "LEGA_RC=" _rc

* ---- LEG B: safety(off), same maxvar, same call ----
dis "=== LEG B: safety(off), cluster(unit) ==="
cap noi wooldid y unit t ttre, att cluster(unit) safety(off)
dis "LEGB_RC=" _rc
if _rc == 0 {
  mat bB = e(b)
  mat VB = e(V)
  dis "LEGB_N=" e(N)
  dis "LEGB_ATT_HEX=" %21x bB[1,1]
  dis "LEGB_SE_HEX="  %21x sqrt(VB[1,1])
  scalar bBv = bB[1,1]
  scalar sBv = sqrt(VB[1,1])
}
dis "R1A_DONE"
