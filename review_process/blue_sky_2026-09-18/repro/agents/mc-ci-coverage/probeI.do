clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* dump the actual permute e() names and the printed interval, at a NON-degenerate
* cluster count: i = unit, cluster(unit), 40 units, constant treatment effect.
clear
set seed 55501
set sortseed 55501
set obs 40
gen long unit = _n
gen byte cg = mod(_n-1,4)+1
gen double ai = rnormal(0,1)
expand 8
bys unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 if cg==1
replace ttre = 2004 if cg==2
replace ttre = 2005 if cg==3
gen double gt = 0.4*(t-1999)
gen double truetau = cond(ttre<. & t>=ttre, 1, 0)
gen double base = ai + gt + truetau
bys unit: gen double ui = rnormal(0,0.4) if _n==1
by unit: replace ui = ui[1]
gen double y = base + ui + rnormal(0,0.6)
save "$WD/fixI.dta", replace

foreach L in 90 95 {
  cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(99) ///
       rseed(2024) permtstat permtci level(`L') dots(0)
  dis "RC|I_`L'|" _rc
  if _rc==0 {
    dis "=== e(scalars) for level `L' ==="
    foreach s in `: e(scalars)' {
      dis "ESC|`L'|`s'|" e(`s')
    }
    dis "=== e(macros) ==="
    dis "EMAC|`L'|" "`: e(macros)'"
    dis "=== e(matrices) ==="
    dis "EMAT|`L'|" "`: e(matrices)'"
    mat list e(wdidmainresults), format(%14.9f)
    cap confirm matrix e(perm_ci)
    if _rc==0 mat list e(perm_ci), format(%14.9f)
  }
}
* plain analytic on the same data, for scale
cap noi wooldid y unit t ttre, att cluster(unit) level(95)
if _rc==0 {
  mat list e(wdidmainresults), format(%14.9f)
}
dis "PROBEI_DONE"

* ---- does esprelength()/espostlength() move the att point estimate or its SE? ----
use "$WD/fixI.dta", clear
cap noi wooldid y cg t ttre, att att_i cluster(unit) level(95)
if _rc==0 {
  mat Z = e(wdidmainresults)
  dis "ESOPT|none|att b=" %21x Z[1,1] "|se=" %21x Z[1,2] "|N=" e(N) "|df=" e(df_r)
  dis "ESOPT|none|att b=" %14.9f Z[1,1] "|se=" %14.9f Z[1,2]
}
foreach k in 1 2 3 {
  cap noi wooldid y cg t ttre, att att_i cluster(unit) esprelength(`k') espostlength(`k') level(95)
  if _rc==0 {
    mat Z = e(wdidmainresults)
    dis "ESOPT|`k'|att b=" %21x Z[1,1] "|se=" %21x Z[1,2] "|N=" e(N) "|df=" e(df_r)
    dis "ESOPT|`k'|att b=" %14.9f Z[1,1] "|se=" %14.9f Z[1,2]
  }
}
dis "PROBEI2_DONE"
