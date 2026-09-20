* R1b -- SAFETY-01: (i) can this Stata raise maxvar above the default at all?
*                   (ii) is the safety(off) answer at a TIGHT maxvar identical to the
*                        default-safety answer with headroom?
*                   (iii) engine split: is the demand real on the reg path?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/safety-guards-refute"

dis "EDITION=" c(edition)
dis "MAXVAR_START=" c(maxvar)
cap set maxvar 20000
dis "SETMAXVAR20000_RC=" _rc
dis "MAXVAR_NOW=" c(maxvar)

* ---- (ii) headroom control on the SAME 4800-unit data as R1a ----
if c(maxvar) > 5000 {
  use "r1a.dta", clear
  dis "=== CONTROL: default safety, maxvar " c(maxvar) " ==="
  cap noi wooldid y unit t ttre, att cluster(unit)
  dis "CTRL_RC=" _rc
  if _rc == 0 {
    mat bC = e(b)
    mat VC = e(V)
    dis "CTRL_N=" e(N)
    dis "CTRL_ATT_HEX=" %21x bC[1,1]
    dis "CTRL_SE_HEX="  %21x sqrt(VC[1,1])
    dis "R1A_LEGB_ATT_HEX=+1.86f429091af9eX+000"
    dis "R1A_LEGB_SE_HEX=+1.3c63b9e9bc101X-007"
  }
}

* ---- (iii) engine split at a tight ceiling ----
clear
cap set maxvar 2048
dis "SETMAXVAR2048_RC=" _rc
dis "MAXVAR_TIGHT=" c(maxvar)
set obs 1900
gen long unit = _n
gen byte tg = (unit <= 15)
expand 8
bysort unit: gen int t = 1999 + _n
gen double ttre = .
replace ttre = 2002 + mod(unit,3) if tg == 1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
qui count
dis "T_NROWS=" r(N)

dis "=== T-LEG A: default safety ==="
cap noi wooldid y unit t ttre, att cluster(unit)
dis "TLEGA_RC=" _rc

dis "=== T-LEG B: safety(off), reghdfe path ==="
cap noi wooldid y unit t ttre, att cluster(unit) safety(off)
dis "TLEGB_RC=" _rc
if _rc == 0 {
  mat bTB = e(b)
  dis "TLEGB_N=" e(N)
  dis "TLEGB_ATT_HEX=" %21x bTB[1,1]
}

dis "=== T-LEG C: safety(off), customvce -> reg path (demand should be REAL) ==="
cap noi wooldid y unit t ttre, att customvce(cluster unit) safety(off)
dis "TLEGC_RC=" _rc

dis "=== T-LEG D: safety(off) + fe(unit) + coarse i (the guard's own advice) ==="
gen long cohortid = cond(tg==1, 1 + mod(unit,3), 999)
cap noi wooldid y cohortid t ttre, att cluster(unit) fe(unit)
dis "TLEGD_RC=" _rc
if _rc == 0 dis "TLEGD_N=" e(N)

dis "R1B_DONE"
