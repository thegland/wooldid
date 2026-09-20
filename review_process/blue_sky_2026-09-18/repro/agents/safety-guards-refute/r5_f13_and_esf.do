* R5 -- (G) does contreatwithin relax the F13 onset-constancy screen, as it does the
*           i-positivity and float screens?  The help says F13 survives.  Measure it,
*           and measure the damage when safety(off) does bypass it.
*       (H) characterise the esfixedbaseperiod degradation under safety(off).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2003 + mod(unit,3) if tg==1
gen double dose = runiform()*2
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
* the classic merge artifact: onset recorded only from the onset period onward
gen double ttre_brk = ttre if t >= ttre & ttre < .
qui count if missing(ttre_brk) & !missing(ttre)
dis "G_BROKEN_ROWS=" r(N)
save "r5.dta", replace

dis "BLOCKG0_START baseline: CORRECT ttre, default safety"
cap noi wooldid y unit t ttre, att cluster(unit)
dis "G0_RC=" _rc
if _rc==0 {
  mat g0=e(b)
  dis "G0_N=" e(N)
  dis "G0_ATT_HEX=" %21x g0[1,1]
  dis "G0_ATT=" g0[1,1]
}
dis "BLOCKG0_END"

dis "BLOCKG1_START broken ttre, default safety (expect the named F13 refusal)"
use "r5.dta", clear
cap noi wooldid y unit t ttre_brk, att cluster(unit)
dis "G1_RC=" _rc
dis "BLOCKG1_END"

dis "BLOCKG2_START broken ttre, safety(off) -- measure the damage"
use "r5.dta", clear
cap noi wooldid y unit t ttre_brk, att cluster(unit) safety(off)
dis "G2_RC=" _rc
if _rc==0 {
  mat g2=e(b)
  dis "G2_N=" e(N)
  dis "G2_ATT_HEX=" %21x g2[1,1]
  dis "G2_ATT=" g2[1,1]
}
dis "BLOCKG2_END"

dis "BLOCKG3_START broken ttre, contreatwithin+contreat, NO safety() typed"
use "r5.dta", clear
cap noi wooldid y unit t ttre_brk, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "G3_RC=" _rc
if _rc==0 {
  mat g3=e(b)
  dis "G3_N=" e(N)
  dis "G3_ATT_HEX=" %21x g3[1,1]
}
dis "BLOCKG3_END"

dis "BLOCKG4_START control: CORRECT ttre, contreatwithin+contreat"
use "r5.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) contreat(dose) contreatcontrols(basic) contreatwithin
dis "G4_RC=" _rc
if _rc==0 {
  mat g4=e(b)
  dis "G4_N=" e(N)
  dis "G4_ATT_HEX=" %21x g4[1,1]
}
dis "BLOCKG4_END"

*=====================================================================
* BLOCK H -- esfixedbaseperiod pre/post refusal under safety(off)
*=====================================================================
clear
set obs 60
gen long unit = _n
gen byte tg = mod(unit,2)==1
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double ttre = .
replace ttre = 2000 if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
save "r5h.dta", replace

dis "BLOCKH1_START esf default safety, no ES window options"
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod
dis "H1_RC=" _rc
dis "BLOCKH1_END"

dis "BLOCKH2_START esf safety(off), no ES window options"
use "r5h.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod safety(off)
dis "H2_RC=" _rc
dis "BLOCKH2_END"

dis "BLOCKH3_START esf safety(off), esprelength(1) espostlength(1)"
use "r5h.dta", clear
cap noi wooldid y unit t ttre, att cluster(unit) esfixedbaseperiod esprelength(1) espostlength(1) safety(off)
dis "H3_RC=" _rc
dis "BLOCKH3_END"

dis "BLOCKH4_START esf safety(off), esprelength(2) espostlength(2) att_i too"
use "r5h.dta", clear
cap noi wooldid y unit t ttre, att att_i cluster(unit) esfixedbaseperiod esprelength(2) espostlength(2) safety(off)
dis "H4_RC=" _rc
dis "BLOCKH4_END"

dis "R5_DONE"
