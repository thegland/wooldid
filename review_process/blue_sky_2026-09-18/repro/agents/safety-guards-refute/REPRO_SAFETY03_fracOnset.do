*==============================================================================
* SELF-CONTAINED REPRO -- SAFETY-03 (fractional onset under safety(off))
* CLAIM UNDER TEST: "wooldid SILENTLY adopts the ROUND-UP convention".
* What this file shows: the numbers do coincide with an honest ttre+1 recode,
* but there is NO rounding and nothing is silent -- the event-study rows are
* LABELLED with the exact fractional horizons (pre1.5 / pre.5 / post.5, and
* pre1.7 / pre.7 / post.3 for a +0.7 offset), in the printed table, in
* e(b_es_att)'s column names and in the relyr column of e(wdidmainresults).
*==============================================================================
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
gen double ttre0 = .
replace ttre0 = 2003 + mod(unit,3) if tg==1
gen double ui = rnormal()
bysort unit (t): replace ui = ui[1]
* DGP: effect switches on AT the integer onset year; pre-trend exactly zero
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t >= ttre0 & ttre0 < .)
gen double ttre_p5 = ttre0 + 0.5
gen double ttre_p3 = ttre0 + 0.3
gen double ttre_p7 = ttre0 + 0.7
gen double ttre_up = ttre0 + 1
save "_repro_s03.dta", replace

cap prog drop shownums
prog def shownums
  args tag
  mat bb = e(b)
  local nm : colnames bb
  forvalues j = 1/`=colsof(bb)' {
    local n1 : word `j' of `nm'
    dis "`tag'_COEF`j'=" %21x bb[1,`j'] " name=`n1'"
  }
end

dis "=== L1 correct integer ttre, default safety (the truth) ==="
cap noi wooldid y unit t ttre0, att cluster(unit) esprelength(2) espostlength(2)
dis "L1_RC=" _rc
if _rc==0 shownums L1

dis "=== L2 ttre+0.5, DEFAULT safety: EXPECT the named r(459) ==="
use "_repro_s03.dta", clear
cap noi wooldid y unit t ttre_p5, att cluster(unit) esprelength(2) espostlength(2)
dis "L2_RC=" _rc

dis "=== L3 ttre+0.5, safety(off): note the pre1.5 / pre.5 / post.5 ROW LABELS ==="
use "_repro_s03.dta", clear
cap noi wooldid y unit t ttre_p5, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "L3_RC=" _rc
if _rc==0 shownums L3

dis "=== L4 honest ttre+1 recode, default safety ==="
use "_repro_s03.dta", clear
cap noi wooldid y unit t ttre_up, att cluster(unit) esprelength(2) espostlength(2)
dis "L4_RC=" _rc
if _rc==0 shownums L4

dis "=== L5 ttre+0.7, safety(off): labels are pre1.7 / pre.7 / post.3 -- NO rounding ==="
use "_repro_s03.dta", clear
cap noi wooldid y unit t ttre_p7, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "L5_RC=" _rc
if _rc==0 shownums L5

dis "=== L6 ttre+0.3, safety(off): dies r(452) -- the bypass is a float lottery ==="
use "_repro_s03.dta", clear
cap noi wooldid y unit t ttre_p3, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "L6_RC=" _rc

dis "=== L7 ttre+0.5, safety(off), DEFAULT esprelength(0): dies r(452) ==="
use "_repro_s03.dta", clear
cap noi wooldid y unit t ttre_p5, att cluster(unit) safety(off)
dis "L7_RC=" _rc
dis "REPRO_SAFETY03_DONE"
