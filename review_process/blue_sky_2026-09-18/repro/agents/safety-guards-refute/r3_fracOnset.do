* R3 -- SAFETY-03: fractional onset under safety(off).
* Tests: (a) does it reproduce on my own data? (b) is the convention really "round UP",
* or is it a by-construction consequence (test +0.3 and +0.7)? (c) are the horizon
* labels fractional (a loud signal) or integer (silent)? (d) the untested esprelength(0) leg.
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
* DGP: effect switches on AT ttre0 (integer year), NO pre-trend whatsoever
gen double y = 1 + 0.2*(t-2000) + ui + rnormal() + 1.5*(t >= ttre0 & ttre0 < .)
gen double ttre_p5 = ttre0 + 0.5
gen double ttre_p3 = ttre0 + 0.3
gen double ttre_p7 = ttre0 + 0.7
gen double ttre_up = ttre0 + 1
save "r3.dta", replace

cap prog drop shownums
prog def shownums
  args tag
  mat bb = e(b)
  local nm : colnames bb
  local k = colsof(bb)
  forvalues j = 1/`k' {
    local n1 : word `j' of `nm'
    dis "`tag'_COEF`j'=" %21x bb[1,`j'] " name=`n1'"
  }
  cap mat rr = e(wdidmainresults)
  if _rc == 0 {
    local rn : rownames rr
    dis "`tag'_MAINROWS=`rn'"
    local cn : colnames rr
    dis "`tag'_MAINCOLS=`cn'"
    local R = rowsof(rr)
    local relc = 0
    local jj = 0
    foreach c of local cn {
      local jj = `jj' + 1
      if "`c'" == "relyr" local relc = `jj'
    }
    if `relc' > 0 {
      forvalues r = 1/`R' {
        local rname : word `r' of `rn'
        dis "`tag'_ROW`r' name=`rname' est=" %12.7f rr[`r',1] " relyr=" %12.7f rr[`r',`relc']
      }
    }
  }
end

dis "=== L1: correct integer ttre, default safety ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre0, att cluster(unit) esprelength(2) espostlength(2)
dis "L1_RC=" _rc
if _rc==0 shownums L1

dis "=== L2: ttre+0.5, DEFAULT safety (expect refusal) ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_p5, att cluster(unit) esprelength(2) espostlength(2)
dis "L2_RC=" _rc

dis "=== L3: ttre+0.5, safety(off) ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_p5, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "L3_RC=" _rc
if _rc==0 shownums L3

dis "=== L4: honest ttre+1 recode, default safety (the round-UP convention) ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_up, att cluster(unit) esprelength(2) espostlength(2)
dis "L4_RC=" _rc
if _rc==0 shownums L4

dis "=== L5: ttre+0.3, safety(off) -- round-to-NEAREST would give the ttre+0 answer ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_p3, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "L5_RC=" _rc
if _rc==0 shownums L5

dis "=== L6: ttre+0.7, safety(off) ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_p7, att cluster(unit) esprelength(2) espostlength(2) safety(off)
dis "L6_RC=" _rc
if _rc==0 shownums L6

dis "=== L7: ttre+0.5, safety(off), DEFAULT esprelength(0) (their untested caveat) ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_p5, att cluster(unit) safety(off)
dis "L7_RC=" _rc
if _rc==0 shownums L7

dis "=== L8: ttre+0.5, safety(off), no ES options at all, att only ==="
use "r3.dta", clear
cap noi wooldid y unit t ttre_p5, att att_i cluster(unit) safety(off) espostlength(2)
dis "L8_RC=" _rc
if _rc==0 shownums L8

dis "R3_DONE"
