*==============================================================================
* SKEPTIC PROBE D -- the untreatedZero "ringer" mechanism, and two side probes.
*
* D1  MECHANISM, wooldid-free.  Under contreatcontrols(untreatedZero) the
*     explicit-dummy engines fit an EXTRA pooled dose term (the DF-1 "ringer",
*     ado:1233, which fires on basic | unconditionalse | customvce | poisson and
*     which the untreatedZero branch at ado:1313 never clears).  margins switches
*     off only targetflag, so the pooled dose slope stays in BOTH predictions and
*     is subtracted out of the reported effect.
* D2  Same trigger on the DEFAULT engine via poisson (which also satisfies the
*     ado:1233 condition, and whose ringer lands inside absorb()).
* D3  semielasticity escalation: strictly positive y, strictly positive true
*     fitted values, but reghdfe's _cons+Xb goes negative.
* D4  customweightmultiplier() x att_it arithmetic check.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do pD_ringer.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 90210
set sortseed 90210
which wooldid

*==============================================================================
* D1 -- wooldid-free mechanism
*==============================================================================
clear
set obs 9
gen int i = _n
gen int ttre = .
replace ttre = 3 if i <= 3
replace ttre = 5 if i > 3 & i <= 6
gen double ui_unit = 100 + 20*runiform()
expand 7
bysort i: gen int t = _n
expand 14
gen double dose = 0.2 + 2.5*runiform()
gen byte post = (t >= ttre) & ttre < .
gen double d2 = dose*post
gen double y = ui_unit + 0.5*t + 3.0*d2 + rnormal()*0.3
gen byte TF = (ttre < .)
gen int relyr = t - ttre
egen long cell = group(i relyr) if TF==1 & relyr != -1
replace cell = 0 if missing(cell)
qui sum d2 if post==1
local truthatt = 3*r(mean)
dis "D1 TRUTH att = " %10.6f `truthatt' "   TRUTH ame = 3.000000"
tempfile d1
save `d1'

* --- model A: the model untreatedZero says it fits (no g() at all)
qui reg y i.i i.t i.cell i.cell#c.d2
predict double pA1, xb
preserve
  replace cell = 0
  predict double pA0, xb
  gen double effA = pA1 - pA0
  qui sum effA if post==1
  dis "D1 model A (no pooled dose term)  mean effect on treated-post = " %10.6f r(mean)
restore

* --- model B: what the explicit-dummy path actually fits (pooled ringer first)
use `d1', clear
qui reg y i.i i.t c.d2 i.cell i.cell#c.d2
predict double pB1, xb
preserve
  replace cell = 0
  predict double pB0, xb
  gen double effB = pB1 - pB0
  qui sum effB if post==1
  dis "D1 model B (pooled dose ringer)   mean effect on treated-post = " %10.6f r(mean)
  dis "D1 pooled dose coefficient c.d2   = " %10.6f _b[d2]
restore

*==============================================================================
* D2 -- poisson also trips ado:1233; its ringer goes inside absorb()
*==============================================================================
use `d1', clear
gen double mu = exp(1.0 + 0.10*t + 0.40*d2)
gen int ycount = rpoisson(mu)
dis "D2 TRUTH: dydx(xb) ame = 0.400000 ; att = 0.40*E[dose|post]"
qui sum d2 if post==1
dis "D2 TRUTH att (xb scale) = " %10.6f 0.40*r(mean)
tempfile d2f
save `d2f'

dis _n "===== D2a poisson + untreatedZero (DEFAULT engine) ====="
cap noi wooldid ycount i t ttre, att contreat(d2) contreatcontrols(untreatedZero) poisson verbose
dis "D2a_RC=" _rc
if _rc==0 dis "D2a att=" %13.7f e(b)[1,1] " ame=" %13.7f e(b)[1,2]

use `d2f', clear
dis _n "===== D2b linear + untreatedZero (DEFAULT engine, reference) ====="
cap noi wooldid ycount i t ttre, att contreat(d2) contreatcontrols(untreatedZero)
dis "D2b_RC=" _rc
if _rc==0 dis "D2b att=" %13.7f e(b)[1,1] " ame=" %13.7f e(b)[1,2]

*==============================================================================
* D3 -- semielasticity: strictly positive y, reghdfe xb negative
*==============================================================================
clear
set obs 40
gen int i = _n
gen int ttre = .
replace ttre = 6 if i <= 16
gen double alpha = cond(i<=16, 26 + 0.5*i, 6 + 0.1*(i-16))
expand 10
bysort i: gen int t = _n
gen byte post = (t >= ttre) & ttre < .
gen double y = alpha - 20*post + rnormal()*0.3
qui sum y
dis "D3 min(y)=" %10.5f r(min) "  mean(y)=" %10.5f r(mean)
qui count if y <= 0
dis "D3 n(y<=0) = " r(N)
tempfile d3
save `d3'

qui reghdfe y post, absorb(i.i i.t)
predict double xbh, xb
qui sum xbh
dis "D3 reghdfe xb: min=" %10.5f r(min) " max=" %10.5f r(max)
qui count if xbh <= 0
dis "D3 n(reghdfe xb <= 0) = " r(N)
qui reg y i.i i.t post
predict double xbf, xb
qui sum xbf
dis "D3 explicit xb: min=" %10.5f r(min) " max=" %10.5f r(max)
qui count if xbf <= 0
dis "D3 n(explicit xb <= 0) = " r(N)

use `d3', clear
dis _n "===== D3a semielasticity, default reghdfe ====="
cap noi wooldid y i t ttre, att semielasticity
dis "D3a_RC=" _rc
if _rc==0 dis "D3a semiel=" %14.9f e(b)[1,1]

use `d3', clear
dis _n "===== D3b semielasticity, customvce(robust) ====="
cap noi wooldid y i t ttre, att semielasticity customvce(robust)
dis "D3b_RC=" _rc
if _rc==0 dis "D3b semiel=" %14.9f e(b)[1,1]

*==============================================================================
* D4 -- customweightmultiplier() x att_it on a one-row-per-cell panel
*       help (sthlp:222-227) calls att_it "the simple, unweighted mean" of the
*       cell effects and spells the weight as w/W_it with no cwm factor.
*==============================================================================
clear
set obs 30
gen int i = _n
gen int ttre = .
replace ttre = 4 if i <= 10
replace ttre = 7 if i > 10 & i <= 20
gen double alpha = 40 + i
gen double cwm = 1 + mod(i,5)              // unit-level multiplier
expand 9
bysort i: gen int t = _n
gen byte post = (t >= ttre) & ttre < .
gen double effi = 1 + 0.4*i                // heterogeneous effect by unit
gen double y = alpha + 0.3*t + effi*post + rnormal()*0.2
tempfile d4
save `d4'

dis _n "===== D4a att_it, no cwm ====="
cap noi wooldid y i t ttre, att_it
dis "D4a_RC=" _rc
if _rc==0 dis "D4a att_it=" %14.9f e(b)[1,1]

use `d4', clear
dis _n "===== D4b att_it, with cwm ====="
cap noi wooldid y i t ttre, att_it customweightmultiplier(cwm)
dis "D4b_RC=" _rc
if _rc==0 dis "D4b att_it=" %14.9f e(b)[1,1]

dis _n "PROBE_D_DONE"
