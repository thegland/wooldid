*=====================================================================
* repro_cet_mislabel.do -- SELF-CONTAINED.
* Claim: contreatelasticitytype() given WITHOUT contreat() is not inert.
* It leaves every number unchanged but RELABELS the results table as
* elasticities / semi-elasticities and sets e(wdid_cetype), so a binary
* difference-in-differences ATT in levels is printed under the heading
* "Marginal Treatment Effects Presented as Elasticities".
* wooldid.ado:5029 asserts in a comment that the option "has no
* requires-parent screen and is inert" -- the second half is false.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

clear
set obs 60
gen long unit = _n
expand 10
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2003 if unit <= 18
qui replace ttre = 2005 if unit > 18 & unit <= 36
set seed 111
gen double y = 0.07*unit + 0.03*(t-1999) + 0.8*(ttre < . & t >= ttre) + rnormal()*0.3
qui save "rpcet.dta", replace

dis "================ A. BASELINE: no elasticity option ================"
use "rpcet.dta", clear
cap noi wooldid y unit t ttre, att robust
dis "A_rc=" _rc
matrix A = e(wdidmainresults)
dis "A_att = " string(A[1,1], "%21x")
dis "A_cetype = |" e(wdid_cetype) "|"

dis "================ B. + contreatelasticitytype(eyex), NO contreat ================"
use "rpcet.dta", clear
cap noi wooldid y unit t ttre, att robust contreatelasticitytype(eyex)
dis "B_rc=" _rc
matrix B = e(wdidmainresults)
dis "B_att = " string(B[1,1], "%21x")
dis "B_cetype = |" e(wdid_cetype) "|"

dis "================ C. + contreatelasticitytype(eydx), NO contreat ================"
use "rpcet.dta", clear
cap noi wooldid y unit t ttre, att robust contreatelasticitytype(eydx)
dis "C_rc=" _rc
matrix C = e(wdidmainresults)
dis "C_att = " string(C[1,1], "%21x")

dis "================ D. + contreatelasticitytype(dyex), NO contreat ================"
use "rpcet.dta", clear
cap noi wooldid y unit t ttre, att robust contreatelasticitytype(dyex)
dis "D_rc=" _rc
matrix D = e(wdidmainresults)
dis "D_att = " string(D[1,1], "%21x")

dis "================ E. the relabel also reaches the event-study tables and plot titles ================"
use "rpcet.dta", clear
cap noi wooldid y unit t ttre, att robust esprelength(2) espostlength(2) contreatelasticitytype(eyex)
dis "E_rc=" _rc

dis "===VERDICT numbers move A->B : " string(mreldif(A,B), "%12.4e")
dis "===VERDICT numbers move A->C : " string(mreldif(A,C), "%12.4e")
dis "===VERDICT numbers move A->D : " string(mreldif(A,D), "%12.4e")
dis "===VERDICT table header changes: see the 'Main Results' line in blocks A and B above."
dis "===VERDICT there is no continuous treatment in ANY of these runs: contreat() was never given."
dis "REPRO_CET_COMPLETE"
