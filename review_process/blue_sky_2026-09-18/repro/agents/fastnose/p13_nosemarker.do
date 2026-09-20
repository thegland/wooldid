* p13 -- mechanism check, no wooldid involved: does `margins, nose' suppress the
* not-estimable omission marker in r(table)'s column names?  That premise is what
* ado:1974-1981 (BS-26) is built on, and it is the reason the binary subgroup-by-horizon
* call keeps its variance in fast mode while the CONTINUOUS twin (ado:3290) does not.
clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916

set obs 200
gen byte g = 1 + mod(_n,4)
gen byte s = (_n <= 100)
gen double x = rnormal()
gen double y = 1 + 0.5*g + x + rnormal()
* keep g==4 in the over() levels but remove it from the ESTIMATION sample, so its
* subpop margin is not estimable and margins must mark it (noestimcheck computes anyway)
replace y = . if g == 4 & s == 1

qui reg y i.g##i.s x
dis "MODEL COLNAMES: `: colnames e(b)'"

qui margins , over(g) subpop(if s == 1) noestimcheck df(`e(df_r)')
dis "WITH_VCE_rows=" rowsof(r(table)) " cols=" colsof(r(table))
dis "WITH_VCE_COLNAMES: `: colnames r(table)'"
dis "WITH_VCE_B_COLNAMES: `: colnames r(b)'"

qui margins , over(g) subpop(if s == 1) noestimcheck nose
dis "NOSE_rows=" rowsof(r(table)) " cols=" colsof(r(table))
dis "NOSE_COLNAMES: `: colnames r(table)'"
dis "NOSE_B_COLNAMES: `: colnames r(b)'"
cap confirm matrix r(V)
dis "NOSE_hasV_rc=" _rc
dis "DONE_P13"
