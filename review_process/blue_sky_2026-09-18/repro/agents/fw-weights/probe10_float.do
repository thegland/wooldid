* probe10: PRECISION OF THE CELL-WEIGHT TOTALS.
* ado:1619  qui gegen `w_per_it' = total(`esamp' * `w')  , by(`i' `t')
* ado:1620  qui gegen `w_per_i'  = total(`w'*`esamp')    , by(`posttreatment' `i')
* neither asks for `double', while the escomp twin at ado:3833 does:
* ado:3833  qui gegen double `f8wpit' = total(`esamp' * `w'), by(`i' `t')
* If gegen defaults to float there, the att_it / att_i / att_itime cell weights carry a
* ~1e-7 relative error whenever the weights are large - which is what population or
* survey-total weights are.  That is exactly the B21 float-cell-weight class of defect.
* TEST: aweights are scale free, so [aw=pop] and [aw=pop/1e6] must agree exactly.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*--- 0. what type does gegen total() actually produce? --------------------------
clear
set obs 100
gen g = mod(_n,10)
gen double big = 1234567 + _n
gegen tot_default = total(big), by(g)
gegen double tot_double = total(big), by(g)
egen  tot_egen = total(big), by(g)
ds tot_default tot_double tot_egen
describe tot_default tot_double tot_egen
gen double d1 = tot_default - tot_double
qui sum d1, meanonly
dis "GEGEN_DEFAULT_MAXDEV=" %12.4e r(max)
dis "GEGEN_TYPE=" "`: type tot_default'"
dis "EGEN_TYPE=" "`: type tot_egen'"

mata:
void wd_rel(string scalar a, string scalar b, string scalar tag) {
  real matrix A, B, D
  A = st_matrix(a); B = st_matrix(b)
  if (rows(A)==0 | rows(B)==0) {
    printf("REL_%s=EMPTY\n", tag)
  }
  else if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    printf("REL_%s=SHAPEMISMATCH\n", tag)
  }
  else {
    D = abs(A:-B) :/ (abs(A):*(abs(A):>1e-12) :+ (abs(A):<=1e-12))
    printf("REL_%s=%10.3e  maxabs=%10.3e\n", tag, colmax(rowmax(D)), colmax(rowmax(abs(A:-B))))
  }
}
end

*--- 1. population-scale weights: scale invariance -----------------------------
use "p2_A.dta", clear
set seed 4711
gen double pop  = 250000 + int(runiform()*9000000)   /* county/state population scale */
gen double popm = pop/1000000
gen double popk = pop/1000
save "p10.dta", replace
sum pop popm, format

foreach o in "att att_i att_it att_itime" "att att_it esprelength(2) espostlength(3)" "att_i att_itime cluster(i)" {
  tempname b1 b2 b3
  qui use "p10.dta", clear
  qui wooldid y i t ttre [aw=pop], `o'
  matrix `b1' = e(b)
  qui wooldid y i t ttre [aw=popm], `o'
  matrix `b2' = e(b)
  qui wooldid y i t ttre [aw=popk], `o'
  matrix `b3' = e(b)
  dis "SCALE opts=`o'"
  mata: wd_rel("`b1'","`b2'","POP_vs_POPM")
  mata: wd_rel("`b1'","`b3'","POP_vs_POPK")
  mata: wd_rel("`b2'","`b3'","POPM_vs_POPK")
}

*--- 2. att alone must be immune (no w_per_it in its weight) --------------------
qui use "p10.dta", clear
tempname a1 a2
qui wooldid y i t ttre [aw=pop], att
matrix `a1' = e(b)
qui wooldid y i t ttre [aw=popm], att
matrix `a2' = e(b)
mata: wd_rel("`a1'","`a2'","ATTONLY_POP_vs_POPM")

*--- 3. an even larger scale: national survey weights (sum ~ 3e8) ---------------
use "p10.dta", clear
gen double swt = pop*30
gen double swts = swt/1e9
save "p10b.dta", replace
tempname c1 c2
qui use "p10b.dta", clear
qui wooldid y i t ttre [aw=swt], att att_i att_it att_itime
matrix `c1' = e(b)
qui wooldid y i t ttre [aw=swts], att att_i att_it att_itime
matrix `c2' = e(b)
mata: wd_rel("`c1'","`c2'","BIGSURVEY")

dis "PROBE10_DONE"
