*==========================================================================
* r4.do — WHICH cell moves, how big is it, and is either arm simply WRONG?
*
* r3 reproduced the default-vs-nofastcollapse gap on my own cohort-shaped panel
* (6.6e-06 on e(V), diagonal included).  Two questions decide severity:
*   Q1 is the moving diagonal element an ORDINARY-SIZED standard error, or a
*      near-zero one whose relative difference is an artifact?  -> name the cells.
*   Q2 adding a CONSTANT to a control changes nothing in exact arithmetic (the
*      fixed effects absorb it).  So compare each arm against ITSELF at c=0 and
*      at c=50: whichever arm moves is the arm carrying numerical error, and the
*      size of that move is the accuracy of wooldid's reported standard errors.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 515151
set sortseed 515151

*------------------------------------------------------------------ data (same DGP as r3)
set obs 150
gen long unit = _n
gen double _u = runiform()
sort _u unit
gen long rk = _n
gen long cohort = 4
gen double ttre = .
replace cohort = 1 in 1/33
replace ttre = 2004 in 1/33
replace cohort = 2 in 34/66
replace ttre = 2006 in 34/66
replace cohort = 3 in 67/98
replace ttre = 2008 in 67/98
gen double a_i = rnormal()
gen long region = mod(unit,5)+1
drop _u rk
expand 10
bysort unit: gen int t = 1999 + _n
sort t unit
by t: gen double _f = rnormal() if _n==1
by t: replace _f = _f[1]
gen double ft = 0.5*_f
drop _f
gen double x1 = rnormal()
gen double x2 = 0.5*a_i + 0.5*rnormal()
gen byte post = (ttre < . & t >= ttre)
gen double tau = 0
replace tau = 0.50 if cohort==1 & post
replace tau = 0.80 if cohort==2 & post
replace tau = 0.30 if cohort==3 & post
gen double y = a_i + ft + tau + 0.5*x1 - 0.25*x2 + 0.30*rnormal()
compress
save "_r4.dta", replace

mata:
void dumpdiff(string scalar an, string scalar bn, real scalar tol, string scalar tag)
{
  real matrix A, B
  string matrix rn, cn
  real scalar i, j, d, s, nshown
  A = st_matrix(an); B = st_matrix(bn)
  rn = st_matrixrowstripe(an)
  cn = st_matrixcolstripe(an)
  nshown = 0
  for (i=1; i<=rows(A); i++) {
    for (j=1; j<=cols(A); j++) {
      s = max((abs(A[i,j]), abs(B[i,j])))
      if (s == 0) continue
      d = abs(A[i,j]-B[i,j])/s
      if (d > tol & nshown < 25) {
        printf("%s  [%s , %s]  A=%21.18g  B=%21.18g  rel=%9.2e\n",
               tag, rn[i,2], cn[j,2], A[i,j], B[i,j], d)
        nshown++
      }
    }
  }
  if (nshown == 0) printf("%s  no cell above tol\n", tag)
}
void worstdiag(string scalar an, string scalar bn, string scalar tag)
{
  real matrix A, B
  string matrix rn
  real scalar i, d, s, bd, bi
  A = st_matrix(an); B = st_matrix(bn)
  rn = st_matrixrowstripe(an)
  bd = 0; bi = 1
  for (i=1; i<=rows(A); i++) {
    s = max((abs(A[i,i]), abs(B[i,i])))
    if (s == 0) continue
    d = abs(A[i,i]-B[i,i])/s
    if (d > bd) { bd = d; bi = i }
  }
  printf("%s  worst DIAGONAL of e(V): [%s]  V_A=%21.18g  V_B=%21.18g  rel=%9.2e   se_A=%21.18g  se_B=%21.18g\n",
         tag, rn[bi,2], A[bi,bi], B[bi,bi], bd, sqrt(A[bi,bi]), sqrt(B[bi,bi]))
}
end

local ES att att_i att_it esprelength(3) espostlength(3)

*------------------------------------------------------------------ c = 0 and c = 50
foreach c in 0 50 {
  use "_r4.dta", clear
  gen double xc = `c' + 8*x1
  qui wooldid y cohort t ttre, `ES' controls(xc)
  matrix VON`c' = e(V)
  matrix MON`c' = e(wdidmainresults)
  use "_r4.dta", clear
  gen double xc = `c' + 8*x1
  qui wooldid y cohort t ttre, `ES' controls(xc) nofastcollapse
  matrix VOF`c' = e(V)
  matrix MOF`c' = e(wdidmainresults)
}

dis _n as text "{hline 100}"
dis as text "Q1: WHICH cells move between the arms (tol 1e-9), with their VALUES"
dis as text "{hline 100}"
mata: worstdiag("VON0","VOF0","c=0 ")
mata: worstdiag("VON50","VOF50","c=50")
dis as text "-- e(wdidmainresults), c=0, collapse ON vs OFF --"
mata: dumpdiff("MON0","MOF0",1e-9,"c=0 ")
dis as text "-- e(wdidmainresults), c=50, collapse ON vs OFF --"
mata: dumpdiff("MON50","MOF50",1e-9,"c=50")

dis _n as text "{hline 100}"
dis as text "Q2: SELF-CONSISTENCY.  Adding a constant to a control cannot change any"
dis as text "    reported quantity in exact arithmetic.  Each arm against itself, c=0 vs c=50."
dis as text "{hline 100}"
dis as text "-- collapse ON  : c=0 vs c=50 --"
mata: worstdiag("VON0","VON50","ON  ")
mata: dumpdiff("MON0","MON50",1e-9,"ON  ")
dis as text "-- collapse OFF : c=0 vs c=50 --"
mata: worstdiag("VOF0","VOF50","OFF ")
mata: dumpdiff("MOF0","MOF50",1e-9,"OFF ")

dis "R4_DONE"
