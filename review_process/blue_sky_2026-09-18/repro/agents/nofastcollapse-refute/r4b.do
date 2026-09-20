*==========================================================================
* r4b.do — isolate the trigger on the cohort panel, name the moving cells, and
*          ask which arm is the inaccurate one.
*
* Facts so far:
*   r3 (150 units, 3 cohorts + never, i = cohort, 10 periods, ES 3/3):
*        controls(x) moves e(V) by 3.3e-07 (mean-zero control) to 4.7e-05.
*   r1/r5 (6000-row micro panel, i = state OR i = 3-value cohort, ES 2/2):
*        clean at 7e-16 for every control scale up to 5000.
*   Difference in the REQUEST: r3 asks esprelength(3) espostlength(3) on a panel
*   whose last cohort (onset 2008 in a 2000-2009 window) has only horizons 0 and 1.
*
* PART 1  event-study window sweep (att only / 1,1 / 2,2 / 3,3)
* PART 2  estimand sweep
* PART 3  name the moving cells with their values (is it an ordinary-sized se?)
* PART 4  self-consistency: adding a CONSTANT to a control cannot change any
*         reported quantity in exact arithmetic - so compare each arm to itself
*         at c=0 vs c=50.  Whichever arm moves carries the numerical error.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 515151
set sortseed 515151

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
save "_r4b.dta", replace

mata:
void cmpm(string scalar an, string scalar bn)
{
    real matrix A, B, S, Dm
    A = st_matrix(an)
    B = st_matrix(bn)
    if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
        st_numscalar("mxall", .)
        st_numscalar("mxdia", .)
        return
    }
    S  = abs(A)
    S  = S :* (S:>=abs(B)) + abs(B) :* (S:<abs(B))
    Dm = abs(A-B) :/ (S :+ (S:==0))
    st_numscalar("mxall", max(Dm))
    if (rows(A)==cols(A)) {
        st_numscalar("mxdia", max(diagonal(Dm)))
    }
    else {
        st_numscalar("mxdia", .)
    }
}
void dumpdiff(string scalar an, string scalar bn, real scalar tol, string scalar tag)
{
    real matrix A, B
    string matrix rn, cn
    real scalar i, j, d, s, nshown
    A = st_matrix(an)
    B = st_matrix(bn)
    rn = st_matrixrowstripe(an)
    cn = st_matrixcolstripe(an)
    nshown = 0
    for (i=1; i<=rows(A); i++) {
        for (j=1; j<=cols(A); j++) {
            s = abs(A[i,j])
            if (abs(B[i,j]) > s) {
                s = abs(B[i,j])
            }
            if (s > 0) {
                d = abs(A[i,j]-B[i,j])/s
                if (d > tol & nshown < 30) {
                    printf("%s [%s,%s] A=%19.17g B=%19.17g rel=%9.2e\n", tag, rn[i,2], cn[j,2], A[i,j], B[i,j], d)
                    nshown = nshown + 1
                }
            }
        }
    }
    if (nshown == 0) {
        printf("%s no cell above tol\n", tag)
    }
}
void worstdiag(string scalar an, string scalar bn, string scalar tag)
{
    real matrix A, B
    string matrix rn
    real scalar i, d, s, bd, bi
    A = st_matrix(an)
    B = st_matrix(bn)
    rn = st_matrixrowstripe(an)
    bd = 0
    bi = 1
    for (i=1; i<=rows(A); i++) {
        s = abs(A[i,i])
        if (abs(B[i,i]) > s) {
            s = abs(B[i,i])
        }
        if (s > 0) {
            d = abs(A[i,i]-B[i,i])/s
            if (d > bd) {
                bd = d
                bi = i
            }
        }
    }
    printf("%s worst V diagonal [%s] se_A=%19.17g se_B=%19.17g relV=%9.2e\n", tag, rn[bi,2], sqrt(A[bi,bi]), sqrt(B[bi,bi]), bd)
}
end

*------------------------------------------------------------- PART 1 + 2
dis _n as text "{hline 96}"
dis as text "spec                                         maxrel V(all)   maxrel V(diag)   maxrel main"
dis as text "{hline 96}"
local s1 att
local s2 att esprelength(1) espostlength(1)
local s3 att esprelength(2) espostlength(2)
local s4 att esprelength(3) espostlength(3)
local s5 att att_i att_it
local s6 att att_i att_it esprelength(3) espostlength(3)
local s7 att esprelength(5) espostlength(5)
forval k = 1/7 {
  local spec `s`k''
  use "_r4b.dta", clear
  gen double xc = 50 + 8*x1
  cap noi qui wooldid y cohort t ttre, `spec' controls(xc)
  local rcA = _rc
  if `rcA'==0 {
    matrix VA = e(V)
    matrix MA = e(wdidmainresults)
  }
  use "_r4b.dta", clear
  gen double xc = 50 + 8*x1
  cap noi qui wooldid y cohort t ttre, `spec' controls(xc) nofastcollapse
  local rcB = _rc
  if `rcB'==0 {
    matrix VB = e(V)
    matrix MB = e(wdidmainresults)
  }
  if `rcA'!=0 | `rcB'!=0 {
    dis as error "spec `k' (`spec') : rcA=`rcA' rcB=`rcB'"
    continue
  }
  mata: cmpm("VA","VB")
  scalar sVall = mxall
  scalar sVdia = mxdia
  mata: cmpm("MA","MB")
  scalar sMall = mxall
  dis as result %-44s "`spec'" "  " %9.2e sVall "       " %9.2e sVdia "      " %9.2e sMall
}
dis as text "{hline 96}"

*------------------------------------------------------------- PART 3 + 4
foreach c in 0 50 {
  use "_r4b.dta", clear
  gen double xc = `c' + 8*x1
  qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc)
  matrix VON`c' = e(V)
  matrix MON`c' = e(wdidmainresults)
  use "_r4b.dta", clear
  gen double xc = `c' + 8*x1
  qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc) nofastcollapse
  matrix VOF`c' = e(V)
  matrix MOF`c' = e(wdidmainresults)
}
dis _n as text "PART 3: which cells move between the arms, with values (tol 1e-9)"
mata: worstdiag("VON0","VOF0","c=0  ON-vs-OFF")
mata: worstdiag("VON50","VOF50","c=50 ON-vs-OFF")
dis as text "-- e(wdidmainresults) c=0, ON vs OFF --"
mata: dumpdiff("MON0","MOF0",1e-9,"c=0 ")
dis as text "-- e(wdidmainresults) c=50, ON vs OFF --"
mata: dumpdiff("MON50","MOF50",1e-9,"c=50")

dis _n as text "PART 4: self-consistency.  c=0 vs c=50 WITHIN an arm (must be exact)"
mata: worstdiag("VON0","VON50","ON  c0-vs-c50")
mata: worstdiag("VOF0","VOF50","OFF c0-vs-c50")
dis as text "-- collapse ON, c=0 vs c=50 --"
mata: dumpdiff("MON0","MON50",1e-9,"ON  ")
dis as text "-- collapse OFF, c=0 vs c=50 --"
mata: dumpdiff("MOF0","MOF50",1e-9,"OFF ")

dis "R4B_DONE"
