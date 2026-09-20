*==========================================================================
* repro_skeptic_NFC01.do
*   RUN: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_skeptic_NFC01.do
*   Self-contained: generates its own data, touches no repo file and no other
*   probe's script.  Uses only the pinned build.
*
* WHAT IT SETTLES (skeptic's verdict on NFC-01)
*
*  1 CONFIRMED  On a 150-unit / 3-cohort / 10-period panel with i = cohort and a
*    percentage-scale control (50 + 8*x1) inside controls(), the default cell
*    collapse and nofastcollapse report different VARIANCES: 6.6e-06 relative on
*    e(V) and on its DIAGONAL, and the difference reaches the PRINTED standard
*    errors and confidence bounds of the pretreatment rows.  wooldid.sthlp:338
*    says variances "were measured at about 1e-16 relative regardless of the
*    scale of the predictions" and sthlp:350 says "no digit wooldid displays is
*    affected".  Both are contradicted.
*
*  2 TRIGGER CORRECTED  The claim's "whenever a continuous, coefficient-bearing
*    control varies within the cell" is NOT sufficient.  On the same panel and the
*    same control, a bare att (no event-study window) is clean at 1.9e-16 and
*    esprelength(1) espostlength(1) is clean at 6.5e-16; the gap needs an
*    event-study window of two or more horizons, and the cells that move are the
*    PRETREATMENT rows.  On a 6,000-row micro panel (30 states x 25 people x 8
*    years, i = state) the same control at the same scale, with the same window,
*    is clean at 6.7e-16.
*
*  3 ATTRIBUTION CORRECTED  The claim says the collapse is only "the larger of two
*    perturbations of the same underlying instability", citing a row permutation
*    under nofastcollapse that moved e(V) by 2.92e-07.  Here the UNCOLLAPSED path
*    is stable: a row permutation with the collapse OFF moves e(V) by 1.3e-15,
*    while the same permutation with the collapse ON moves it by 6.6e-06.  The
*    instability is the collapsed path's alone.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 515151
set sortseed 515151

*---------------- panel A: 150 units, 3 treated cohorts + never-treated ----------------
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
gen double xc = 50 + 8*x1
bysort cohort t (unit): gen double xcellc = xc[1]
compress
save "_sk_cohort.dta", replace
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "_sk_cohort_shuf.dta", replace

*---------------- panel B: 6,000-row micro panel, i = state ----------------
clear
set seed 71113
set sortseed 71113
set obs 30
gen int state = _n
gen double u_s = rnormal()
gen int ttre = .
replace ttre = 2013 if state <= 8
replace ttre = 2015 if state > 8 & state <= 16
replace ttre = 2017 if state > 16 & state <= 22
expand 25
bysort state: gen int pid = _n
gen double u_i = rnormal()
expand 8
bysort state pid: gen int year = 2010 + _n
gen double z = rnormal()
gen byte D = (ttre < . & year >= ttre)
gen double y = 2 + u_s + u_i + 0.05*(year-2010) + 0.4*D + 0.10*z + 0.5*rnormal()
gen double xc = 50 + 8*z
compress
save "_sk_micro.dta", replace

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
end

dis _n as text "{hline 100}"
dis as text "case                                                   maxrel e(V)   maxrel diag(e(V))   verdict"
dis as text "{hline 100}"

* ---- 1/2: window sweep on panel A, identical data and identical control ----
foreach k in 1 2 3 4 {
  if `k'==1 local spec att
  if `k'==2 local spec att esprelength(1) espostlength(1)
  if `k'==3 local spec att esprelength(2) espostlength(2)
  if `k'==4 local spec att att_i att_it esprelength(3) espostlength(3)
  use "_sk_cohort.dta", clear
  qui wooldid y cohort t ttre, `spec' controls(xc)
  matrix VA = e(V)
  use "_sk_cohort.dta", clear
  qui wooldid y cohort t ttre, `spec' controls(xc) nofastcollapse
  matrix VB = e(V)
  mata: cmpm("VA","VB")
  scalar s1 = mxall
  scalar s2 = mxdia
  local vd "clean"
  if s1 > 1e-10 local vd "ABOVE rtol 1e-10"
  dis as result %-52s "A  `spec'" "  " %9.2e s1 "        " %9.2e s2 "      `vd'"
}

* ---- cell-constant control: the within-cell variation is necessary ----
use "_sk_cohort.dta", clear
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xcellc)
matrix VA = e(V)
use "_sk_cohort.dta", clear
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xcellc) nofastcollapse
matrix VB = e(V)
mata: cmpm("VA","VB")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-52s "A  cell-CONSTANT control, ES 3/3" "  " %9.2e s1 "        " %9.2e s2 "      clean"

* ---- 3: attribution.  row permutation, collapse OFF in both, then ON in both ----
use "_sk_cohort.dta", clear
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc) nofastcollapse
matrix VA = e(V)
use "_sk_cohort_shuf.dta", clear
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc) nofastcollapse
matrix VB = e(V)
mata: cmpm("VA","VB")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-52s "A  ROW PERMUTATION, collapse OFF in both" "  " %9.2e s1 "        " %9.2e s2 "      stable"
use "_sk_cohort.dta", clear
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc)
matrix VA = e(V)
use "_sk_cohort_shuf.dta", clear
qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc)
matrix VB = e(V)
mata: cmpm("VA","VB")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-52s "A  ROW PERMUTATION, collapse ON  in both" "  " %9.2e s1 "        " %9.2e s2 "      ABOVE rtol"

* ---- 2b: the other realistic shape stays clean ----
use "_sk_micro.dta", clear
qui wooldid y state year ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc)
matrix VA = e(V)
use "_sk_micro.dta", clear
qui wooldid y state year ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc) nofastcollapse
matrix VB = e(V)
mata: cmpm("VA","VB")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-52s "B  micro panel i=state, same control, ES 3/3" "  " %9.2e s1 "        " %9.2e s2 "      clean"
dis as text "{hline 100}"

* ---- the printed table: diff these two blocks ----
dis "<<<<< PRINTED_ON  (default collapse) >>>>>"
use "_sk_cohort.dta", clear
wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc)
dis "<<<<< PRINTED_OFF (nofastcollapse) >>>>>"
use "_sk_cohort.dta", clear
wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc) nofastcollapse
dis "<<<<< PRINTED_END >>>>>"
dis as text "Expect, in the pretreatment block: se .11036045 vs .11036009, se .11064408 vs"
dis as text ".11064417, lb -.24568926 vs -.24568854, ub .18727346 vs .18727275."

erase "_sk_cohort.dta"
erase "_sk_cohort_shuf.dta"
erase "_sk_micro.dta"
dis "REPRO_SKEPTIC_NFC01_DONE"
