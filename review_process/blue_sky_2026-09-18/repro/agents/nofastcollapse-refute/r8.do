*==========================================================================
* r8.do — MECHANISM + REMEDY, with a one-line patched COPY of wooldid.ado.
*
* The copy at ./patchbuild/wooldid.ado differs from the pinned build by three
* lines inserted after "local fccont" in the collapse gate:
*     if "$WDPATCH_AVGCTRL" != "" & "`controls'" != "" local fccont contreat(`controls')
* i.e. it hands the user's controls() variables to wooldid_fcol's dose-averaging
* branch, so the retained row carries the cell's WEIGHTED MEAN of each control
* instead of the first row's raw value.  That is exactly what the code already
* does for the contreat() dose (ado ~1894-1901).
*
* ARM P0  patch INERT (global unset): must reproduce the 6.59e-06 gap
* ARM P1  patch ARMED               : if the gap collapses to machine precision,
*         the retained row's control VALUE is the mechanism and averaging is the
*         remedy.
* ARM P2  patch ARMED, row-shuffled data, collapse ON in both: the storage-order
*         dependence must vanish too.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse-refute/patchbuild"
which wooldid
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
compress
save "_r8.dta", replace
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "_r8_shuf.dta", replace

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

local ES att att_i att_it esprelength(3) espostlength(3)

dis _n as text "{hline 96}"
dis as text "arm                                              maxrel V(all)   maxrel V(diag)"
dis as text "{hline 96}"

* P0: patch inert
global WDPATCH_AVGCTRL
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix VA = e(V)
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc) nofastcollapse
matrix VB = e(V)
mata: cmpm("VA","VB")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-48s "P0 patch INERT, ON vs OFF" "  " %9.2e s1 "       " %9.2e s2

* P1: patch armed
global WDPATCH_AVGCTRL 1
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix VC = e(V)
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc) nofastcollapse
matrix VD = e(V)
mata: cmpm("VC","VD")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-48s "P1 patch ARMED, ON vs OFF" "  " %9.2e s1 "       " %9.2e s2

* P2: patch armed, storage-order dependence
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix VE = e(V)
use "_r8_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix VF = e(V)
mata: cmpm("VE","VF")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-48s "P2 patch ARMED, shuffle, ON in both" "  " %9.2e s1 "       " %9.2e s2

* P3: control: patch inert, storage-order dependence (must be ~6.6e-06)
global WDPATCH_AVGCTRL
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix VG = e(V)
use "_r8_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix VH = e(V)
mata: cmpm("VG","VH")
scalar s1 = mxall
scalar s2 = mxdia
dis as result %-48s "P3 patch INERT, shuffle, ON in both" "  " %9.2e s1 "       " %9.2e s2

* P4: does the patch move the POINT ESTIMATES?  (it must not, materially)
global WDPATCH_AVGCTRL 1
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc)
matrix BC = e(b)
global WDPATCH_AVGCTRL
use "_r8.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xc) nofastcollapse
matrix BD = e(b)
mata: cmpm("BC","BD")
scalar s1 = mxall
dis as result %-48s "P4 patched-ON vs OFF, max rel on e(b)" "  " %9.2e s1

dis as text "{hline 96}"
dis "R8_DONE"
