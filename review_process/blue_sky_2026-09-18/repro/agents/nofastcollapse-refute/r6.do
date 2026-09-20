*==========================================================================
* r6.do — three loose ends.
*  A  ROW SHUFFLE arms on the FAILING panel: is the uncollapsed path also
*     storage-order dependent there (the prober's arms F/G), or is the
*     instability confined to the collapsed path?
*  B  cell-CONSTANT control with the same event-study window: is the within-cell
*     variation necessary on my data too?
*  C  the 6000-row micro panel that stayed clean in r1/r5, now with THREE onset
*     cohorts and esprelength(3) espostlength(3): does the failing ingredient
*     travel to a different panel shape?
*  D  the PRINTED table, both arms, on the failing panel.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 515151
set sortseed 515151

*---------------- failing panel (same DGP as r3/r4b) ----------------
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
save "_r6_fail.dta", replace
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "_r6_fail_shuf.dta", replace

*---------------- micro panel with THREE onsets ----------------
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
save "_r6_micro3.dta", replace

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

* A: shuffle, collapse OFF in both / ON in both, on the FAILING panel
foreach a in SHUFOFF SHUFON CELLC MICRO3 {
  local dA "_r6_fail.dta"
  local dB "_r6_fail.dta"
  local ctrl controls(xc)
  local oA
  local oB nofastcollapse
  local lbl
  if "`a'"=="SHUFOFF" {
    local lbl "A1 row shuffle, collapse OFF in BOTH"
    local oA nofastcollapse
    local dB "_r6_fail_shuf.dta"
  }
  if "`a'"=="SHUFON" {
    local lbl "A2 row shuffle, collapse ON in BOTH"
    local oB
    local dB "_r6_fail_shuf.dta"
  }
  if "`a'"=="CELLC" {
    local lbl "B  cell-CONSTANT control, ON vs OFF"
    local ctrl controls(xcellc)
  }
  if "`a'"=="MICRO3" {
    local lbl "C  micro panel, 3 onsets, ES 3/3, ON vs OFF"
    local dA "_r6_micro3.dta"
    local dB "_r6_micro3.dta"
  }
  local iv cohort
  local tv t
  if "`a'"=="MICRO3" {
    local iv state
    local tv year
  }
  use "`dA'", clear
  cap noi qui wooldid y `iv' `tv' ttre, `ES' `ctrl' `oA'
  local rcA = _rc
  if `rcA'==0 {
    matrix VA = e(V)
  }
  use "`dB'", clear
  cap noi qui wooldid y `iv' `tv' ttre, `ES' `ctrl' `oB'
  local rcB = _rc
  if `rcB'==0 {
    matrix VB = e(V)
  }
  if `rcA'!=0 | `rcB'!=0 {
    dis as error "`lbl' : rcA=`rcA' rcB=`rcB'"
    continue
  }
  mata: cmpm("VA","VB")
  scalar sVall = mxall
  scalar sVdia = mxdia
  dis as result %-48s "`lbl'" "  " %9.2e sVall "       " %9.2e sVdia
}
dis as text "{hline 96}"

*---------------- D: the printed table ----------------
dis "<<<<< PRINTED_ON  (default collapse) >>>>>"
use "_r6_fail.dta", clear
cap noi wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc)
dis "<<<<< PRINTED_OFF (nofastcollapse) >>>>>"
use "_r6_fail.dta", clear
cap noi wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) controls(xc) nofastcollapse
dis "<<<<< PRINTED_END >>>>>"

dis "R6_DONE"
