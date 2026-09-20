* ============================================================================
* p2_cwm_adjacent.do -- novel probes on the customweightmultiplier() surface
*   F  : cwm that GROWS over time -> do the pretreatment and main att_i rows
*        weight the cohorts differently? (they are both equal-cohort with no cwm)
*   G  : cwm as a 0/1 subpopulation indicator that is MIXED within a cell
*   H  : contreat() + within-cell-varying cwm, default cell collapse vs
*        nofastcollapse  (the collapse's own comment claims they agree)
*   K  : binary + a cwm that is exactly ZERO on the cell's first row,
*        default collapse vs nofastcollapse
* Self-contained; generates its own data.
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

mata:
void wd_cmp(string scalar tag, string scalar m1, string scalar m2)
{
    real matrix A, B
    string colvector an, bn
    real scalar j, k, mx, d
    A = st_matrix(m1); an = st_matrixcolstripe(m1)[.,2]
    B = st_matrix(m2); bn = st_matrixcolstripe(m2)[.,2]
    mx = 0
    for (j=1; j<=length(an); j++) {
        k = 0
        if (anyof(bn, an[j])) k = selectindex(bn:==an[j])[1]
        if (k>0) {
            d = abs(A[1,j]-B[1,k])/max((abs(A[1,j]),1e-12))
            if (d>mx) mx = d
        }
    }
    printf("CMP %s maxreldif=%g\n", tag, mx)
}
end

* ---------------------------------------------------------------- data ----
clear
set obs 8
gen int i = _n
gen int ttre = .
replace ttre = 4 if i==1
replace ttre = 5 if i==2
replace ttre = 6 if i==3
replace ttre = 7 if i==4
gen double pop0 = cond(inlist(i,3,4), 50, 1)
expand 10
bys i: gen int t = _n
expand 5
bys i t: gen byte k = _n
gen long obsid = _n
gen double w = 0.5 + runiform()
* F: population-like multiplier that GROWS over time (cell-constant)
gen double popgrow = pop0 * (1.15^t)
* G: 0/1 subpopulation indicator, mixed within the cell
gen byte female = (runiform() < 0.5)
* H/K: cwm varying within cell, and one that is exactly zero on the cell's 1st row
gen double cwmind = 0.2 + 2*runiform()
gen double cwmz   = cond(k==1, 0, cwmind)
gen byte D = (ttre < .) & (t >= ttre)
gen double dose = 0
replace dose = 0.5 + runiform() if D==1            // varies WITHIN the cell
gen double tau = cond(i==1,1.0,cond(i==2,2.0,cond(i==3,3.0,cond(i==4,4.0,0))))
replace tau = tau + 0.2*(t-ttre) if D==1
replace tau = 0 if D==0
gen double y  = 0.3*i + 0.1*t + tau*D + 0.05*rnormal()
gen double yc = 0.3*i + 0.1*t + 1.7*dose + 0.05*rnormal()
compress
save "panel2.dta", replace

* ================================ ARM F ====================================
dis _n(2) "======================== ARM F (growing cwm) ========================"
use "panel2.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects customweightmultiplier(popgrow)
dis "ARMF_RC=" _rc
if _rc == 0 {
  mat G = e(gteffects)
  preserve
    clear
    svmat double G, names(col)
    gen byte post = relyr >= 0
    collapse (sum) w_att_i w_att_itime w_att_it (count) n=relyr, by(i post)
    dis as result "ARM F: per-cohort TOTAL weight, pre block (post=0) vs post block (post=1)"
    list, noobs sepby(i)
  restore
}
* the same table with NO cwm, for contrast
use "panel2.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime gteffects
dis "ARMF0_RC=" _rc
if _rc == 0 {
  mat G = e(gteffects)
  preserve
    clear
    svmat double G, names(col)
    gen byte post = relyr >= 0
    collapse (sum) w_att_i w_att_itime w_att_it (count) n=relyr, by(i post)
    dis as result "ARM F0 (no cwm): per-cohort TOTAL weight, pre vs post"
    list, noobs sepby(i)
  restore
}

* ================================ ARM G ====================================
dis _n(2) "======================== ARM G (0/1 cwm, mixed cells) ==============="
use "panel2.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it gteffects customweightmultiplier(female)
dis "ARMG_RC=" _rc
if _rc == 0 {
  mat G = e(gteffects)
  preserve
    gen double _w = w
    gen double _wc = w*female
    collapse (sum) _w _wc, by(i t)
    gen double share = _wc/_w
    tempfile ORC
    save `ORC'
    clear
    svmat double G, names(col)
    rename caltime t
    merge m:1 i t using `ORC', keep(match master) nogen
    gen double reld = abs(w_att_it-share)/max(abs(share),1e-30)
    qui su reld
    dis "ARMG_MAXRELD=" r(max)
    qui su w_att_it
    dis "ARMG_WIT min=" r(min) " max=" r(max)
    list i t relyr estimate w_att_it share if relyr>=0 & i<=2, noobs
  restore
}

* ============================ ARMS H and K =================================
foreach case in H K {
  foreach cw in cwmind cwmz {
    dis _n(2) "==================== ARM `case' / `cw' ===================="
    local yv = cond("`case'"=="H", "yc", "y")
    local ct = cond("`case'"=="H", "contreat(dose)", "")
    use "panel2.dta", clear
    cap noi wooldid `yv' i t ttre [aw=w], att att_it att_i att_itime `ct' customweightmultiplier(`cw')
    dis "ARM`case'`cw'_FAST_RC=" _rc
    if _rc == 0 mat bF = e(b)
    use "panel2.dta", clear
    cap noi wooldid `yv' i t ttre [aw=w], att att_it att_i att_itime `ct' customweightmultiplier(`cw') nofastcollapse
    dis "ARM`case'`cw'_SLOW_RC=" _rc
    if _rc == 0 mat bS = e(b)
    cap mata: wd_cmp("`case'_`cw'", "bF", "bS")
    dis "CMPRC_`case'_`cw'=" _rc
    cap mat drop bF
    cap mat drop bS
  }
}
dis "DONE_P2"
