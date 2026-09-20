* ============================================================================
* p3_cwm_collapse.do
*   H : contreat() + within-cell-varying cwm; default cell collapse vs
*       nofastcollapse.  The collapse's own comment (ado:5627-5644) claims the
*       two agree because every estimand weight is w*cwm over a cell-constant
*       denominator.  cwmz is exactly ZERO on the cell's FIRST row, which the
*       same comment names as the case that once broke.
*   K : the same comparison on a plain binary run.
*   L : cross-path consistency: e(escomp)'s wshare_att_it vs the same share
*       rebuilt from e(gteffects)'s w_att_it column, under a cwm.
* Self-contained.
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
    real scalar j, k, mx, d, nm
    A = st_matrix(m1); an = st_matrixcolstripe(m1)[.,2]
    B = st_matrix(m2); bn = st_matrixcolstripe(m2)[.,2]
    mx = 0 ; nm = 0
    for (j=1; j<=length(an); j++) {
        k = 0
        if (anyof(bn, an[j])) k = selectindex(bn:==an[j])[1]
        if (k>0) {
            nm++
            d = abs(A[1,j]-B[1,k])/max((abs(A[1,j]),1e-12))
            if (d>mx) mx = d
        }
    }
    printf("CMP %s ncommon=%g maxreldif=%g\n", tag, nm, mx)
}
end

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
gen double w = 0.5 + runiform()
gen double cwmind = 0.2 + 2*runiform()
gen double cwmz   = cond(k==1, 0, cwmind)
gen byte D = (ttre < .) & (t >= ttre)
gen double dose = 0
replace dose = 0.5 + runiform() if D==1
gen double tau = cond(i==1,1.0,cond(i==2,2.0,cond(i==3,3.0,cond(i==4,4.0,0))))
replace tau = tau + 0.2*(t-ttre) if D==1
replace tau = 0 if D==0
gen double y  = 0.3*i + 0.1*t + tau*D + 0.05*rnormal()
gen double yc = 0.3*i + 0.1*t + 1.7*dose + 0.05*rnormal()
compress
save "panel3.dta", replace

foreach case in H K {
  foreach cw in cwmind cwmz {
    dis _n(2) "==================== ARM `case' / `cw' ===================="
    local yv = cond("`case'"=="H", "yc", "y")
    local ct = cond("`case'"=="H", "contreat(dose) contreatcontrols(basic)", "")
    use "panel3.dta", clear
    cap noi wooldid `yv' i t ttre [aw=w], att att_it att_i att_itime `ct' customweightmultiplier(`cw')
    dis "ARM`case'`cw'_FAST_RC=" _rc
    if _rc == 0 {
      mat bF = e(b)
      dis "FASTb:"
      mat list bF, format(%21x)
    }
    use "panel3.dta", clear
    cap noi wooldid `yv' i t ttre [aw=w], att att_it att_i att_itime `ct' customweightmultiplier(`cw') nofastcollapse
    dis "ARM`case'`cw'_SLOW_RC=" _rc
    if _rc == 0 {
      mat bS = e(b)
      dis "SLOWb:"
      mat list bS, format(%21x)
    }
    cap noi mata: wd_cmp("`case'_`cw'", "bF", "bS")
    dis "CMPRC_`case'_`cw'=" _rc
    cap mat drop bF
    cap mat drop bS
  }
}

* ================================ ARM L ====================================
dis _n(2) "==================== ARM L (escomp vs gteffects) ===================="
use "panel3.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it gteffects customweightmultiplier(pop0) esprelength(3) espostlength(4)
dis "ARML_RC=" _rc
if _rc == 0 {
  mat G = e(gteffects)
  mat EC = e(escomp)
  dis "ESCOMP:"
  mat list EC
  preserve
    clear
    svmat double G, names(col)
    keep if inwindow == 1 & relyr >= 0
    collapse (sum) w_att_it w_att, by(relyr)
    qui su w_att_it
    gen double share_it = w_att_it/r(sum)
    qui su w_att
    gen double share_a = w_att/r(sum)
    dis as result "ARM L: horizon shares rebuilt from gteffects"
    list relyr w_att_it share_it w_att share_a, noobs
  restore
}
dis "DONE_P3"
