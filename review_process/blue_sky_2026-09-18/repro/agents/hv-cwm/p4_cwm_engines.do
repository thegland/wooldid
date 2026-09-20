* ============================================================================
* p4_cwm_engines.do -- does customweightmultiplier() survive the replication
* engines?  sthlp:532 promises "point estimates are unchanged" under
* bootstrap(); the same holds for permute().  cwm is forwarded into the
* replication option list at ado:134 (permute) and ado:329 (bootstrap).
* Also confirms the bayes x cell-weighted gate actually fires under cwm.
* Self-contained.
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242
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
gen int county = _n
gen int adopt = .
replace adopt = 4 if county==1
replace adopt = 5 if county==2
replace adopt = 6 if county==3
replace adopt = 7 if county==4
gen double pop = cond(inlist(county,3,4), 50, 1)
expand 10
bys county: gen int year = _n
expand 5
bys county year: gen byte estab = _n
gen double sw = 0.5 + runiform()
gen byte post = (adopt < .) & (year >= adopt)
gen double tau = cond(county==1,1,cond(county==2,2,cond(county==3,3,cond(county==4,4,0))))
replace tau = tau + 0.2*(year-adopt) if post==1
replace tau = 0 if post==0
gen double y = 0.3*county + 0.1*year + tau*post + 0.05*rnormal()
compress
save "panel4.dta", replace

dis _n "=== N0: base (cwm, cluster) ==="
use "panel4.dta", clear
cap noi wooldid y county year adopt [aw=sw], att att_it cluster(county) customweightmultiplier(pop)
dis "N0_RC=" _rc
if _rc == 0 mat b0 = e(b)

dis _n "=== N1: + bootstrap(standard) ==="
use "panel4.dta", clear
cap noi wooldid y county year adopt [aw=sw], att att_it cluster(county) customweightmultiplier(pop) bootstrap(standard) reps(25) rseed(7) dots(0)
dis "N1_RC=" _rc
if _rc == 0 {
  mat b1 = e(b)
  cap noi mata: wd_cmp("N1_boot_vs_base", "b0", "b1")
}

dis _n "=== N2: + permute(shuffle) ==="
use "panel4.dta", clear
cap noi wooldid y county year adopt [aw=sw], att att_it cluster(county) customweightmultiplier(pop) permute(shuffle) permreps(25) rseed(7) dots(0)
dis "N2_RC=" _rc
if _rc == 0 {
  mat b2 = e(b)
  cap noi mata: wd_cmp("N2_perm_vs_base", "b0", "b2")
}

dis _n "=== N3: bootstrap(bayes) + att_it (gate should refuse) ==="
use "panel4.dta", clear
cap noi wooldid y county year adopt [aw=sw], att att_it cluster(county) customweightmultiplier(pop) bootstrap(bayes) reps(25) rseed(7) dots(0)
dis "N3_RC=" _rc

dis _n "=== N4: bootstrap(bayes) + att only ==="
use "panel4.dta", clear
cap noi wooldid y county year adopt [aw=sw], att cluster(county) customweightmultiplier(pop) bootstrap(bayes) reps(25) rseed(7) dots(0)
dis "N4_RC=" _rc
if _rc == 0 {
  mat b4 = e(b)
  cap noi mata: wd_cmp("N4_bayes_vs_base", "b0", "b4")
}
dis "DONE_P4"
