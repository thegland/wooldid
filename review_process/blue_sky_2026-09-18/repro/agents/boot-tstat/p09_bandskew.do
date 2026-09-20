clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid

*---- DGP: 12 states, ONE treated state with a much larger effect.
*     Heterogeneous treatment effects with a long right tail is an ordinary
*     applied situation; it makes the cluster-bootstrap draw distribution of
*     the ES coefficients strongly right-skewed.
clear
set obs 12
gen int state = _n
gen int onset = .
replace onset = 2002 if inrange(state,1,3)
replace onset = 2003 if inrange(state,4,6)
replace onset = 2004 if inrange(state,7,9)
gen double taus = 1
replace taus = 14 if state == 1
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.3
expand 6
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal()*0.2 if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = taus
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.15
gen int ttre = onset
dis "NOBS=" _N
save skewclus, replace

dis "=== K1: percentile display + simultaneousbands (m = 2 horizons) ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(999) rseed(4) att esprelength(1) espostlength(1) simultaneousbands dots(0)
dis "K1_RC=" _rc
dis "K1_ENGINE=" "`e(esband_engine)'"
dis "K1_CRIT=" e(esband_crit_att)
dis "K1_BANDN=" e(esband_n_att) " COMPLETED=" e(boot_completed)
local BFAM att
local BTAG K1
include "`WD'/bandcmp.do"
* skewness of the draws
mat BDK = e(bootdist)
mat BOK = e(b)
mata:
  A = st_matrix("BDK"); b = st_matrix("BOK")
  cs = st_matrixcolstripe("BDK")
  for (k=1;k<=cols(A);k++) {
    v = select(A[.,k], A[.,k] :< .)
    n = rows(v)
    m = mean(v); s = sqrt(variance(v))
    sk = mean((v:-m):^3)/s^3
    sv = sort(v,1)
    printf("K1 draws %s  n=%f sd=%9.5f skew=%9.4f  (q975-b)/sd=%8.4f  (b-q025)/sd=%8.4f\n",
      cs[k,2], n, s, sk, (sv[ceil(0.975*n)]-b[1,k])/s, (b[1,k]-sv[ceil(0.025*n)])/s)
  }
end

dis "=== K2: same run with boottstat (studentized band) ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(999) rseed(4) boottstat att esprelength(1) espostlength(1) simultaneousbands dots(0)
dis "K2_RC=" _rc
dis "K2_ENGINE=" "`e(esband_engine)'"
dis "K2_CRIT=" e(esband_crit_att)
local BFAM att
local BTAG K2
include "`WD'/bandcmp.do"

dis "=== K3: bootdisplay(wald) + band (the family the band matches) ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(999) rseed(4) bootdisplay(wald) att esprelength(1) espostlength(1) simultaneousbands dots(0)
dis "K3_RC=" _rc
local BFAM att
local BTAG K3
include "`WD'/bandcmp.do"
dis "P09_DONE"
