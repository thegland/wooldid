* ---------------------------------------------------------------------------------------
* probe_determ.do — deterministic (non-MC) checks around the finer-than-i cluster machinery.
* One DGP-S draw; everything wrapped in capture noisily with sentinels.
* ---------------------------------------------------------------------------------------
clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 424242
set sortseed 424242

qui set obs 40
gen double state = _n
gen double s_g   = rnormal()*0.8
gen double ttre  = cond(state<=10, 5, cond(state<=20, 8, .))
qui expand 5
bysort state: gen double upos = _n
gen double unit = (state-1)*5 + upos
gen double a_u  = rnormal()*0.5
qui expand 12
bysort unit: gen double t = _n
tempvar first vv
bysort state t: gen byte `first' = (_n==1)
gen double `vv' = rnormal()*0.4 if `first'
bysort state t (`first'): replace `vv' = `vv'[_N]
gen double v_gt = `vv'
tempvar eps
gen double `eps' = rnormal()*0.4
bysort unit (t): gen double e = `eps'[1]
by unit (t): replace e = 0.5*e[_n-1] + sqrt(1-0.25)*`eps'[_n] if _n > 1
gen double rel = t - ttre
gen double tau = 0
qui replace tau = cond(ttre==5, 1.0, 0.4) + 0.05*rel if ttre<. & t>=ttre
gen double y = s_g + v_gt + a_u + 0.2*t + 0.5*cos(t) + tau + e
gen double coh = cond(ttre<., ttre, 99)
gen double g   = cond(ttre<., ttre, 0)
* a cluster variable that CROSS-CUTS coh: region cuts across states inside each cohort
gen double region = mod(state, 8) + 1
* a cluster variable with missing values on 3 states (realistic: unmatched administrative ids)
gen double stmiss = state
qui replace stmiss = . if inlist(state, 4, 15, 33)
tempfile base
qui save `base'

* ---- D1: nofastcollapse must not move any SE when cluster() is finer than i -------------
use `base', clear
cap noisily qui wooldid y coh t ttre, att att_i att_it cluster(state) esprelength(3) espostlength(4) esfixedbaseperiod
dis "D1A_RC=" _rc
matrix D1A_b = e(b)
matrix D1A_V = e(V)
use `base', clear
cap noisily qui wooldid y coh t ttre, att att_i att_it cluster(state) esprelength(3) espostlength(4) esfixedbaseperiod nofastcollapse
dis "D1B_RC=" _rc
matrix D1B_b = e(b)
matrix D1B_V = e(V)
dis "D1_BRELDIF=" %21x mreldif(D1A_b, D1B_b)
dis "D1_VRELDIF=" %21x mreldif(D1A_V, D1B_V)
dis "D1_BRELDIF_dec=" mreldif(D1A_b, D1B_b)
dis "D1_VRELDIF_dec=" mreldif(D1A_V, D1B_V)

* ---- D2: two-way cluster() with a coarse i ---------------------------------------------
use `base', clear
cap noisily wooldid y coh t ttre, att cluster(state region) esprelength(3) espostlength(4) esfixedbaseperiod
dis "D2_RC=" _rc
if _rc == 0 {
  dis "D2_NCLUST=" e(N_clust)
  dis "D2_DFR=" e(df_r)
  matrix list e(b)
}

* ---- D3: the jointtests caution across the three geometries ----------------------------
* i nested in cluster (at/above i)  -> caution EXPECTED
use `base', clear
cap noisily wooldid y state t ttre, att cluster(state) jointtests esprelength(3) espostlength(4) esfixedbaseperiod
dis "D3_ATI_RC=" _rc
* cluster strictly finer than i     -> caution must NOT fire
use `base', clear
cap noisily wooldid y coh t ttre, att cluster(state) jointtests esprelength(3) espostlength(4) esfixedbaseperiod
dis "D3_FINER_RC=" _rc
* cross-cutting cluster             -> caution must NOT fire
use `base', clear
cap noisily wooldid y coh t ttre, att cluster(region) jointtests esprelength(3) espostlength(4) esfixedbaseperiod
dis "D3_CROSS_RC=" _rc

* ---- D4: cluster variable with missing values ------------------------------------------
use `base', clear
cap noisily wooldid y coh t ttre, att cluster(stmiss) esprelength(3) espostlength(4) esfixedbaseperiod
dis "D4_RC=" _rc
if _rc == 0 {
  dis "D4_NCLUST=" e(N_clust)
  dis "D4_N=" e(N)
  dis "D4_NDROP=" e(N_dropped)
}

* ---- D5: b must not move with the VCE choice -------------------------------------------
use `base', clear
cap noisily qui wooldid y coh t ttre, att cluster(state) esprelength(3) espostlength(4) esfixedbaseperiod
matrix D5A = e(b)
use `base', clear
cap noisily qui wooldid y coh t ttre, att robust esprelength(3) espostlength(4) esfixedbaseperiod
matrix D5B = e(b)
use `base', clear
cap noisily qui wooldid y coh t ttre, att cluster(unit) esprelength(3) espostlength(4) esfixedbaseperiod
matrix D5C = e(b)
dis "D5_cluststate_vs_robust=" mreldif(D5A, D5B)
dis "D5_cluststate_vs_clustunit=" mreldif(D5A, D5C)

* ---- D6: e(N_clust) against the true count ---------------------------------------------
use `base', clear
qui levelsof state, local(LS)
dis "D6_TRUE_NSTATE=" wordcount("`LS'")
qui levelsof unit, local(LU)
dis "D6_TRUE_NUNIT=" wordcount("`LU'")
cap noisily qui wooldid y coh t ttre, att cluster(unit) esprelength(3) espostlength(4) esfixedbaseperiod
dis "D6_NCLUST_unit=" e(N_clust)
dis "D6_DFR_unit=" e(df_r)

* ---- D7: does a cluster FINER than i change the point estimate at all? -----------------
* (already covered per-draw in the MC; this prints the hex so a skeptic can eyeball it)
dis "D7_att_clustunit=" %21x D5C[1,1]
dis "D7_att_cluststate=" %21x D5A[1,1]
dis "D7_att_robust=" %21x D5B[1,1]
dis "PROBE_DETERM_DONE"
