* r13_disc.do -- SELF-CONTAINED.  Two discriminators.
* (A) BX-N3: is the exact rank-1 deficiency of the analytic pre block on design A the
*     SINGLETON-TREATED-CELL degeneracy (i = county means one cohort per unit) reaching
*     `robust', or something else?  Coarsen i to the onset cohort with fe(county) so the
*     treated cells are no longer singletons, and look at the same eigenvalues.
* (B) BOOTX-4: the prober's exact geometry (40 states, 12 periods, onsets 5/7/9 with ten
*     states each and ten never-treated, timetrends, cluster(state)) rebuilt with MY
*     generator and MY seed.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

dis "======== (A) singleton-cell discriminator"
clear
set obs 30
gen state = _n
gen onset = .
replace onset = 7 + mod(_n, 6) if mod(_n,3)==0
expand 2
bysort state: gen cty = _n
gen county = state*10 + cty
expand 16
bysort county: gen t = _n
bysort county (t): gen double ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
drop rel te
gen coh = cond(onset<., onset, 0)
tab coh, m

local OPT "esfixedbaseperiod esprelength(6) espostlength(4)"
qui wooldid y county t onset, robust `OPT'
matrix A_rob_i_is_county = e(V_es_att)[1..6,1..6]
qui wooldid y coh t onset, fe(county) robust `OPT'
dis "COH_ROB_RC=" _rc
matrix A_rob_i_is_coh = e(V_es_att)[1..6,1..6]
qui wooldid y coh t onset, fe(county) cluster(state) `OPT'
dis "COH_CLS_RC=" _rc
matrix A_cls_i_is_coh = e(V_es_att)[1..6,1..6]
mata:
  nmv = ("A_rob_i_is_county","A_rob_i_is_coh","A_cls_i_is_coh")
  for (q=1; q<=3; q++) {
    S = st_matrix(nmv[q]); d = sqrt(diagonal(S))
    printf("\nDISC| %s\n", nmv[q])
    printf("DISC|   eigenvalues (ascending): ")
    ev = sort(symeigenvalues(makesymmetric(S))', 1)'
    for (z=1; z<=cols(ev); z++) printf("%12.4e ", ev[z])
    printf("\nDISC|   max off-diagonal |corr| = %12.8f\n", max(abs(select(vec(S:/(d*d')), vec(I(rows(S))):==0))))
  }
end

dis "======== (B) BOOTX-4 geometry, my generator"
clear
set obs 40
gen state = _n
gen onset = .
replace onset = 5 if state <= 10
replace onset = 7 if state > 10 & state <= 20
replace onset = 9 if state > 20 & state <= 30
expand 12
bysort state: gen t = _n
bysort state (t): gen double sfe = rnormal(0,1) if _n==1
bysort state (t): replace sfe = sfe[_n-1] if _n>1
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .30 + .05*rel, 0)
gen double y = 1 + sfe + .05*t + te + rnormal(0,.8)
drop rel te
local O "cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) timetrends"
cap noi wooldid y state t onset, `O'
dis "TT4_AN_RC=" _rc " VCV=[`e(wdid_vcvstate)']"
cap confirm matrix e(V)
dis "TT4_AN_HASV_RC=" _rc
cap noi wooldid y state t onset, `O' bootstrap(standard) reps(40) rseed(4242)
dis "TT4_BO_RC=" _rc " VCV=[`e(wdid_vcvstate)'] VPOSTED=" e(boot_v_posted) " COMPLETED=" e(boot_completed)
cap confirm matrix e(V)
dis "TT4_BO_HASV_RC=" _rc
cap confirm matrix e(V_boot)
dis "TT4_BO_HASVBOOT_RC=" _rc
cap confirm matrix e(V_analytic)
dis "TT4_BO_HASVANALYTIC_RC=" _rc
cap noi test ES_att_post0 = 0
dis "TT4_BO_TEST_RC=" _rc
dis "DONE"
