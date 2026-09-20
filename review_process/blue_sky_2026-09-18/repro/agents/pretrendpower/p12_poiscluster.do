* ================================================================
* p12_poiscluster.do -- does the poisson path fall into the
*   chi-squared shortcut EVEN WHEN cluster() leaves few degrees of
*   freedom?  wooldid_ptp takes the chi2 branch whenever e(df_r) is
*   missing (ado:10357-10366, 10545-10548); the code's own comment
*   says "a chi-squared shortcut is optimistic by ~25 percent at 10
*   clusters".  Compare the same design under reghdfe and ppmlhdfe.
* ================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

clear
set obs 16
gen long st = _n
gen byte coh = 4
replace coh = 1 if st <= 4
replace coh = 2 if st > 4 & st <= 8
replace coh = 3 if st > 8 & st <= 12
gen int ttre = .
replace ttre = 5 if coh == 1
replace ttre = 7 if coh == 2
replace ttre = 9 if coh == 3
expand 10
bysort st: gen int u = _n
expand 12
bysort st u: gen int pid = _n
gen int t = u
sort st t pid
gen byte D = (ttre<. & t>=ttre)
gen double mu = exp(2.5 + 0.03*t + 0.3*D)
gen double y = rpoisson(mu)
gen double ylin = 2.5 + 0.03*t + 0.3*D + rnormal()

dis "TAG_P1 gaussian + cluster(st)"
cap noi wooldid ylin coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "P1_RC=" _rc
dis "P1 df_r = " e(df_r) "  keff = " e(roth_keff_att) "  prep = " e(roth_prep_att)
mat P1 = e(roth_power_att)
dis "P1 mdt_wald = " P1[1,2]

dis "TAG_P2 poisson + cluster(st)"
cap noi wooldid y coh t ttre, poisson fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
dis "P2_RC=" _rc
dis "P2 df_r = " e(df_r) "  N_clust = " e(N_clust) "  keff = " e(roth_keff_att) "  prep = " e(roth_prep_att)
mat P2 = e(roth_power_att)
dis "P2 mdt_wald = " P2[1,2]
dis "P2 chi2-vs-F on the SAME information: F-branch mdt would be " ///
    P2[1,2]*sqrt(npnF(3,15,invF(3,15,.95),.2)/npnchi2(3,invchi2(3,.95),.2))
dis "P2 ratio chi2/F = " sqrt(npnchi2(3,invchi2(3,.95),.2)/npnF(3,15,invF(3,15,.95),.2))
dis "P2 true size of a chi2(3) 5% cutoff when the statistic is F(3,15): " ///
    Ftail(3,15,invchi2(3,.95)/3)
dis "TAG_PEND"
