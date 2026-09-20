clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 71113
set sortseed 71113

which wooldid

* ---------------------------------------------------------------
* L0 DGP: canonical Wooldridge ETWFE shape.
*   i = cg (cohort group) -- the fixed effects are at cohort level
*   cluster(unit) -- BELOW the level of i, which sthlp:251 calls
*   "the intended practice".
*   60 units: 4 groups of 15. cg 1/2/3 treated at t = 5/7/9,
*   cg 4 never treated. T = 1..12.
* ---------------------------------------------------------------
clear
set obs 60
gen unit = _n
gen cg = ceil(_n/15)
gen ttre = .
replace ttre = 5 if cg==1
replace ttre = 7 if cg==2
replace ttre = 9 if cg==3
* unit random intercept (NOT absorbed: the FE is at cg level)
gen double ui = rnormal(0,0.80)
* unit-level effect heterogeneity within cohort
gen double uhet = rnormal(0,0.30)
expand 12
bysort unit: gen t = _n
gen double te = 0.05*t
gen double eps = rnormal(0,0.60)
* true cohort effects
gen double effc = .
replace effc = 1.00 if cg==1
replace effc = 1.20 if cg==2
replace effc = 0.80 if cg==3
replace effc = 0    if cg==4
gen byte post = !missing(ttre) & t>=ttre
gen double y = 2 + ui + te + post*(effc + uhet) + eps

sum y
tab cg, sum(ttre)

dis "=== ARM A: impvalues, cluster(unit) [cluster FINER than i] ==="
capture noisily wooldid y cg t ttre, att cluster(unit) impvalues
dis "A_RC=" _rc
capture scalar bA = e(b)[1,1]
scalar vA = e(V)[1,1]
scalar seA = sqrt(vA)
scalar dfminA = e(imdf_min)
scalar dfmaxA = e(imdf_max)
scalar dfrA = e(df_r)
ereturn list
dis "A_b = " %21x bA
dis "A_se = " %21x seA
dis "A_dfmin/dfmax/dfr = " dfminA " " dfmaxA " " dfrA
capture noisily matrix list e(wdidmainresults)
dis "A_mainres_rc=" _rc

dis "=== ARM B: impvalues, cluster(cg) [cluster == i : textbook IM] ==="
capture noisily wooldid y cg t ttre, att cluster(cg) impvalues
dis "B_RC=" _rc
scalar bB = e(b)[1,1]
scalar seB = sqrt(e(V)[1,1])
scalar dfminB = e(imdf_min)
scalar dfmaxB = e(imdf_max)
dis "B_b = " %21x bB
dis "B_se = " %21x seB
dis "B_dfmin/dfmax = " dfminB " " dfmaxB

dis "=== ARM C: analytic, cluster(unit) ==="
capture noisily wooldid y cg t ttre, att cluster(unit)
dis "C_RC=" _rc
scalar bC = e(b)[1,1]
scalar seC = sqrt(e(V)[1,1])
dis "C_b = " %21x bC
dis "C_se = " %21x seC

* ---------------------------------------------------------------
* DETERMINISTIC IDENTITY TEST (no Monte Carlo needed).
* If the 45 per-unit reads are the 3 per-cohort values each repeated
* m=15 times, then
*    se_unit = se_cg * sqrt( k*m*(k-1) / ((G-1)*G) ),  k=3,m=15,G=45
* ---------------------------------------------------------------
scalar k = 3
scalar m = 15
scalar G = k*m
scalar ratio_pred = sqrt(k*m*(k-1)/((G-1)*G))
scalar se_pred = seB*ratio_pred
dis "PRED_ratio  = " %21x ratio_pred "   (" ratio_pred ")"
dis "OBS_ratio   = " %21x (seA/seB) "   (" seA/seB ")"
dis "IDENT_relerr= " abs(seA-se_pred)/se_pred
dis "B_EQUAL_relerr= " abs(bA-bB)/abs(bB)
dis "SE_DEFLATION_FACTOR= " seB/seA
dis "ANALYTIC_over_IMunit= " seC/seA
dis "DONE_EXPLORE"
