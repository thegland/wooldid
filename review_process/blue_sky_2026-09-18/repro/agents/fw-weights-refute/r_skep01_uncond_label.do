* SELF-CONTAINED REPRO for SKEP-01.
* CLAIM: on the unconditionalse path e(vce)/e(vcetype) name the STAGE-1 regression's vce,
* not the covariance e(V) actually is. ado:4291's ladder tests `robust'/`cluster' BEFORE
* `unconditionalse', so:
*   (i) `unconditionalse robust'  -> e(vce)="robust",  e(vcetype)="Robust"
*   (ii) plain `robust'           -> e(vce)="robust",  e(vcetype)="Robust"
* yet (i) and (ii) are DIFFERENT matrices; and
*   (iii) `unconditionalse' with a bare [pw=] (the one call that reaches the fallthrough)
*         -> e(vce)="unconditional"
* with (i) and (iii) being the SAME matrix bitwise. One covariance, two names; two
* covariances, one name. ado:4283 says "e(vce) and e(vcetype) name what e(V) IS".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 24
gen long state = _n
gen ttre = .
replace ttre = 5  if mod(_n,4)==1
replace ttre = 8  if mod(_n,4)==2
replace ttre = 10 if mod(_n,4)==3
gen double as = rnormal()
expand 12
bys state: gen int year = _n
gen double y = as + 0.2*year + rnormal() + 0.4*(year>=ttre & ttre<.)
gen double svywt = 0.5 + runiform()

dis _n "=== (i) [pw=svywt], att, unconditionalse robust ==="
capture noisily wooldid y state year ttre [pw=svywt], att unconditionalse robust
dis "I_RC=" _rc
dis "I  se=" %21x (sqrt(e(V)[1,1])) "  e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"

dis _n "=== (ii) [pw=svywt], att, plain robust (NO unconditionalse) ==="
capture noisily wooldid y state year ttre [pw=svywt], att robust
dis "II_RC=" _rc
dis "II se=" %21x (sqrt(e(V)[1,1])) "  e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"

dis _n "=== (iii) [pw=svywt], att, unconditionalse, no robust/cluster ==="
capture noisily wooldid y state year ttre [pw=svywt], att unconditionalse
dis "III_RC=" _rc
dis "III se=" %21x (sqrt(e(V)[1,1])) "  e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"

dis _n "=== (iv) [aw=svywt], att, unconditionalse robust  (the ordinary aweight user) ==="
capture noisily wooldid y state year ttre [aw=svywt], att unconditionalse robust
dis "IV_RC=" _rc
dis "IV se=" %21x (sqrt(e(V)[1,1])) "  e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"

dis _n "=== (v) [aw=svywt], att, plain robust ==="
capture noisily wooldid y state year ttre [aw=svywt], att robust
dis "V_RC=" _rc
dis "V  se=" %21x (sqrt(e(V)[1,1])) "  e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"
dis "SKEP01_DONE"
