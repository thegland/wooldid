* SELF-CONTAINED REPRO — fw-weights/FW-03
* With [pw=...] and no vce option, wooldid's e(V) IS the robust covariance (bitwise equal
* to the same run with `robust'), but the run is LABELLED as if it were OLS:
*   header    "Standard Errors: "   (blank, the unweighted-iid look)
*   e(vce)    "ols"
*   e(vcetype) not posted
* ado:4283 states the contract this breaks in its own words: "e(vce) and e(vcetype) name
* what e(V) IS".  ado:4301's fallthrough,
*     else ereturn local vce = cond("`poisson'" == "poisson", "oim", "ols")
* never looks at the weight type.  reghdfe, the estimator underneath, gets it right on the
* identical fit and posts e(vce)="robust", e(vcetype)="Robust".
* Downstream reporting tools (estout/esttab, outreg2) read e(vcetype) to decide whether to
* print "Robust standard errors", so a pweighted wooldid table is labelled OLS.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 6
replace ttre = 5 if i > 6  & i <= 12
replace ttre = 6 if i > 12 & i <= 18
gen double ufe = rnormal()
expand 8
bysort i: gen t = _n
expand 3
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ufe + 0.15*t + 1.25*post + rnormal()*0.5
gen double svywt = 500 + int(runiform()*4500)     /* an ordinary survey weight */

dis "===== A. [pw=svywt], no vce option ====="
wooldid y i t ttre [pw=svywt], att
scalar se_pw = sqrt(e(V)[1,1])
dis "A: se = " %14.10f se_pw
dis "A: e(vce) = [" e(vce) "]   e(vcetype) = [" e(vcetype) "]"

dis "===== B. [aw=svywt], robust  (the covariance A actually computed) ====="
wooldid y i t ttre [aw=svywt], att robust
scalar se_awr = sqrt(e(V)[1,1])
dis "B: se = " %14.10f se_awr
dis "B: e(vce) = [" e(vce) "]   e(vcetype) = [" e(vcetype) "]"

dis "===== C. [aw=svywt], default iid (what A's label claims) ====="
wooldid y i t ttre [aw=svywt], att
scalar se_aw = sqrt(e(V)[1,1])
dis "C: se = " %14.10f se_aw
dis "C: e(vce) = [" e(vce) "]   e(vcetype) = [" e(vcetype) "]"

dis "===== verdict ====="
dis "se(pw) == se(aw,robust) ?  " (se_pw == se_awr)   "   |diff| = " %12.4e abs(se_pw-se_awr)
dis "se(pw) == se(aw,iid)    ?  " (se_pw == se_aw)    "   |diff| = " %12.4e abs(se_pw-se_aw)
dis "so e(V) is ROBUST, and e(vce) says [" "ols" "]"

dis "===== attribution: reghdfe on the identical pweighted fit ====="
qui reghdfe y post [pw=svywt], absorb(i t)
dis "reghdfe pw : e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"
qui reghdfe y post [aw=svywt], absorb(i t)
dis "reghdfe aw : e(vce)=[" e(vce) "]  e(vcetype)=[" e(vcetype) "]"
dis "REPRO_DONE"
