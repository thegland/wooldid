*=========================================================================
* p1_quick.do -- SKEPTIC probes, part 1
*   (A) mcci-03: does requesting esprelength()/espostlength() move the att?
*   (B) mcci-01: one-dataset discriminator for the impvalues SE.
*       impvalues cluster(unit)  [cluster FINER than i]  vs
*       impvalues cluster(cg)    [cluster == i, textbook IM] vs
*       analytic   cluster(unit)
* Own DGP: unbalanced panel, AR(1) errors, within-cohort effect
* heterogeneity, unit random intercepts NOT absorbed by i=cg.
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p1_quick.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 771103
set sortseed 771103

*---------------- own generator ------------------------------------------
cap prog drop mkdat
prog def mkdat
  * args: nunit  seed  unbal(0/1)
  args nunit sd unbal
  clear
  set seed `sd'
  set sortseed `sd'
  set obs `nunit'
  gen long unit = _n
  * 3 treated cohorts + never-treated, roughly equal
  gen byte cg = mod(_n-1,4)+1
  gen double ttre = .
  replace ttre = 2004 if cg==1
  replace ttre = 2006 if cg==2
  replace ttre = 2008 if cg==3
  gen double ui   = rnormal(0,0.70)          // unit random intercept
  gen double hmul = 1 + rnormal(0,0.25)      // WITHIN-cohort effect heterogeneity
  expand 12
  bys unit: gen int t = 2000 + _n - 1
  * AR(1) idiosyncratic error, rho = 0.4
  gen double e = .
  bys unit (t): replace e = rnormal(0,0.55) if _n==1
  bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
  * true effects: 1.00 / 1.15 / 0.85 times a unit multiplier
  gen double truetau = 0
  replace truetau = 1.00*hmul if cg==1 & t>=ttre
  replace truetau = 1.15*hmul if cg==2 & t>=ttre
  replace truetau = 0.85*hmul if cg==3 & t>=ttre
  gen double y = ui + 0.25*(t-2000) + truetau + e
  if `unbal' == 1 {
    * administrative-shape unbalancedness: late entry / early exit / holes
    gen double r1 = runiform()
    gen double r2 = runiform()
    bys unit (t): gen byte entrylag = floor(3*r1[1]) if r1[1] < .35
    bys unit (t): replace entrylag = 0 if missing(entrylag)
    bys unit (t): replace entrylag = entrylag[1]
    drop if (t - 2000) < entrylag
    drop if r2 < 0.03
  }
end

*========================= (A) mcci-03 ===================================
di as txt _n "########## BLOCK A : esprelength moves the att ? ##########"
mkdat 60 4411 0
tempfile dA
qui save "`dA'"

foreach spec in none p1 p2 p3 {
  use "`dA'", clear
  local opt ""
  if "`spec'"=="p1" local opt "esprelength(1) espostlength(1)"
  if "`spec'"=="p2" local opt "esprelength(2) espostlength(2)"
  if "`spec'"=="p3" local opt "esprelength(3) espostlength(3)"
  di as txt _n "=== ESOPT spec=`spec'  opts: `opt'"
  cap noi wooldid y cg t ttre, att cluster(unit) `opt'
  di "ESOPT_RC_`spec'=" _rc
  if _rc==0 {
    mat T = e(wdidmainresults)
    di "ESOPT_`spec'_b="       %21x T[1,1]
    di "ESOPT_`spec'_se="      %21x T[1,2]
    di "ESOPT_`spec'_b_d="     %12.9f T[1,1]
    di "ESOPT_`spec'_se_d="    %12.9f T[1,2]
    di "ESOPT_`spec'_N="       e(N)
    di "ESOPT_`spec'_dfr="     e(df_r)
    cap di "ESOPT_`spec'_refN=" e(escomp_refN)
    cap di "ESOPT_`spec'_refW=" e(escomp_refW)
  }
}

*========================= (B) mcci-01 ===================================
di as txt _n "########## BLOCK B : impvalues SE, one dataset ##########"
* balanced version first (identical within-cohort weight profiles)
mkdat 60 9021 0
tempfile dB
qui save "`dB'"

di as txt _n "=== B1 impvalues, cluster(unit)  [cluster FINER than i=cg] ==="
use "`dB'", clear
cap noi wooldid y cg t ttre, att cluster(unit) impvalues level(95)
di "B1_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B1_b="   %12.9f T[1,1]
  di "B1_se="  %12.9f T[1,2]
  di "B1_dfr=" e(df_r)
}

di as txt _n "=== B2 impvalues, cluster(cg)  [cluster == i, textbook IM] ==="
use "`dB'", clear
cap noi wooldid y cg t ttre, att cluster(cg) impvalues level(95)
di "B2_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B2_b="   %12.9f T[1,1]
  di "B2_se="  %12.9f T[1,2]
  di "B2_dfr=" e(df_r)
}

di as txt _n "=== B3 analytic, cluster(unit) ==="
use "`dB'", clear
cap noi wooldid y cg t ttre, att cluster(unit) level(95)
di "B3_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B3_b="   %12.9f T[1,1]
  di "B3_se="  %12.9f T[1,2]
  di "B3_dfr=" e(df_r)
}

di as txt _n "=== B4 analytic, robust ==="
use "`dB'", clear
cap noi wooldid y cg t ttre, att robust level(95)
di "B4_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B4_b="   %12.9f T[1,1]
  di "B4_se="  %12.9f T[1,2]
}

*--- UNBALANCED: does the degeneracy scrub have any chance of firing? -----
di as txt _n "=== B5 impvalues cluster(unit), UNBALANCED admin-shape panel ==="
mkdat 60 3377 1
qui count
di "B5_rows=" r(N)
cap noi wooldid y cg t ttre, att cluster(unit) impvalues level(95)
di "B5_RC=" _rc
if _rc==0 {
  mat T = e(wdidmainresults)
  di "B5_b="   %12.9f T[1,1]
  di "B5_se="  %12.9f T[1,2]
  di "B5_dfr=" e(df_r)
}

*--- direct look at the per-cluster reads margins produces -----------------
di as txt _n "=== B6 the per-cluster reads themselves (balanced data) ==="
use "`dB'", clear
qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
di "B6 done, RC=" _rc

di as txt _n "ALLDONE"
