* =================================================================
* probe_br4.do -- SKEPTIC repro of BR-4.
* Claim under test: bootstrap() builds its resampling universe from the
* user-if sample, NOT e(sample), so whole units that contribute nothing
* to the estimation still count as resampling clusters -- inflating
* e(boot_df) (= G-1) and silencing the <10-cluster and <6-treated-cluster
* cautions, on two runs with a BIT-IDENTICAL estimation sample.
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_br4.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909090
set sortseed 909090

* ---------- BLOCK 1: all-missing-outcome units ----------
* 9 REAL units (2 treated 2004, 2 treated 2006, 5 never) + 6 PHANTOM units
* (4 treated 2004, 2 never) whose outcome is missing on every row.
set obs 15
gen long unit = _n
gen byte phantom = unit > 9
gen int ttre = .
replace ttre = 2004 if inlist(unit,1,2)
replace ttre = 2006 if inlist(unit,3,4)
replace ttre = 2004 if inlist(unit,10,11,12,13)
gen double ufe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.6, 0))
gen double y  = 1 + ufe + 0.1*(year-1999) + te + rnormal()*0.6
replace y = .  if phantom            /* outcome entirely missing for phantoms */
save "sk_b4_A.dta", replace
drop if phantom
save "sk_b4_B.dta", replace

foreach arm in A B {
  use "sk_b4_`arm'.dta", clear
  dis "########## BLOCK1 ARM `arm' ##########"
  cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(40) rseed(31337) dots(0)
  dis "B4_`arm'_RC=" _rc
  dis "B4_`arm'_N="       e(N)
  dis "B4_`arm'_NCLUST="  e(N_clust)
  dis "B4_`arm'_BOOTDF="  e(boot_df)
  dis "B4_`arm'_DFR="     e(df_r)
  dis "B4_`arm'_ATT="     %21x el(e(b),1,1)
  dis "B4_`arm'_SE="      %21x sqrt(el(e(V),1,1))
}

* ---------- BLOCK 2: which clusters actually get drawn? ----------
* 30 REAL + 10 PHANTOM, bootsave() the multiplicity matrix and read it back.
clear
set seed 171717
set obs 40
gen long unit = _n
gen byte phantom = unit > 30
gen int ttre = .
replace ttre = 2004 if unit <= 10
replace ttre = 2006 if unit > 10 & unit <= 20
replace ttre = 2004 if unit > 30 & unit <= 35
gen double ufe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double te = cond(missing(ttre), 0, cond(year >= ttre, 0.6, 0))
gen double y  = 1 + ufe + 0.1*(year-1999) + te + rnormal()*0.6
replace y = . if phantom
save "sk_b4_W.dta", replace

cap erase "skw.dta"
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(200) rseed(24680) bootsave(skw) dots(0)
dis "B4_W_RC=" _rc
dis "B4_W_N="      e(N)
dis "B4_W_NCLUST=" e(N_clust)
dis "B4_W_BOOTDF=" e(boot_df)

use "skw.dta", clear
dis "B4_W_FILEOBS=" _N
qui ds w_*
local wcols `r(varlist)'
local nw : word count `wcols'
dis "B4_W_NWCOLS=" `nw'
* real units are 1..30, phantoms 31..40
qui egen double realmass = rowtotal(w_1-w_30) if _n > 1
qui egen double phanmass = rowtotal(w_31-w_40) if _n > 1
qui sum realmass
dis "B4_W_REALMASS_MEAN=" r(mean) "  SD=" r(sd) "  MIN=" r(min) "  MAX=" r(max)
qui sum phanmass
dis "B4_W_PHANMASS_MEAN=" r(mean) "  SD=" r(sd) "  MIN=" r(min) "  MAX=" r(max)
qui sum w_35
dis "B4_W_PHANTOMUNIT35_MEAN=" r(mean) "  MAX=" r(max)

dis "PROBE_BR4_DONE"
