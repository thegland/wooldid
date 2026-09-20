*** LENS 0 (independent reproduction) of BR-4.
*** Built from the CLAIM ALONE. Own DGP, own geometry, own routes.
*** Falsification stated before running:
***   REFUTED  if e(boot_df) == e(N_clust)-1 in BOTH arms, or both cautions fire in both arms.
***   CONFIRMED if e(N_clust) is EQUAL across arms and e(boot_df) is NOT.
***   PARTIAL   if e(N_clust) itself moves.

clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

dis "BUILDCHECK"
which wooldid

cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-reproduce-lens0"

*** ------------------------------------------------------------------
*** DGP: 9 real counties, 2010-2017. 3 treated 2013, 2 treated 2015, 4 never.
***      -> 9 resampling clusters, 5 treated clusters: BOTH documented
***         cautions (<10 total, <6 treated) must fire on the clean file.
*** ------------------------------------------------------------------
prog def mkreal
  clear
  set obs 9
  gen int id = _n
  gen int onset = .
  replace onset = 2013 if inrange(id,1,3)
  replace onset = 2015 if inrange(id,4,5)
  gen double ufe = rnormal(0,1)
  expand 8
  bys id: gen int year = 2009 + _n
  gen double x = rnormal(0,1) + 0.2*year/1000
  gen byte post = (year >= onset) & onset < .
  gen double y = 2 + ufe + 0.10*(year-2009) + 0.35*x + 0.50*post + rnormal(0,0.5)
  gen byte phantom = 0
end

*** ------------------------------------------------------------------
*** Reporting helper: dump the deterministic discriminators.
*** ------------------------------------------------------------------
prog def dumpit
  syntax , tag(string)
  dis "RESULT `tag' N="        e(N)
  dis "RESULT `tag' NCLUST="   e(N_clust)
  dis "RESULT `tag' BOOTDF="   e(boot_df)
  dis "RESULT `tag' DFR="      e(df_r)
  dis "RESULT `tag' DFRAN="    e(df_r_analytic)
  dis "RESULT `tag' BFAIL="    e(boot_failed)
  dis "RESULT `tag' BCOMP="    e(boot_completed)
  cap mat wdb = e(b)
  if _rc == 0 dis "RESULT `tag' ATTHEX=" %21x wdb[1,1]
  cap mat wdv = e(V)
  if _rc == 0 dis "RESULT `tag' VHEX="   %21x wdv[1,1]
end

*** ==================================================================
*** BLOCK 1 -- SINGLETON ROUTE (headline).
***   5 extra counties observed in exactly ONE year each (2 treated 2013,
***   3 never). reghdfe drops singletons, so they cannot be in e(sample).
*** ==================================================================
dis "===== BLOCK1 SINGLETON ====="
mkreal
* five singleton units, one row each, in well-populated years
local k = _N
set obs `=_N+5'
replace id      = 100 + (_n - `k')   in `=`k'+1'/`=`k'+5'
replace year    = 2011               in `=`k'+1'/`=`k'+2'
replace year    = 2012               in `=`k'+3'/`=`k'+3'
replace year    = 2014               in `=`k'+4'/`=`k'+5'
replace onset   = 2013               in `=`k'+1'/`=`k'+2'
replace onset   = .                  in `=`k'+3'/`=`k'+5'
replace x       = rnormal(0,1)       in `=`k'+1'/`=`k'+5'
replace y       = 2 + rnormal(0,1)   in `=`k'+1'/`=`k'+5'
replace phantom = 1                  in `=`k'+1'/`=`k'+5'
save "b1.dta", replace

dis "=== B1 ARM A (singletons LEFT IN) ==="
use "b1.dta", clear
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(60) dots(0) rseed(777)
dis "B1A_RC=" _rc
dumpit, tag(B1A)

dis "=== B1 ARM B (singletons DROPPED) ==="
use "b1.dta", clear
drop if phantom == 1
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(60) dots(0) rseed(777)
dis "B1B_RC=" _rc
dumpit, tag(B1B)

*** ==================================================================
*** BLOCK 2 -- ALL-MISSING-OUTCOME ROUTE (the original claim's route).
***   5 extra counties, full 8-year panel, y missing on every row.
*** ==================================================================
dis "===== BLOCK2 MISSY ====="
mkreal
local k = _N
expand 1
set obs `=_N+40'
local j = `k'
forvalues u = 1/5 {
  forvalues yy = 2010/2017 {
    local ++j
    replace id      = 100 + `u' in `j'
    replace year    = `yy'      in `j'
    replace onset   = cond(`u' <= 2, 2013, .) in `j'
    replace x       = rnormal(0,1) in `j'
    replace y       = .         in `j'
    replace phantom = 1         in `j'
  }
}
save "b2.dta", replace

dis "=== B2 ARM A (missing-y units LEFT IN) ==="
use "b2.dta", clear
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(60) dots(0) rseed(777)
dis "B2A_RC=" _rc
dumpit, tag(B2A)

dis "=== B2 ARM B (missing-y units DROPPED) ==="
use "b2.dta", clear
drop if phantom == 1
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(60) dots(0) rseed(777)
dis "B2B_RC=" _rc
dumpit, tag(B2B)

*** ==================================================================
*** BLOCK 3 -- MISSING-CONTROL ROUTE.
***   Same 5 extra counties, y fully present, controls(x) missing on every row.
*** ==================================================================
dis "===== BLOCK3 MISSX ====="
use "b2.dta", clear
replace y = 2 + rnormal(0,1) + 0.1*(year-2009) if phantom == 1
replace x = .                                   if phantom == 1
save "b3.dta", replace

dis "=== B3 ARM A (missing-x units LEFT IN) ==="
use "b3.dta", clear
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(60) dots(0) rseed(777)
dis "B3A_RC=" _rc
dumpit, tag(B3A)

dis "=== B3 ARM B (missing-x units DROPPED) ==="
use "b3.dta", clear
drop if phantom == 1
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(60) dots(0) rseed(777)
dis "B3B_RC=" _rc
dumpit, tag(B3B)

*** ==================================================================
*** BLOCK 4 -- MECHANISM: bootsave() on an arm-A run.
***   A correct cluster bootstrap over G estimation clusters writes exactly
***   G w_ columns and holds the total drawn mass over those G at exactly G
***   in every replication.
*** ==================================================================
dis "===== BLOCK4 BOOTSAVE ====="
use "b1.dta", clear
cap noi wooldid y id year onset, att cluster(id) controls(x) bootstrap(standard) reps(120) dots(0) rseed(4242) bootsave(bsv_a) replace
dis "B4A_RC=" _rc
dis "RESULT B4A NCLUST=" e(N_clust)
dis "RESULT B4A BOOTDF=" e(boot_df)

use "bsv_a.dta", clear
desc, varlist
local wl : list clean "`r(varlist)'"
local nw 0
foreach v of local wl {
  if substr("`v'",1,2) == "w_" local ++nw
}
dis "RESULT B4A WCOLS=`nw'"
dis "RESULT B4A WNAMES: `wl'"
* drawn mass over the NINE real clusters (w_1..w_9), per replication
cap confirm variable w_1
if _rc == 0 {
  egen double realmass = rowtotal(w_1 w_2 w_3 w_4 w_5 w_6 w_7 w_8 w_9)
  qui sum realmass
  dis "RESULT B4A REALMASS mean=" r(mean) " sd=" r(sd) " min=" r(min) " max=" r(max)
  foreach c of local wl {
    if substr("`c'",1,2) == "w_" {
      qui sum `c'
      dis "RESULT B4A COL `c' mean=" r(mean) " max=" r(max) " sd=" r(sd)
    }
  }
}

dis "===== LENS0 DONE ====="
