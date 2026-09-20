*** LENS 0 follow-up: mechanism check on the bootsave() draw matrix from BLOCK 4
*** (arm A of the SINGLETON route: 9 estimation clusters, 5 singleton phantoms).
clear all
set more off
set processors 2
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-reproduce-lens0"

use "bsv_a.dta", clear
drop if _rep == 0
dis "RESULT MASS NREPS=" _N

egen double realmass = rowtotal(w_1 w_2 w_3 w_4 w_5 w_6 w_7 w_8 w_9)
egen double phanmass = rowtotal(w_101 w_102 w_103 w_104 w_105)
qui sum realmass
dis "RESULT MASS REAL mean=" r(mean) " sd=" r(sd) " min=" r(min) " max=" r(max)
qui sum phanmass
dis "RESULT MASS PHANTOM mean=" r(mean) " sd=" r(sd) " min=" r(min) " max=" r(max)
qui count if phanmass > 0
dis "RESULT MASS REPSWITHPHANTOM=" r(N) " of " _N

foreach v of varlist w_* {
  qui sum `v'
  dis "RESULT MASS COL `v' mean=" r(mean) " min=" r(min) " max=" r(max)
}

* how many of the NINE real clusters a replication actually holds
gen byte nreal = 0
foreach v of varlist w_1 w_2 w_3 w_4 w_5 w_6 w_7 w_8 w_9 {
  qui replace nreal = nreal + (`v' > 0)
}
qui sum nreal
dis "RESULT MASS NREALDRAWN mean=" r(mean) " min=" r(min) " max=" r(max)
qui tab nreal
dis "===== MASS DONE ====="
