* r2_dots.do - SKEPTIC's independent repro of hv-permsave-comma-03 (dots() validation)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909090
set sortseed 909090
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w2"
which wooldid

clear
set obs 90
gen long unit = _n
gen byte grp = mod(_n, 5)
gen ttre = .
replace ttre = 2005 if grp == 1
replace ttre = 2008 if grp == 2
replace ttre = 2011 if grp == 3
gen double ue = rnormal()*0.8
expand 14
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + 0.4*(year - ttre + 1)/3 if year >= ttre & ttre < .

* A: dots(-5) under bootstrap()
capture noisily wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(20) dots(-5)
dis "A_boot_dotsneg5_RC=" _rc

* B: dots(-5) under permute()
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) dots(-5)
dis "B_perm_dotsneg5_RC=" _rc

* C: dots(0) under permute()  (legal: suppress)
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) dots(0)
dis "C_perm_dots0_RC=" _rc

* D: dots(5) under permute()  (legal: show)
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) dots(5)
dis "D_perm_dots5_RC=" _rc

* E [NEW]: dots(-1) is the internal "not given" sentinel. Under bootstrap() does the
*          user's explicit dots(-1) get treated as "not given" (i.e. default dots(1))?
capture noisily wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(20) dots(-1)
dis "E_boot_dotsneg1_RC=" _rc

* F [NEW]: no engine at all + dots(0) -> stray-companion screen should refuse
capture noisily wooldid y unit year ttre, att cluster(unit) dots(0)
dis "F_noengine_dots0_RC=" _rc

* G [NEW]: no engine at all + dots(-1) -> sentinel means the stray screen cannot see it
capture noisily wooldid y unit year ttre, att cluster(unit) dots(-1)
dis "G_noengine_dotsneg1_RC=" _rc

* H [NEW]: no engine + dots(7)
capture noisily wooldid y unit year ttre, att cluster(unit) dots(7)
dis "H_noengine_dots7_RC=" _rc

* I [NEW]: rseed() range check IS hoisted - confirm it fires under permute()
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(-3)
dis "I_perm_badseed_RC=" _rc

dis "ALLDONE"
