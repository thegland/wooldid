* p08_more.do -- further permexact probes: multi-way cluster(), the design sample under
* `if', level() independence, the few-treated-cluster caution, and permtci.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

*================================================================ F1-F3 two-way cluster()
* i = unit, state nests unit, onset constant within state.  cluster(state t) is Stata's
* ordinary two-way clustering and wooldid's parser accepts cluster(varlist numeric).
clear
set obs 12
gen long unit = _n
gen byte state = ceil(unit/3)
gen double ons = .
replace ons = 5 if state == 1
replace ons = 7 if state == 2
expand 10
bysort unit: gen int t = _n
expand 2
gen double y = 0.6*unit + 0.15*t + 0.9*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p08_data, replace

cap noi wooldid y unit t ons, cluster(state t)
dis "F1_RC=" _rc
cap noi wooldid y unit t ons, cluster(state t) permute(shuffle) permexact
dis "F2_RC=" _rc
cap noi wooldid y unit t ons, cluster(state t) permute(shuffle) permreps(30)
dis "F3_RC=" _rc
cap noi wooldid y unit t ons, cluster(state t) bootstrap(standard) reps(20)
dis "F3b_boot_RC=" _rc

*================================================================ F4 design sample under `if'
clear
set obs 8
gen long unit = _n
gen double ons = .
replace ons = 4 if unit == 1
replace ons = 7 if unit == 2
expand 10
bysort unit: gen int t = _n
expand 3
gen double y = 0.6*unit + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p08_data2, replace
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact
dis "F4a_RC=" _rc " NPOSS=" e(perm_npossible) "  (8 clusters -> expect 56)"
cap noi wooldid y unit t ons if unit <= 6, cluster(unit) permute(shuffle) permexact
dis "F4b_RC=" _rc " NPOSS=" e(perm_npossible) "  (if unit<=6 -> expect 30)"

*================================================================ F5 level() independence
use p08_data2, clear
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact level(90)
dis "F5a_RC=" _rc
mat L90p = e(perm_p)
mat L90c = e(perm_ci)
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact level(99)
dis "F5b_RC=" _rc
mat L99p = e(perm_p)
mat L99c = e(perm_ci)
mat li L90p
mat li L99p
mat li L90c
mat li L99c

*================================================================ F6 few-treated caution
use p08_data2, clear
replace ons = . if unit == 2
cap noi wooldid y unit t ons, cluster(unit) permute(shuffle) permexact
dis "F6_RC=" _rc " NPOSS=" e(perm_npossible) " p="
mat li e(perm_p)

*================================================================ F7 permtci on the p01 design
use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permtstat permtci
dis "F7_RC=" _rc
dis "F7 perm_tlb_att=" %21x e(perm_tlb_att) "  perm_tub_att=" %21x e(perm_tub_att)
dis "F7 perm_tlb_att=" e(perm_tlb_att) "  perm_tub_att=" e(perm_tub_att)
dis "F7 b=" %21x _b[main_att] " se=" %21x _se[main_att]
mat li e(perm_p)
mat li e(perm_ci)

dis "ALLDONE"
