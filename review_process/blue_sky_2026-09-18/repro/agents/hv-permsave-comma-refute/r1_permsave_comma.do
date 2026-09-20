* r1_permsave_comma.do - SKEPTIC's independent repro of hv-permsave-comma-01
* Self-contained: generates its own panel, works in its own empty subdirectory.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242424
set sortseed 4242424

cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w1"
* start from an EMPTY directory so every file below is ours
cap erase "pd.dta"
cap shell rm -f *.dta

which wooldid
dis "BUILDPATH_CHECKED"

* ---------------- my own panel: 90 firms, 14 years, 3 onset cohorts + never-treated
clear
set obs 90
gen long unit = _n
gen byte grp = mod(_n, 5)
gen ttre = .
replace ttre = 2005 if grp == 1
replace ttre = 2008 if grp == 2
replace ttre = 2011 if grp == 3
* grp 0 and 4 never treated
gen double ue = rnormal()*0.8
expand 14
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + 0.4*(year - ttre + 1) / 3 if year >= ttre & ttre < .
compress
count
save "panel.dta", replace

* ============ T-A: permsave with a comma, Stata-idiomatic saving() form
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(pd, replace)
dis "TA_RC=" _rc
shell ls -1 > ls_A.txt
type ls_A.txt

* ============ T-B: the twin, bootsave with a comma
capture noisily wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(bd, replace)
dis "TB_RC=" _rc
shell ls -1 > ls_B.txt
type ls_B.txt

* ============ T-C: repeat T-A verbatim (does it 602 telling me to add replace?)
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(pd, replace)
dis "TC_RC=" _rc

* ============ T-D: the CORRECT wooldid form
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(pd) replace
dis "TD_RC=" _rc
capture confirm file "pd.dta"
dis "TD_pd_dta_exists_rc=" _rc

* ============ T-E: can the weirdly named file be opened, and is it real?
capture noisily use "pd, replace.dta", clear
dis "TE_RC=" _rc
if _rc == 0 {
  describe, short
  dis "TE_N=" _N
  ds
  dis "TE_VARS=" r(varlist)
}

* ============ T-F: what does a plain `use pd` do after T-A alone?
capture noisily use pd, clear
dis "TF_RC=" _rc

shell ls -1 > ls_final.txt
type ls_final.txt
dis "ALLDONE"
