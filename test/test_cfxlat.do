* Targeted suite for the cfxlattice() screen domain narrowing (VAL-G12A/VAL-G12B):
* both content screens now tolerate rows with a missing recorded onset (never-treated
* rows feed no executed CFX subpop), curing false rejections / false cautions, while
* keeping full protection on nonmissing-onset rows.
* Run from test/: "<stata>" -e do test_cfxlat.do -> CFXLAT_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 13579
set sortseed 13579
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop mkcfx
program define mkcfx
  use "data/dta_cont", clear
  qui gen byte clat = mod(i, 3)
end

* ---- L1: bad value confined to one never-treated row -> rc 0 with the CFX block
* actually executing (was r411 at parse) ----
mkcfx
qui replace clat = 0.5 if ttre >= . & _n == _N
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
if _rc != 0 {
  dis as error "L1 FAIL: bad-value-on-never-treated rc=" _rc " (expected 0)"
  local ++fails
}
else {
  cap confirm matrix e(cfxatt_cfxpost)
  if _rc != 0 {
    dis as error "L1 FAIL: rc 0 but no e(cfxatt_cfxpost) - CFX block did not execute (vacuous pass)"
    local ++fails
  }
}

* ---- L2: twin parity - the poisoned never-treated row vs a benign integer there;
* returned cfx matrices must be bit-identical ----
tempname T1 T2
mkcfx
qui replace clat = 0.5 if ttre >= . & _n == _N
qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
matrix `T1' = e(cfxatt_cfxpost)
mkcfx
qui replace clat = 1 if ttre >= . & _n == _N
qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
matrix `T2' = e(cfxatt_cfxpost)
if mreldif(`T1', `T2') != 0 {
  dis as error "L2 FAIL: never-treated-row value moved cfx results, mreldif=" mreldif(`T1', `T2')
  local ++fails
}

* ---- L3: bad value on a treated-unit row still errors r(411) ----
mkcfx
qui replace clat = 0.5 if ttre < . & _n == 1
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
if _rc != 411 {
  dis as error "L3 FAIL: bad-value-on-treated rc=" _rc " (expected 411)"
  local ++fails
}

* ---- L4: missing on never-treated only -> NO caution; missing on a treated row -> caution ----
mkcfx
qui replace clat = . if ttre >= . & _n == _N
qui log using cfxlat_l4a.txt, text replace name(cfxl)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
qui log close cfxl
if _rc != 0 {
  dis as error "L4a FAIL: rc=" _rc
  local ++fails
}
else if strpos(fileread("cfxlat_l4a.txt"), "cfxlattice() variable is missing") != 0 {
  dis as error "L4a FAIL: caution fired for a never-treated-only missing"
  local ++fails
}
mkcfx
qui replace clat = . if ttre < . & t >= ttre & _n <= 400
qui count if clat >= .
local nmiss = r(N)
qui log using cfxlat_l4b.txt, text replace name(cfxl)
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
qui log close cfxl
if `nmiss' == 0 {
  dis as error "L4b FAIL: fixture built no missing treated rows"
  local ++fails
}
else if strpos(fileread("cfxlat_l4b.txt"), "observation(s) with a nonmissing treatment onset") == 0 {
  dis as error "L4b FAIL: onset-scoped caution missing for treated-row missings"
  local ++fails
}

* ---- L5: safety(off) still bypasses the value screen ----
mkcfx
qui replace clat = 0.5 if ttre >= . & _n == _N
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat) safety(off)
if _rc != 0 {
  dis as error "L5 FAIL: safety(off) leg rc=" _rc
  local ++fails
}

* ---- L6 P7 accident 2 (author ruling 2026-09-08): the <saveplots>_cfx.dta plot-data file is
* written on the ANALYTIC path too. It used to be produced only under bootstrap() or permute():
* analytically the figures come from marginsplot and saveplots() wrote only the .gph files, so
* the one machine-readable record of what was plotted was an engine-only artifact. The contract:
* the file exists, carries the plotted points, and the .gph files are still written beside it. ----
use "`DATA'/dta_cont", clear
cap erase l6probe_cfx.dta
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots saveplots(l6probe) replace
local l6rc = _rc
local l6stem "`e(wdid_cfxsave)'"
cap confirm file "l6probe_cfx.dta"
local l6dta = (_rc == 0)
cap confirm file "l6probe_cfxatt_cfxpost.gph"
local l6gph = (_rc == 0)
local l6rows 0
local l6cols ""
if `l6dta' {
  preserve
  use "l6probe_cfx.dta", clear
  qui ds
  local l6cols "`r(varlist)'"
  qui count
  local l6rows = r(N)
  restore
}
if `l6rc' != 0 | `l6dta' != 1 | `l6gph' != 1 | `l6rows' < 1 | strpos("`l6cols'", "matname") == 0 | strpos("`l6cols'", "b") == 0 {
  dis as error "L6 FAIL: analytic cfx plot data rc=`l6rc' stem=[`l6stem'] dta=`l6dta' gph=`l6gph' rows=`l6rows' cols=[`l6cols']"
  local ++fails
}
else dis as result "L6 PASS: the analytic run wrote l6probe_cfx.dta (`l6rows' plotted points) beside its .gph files"
cap erase l6probe_cfx.dta

* ---- L7 cfxplots-userlattice-axis (2026-09-23): a one-word lattice(var) with no cfxlattice()
* plots over the user's own variable, so the engine-path figures (drawn by wooldid_cfxplots from
* e() under bootstrap()/permute()) must carry its name or label and its value labels on the
* x-axis the way the analytic marginsplot figure does. The worker posts that variable in
* e(wdid_cfxlatvar) (it used to post only cfxlattice(), leaving the macro empty here) and
* wooldid_cfxplots defaults the x-title to the variable label, else its name, when the worker
* posted no xtitle. The analytic figure is untouched, an auto lattice keeps the macro empty
* (its tempvar does not outlive the worker, and a posted tempvar name could later confirm
* against an unrelated recycled tempvar), and cfxlattice() still wins when both are given.
* Figures are read back from the saveplots() .gph files: under set graphics off the analytic
* marginsplot figure cannot be exported by name ("could not find Graph window"). ----
cap program drop mkcfxlat3
program define mkcfxlat3
  use "data/dta_cont", clear
  qui gen byte lat3 = 1 + mod(i, 3)
  label define l7latl 1 "low dose" 2 "mid dose" 3 "high dose", replace
  label values lat3 l7latl
  label variable lat3 "Dose tercile"
end
* svg text probe: r(hit) = 1 when the exported figure's text holds the needle, 0 otherwise
cap program drop l7svghas
program define l7svghas, rclass
  syntax , file(string) needle(string)
  local hit 0
  cap confirm file "`file'"
  if _rc == 0 mata: st_local("hit", strofreal(sum(strpos(cat("`file'"), "`needle'") :> 0)))
  return scalar hit = `hit' > 0
end
local l7files l7probe_cfxatt_cfxpost.gph l7probe_cfxame_cfxpost.gph l7probe_cfx.dta l7probe_an.svg l7probe_bt.svg l7probe_pm.svg
foreach l7f of local l7files {
  cap erase `l7f'
}
* L7a analytic: the macro names the variable; the figure (marginsplot's own) is unchanged
mkcfxlat3
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(lat3) esprelength(1) espostlength(1) makecfxplots saveplots(l7probe) replace
local l7arc = _rc
local l7alat "`e(wdid_cfxlatvar)'"
local l7axt `"`e(wdid_cfxxt)'"'
cap graph use "l7probe_cfxatt_cfxpost.gph", name(l7g, replace)
cap graph export "l7probe_an.svg", name(l7g) replace
local l7aexp = _rc
l7svghas , file(l7probe_an.svg) needle(Dose tercile)
local l7a_xt = r(hit)
l7svghas , file(l7probe_an.svg) needle(low dose)
local l7a_vl = r(hit)
if `l7arc' != 0 | "`l7alat'" != "lat3" | `"`l7axt'"' != "" | `l7aexp' != 0 | `l7a_xt' != 1 | `l7a_vl' != 1 {
  dis as error `"L7a FAIL: analytic lattice(lat3) rc=`l7arc' latvar=[`l7alat'] xt=[`l7axt'] exprc=`l7aexp' svg_varlabel=`l7a_xt' svg_valuelabel=`l7a_vl'"'
  local ++fails
}
* L7b bootstrap: the same spec drawn by wooldid_cfxplots must show the same axis. the analytic
* leg's .gph is erased first, so a run that draws nothing cannot hand this leg that figure
cap erase l7probe_cfxatt_cfxpost.gph
mkcfxlat3
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(lat3) esprelength(1) espostlength(1) makecfxplots saveplots(l7probe) replace bootstrap(standard) reps(10) rseed(7) dots(0)
local l7brc = _rc
local l7blat "`e(wdid_cfxlatvar)'"
cap graph use "l7probe_cfxatt_cfxpost.gph", name(l7g, replace)
cap graph export "l7probe_bt.svg", name(l7g) replace
local l7bexp = _rc
l7svghas , file(l7probe_bt.svg) needle(Dose tercile)
local l7b_xt = r(hit)
l7svghas , file(l7probe_bt.svg) needle(low dose)
local l7b_vl = r(hit)
if `l7brc' != 0 | "`l7blat'" != "lat3" | `l7bexp' != 0 | `l7b_xt' != 1 | `l7b_vl' != 1 {
  dis as error `"L7b FAIL: bootstrap lattice(lat3) rc=`l7brc' latvar=[`l7blat'] exprc=`l7bexp' svg_varlabel=`l7b_xt' svg_valuelabel=`l7b_vl'"'
  local ++fails
}
* L7c permute: the shell's other wooldid_cfxplots call site, same contract
cap erase l7probe_cfxatt_cfxpost.gph
mkcfxlat3
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(lat3) esprelength(1) espostlength(1) makecfxplots saveplots(l7probe) replace permute(shuffle) permreps(20) rseed(7)
local l7crc = _rc
local l7clat "`e(wdid_cfxlatvar)'"
cap graph use "l7probe_cfxatt_cfxpost.gph", name(l7g, replace)
cap graph export "l7probe_pm.svg", name(l7g) replace
local l7cexp = _rc
l7svghas , file(l7probe_pm.svg) needle(Dose tercile)
local l7c_xt = r(hit)
l7svghas , file(l7probe_pm.svg) needle(low dose)
local l7c_vl = r(hit)
if `l7crc' != 0 | "`l7clat'" != "lat3" | `l7cexp' != 0 | `l7c_xt' != 1 | `l7c_vl' != 1 {
  dis as error `"L7c FAIL: permute lattice(lat3) rc=`l7crc' latvar=[`l7clat'] exprc=`l7cexp' svg_varlabel=`l7c_xt' svg_valuelabel=`l7c_vl'"'
  local ++fails
}
* L7d auto lattice: the tempvar must not be posted
mkcfxlat3
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots
if _rc != 0 | "`e(wdid_cfxlatvar)'" != "" {
  dis as error "L7d FAIL: auto lattice rc=" _rc " latvar=[`e(wdid_cfxlatvar)'] (expected empty)"
  local ++fails
}
* L7e cfxlattice() beside a one-word lattice(): cfxlattice() is the plot variable and wins
mkcfxlat3
qui gen byte clat = mod(i, 3)
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(lat3) cfxlattice(clat) esprelength(1) espostlength(1) makecfxplots
if _rc != 0 | "`e(wdid_cfxlatvar)'" != "clat" | strpos(`"`e(wdid_cfxxt)'"', "estimation lattice: clat") == 0 {
  dis as error `"L7e FAIL: cfxlattice(clat)+lattice(lat3) rc=`=_rc' latvar=[`e(wdid_cfxlatvar)'] xt=[`e(wdid_cfxxt)']"'
  local ++fails
}
if `fails' == 0 dis as result "L7 PASS: engine-path cfx figures carry the user lattice variable's label and value labels"
cap graph drop l7g
foreach l7f of local l7files {
  cap erase `l7f'
}

if `fails' == 0 dis as result "CFXLAT_ALL_PASS"
else dis as error "CFXLAT_FAILS=`fails'"
