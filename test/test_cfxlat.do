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

* ---- L8 shell-cfxtwin-frame (change audit 2026-09-22): when the analytic plot-data write is
* refused (a <stem>_cfx.dta exists and replace was not given), wooldid_cfxplots dies inside its
* tempname frame and the shell prints the caution with rc 0 - and must hand the caller back the
* frame they called from, as the engine branches do after the same failure. It used to leave the
* session in default, silently, rc 0. The caution must actually have printed (checked through a
* named log), else a twin that never fires would pass this vacuously ----
foreach f in l8probe_cfx.dta l8probe_cfxatt_cfxpost.gph l8probe_cfxame_cfxpost.gph {
  cap erase "`f'"
}
cap frame drop l8fr
frame create l8fr
frame change l8fr
use "`DATA'/dta_cont", clear
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots saveplots(l8probe) replace
local l8rc1 = _rc
cap confirm file "l8probe_cfx.dta"
local l8dta = (_rc == 0)
* the .gph files go, so marginsplot's own save succeeds and only the plot-data write collides
foreach f in l8probe_cfxatt_cfxpost.gph l8probe_cfxame_cfxpost.gph {
  cap erase "`f'"
}
cap frame change l8fr
qui log using l8probe_log.txt, text replace name(l8)
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots saveplots(l8probe)
local l8rc2 = _rc
qui log close l8
qui frame
local l8frame "`r(currentframe)'"
* a short fragment: the log line may wrap
local l8caut = strpos(fileread("l8probe_log.txt"), "plot-data file l8probe_cfx.dta") > 0
if `l8rc1' != 0 | `l8dta' != 1 | `l8rc2' != 0 | `l8caut' != 1 | "`l8frame'" != "l8fr" | "`e(cmd)'" != "wooldid" {
  dis as error "L8 FAIL: refused plot-data write rc1=`l8rc1' dta=`l8dta' rc2=`l8rc2' caution=`l8caut' frame=[`l8frame'] (expected l8fr) ecmd=[`e(cmd)']"
  local ++fails
}
else dis as result "L8 PASS: the refused plot-data write printed its caution with rc 0 and handed back frame l8fr"
frame change default
cap frame drop l8fr
foreach f in l8probe_cfx.dta l8probe_cfxatt_cfxpost.gph l8probe_cfxame_cfxpost.gph l8probe_log.txt {
  cap erase "`f'"
}

* ---- L9 shell-update-stale-cfxsave (change audit 2026-09-22): wooldid, update exits without
* posting e(), so after an analytic makecfxplots saveplots() run the previous run's
* e(wdid_cfxsave) is still set - and the shell's plot-data twin, gated on it, used to re-write
* that run's <stem>_cfx.dta on update (silently under replace, clobbering anything the user had
* done to the file; with a spurious caution without it). The twin must fire only on a stem THIS
* call's worker handed back. Run under a copy of the ado whose one net install line is commented
* out, so the leg can never reach the network or the installed build: the copy must differ in
* exactly that line, else update is not called at all. Non-vacuous by construction: the stale
* stem must still be in e() after update and the update message must have printed ----
foreach f in l9probe_cfx.dta l9probe_cfxatt_cfxpost.gph l9probe_cfxame_cfxpost.gph l9probe_log.txt l9probe_stub.do {
  cap erase "`f'"
}
local l9rc2 .
local l9stale ""
filefilter "../wooldidinstall/wooldid.ado" "l9probe_stub.do", from("net install wooldid,") to("*net install wooldid,") replace
local l9occ = r(occurrences)
* loaded whatever the count - a copy with no line changed is the real build, harmless to load -
* and only the update call is gated on the count. a second do of the ado in one session dies
* r(3000) (measured: wooldid_sg1live() has no mata drop guard, and the earlier legs' reghdfe
* calls leave matastrict on, which the ado's mata does not compile under), hence the clear all
* (macros survive it) and matastrict off before each load
clear all
mata: mata set matastrict off
qui do "l9probe_stub.do"
use "`DATA'/dta_cont", clear
cap noi wooldid y i t ttre, att cluster(i) contreat(intensity) contreatcontrols(basic) lattice(all 3) esprelength(1) espostlength(1) makecfxplots saveplots(l9probe) replace
local l9rc1 = _rc
cap confirm file "l9probe_cfx.dta"
local l9dta = (_rc == 0)
* the file goes, so a re-write by update would be visible as its return
cap erase "l9probe_cfx.dta"
qui log using l9probe_log.txt, text replace name(l9)
if `l9occ' == 1 {
  cap noi wooldid, update
  local l9rc2 = _rc
}
qui log close l9
local l9stale "`e(wdid_cfxsave)'"
cap confirm file "l9probe_cfx.dta"
local l9back = (_rc == 0)
local l9upd = strpos(fileread("l9probe_log.txt"), "Update complete") > 0
local l9caut = strpos(fileread("l9probe_log.txt"), "plot-data file") > 0
if `l9occ' != 1 | `l9rc1' != 0 | `l9dta' != 1 | `l9rc2' != 0 | `l9upd' != 1 | "`l9stale'" != "l9probe" | "`e(cmd)'" != "wooldid" | `l9back' != 0 | `l9caut' != 0 {
  dis as error "L9 FAIL: update after a saveplots() cfx run stub=`l9occ' rc1=`l9rc1' dta=`l9dta' rc2=`l9rc2' updmsg=`l9upd' stale=[`l9stale'] ecmd=[`e(cmd)'] rewritten=`l9back' (expected 0) caution=`l9caut' (expected 0)"
  local ++fails
}
else dis as result "L9 PASS: update left the previous run's plot-data file alone and printed no plot-data caution"
* back to the real build
clear all
mata: mata set matastrict off
qui do "../wooldidinstall/wooldid.ado"
foreach f in l9probe_cfx.dta l9probe_cfxatt_cfxpost.gph l9probe_cfxame_cfxpost.gph l9probe_log.txt l9probe_stub.do {
  cap erase "`f'"
}

* ---- L10 cfx-omitted-rows (author ruling 2026-09-23): on a run that computes the plot cells'
* standard errors, margins marks a lattice level whose every dose-slope term was dropped as
* omitted in the row name (1bno. on the first level, 2o. on a later one) and carries a structural
* 0 as its estimate. wooldid_cfxstripe keeps such a row under its lattice_# key with every column
* missing. L10a: i = unit with one observation per cell, so every treated dose-slope term drops
* and the all-slopes-dropped caution prints - the run completes, both figures are saved, the ame
* matrix is blank, the att matrix keeps its estimates and standard errors, and the plot-data file
* carries the blank ame rows. L10b/L10c: the level-1 (L10b) or level-2 (L10c) units each form
* their own i, one observation per cell, so only that lattice level is omitted - its ame row is
* blank and the other two keep their estimates ----
cap program drop mkq33c
program define mkq33c
  clear
  qui set obs 20
  qui gen i = _n
  qui gen ttre = cond(i <= 8, 5, cond(i <= 15, 7, .))
  qui expand 10
  qui bysort i: gen t = _n
  set seed 4411
  qui gen dose = runiform() + 0.5
  qui gen y = 1.1*(t >= ttre & ttre < .)*dose + 0.08*i + 0.04*t + rnormal()
end
* l10chk: bad = number of failed shape checks on e(cfxame_cfxpost)/e(cfxatt_cfxpost); blank =
* the ame lattice levels expected all-missing, every other row must carry b and a positive se
cap program drop l10chk
program define l10chk
  syntax , [blank(numlist)]
  local bad 0
  if "`e(wdid_cfxmats)'" != "cfxame_cfxpost cfxatt_cfxpost" local ++bad
  foreach nm in cfxame_cfxpost cfxatt_cfxpost {
    cap confirm matrix e(`nm')
    if _rc != 0 {
      local ++bad
      continue
    }
    tempname M
    matrix `M' = e(`nm')
    if rowsof(`M') != 3 | colsof(`M') != 6 | "`: rowfullnames `M''" != "lattice_1 lattice_2 lattice_3" local ++bad
    else {
      forval r = 1/3 {
        local isblank = ("`nm'" == "cfxame_cfxpost" & `: list r in blank')
        forval c = 1/6 {
          if `isblank' & `M'[`r', `c'] < . local ++bad
        }
        if !`isblank' & (`M'[`r', 1] >= . | `M'[`r', 2] >= . | `M'[`r', 2] <= 0) local ++bad
      }
    }
  }
  c_local l10bad `bad'
end
foreach f in l10probe_cfx.dta l10probe_cfxatt_cfxpost.gph l10probe_cfxame_cfxpost.gph l10probe_log.txt {
  cap erase "`f'"
}
mkq33c
qui log using l10probe_log.txt, text replace name(l10)
cap noi wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) makecfxplots lattice(all 3) esprelength(1) espostlength(1) saveplots(l10probe) replace
local l10arc = _rc
qui log close l10
local l10bad 1
if `l10arc' == 0 l10chk , blank(1 2 3)
local l10abad = `l10bad'
local l10caut = strpos(fileread("l10probe_log.txt"), "every treated cohort-period dose-slope term") > 0
local l10gph 0
foreach f in l10probe_cfxatt_cfxpost.gph l10probe_cfxame_cfxpost.gph {
  cap confirm file "`f'"
  if _rc == 0 local ++l10gph
}
local l10dame .
local l10datt .
cap confirm file "l10probe_cfx.dta"
if _rc == 0 {
  preserve
  qui use "l10probe_cfx.dta", clear
  qui count if matname == "cfxame_cfxpost" & b >= .
  local l10dame = r(N)
  qui count if matname == "cfxatt_cfxpost" & b < .
  local l10datt = r(N)
  restore
}
if `l10arc' != 0 | `l10abad' != 0 | `l10caut' != 1 | `l10gph' != 2 | "`l10dame'" != "3" | "`l10datt'" != "3" {
  dis as error "L10a FAIL: all dose slopes dropped, makecfxplots with standard errors rc=`l10arc' (r(498) = the row-name parse) shape_bad=`l10abad' caution=`l10caut' gph=`l10gph' dta_blank_ame=`l10dame' dta_att=`l10datt'"
  local ++fails
}
else dis as result "L10a PASS: every dose slope dropped - the run completes, the ame cells are blank under lattice_1-3, the att cells and both figures stand"
cap program drop mkpartl
program define mkpartl
  args lconst
  clear
  set seed 9031
  set sortseed 9031
  qui set obs 45
  qui gen u = _n
  qui gen ttre = cond(u <= 15, 5, cond(u <= 30, 7, .))
  qui gen L = 1 + mod(u, 3)
  qui expand 10
  qui bysort u: gen t = _n
  qui gen dose = runiform() + 1
  qui gen i = cond(L == `lconst', u, cond(ttre == 5, 1001, cond(ttre == 7, 1002, u)))
  qui gen y = 1.1*(t >= ttre & ttre < .)*dose + 0.08*u + 0.04*t + rnormal()
end
foreach lc in 1 2 {
  local l10lg = cond(`lc' == 1, "L10b", "L10c")
  mkpartl `lc'
  cap noi wooldid y i t ttre, att cluster(u) contreat(dose) contreatcontrols(basic) makecfxplots lattice(L) esprelength(1) espostlength(1)
  local l10rc = _rc
  local l10bad 1
  if `l10rc' == 0 l10chk , blank(`lc')
  if `l10rc' != 0 | `l10bad' != 0 {
    dis as error "`l10lg' FAIL: lattice level `lc' omitted rc=`l10rc' (r(498) = the row-name parse) shape_bad=`l10bad' (its ame row must be blank, every other row carry b and se)"
    local ++fails
  }
  else dis as result "`l10lg' PASS: lattice level `lc' omitted - its ame row is blank, the other levels keep their estimates"
}
foreach f in l10probe_cfx.dta l10probe_cfxatt_cfxpost.gph l10probe_cfxame_cfxpost.gph l10probe_log.txt {
  cap erase "`f'"
}

if `fails' == 0 dis as result "CFXLAT_ALL_PASS"
else dis as error "CFXLAT_FAILS=`fails'"
