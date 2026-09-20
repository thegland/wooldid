* r6_onsetcoh_both.do - NEW probes (not seeded), continuation of r5.
*  S2: permsave() under permonsetcohorts with i = the onset grouping and cluster() finer
*      (the T1/T2/T3 finer-than-i machinery is one day old on the pinned build)
*  S3: permute() and bootstrap() given TOGETHER - refused, or does one silently win?
*  S4: bootsave()/permsave() name-rule symmetry on a directory part and an empty name
*  S5: permsave() when a cluster id would collide with a coefficient column name
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771177
set sortseed 771177
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w6"
which wooldid
cap shell rm -f *.dta

clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n, 5)
gen ttre = .
replace ttre = 2005 if grp == 1
replace ttre = 2008 if grp == 2
replace ttre = 2011 if grp == 3
* i = the onset grouping itself (4 levels incl. never-treated)
gen byte coh = cond(ttre >= ., 0, grp)
gen double ue = rnormal()*0.8
expand 14
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + 0.4*(year - ttre + 1)/3 if year >= ttre & ttre < .
save "panel.dta", replace

preserve
  bys unit: keep if _n == 1
  keep unit ttre
  rename ttre obsttre
  save "obs_onsets.dta", replace
restore

prog def chkfile
  syntax , file(string) tag(string)
  preserve
  use "`file'", clear
  qui ds onset_*
  local nons : word count `r(varlist)'
  qui su _rep
  local nreps = r(max)
  keep _rep onset_*
  qui reshape long onset_, i(_rep) j(unit)
  rename onset_ drawn
  qui merge m:1 unit using "obs_onsets.dta", keep(match master)
  qui count if _merge != 3
  dis "`tag'_UNMATCHED_UNITS=" r(N) "  ncols=`nons'  nreps=`nreps'"
  qui count if _rep == 0 & !(drawn == obsttre | (drawn >= . & obsttre >= .))
  dis "`tag'_ROW0_BAD=" r(N)
  qui gen byte lev = 0
  qui replace lev = 1 if drawn == 2005
  qui replace lev = 2 if drawn == 2008
  qui replace lev = 3 if drawn == 2011
  qui gen byte olev = 0
  qui replace olev = 1 if obsttre == 2005
  qui replace olev = 2 if obsttre == 2008
  qui replace olev = 3 if obsttre == 2011
  local badg 0
  forval L = 0/3 {
    qui count if _rep == 0 & olev == `L'
    local want = r(N)
    forval r = 0/`nreps' {
      qui count if _rep == `r' & lev == `L'
      if r(N) != `want' {
        local ++badg
        if `badg' <= 5 dis "`tag'_MULTISET_MISMATCH rep=`r' lev=`L' got=" r(N) " want=`want'"
      }
    }
  }
  dis "`tag'_MULTISET_BAD=" `badg'
  restore
end

* ---------------- S2: permonsetcohorts, i = coh, cluster = unit (finer than i)
use "panel.dta", clear
capture noisily wooldid y coh year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    permonsetcohorts permsave(po) replace dots(0)
dis "S2_RC=" _rc
dis "S2_irebuilt=" e(perm_irebuilt)
capture confirm file "po.dta"
dis "S2_FILE_RC=" _rc
capture noisily chkfile, file("po.dta") tag(S2)
dis "S2_CHK_RC=" _rc

* ---------------- S2b: same design WITHOUT permonsetcohorts, for contrast
use "panel.dta", clear
capture noisily wooldid y coh year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    permsave(pn) replace dots(0)
dis "S2b_RC=" _rc
capture noisily chkfile, file("pn.dta") tag(S2b)
dis "S2b_CHK_RC=" _rc

* ---------------- S3: both engines at once
use "panel.dta", clear
capture noisily wooldid y coh year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    bootstrap(standard) reps(20) dots(0)
dis "S3_RC=" _rc
dis "S3_perm_reps=" e(perm_reps) " S3_boot_reps=" e(boot_reps)
dis "S3_perm_method=`e(perm_method)' S3_boot_method=`e(boot_method)'"
dis "S3_cmdline=`e(cmdline)'"

* ---------------- S4: name-rule symmetry
use "panel.dta", clear
capture noisily wooldid y coh year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(sub/dir) dots(0)
dis "S4a_perm_dir_RC=" _rc
capture noisily wooldid y coh year ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(sub/dir) dots(0)
dis "S4b_boot_dir_RC=" _rc
capture noisily wooldid y coh year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(.dta) dots(0)
dis "S4c_perm_empty_RC=" _rc
capture noisily wooldid y coh year ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(.dta) dots(0)
dis "S4d_boot_empty_RC=" _rc

shell ls -1
dis "ALLDONE"
