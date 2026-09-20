* r5_strata_both.do - NEW probes (not seeded), same surface as the seeded leads.
*  S1: permsave() under permstrata() - is the onset block still keyed to the right cluster,
*      and is the onset multiset preserved WITHIN each stratum in every replication?
*      (ado:7500 documents a frlink reorder trap on exactly this frame; the region labels
*       here deliberately CROSS ascending cluster order.)
*  S2: permsave() under permonsetcohorts (one day old on the pinned build)
*  S3: permute() and bootstrap() given TOGETHER - refused, or does one silently win?
*  S4: bootsave()/permsave() name-rule symmetry on a directory part and an empty name
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 60660
set sortseed 60660
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w5"
which wooldid
cap shell rm -f *.dta

* 80 units in 4 regions; onset cohort assigned so region order CROSSES unit order
clear
set obs 80
gen long unit = _n
gen byte region = mod(_n, 4) + 1
gen byte grp = mod(floor((_n-1)/4), 5)
gen ttre = .
replace ttre = 2005 if grp == 1
replace ttre = 2008 if grp == 2
replace ttre = 2011 if grp == 3
gen double ue = rnormal()*0.8
expand 14
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + 0.4*(year - ttre + 1)/3 if year >= ttre & ttre < .
save "panel.dta", replace

preserve
  bys unit: keep if _n == 1
  keep unit region ttre
  rename ttre obsttre
  save "obs_onsets.dta", replace
  tab region obsttre, missing
restore

* ---- checker: reshape the onset block long and compare to the data
prog def chkfile
  syntax , file(string) tag(string) [strata]
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
  * (1) observed row
  qui count if _rep == 0 & !(drawn == obsttre | (drawn >= . & obsttre >= .))
  dis "`tag'_ROW0_BAD=" r(N)
  * (2) global multiset preservation in every row
  tempvar n1
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
  * (3) within-stratum multiset preservation
  if "`strata'" != "" {
    local bads 0
    forval g = 1/4 {
      forval L = 0/3 {
        qui count if _rep == 0 & olev == `L' & region == `g'
        local want = r(N)
        forval r = 0/`nreps' {
          qui count if _rep == `r' & lev == `L' & region == `g'
          if r(N) != `want' {
            local ++bads
            if `bads' <= 8 dis "`tag'_STRATUM_LEAK rep=`r' region=`g' lev=`L' got=" r(N) " want=`want'"
          }
        }
      }
    }
    dis "`tag'_STRATA_BAD=" `bads'
  }
  restore
end

* ---------------- S1: permstrata(region)
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    permstrata(region) permsave(ps) replace dots(0)
dis "S1_RC=" _rc
capture confirm file "ps.dta"
dis "S1_FILE_RC=" _rc
chkfile, file("ps.dta") tag(S1) strata

* ---------------- S2: permonsetcohorts
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    permonsetcohorts permsave(po) replace dots(0)
dis "S2_RC=" _rc
dis "S2_irebuilt=" e(perm_irebuilt)
capture confirm file "po.dta"
dis "S2_FILE_RC=" _rc
chkfile, file("po.dta") tag(S2)

* ---------------- S3: both engines at once
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    bootstrap(standard) reps(20) dots(0)
dis "S3_RC=" _rc
dis "S3_perm_reps=" e(perm_reps) " S3_boot_reps=" e(boot_reps)
dis "S3_perm_method=`e(perm_method)' S3_boot_method=`e(boot_method)'"

* ---------------- S4: name-rule symmetry
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(sub/dir) dots(0)
dis "S4a_perm_dir_RC=" _rc
capture noisily wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(sub/dir) dots(0)
dis "S4b_boot_dir_RC=" _rc
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) permsave(.dta) dots(0)
dis "S4c_perm_empty_RC=" _rc
capture noisily wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(20) bootsave(.dta) dots(0)
dis "S4d_boot_empty_RC=" _rc

shell ls -1
dis "ALLDONE"
