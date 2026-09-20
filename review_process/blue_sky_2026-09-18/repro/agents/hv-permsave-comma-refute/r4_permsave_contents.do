* r4_permsave_contents.do - NEW probe (not seeded): is the permsave() FILE correct?
* Invariants checked against e() and against the data:
*   (1) row _rep==0 of onset_<c> equals cluster c's observed onset
*   (2) every row's onset multiset equals the observed onset multiset (shuffle = relabel)
*   (3) SD of each coefficient column over reps 1..R equals e(perm_sd) for that column
*   (4) non-missing replication count equals e(perm_n)
*   (5) row _rep==0 of each coefficient column equals e(perm_comparator_b)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w4"
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
gen double ue = rnormal()*0.8
expand 14
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + 0.4*(year - ttre + 1)/3 if year >= ttre & ttre < .
save "panel.dta", replace

preserve
  bys unit: keep if _n == 1
  keep unit ttre
  save "obs_onsets.dta", replace
  qui count
  dis "NCLUS=" r(N)
  qui count if ttre == 2005
  local obsc1 = r(N)
  qui count if ttre == 2008
  local obsc2 = r(N)
  qui count if ttre == 2011
  local obsc3 = r(N)
  qui count if ttre >= .
  local obscm = r(N)
restore
dis "OBSCOUNTS 2005=`obsc1' 2008=`obsc2' 2011=`obsc3' miss=`obscm'"

wooldid y unit year ttre, att cluster(unit) esprelength(2) espostlength(2) ///
    permute(shuffle) permreps(40) rseed(12345) permsave(pv) replace dots(0)
dis "RUN_RC=" _rc

tempname PSD PN PCB
matrix `PSD' = e(perm_sd)
matrix `PN'  = e(perm_n)
matrix `PCB' = e(perm_comparator_b)
local cn : colnames `PSD'
local cnb : colnames `PCB'
dis "PERMSD_COLS=`cn'"

use "pv.dta", clear
dis "FILE_N=" _N
qui ds onset_*
local ons `r(varlist)'
local nons : word count `ons'
dis "FILE_NONSETCOLS=" `nons'

* ---- (1) observed row vs the data
local bad1 0
local ck 0
foreach v of local ons {
  local u = subinstr("`v'", "onset_", "", 1)
  qui su `v' if _rep == 0, meanonly
  local got = r(mean)
  local gotn = r(N)
  preserve
    qui use "obs_onsets.dta", clear
    qui su ttre if unit == `u', meanonly
    local want = r(mean)
    local wantn = r(N)
  restore
  local ++ck
  local mism 0
  if `wantn' == 0 & `gotn' != 0 local mism 1
  if `wantn' == 1 & `gotn' == 0 local mism 1
  if `wantn' == 1 & `gotn' == 1 {
    if abs(`got' - `want') > 1e-10 local mism 1
  }
  if `mism' {
    local ++bad1
    if `bad1' <= 5 dis "ROW0_MISMATCH unit=`u' file=`got'(n=`gotn') data=`want'(n=`wantn')"
  }
}
dis "CHK1_checked=" `ck' " CHK1_bad=" `bad1'

* ---- (2) multiset preservation, every row including rep 0
local bad2 0
forval r = 0/40 {
  local c1 0
  local c2 0
  local c3 0
  local cm 0
  foreach v of local ons {
    qui count if _rep == `r' & `v' == 2005
    local c1 = `c1' + r(N)
    qui count if _rep == `r' & `v' == 2008
    local c2 = `c2' + r(N)
    qui count if _rep == `r' & `v' == 2011
    local c3 = `c3' + r(N)
    qui count if _rep == `r' & `v' >= .
    local cm = `cm' + r(N)
  }
  if `c1' != `obsc1' | `c2' != `obsc2' | `c3' != `obsc3' | `cm' != `obscm' {
    local ++bad2
    if `bad2' <= 6 dis "MULTISET_MISMATCH rep=`r' got(`c1',`c2',`c3',`cm') want(`obsc1',`obsc2',`obsc3',`obscm')"
  }
}
dis "CHK2_bad=" `bad2'

* ---- (3)(4)(5) coefficient columns against e()
qui ds onset_* _rep, not
local coefs `r(varlist)'
dis "FILE_COEFCOLS=`coefs'"
local bad3 0
local bad4 0
local bad5 0
foreach v of local coefs {
  local pos : list posof "`v'" in cn
  if `pos' == 0 {
    dis "NO_E_COLUMN_FOR `v'"
    continue
  }
  qui su `v' if _rep > 0
  local fsd = r(sd)
  local fn  = r(N)
  qui su `v' if _rep == 0, meanonly
  local f0 = r(mean)
  local esd = `PSD'[1, `pos']
  local en  = `PN'[1, `pos']
  local posb : list posof "`v'" in cnb
  local eb = cond(`posb' > 0, `PCB'[1, `posb'], .)
  local d3 = abs(`fsd' - `esd')/max(abs(`esd'), 1e-12)
  local d5 = abs(`f0' - `eb')/max(abs(`eb'), 1e-12)
  dis "COL `v' fsd=" %21x `fsd' " esd=" %21x `esd' " rd=" `d3'
  dis "COL `v' fn=" `fn' " en=" `en' " f0=" %21x `f0' " ecmp=" %21x `eb' " rd=" `d5'
  if `d3' > 1e-10 local ++bad3
  if `fn' != `en' local ++bad4
  if `posb' > 0 & `d5' > 1e-10 local ++bad5
}
dis "CHK3_bad=" `bad3' " CHK4_bad=" `bad4' " CHK5_bad=" `bad5'
dis "ALLDONE"
