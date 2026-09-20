clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 121212
set sortseed 121212
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat"

use "`D'/p1like.dta", clear
qui gquantiles latu = dose_it , xtile nquantiles(3)
qui save "`D'/p1like_u.dta", replace

foreach spec in ///
  "lattice(all 3) robust" ///
  "lattice(all 3) robust regtol(14)" ///
  "lattice(all 3) robust regtol(6)" ///
  "lattice(all 3) ccnoabsorb robust" ///
  "lattice(all 3) customvce(robust)" ///
  "lattice(all 3) unconditionalse robust" ///
  "lattice(latu) robust" ///
  "lattice(all 3) cluster(st)" {
  use "`D'/p1like_u.dta", clear
  dis "########## LEG: `spec'"
  cap noisily qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) `spec' esprelength(2) espostlength(2)
  local rc = _rc
  if `rc' == 0 {
    mat M = e(wdidmainresults)
    dis "DIAG|`spec'|N=" e(N) "|df_r=" e(df_r) "|rss=" %10.4e e(rss) "|r2=" %8.6f e(r2)
    dis "DIAGb|ame=" %14.8f M[2,1] "|se_ame=" %14.8g M[2,2] "|att=" %14.8f M[1,1] "|se_att=" %14.8g M[1,2]
  }
  else dis "DIAG|`spec'|RC=`rc'"
}
dis "DONE"
