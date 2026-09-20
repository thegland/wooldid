* p07_subgroup.do -- (a) can the SUBGROUP estimands be rebuilt from the single
* full-sample e(gteffects) table? (the help returns one table for all subgroups)
* (b) which weight columns are populated when only SOME estimands are requested?
* (c) does the identity survive bootstrap()?
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime

mkpanel , n(90) t(9) onsets("2004 2006") ysd(0.4) seed(6161)
* subgroup must be constant within the (i,t) cell -> define it at the cohort level
gen byte sg = mod(cohp, 2)
tempfile base
qui save `base'

*================= S1: subgroup rebuild ===============================
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) subgroups(sg)
dis "S1_RC=" _rc
gtcheck , tag(S1_sg_full) pre es
dis "S1 e(b) colnames:"
dis "`: colnames e(b)'"
tempname GT
matrix `GT' = e(gteffects)
* map each table row to its subgroup via (i, caltime)
preserve
qui keep cohp year sg
qui duplicates drop
tempfile sgmap
qui save `sgmap'
restore
preserve
clear
qui svmat double `GT', names(col)
rename i cohp
rename caltime year
qui merge m:1 cohp year using `sgmap', keep(master match) nogen
qui keep if relyr >= 0
foreach v in att att_it att_i att_itime {
  tempvar num`v'
  qui gen double `num`v'' = estimate*w_`v'
}
levelsof sg, local(sgl)
foreach s of local sgl {
  foreach v in att att_it att_i att_itime {
    qui sum `num`v'' if sg == `s', meanonly
    local nu = r(sum)
    qui sum w_`v' if sg == `s', meanonly
    local de = r(sum)
    dis as result "SG|rebuilt|sg=`s'|`v'|" %21x `nu'/`de' "|sumw=" %12.6g `de'
  }
}
restore
dis "S1 subgroup e() matrices:"
cap noi matrix list e(wdidsgresults_0), format(%21x)
cap noi matrix list e(wdidsgresults_1), format(%21x)

*================= S2: which weight columns are populated? ============
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, att_i gteffects esprelength(2) espostlength(2)
dis "S2_RC=" _rc
dis "S2 e(b): `: colnames e(b)'"
preserve
clear
cap noi svmat double e(gteffects), names(col)
foreach v in w_att w_att_it w_att_i w_att_itime {
  cap confirm variable `v'
  if _rc == 0 {
    qui count if `v' < .
    dis as result "S2|COL|`v'|nonmissing=" r(N)
  }
  else dis as result "S2|COL|`v'|ABSENT"
}
restore

use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, att_it gteffects esprelength(2) espostlength(2)
dis "S2b_RC=" _rc
preserve
clear
cap noi svmat double e(gteffects), names(col)
foreach v in w_att w_att_it w_att_i w_att_itime {
  cap confirm variable `v'
  if _rc == 0 {
    qui count if `v' < .
    dis as result "S2b|COL|`v'|nonmissing=" r(N)
  }
  else dis as result "S2b|COL|`v'|ABSENT"
}
restore

*================= S3: bootstrap() -- points must still rebuild ========
use `base', clear
ereturn clear
cap noi wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) bootstrap(standard) reps(25) cluster(unit) rseed(77)
dis "S3_RC=" _rc
gtcheck , tag(S3_boot) pre es

dis "DONE_P07"
