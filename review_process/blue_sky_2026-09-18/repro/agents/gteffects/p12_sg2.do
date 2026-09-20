* p12_sg2.do -- re-run of p07 with the case bug fixed, plus:
*  - can the SUBGROUP estimands be rebuilt from the single full-sample table?
*  - which weight columns are populated when only SOME estimands are requested?
*  - identity under bootstrap()
*  - unconditionalse (needs robust/cluster)
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
gen byte sg = mod(unit, 2)
tempfile base
qui save `base'

*================= S1: subgroup rebuild (i = unit so sg is cell-constant) ====
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects esprelength(2) espostlength(2) subgroups(sg)
dis "S1_RC=" _rc
gtcheck , tag(S1_sg_full) pre es
tempname GT
matrix `GT' = e(gteffects)
foreach s in 0 1 {
  foreach v in att att_it att_i att_itime {
    local c = colnumb(e(b), "main_SG`s'_`v'")
    if `c' < . dis as result "SG|reported|sg=`s'|`v'|" %21x e(b)[1,`c']
  }
}
preserve
qui keep unit year sg
qui duplicates drop
tempfile sgmap
qui save `sgmap'
restore
preserve
clear
qui svmat double `GT', names(col)
rename i unit
rename caltime year
qui merge m:1 unit year using `sgmap', keep(master match) nogen
qui keep if relyr >= 0
foreach v in att att_it att_i att_itime {
  qui gen double nm_`v' = estimate*w_`v'
}
levelsof sg, local(sgl)
foreach s of local sgl {
  foreach v in att att_it att_i att_itime {
    qui sum nm_`v' if sg == `s', meanonly
    local nu = r(sum)
    qui sum w_`v' if sg == `s', meanonly
    local de = r(sum)
    dis as result "SG|rebuilt |sg=`s'|`v'|" %21x `nu'/`de' " sumw=" %12.6g `de'
  }
}
restore

*================= S2: which weight columns are populated? ============
foreach req in "att_i" "att_it" "att_itime" "att" "att att_i" {
  use `base', clear
  ereturn clear
  cap noi wooldid y unit year ttre, `req' gteffects esprelength(2) espostlength(2)
  dis "S2_[`req']_RC=" _rc
  if _rc == 0 {
    preserve
    clear
    cap noi svmat double e(gteffects), names(col)
    foreach v in w_att w_att_it w_att_i w_att_itime {
      cap confirm variable `v'
      if _rc == 0 {
        qui count if `v' < .
        dis as result "S2|[`req']|`v'|nonmissing=" r(N)
      }
      else dis as result "S2|[`req']|`v'|ABSENT"
    }
    restore
  }
}

*================= S3: bootstrap() -- points must still rebuild ========
use `base', clear
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects esprelength(2) espostlength(2) bootstrap(standard) reps(25) cluster(unit) rseed(77)
dis "S3_RC=" _rc
gtcheck , tag(S3_boot) pre es

*================= S4: unconditionalse with robust ====================
use `base', clear
ereturn clear
cap noi wooldid y unit year ttre, att gteffects esprelength(2) espostlength(2) unconditionalse robust
dis "S4_RC=" _rc
gtcheck , tag(S4_uncond) pre es

dis "DONE_P12"
