* p11 -- FN-01 follow-up: with the RNG pinned to the SAME seed before each call, is
* permute()+permanalyticse+fastnose numerically identical to the permute() fast-mode
* default (so the defect is purely the label), or does it report different inference?
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fastnose"
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
run "`WD'/wdump.do"
which wooldid
cap erase "`WD'/p11dump.txt"

local gN 90
local gT 10
local gonsets "2003 2005 2007"
local gysd 0.30
include "`GEN'/gen_core.do"
qui save "`WD'/p11panel.dta", replace

local base "y cohort t ttre , fe(unit) cluster(cohort) permute(shuffle) permreps(60) esprelength(2) espostlength(2)"

* R1 = fast-mode default ; R2 = permanalyticse ; R3 = permanalyticse + fastnose
local o1 ""
local o2 "permanalyticse"
local o3 "permanalyticse fastnose"
forval k = 1/3 {
  use "`WD'/p11panel.dta", clear
  ereturn clear
  set seed 777777
  set sortseed 777777
  dis "########## R`k' : `o`k''"
  cap noi wooldid `base' `o`k''
  dis "R`k'_RC=" _rc
  if _rc == 0 {
    dis "R`k'_vcvstate=" e(wdid_vcvstate)
    foreach m in perm_p perm_ci perm_sd perm_se_shift perm_n perm_comparator_b {
      tempname MM
      cap matrix `MM' = e(`m')
      if _rc == 0 {
        dis "R`k'_`m':"
        matrix list `MM', format(%21x)
        cap matrix drop `MM'
      }
      else dis "R`k'_`m'_ABSENT"
    }
    wdump , tag("R`k'_X") file("`WD'/p11dump.txt")
  }
}
dis "DONE_P11"
