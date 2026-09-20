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

dis "WHICH:"
which wooldid

* ---------------- build ONE panel used by every spec ----------------
local gN 80
local gT 10
local gonsets "2004 2006 2008"
local gysd 0.30
local gdose unitconst
local gwtype unit
local gpois 1
include "`GEN'/gen_core.do"
dis "NOBS=" _N
qui save "`WD'/p2panel.dta", replace

cap erase "`WD'/p2dump.txt"

* --- e(scalars) sanity: does the empty-member-list behaviour belong to wooldid? ---
qui reg y x1 x2
dis "REG_SCALARS: `e(scalars)'"
dis "REG_MACROS: `e(macros)'"
dis "REG_MATRICES: `e(matrices)'"
ereturn clear
cap noi qui wooldid y unit t ttre , esprelength(2) espostlength(2)
dis "WDID_SCALARS: `e(scalars)'"
dis "WDID_MACROS: `e(macros)'"
dis "WDID_MATRICES: `e(matrices)'"
dis "WDID_properties=" e(properties)
tempname GTX
cap matrix `GTX' = e(gteffects)
dis "WDID_gteffects_rc=" _rc



* ---------------- the spec matrix ----------------
local nspec 20
local sp1  "y unit t ttre ,"
local sp2  "y unit t ttre , att att_it att_i att_itime"
local sp3  "y unit t ttre , esprelength(3) espostlength(4)"
local sp4  "y unit t ttre , esprelength(3) espostlength(4) esfixedbaseperiod esrelativeto(-2)"
local sp5  "y unit t ttre , subgroups(sg2)"
local sp6  "y unit t ttre , subgroups(sg2) esprelength(2) espostlength(3)"
local sp7  "y unit t ttre , subgroups(sg2) esprelength(2) espostlength(3) att att_it att_i att_itime"
local sp8  "y cohort t ttre , fe(unit) esprelength(2) espostlength(3) att att_it att_i att_itime"
local sp9  "y unit t ttre , controls(x1 x2) esprelength(2) espostlength(3)"
local sp10 "y unit t ttre [aw=w] , esprelength(2) espostlength(3) att att_it att_i att_itime"
local sp11 "y unit t ttre if x2 > -0.5 , esprelength(2) espostlength(3)"
local sp12 "ycount unit t ttre , poisson esprelength(2) espostlength(3)"
local sp13 "y unit t ttre , gteffects esprelength(2) espostlength(3)"
local sp14 "y unit t ttre , cluster(region) esprelength(2) espostlength(3)"
local sp15 "y unit t ttre , robust esprelength(2) espostlength(3)"
local sp16 "y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) esprelength(2) espostlength(3)"
local sp17 "y cohort t ttre , contreat(dose) contreatcontrols(untreatedZero) fe(unit) makecfxplots cfxplottypes(0.5 1 2) esprelength(2) espostlength(3)"
local sp18 "y unit t ttre , interactivecontrols(x2) esprelength(2) espostlength(3)"
local sp19 "y unit t ttre , coarsecohortcontrols esprelength(2) espostlength(3)"
local sp20 "y unit t ttre , subgroups(sgtiny) esprelength(2) espostlength(3)"

forval s = 1/`nspec' {
  local cmd "`sp`s''"
  dis "########## SPEC `s' : `cmd'"
  * ---- analytic ----
  use "`WD'/p2panel.dta", clear
  ereturn clear
  dis "--- SPEC`s' ANALYTIC ---"
  cap noi wooldid `cmd'
  local rcA = _rc
  dis "SPEC`s'_AN_RC=`rcA'"
  if `rcA' == 0 wdump , tag("S`s'_AN") file("`WD'/p2dump.txt")
  * ---- fastnose ----
  use "`WD'/p2panel.dta", clear
  ereturn clear
  dis "--- SPEC`s' FASTNOSE ---"
  local hascomma = strpos("`cmd'", ",")
  cap noi wooldid `cmd' fastnose
  local rcF = _rc
  dis "SPEC`s'_FN_RC=`rcF'"
  if `rcF' == 0 wdump , tag("S`s'_FN") file("`WD'/p2dump.txt")
  dis "SPEC`s'_DONE rcA=`rcA' rcF=`rcF'"
}
dis "DONE_P2"
