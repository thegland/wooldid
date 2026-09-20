* p8 -- L3: wooldid_fcdeg (ado:5768) reads r(V).  Under `nose' K==0 and it returns at
* once, so the cell-collapse degeneracy fallback cannot fire in fast mode.  When it fires
* in the ANALYTIC run the code does `restore, preserve' and evaluates every later margins
* call on the FULL data; fast mode stays on the collapsed rows.  Different summation
* order -> possibly different numbers, and (per the fcdeg caution text) possibly a
* different set of reported subgroup horizons.
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

local gN 100
local gT 12
local gonsets "2003 2005 2007 2009"
local gysd 0.30
include "`GEN'/gen_core.do"
qui save "`WD'/p8panel.dta", replace
cap erase "`WD'/p8dump.txt"

* few clusters, many horizons, two subgroups: the geometry the fcdeg comment names
local q1 "y unit t ttre , cluster(cohort) subgroups(sg2) esprelength(5) espostlength(6) att att_it"
local q2 "y unit t ttre , cluster(region) subgroups(sg2) esprelength(5) espostlength(6) att att_it"
local q3 "y cohort t ttre , fe(unit) cluster(cohort) subgroups(sg2) esprelength(5) espostlength(6) att att_it"
* and the same three with the collapse OFF on both sides -- isolates "collapse vs full"
* from "nose vs df()"
local q4 "y unit t ttre , cluster(cohort) subgroups(sg2) esprelength(5) espostlength(6) att att_it nofastcollapse"

forval q = 1/4 {
  foreach mode in AN FN {
    local ex
    if "`mode'" == "FN" local ex fastnose
    dis "########## Q`q' `mode'"
    use "`WD'/p8panel.dta", clear
    ereturn clear
    cap noi wooldid `q`q'' `ex'
    dis "Q`q'_`mode'_RC=" _rc
    if _rc == 0 {
      dis "Q`q'_`mode'_fcol=" e(wdid_fcol)
      dis "Q`q'_`mode'_fcguard=" e(wdid_fcguard)
      wdump , tag("Q`q'_`mode'") file("`WD'/p8dump.txt")
    }
  }
}
dis "DONE_P8"
