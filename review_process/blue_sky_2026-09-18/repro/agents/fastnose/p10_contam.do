* p10 -- does a fastnose run leave STATE that changes the next run?
* fastnose sets locals (suppressjointtests, marginsdof, fcguard...) and wooldid keeps
* globals (WOOLDID_LASTRC) and Mata functions across calls.  Run the SAME analytic spec
* before and after a fastnose run and demand bit identity; then the mirror image.
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
cap erase "`WD'/p10dump.txt"

local gN 80
local gT 11
local gonsets "2003 2005 2008"
local gysd 0.30
include "`GEN'/gen_core.do"
qui save "`WD'/p10panel.dta", replace

local SPEC "y unit t ttre , att att_it att_i att_itime subgroups(sg2) esprelength(3) espostlength(3) jointtests adversarial gteffects"
local SPECF "y unit t ttre , att att_it att_i att_itime subgroups(sg2) esprelength(3) espostlength(3) gteffects"

* order: A1 (analytic) -> F1 (fastnose) -> A2 (analytic, same spec) -> F2
foreach r in A1 F1 A2 F2 {
  local ex
  if substr("`r'",1,1) == "F" local ex fastnose
  local sp "`SPEC'"
  if "`ex'" != "" local sp "`SPECF'"
  use "`WD'/p10panel.dta", clear
  ereturn clear
  dis "########## RUN `r'"
  cap noi qui wooldid `sp' `ex'
  dis "RUN_`r'_RC=" _rc
  if _rc == 0 wdump , tag("C1_`r'") file("`WD'/p10dump.txt")
}
* mirror: fastnose FIRST in a fresh session-order sense (same session, but F before A)
foreach r in F3 A3 {
  local ex
  if substr("`r'",1,1) == "F" local ex fastnose
  local sp "`SPEC'"
  if "`ex'" != "" local sp "`SPECF'"
  use "`WD'/p10panel.dta", clear
  ereturn clear
  dis "########## RUN `r'"
  cap noi qui wooldid `sp' `ex'
  dis "RUN_`r'_RC=" _rc
  if _rc == 0 wdump , tag("C2_`r'") file("`WD'/p10dump.txt")
}
dis "DONE_P10"
