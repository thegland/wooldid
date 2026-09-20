* CRITICAL: is the dCc controls() result run-to-run DETERMINISTIC?
* Four runs, two per arm, on the dataset where the effect fires.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"
local OO "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/out"
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"
local ES att att_i att_it esprelength(3) espostlength(3) jointtests
foreach r in 1 2 {
  qui use "`DD'/dCc.dta", clear
  qui wooldid y cohort t ttre, `ES' controls(xpct)
  wddump , file("`OO'/nfA`r'.txt")
  qui use "`DD'/dCc.dta", clear
  qui wooldid y cohort t ttre, `ES' controls(xpct) nofastcollapse
  wddump , file("`OO'/nfB`r'.txt")
}
* a THIRD arm-A run after a different intervening command, to break any session-order artefact
qui use "`DD'/dCc.dta", clear
qui reg y x1
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xpct)
wddump , file("`OO'/nfA3.txt")
* and the same in a fresh-order run: B first
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xpct) nofastcollapse
wddump , file("`OO'/nfB3.txt")
dis "PNF_DONE"
