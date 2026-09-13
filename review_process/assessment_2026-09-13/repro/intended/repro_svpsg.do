* repro_svpsg.do - adm_svpsg (C-esfn saveplots x subgroups filename scheme), one build per process
*   cd <workdir>/old ; stata-mp -e do <repro>/intended/repro_svpsg.do old <repo> <workdir>
*   cd <workdir>/new ; stata-mp -e do <repro>/intended/repro_svpsg.do new <repo> <workdir>
* ARGUMENTS (positional, all three required):
*   `1' build   - "old" or "new"; adopath ++ <workdir>/<build>
*   `2' repo    - the wooldid_dev checkout root
*   `3' workdir - scratch directory holding old/ and new/
clear all
set more off
set maxvar 20000
set graphics off
local build "`1'"
local REPO "`2'"
local W "`3'"
if "`build'" == "" | "`REPO'" == "" | "`W'" == "" {
  di as error "usage: do repro_svpsg.do <old|new> <repo-root> <workdir>"
  exit 198
}
adopath ++ "`W'/`build'"
which wooldid
local pre : dir . files "svpadm*"
foreach f of local pre {
  cap erase `f'
}
set seed 98765
set sortseed 98765
ereturn clear
use "`REPO'/test/suite/data/dta_s_admin.dta", clear
qui keep if mod(fips - 1000*floor(fips/1000), 7) == 1
cap noi wooldid y fips year ttre, att cluster(state) subgroups(region) esprelength(1) espostlength(1) makeplots saveplots(svpadm)
local rc = _rc
local files : dir . files "svpadm*"
local files : list sort files
di "REPRO_C build=`build' rc=`rc' saveplots wrote the following files (one per line):"
foreach f of local files {
  di "REPRO_C_FILE build=`build' `f'"
  cap erase `f'
}
di "REPRO_C_DONE"
exit, clear
