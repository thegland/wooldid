* Minimise B0062: a RAW `margins' error carrying a wooldid tempvar name, no wooldid Error:
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
local BASE "att robust permute(shuffle) permreps(30) cluster(unit) permstrata(region) customweightmultiplier(w) esprelength(1) espostlength(3) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3) ccnoabsorb semielasticity nofastcollapse safety(off)"

local n 0
foreach o in "`BASE'" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) permstrata(region) customweightmultiplier(w) esprelength(1) espostlength(3) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3) ccnoabsorb nofastcollapse safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) permstrata(region) customweightmultiplier(w) esprelength(1) espostlength(3) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) esprelength(1) espostlength(3) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic)" ///
             "att robust contreat(dose) contreatcontrols(basic) lattice(baseTreatedPost 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic) lattice(all 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic) lattice(baseTreated 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic) lattice(baseTreatedPre 3) safety(off)" ///
             "att robust permute(shuffle) permreps(30) cluster(unit) contreat(dose) contreatcontrols(basic) lattice(byCohort 3) safety(off)" {
  local n = `n' + 1
  use "`D'/ds_dose_spike.dta", clear
  ereturn clear
  dis "@@MBEGIN|`n'|`o'"
  capture noisily wooldid ync unit t ttre if inrange(t,2001,2007) [aw=w], `o'
  dis "@@MEND|`n'|rc=" _rc
}
dis "MARGFAIL_DONE"
