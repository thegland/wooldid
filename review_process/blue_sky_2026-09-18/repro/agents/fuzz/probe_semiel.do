* Pin down the semielasticity / contreat crash.  Analytic path, no engine.
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
local n 0
foreach spec in ///
 "dose_spike|unit|robust|contreat(dose) contreatcontrols(basic) semielasticity safety(off)" ///
 "dose_spike|unit|robust|contreat(dose) contreatcontrols(basic) semielasticity nofastcollapse safety(off)" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity nofastcollapse" ///
 "dose_spike|cohort|robust|contreat(dose) contreatcontrols(basic) semielasticity" ///
 "dose_spike|cohort||contreat(dose) contreatcontrols(basic) semielasticity" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) contreatelasticitytype(eydx)" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) contreatelasticitytype(eyex)" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) contreatelasticitytype(dyex)" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) contreatelasticitytype(dydx)" ///
 "dose_lattice|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity" ///
 "dose_cellconst|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity" ///
 "dose_lattice|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity nofastcollapse" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(byCohort) semielasticity" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(untreatedZero) semielasticity" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity esprelength(2) espostlength(2)" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity lattice(byCohort 3)" ///
 "dose_spike|cohort|cluster(cohort)|contreat(dose) contreatcontrols(basic) semielasticity lattice(byCohort 3) nofastcollapse" {
  local n = `n' + 1
  local ds  = substr("`spec'",1,strpos("`spec'","|")-1)
  local rr  = substr("`spec'",strpos("`spec'","|")+1,.)
  local iv  = substr("`rr'",1,strpos("`rr'","|")-1)
  local rr2 = substr("`rr'",strpos("`rr'","|")+1,.)
  local vc  = substr("`rr2'",1,strpos("`rr2'","|")-1)
  local oo  = substr("`rr2'",strpos("`rr2'","|")+1,.)
  use "`D'/ds_`ds'.dta", clear
  ereturn clear
  dis "@@SBEGIN|`n'|wooldid ync `iv' t ttre, att `vc' `oo'"
  capture noisily wooldid ync `iv' t ttre, att `vc' `oo'
  dis "@@SEND|`n'|rc=" _rc
}
dis "SEMIEL_DONE"
