* Isolate the raw `margins' failure.  Is the ANALYTIC path affected, or only the engine?
clear all
set more off
set processors 2
set graphics off
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

dis "@@MBEGIN|a|analytic i=unit + semielasticity + safety(off)"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync unit t ttre if inrange(t,2001,2007) [aw=w], att robust contreat(dose) contreatcontrols(basic) semielasticity safety(off)
dis "@@MEND|a|rc=" _rc

dis "@@MBEGIN|b|i=cohort, NO safety(off), permute + semielasticity"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync cohort t ttre [aw=w], att cluster(cohort) contreat(dose) contreatcontrols(basic) semielasticity permute(shuffle) permreps(20)
dis "@@MEND|b|rc=" _rc

dis "@@MBEGIN|c|minimal permute i=unit + semielasticity + safety(off)"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync unit t ttre if inrange(t,2001,2007) [aw=w], att robust contreat(dose) contreatcontrols(basic) semielasticity permute(shuffle) permreps(20) cluster(unit) safety(off)
dis "@@MEND|c|rc=" _rc

dis "@@MBEGIN|d|c + lattice"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync unit t ttre if inrange(t,2001,2007) [aw=w], att robust contreat(dose) contreatcontrols(basic) semielasticity permute(shuffle) permreps(20) cluster(unit) lattice(baseTreatedPost 3) safety(off)
dis "@@MEND|d|rc=" _rc

dis "@@MBEGIN|e|d but full time window, no if"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync unit t ttre [aw=w], att robust contreat(dose) contreatcontrols(basic) semielasticity permute(shuffle) permreps(20) cluster(unit) lattice(baseTreatedPost 3) safety(off)
dis "@@MEND|e|rc=" _rc

dis "@@MBEGIN|f|c with NO weights"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync unit t ttre if inrange(t,2001,2007), att robust contreat(dose) contreatcontrols(basic) semielasticity permute(shuffle) permreps(20) cluster(unit) safety(off)
dis "@@MEND|f|rc=" _rc

dis "@@MBEGIN|g|c but contreatelasticitytype(eydx) instead of semielasticity"
use "`D'/ds_dose_spike.dta", clear
ereturn clear
capture noisily wooldid ync unit t ttre if inrange(t,2001,2007) [aw=w], att robust contreat(dose) contreatcontrols(basic) contreatelasticitytype(eydx) permute(shuffle) permreps(20) cluster(unit) safety(off)
dis "@@MEND|g|rc=" _rc
dis "MARGFAIL2_DONE"
