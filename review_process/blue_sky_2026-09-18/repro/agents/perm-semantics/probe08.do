* probe08 -- FAILURE ACCOUNTING, read straight off the permsave dump.
* The help says a failed replication "keeps its row, all missing" and that failures are
* excluded and counted in e(perm_failed).  A replication that drew ZERO treated clusters
* cannot be estimated, so #\{zero-treated draws\} must equal e(perm_failed).
clear all
set more off
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
use f_prob.dta, clear
dis "rows in dump (incl _rep=0): " _N
qui ds onset_*
local ons `r(varlist)'
dis "onset columns: `: word count `ons''"
qui ds onset_* _rep, not
local bc `r(varlist)'
dis "coef columns: `: word count `bc''"
egen byte ntreat = rownonmiss(`ons')
egen byte nbmiss = rowmiss(`bc')
local nb : word count `bc'
gen byte allmiss = (nbmiss == `nb')
gen byte anymiss = (nbmiss > 0 & nbmiss < `nb')
tab ntreat allmiss, m
count if _rep > 0
local NR = r(N)
count if _rep > 0 & allmiss
dis "ALLMISS rows (excluding rep 0) = " r(N) "   of `NR'"
count if _rep > 0 & ntreat == 0
dis "ZERO-TREATED draws            = " r(N)
count if _rep > 0 & ntreat == 0 & !allmiss
dis "ZERO-TREATED but NOT all-missing (would be a defect) = " r(N)
count if _rep > 0 & ntreat > 0 & allmiss
dis "TREATED-draws that still failed = " r(N)
count if _rep > 0 & anymiss
dis "PARTIALLY missing rows = " r(N)
list _rep ntreat nbmiss in 1/12, noobs
* the dataset characteristics the help promises
char list _dta[]
dis "PROBE08_DONE"
