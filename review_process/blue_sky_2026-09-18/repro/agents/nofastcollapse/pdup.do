* The collapse key is
*   key(`_i' `_t' `targetflag' `ifvar' `t' `ttre' `subgroups' `fckeyx')
* with   fckeyx = [cluster if impvalues] [i if contreat]   (ado:1891-1896).
* So a user who clusters on the PERIOD variable, on the ONSET variable, or - with
* contreat() - on the cohort variable `i' itself, hands gcollapse a by() list with
* the SAME variable twice.  Does that error, or silently change the cells?
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

* cluster on the PERIOD variable -> t appears twice in the key
wdcell , dta("`DD'/dC.dta") tag(d01_clt)   base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(2) espostlength(2) impvalues cluster(t))
* cluster on the ONSET variable -> ttre appears twice in the key
wdcell , dta("`DD'/dC.dta") tag(d02_clttre) base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(2) espostlength(2) impvalues cluster(ttre))
* impvalues + contreat + cluster on i -> i appears twice in fckeyx
wdcell , dta("`DD'/dC.dta") tag(d03_cli_ct) base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(2) espostlength(2) impvalues cluster(cohort) contreat(dose) contreatcontrols(basic))
* and the ordinary control: cluster on a crossing variable
wdcell , dta("`DD'/dC.dta") tag(d04_clreg_ct) base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(2) espostlength(2) impvalues cluster(region) contreat(dose) contreatcontrols(basic))
dis "PDUP_DONE"
