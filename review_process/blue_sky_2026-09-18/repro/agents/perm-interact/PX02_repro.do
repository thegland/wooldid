*==============================================================================
* PX-02 self-contained repro.
* customvce() can NEVER reach a permutation replication, so the customvce()
* forwarding line in the permute rep-option string (ado:128) is dead code, and
* the "permute forwards customvce while bootstrap refuses it" asymmetry has no
* behavioural consequence.  Both doors are shut:
*   with cluster():    ado:4667  "customvce() is not compatible with robust and cluster()"
*   without cluster(): wooldid_permctl  "permute() requires cluster() with a single variable"
* and permute() admits exactly one cluster variable, so there is no third case.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260916
set sortseed 20260916
clear
set obs 40
gen long unit = _n
gen double ttre = .
replace ttre = 2003 if unit <= 12
replace ttre = 2005 if unit > 12 & unit <= 24
expand 8
bysort unit: gen int t = 1999 + _n
gen double a_i = rnormal()
gen double g_t = 0.25*(t - 2002)
gen double y = a_i + g_t + 2.0*(ttre < . & t >= ttre) + rnormal()*0.6

dis _n "### A. customvce() WITH cluster() -> refused by the customvce/cluster screen"
cap noi wooldid y unit t ttre , att cluster(unit) customvce(cluster unit) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "PX02|with_cluster|RC=" _rc

dis _n "### B. customvce() WITHOUT cluster() -> refused by wooldid_permctl"
cap noi wooldid y unit t ttre , att customvce(cluster unit) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "PX02|no_cluster|RC=" _rc
cap noi wooldid y unit t ttre , att customvce(robust) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "PX02|no_cluster_robust|RC=" _rc

dis _n "### C. for contrast, the bootstrap twin's explicit refusal"
cap noi wooldid y unit t ttre , att customvce(cluster unit) bootstrap(standard) reps(30) rseed(777) dots(0)
dis "PX02|boot|RC=" _rc
dis "PX02_DONE"
