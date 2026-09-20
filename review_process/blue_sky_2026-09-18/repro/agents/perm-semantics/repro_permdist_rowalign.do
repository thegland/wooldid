*! SELF-CONTAINED REPRO -- e(permdist_b) and e(perm_jointW) are NOT row-aligned
*! when replications fail, and neither is indexed by permsave's _rep.
*! Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_permdist_rowalign.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
* --- a 12-cluster panel; 3 clusters treated at 2004, 9 never treated ---
clear
set obs 12
gen long unit = _n
gen double ttre = cond(inlist(unit,2,5,9), 2004, .)
gen double ai = rnormal()
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ai + 0.2*(t-2000) + 1.5*post + 0.5*rnormal()
* --- a probabilities design in which a draw can leave ZERO treated clusters
*     (p = .85^12 = .142 per draw), so a sizeable share of replications fail ---
local PP (ttre 2004 .)
forvalues u = 1/12 {
  local PP `PP' (`u' 0.15 0.85)
}
wooldid y unit t ttre, cluster(unit) att esprelength(2) espostlength(2) jointtests ///
    adversarial esfixedbaseperiod permute(probabilities) permprobs(`PP') permreps(120) ///
    rseed(31) permsave(repro_rowalign) replace
dis "requested=" e(perm_reps) "  completed=" e(perm_completed) "  failed=" e(perm_failed)
foreach m in permdist_b permdist perm_jointW perm_advW {
  matrix TMP = e(`m')
  dis "e(`m') : " rowsof(TMP) " rows x " colsof(TMP) " cols"
}
preserve
use repro_rowalign.dta, clear
dis "permsave file rows (incl _rep=0) = " _N
restore
dis ""
dis "EXPECTED: all the replication-indexed objects of one run carry the same rows,"
dis "or the help says which do not.  OBSERVED: e(permdist_b) and e(permdist) are"
dis "compacted to the completed replications while e(perm_jointW) and e(perm_advW)"
dis "keep the full requested length, so row r of one is a DIFFERENT replication from"
dis "row r of the other as soon as one replication has failed."
