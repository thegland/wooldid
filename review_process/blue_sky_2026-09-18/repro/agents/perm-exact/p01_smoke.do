* p01_smoke.do -- assert the pinned build, run a fully enumerable permexact design,
*                 dump every permutation e() object for an independent oracle.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"

which wooldid
dis "BUILDCHECK done"

*---------------------------------------------------------------- data
* 6 clusters (= units), 10 periods, onsets {4, 7, ., ., ., .}
*   M = 6!/(1! 1! 4!) = 30 distinct assignments
clear
set obs 6
gen long id = _n
gen double ons = .
replace ons = 4 in 1
replace ons = 7 in 2
expand 10
bysort id: gen int t = _n
gen double y = 0.6*id + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
save p01_data, replace

*---------------------------------------------------------------- run
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permsave(p01_asg) replace
dis "RUN_RC=" _rc

dis "NPOSSIBLE=" e(perm_npossible)
dis "NPOSSTR=" "`e(perm_npossible_str)'"
dis "REPS=" e(perm_reps)
dis "COMPLETED=" e(perm_completed)
dis "FAILED=" e(perm_failed)
dis "EXACTFLAG=" e(perm_exact)
dis "INCOMPLETE=" "`e(perm_incomplete)'"
dis "ADDONE=" "`e(perm_addone)'"
dis "OFFCENTRE=" "`e(perm_offcentre)'"

mat list e(perm_p)
mat list e(perm_n)
mat list e(perm_p_lo)
mat list e(perm_p_hi)
mat list e(perm_p_se)
mat list e(perm_comparator_b)

*---------------------------------------------------------------- dump
capture program drop dumpmat
program define dumpmat
  args mname fname
  mata: wdmp("`mname'", "`fname'")
end
mata:
void wdmp(string scalar mname, string scalar fname)
{
  real matrix X
  real scalar fh, r
  string rowvector cs
  X = st_matrix(mname)
  fh = fopen(fname, "w")
  cs = st_matrixcolstripe(mname)[., 2]'
  fput(fh, invtokens(cs, ","))
  for (r = 1; r <= rows(X); r++) fput(fh, invtokens(strofreal(X[r, .], "%21x"), ","))
  fclose(fh)
}
end

dumpmat "e(permdist_b)" "p01_permdist_b.csv"
dumpmat "e(perm_comparator_b)" "p01_cmp_b.csv"
dumpmat "e(perm_p)" "p01_p.csv"
dumpmat "e(perm_n)" "p01_n.csv"
dumpmat "e(perm_p_lo)" "p01_plo.csv"
dumpmat "e(perm_p_hi)" "p01_phi.csv"
dumpmat "e(perm_sd)" "p01_sd.csv"
dumpmat "e(perm_ci)" "p01_ci.csv"
dumpmat "e(b)" "p01_b.csv"

dis "DUMP done"

* the saved assignment file: onset_ columns, one per cluster
use p01_asg, clear
desc, fullnames
list _rep onset_1 onset_2 onset_3 onset_4 onset_5 onset_6 in 1/10, noobs
export delimited using "p01_asg.csv", replace
dis "ALLDONE"
