* p11_joint.do -- dump the permutation JOINT pre-trend machinery for an independent check.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

use p04_data, clear
cap noi wooldid y id t ons, cluster(id) esprelength(2) espostlength(3) att jointtests ///
    permute(shuffle) permexact
dis "RC=" _rc
dis "M=" e(perm_npossible) " done=" e(perm_completed)
ereturn list
mata:
void wdmp(string scalar mname, string scalar fname)
{
    real matrix X
    real scalar fh, r
    X = st_matrix(mname)
    if (rows(X) == 0) return
    fh = fopen(fname, "w")
    fput(fh, invtokens(st_matrixcolstripe(mname)[., 2]', "|"))
    for (r = 1; r <= rows(X); r++) fput(fh, invtokens(strofreal(X[r, .], "%26.17e"), "|"))
    fclose(fh)
}
end
foreach m in perm_jointW perm_advW permdist_b perm_comparator_b perm_p {
  cap mata: wdmp("e(`m')", "p11_`m'.csv")
}
dis "ALLDONE"
