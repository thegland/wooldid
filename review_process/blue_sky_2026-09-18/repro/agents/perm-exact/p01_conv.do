* p01_conv.do -- fully enumerable permexact design; run every p-value convention and
*                dump the draws, the comparator and every reported p for an independent
*                oracle (see p01_oracle.py).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"

which wooldid
dis "BUILDCHECK"

*---------------------------------------------------------------- data
* 6 clusters (= units), 10 periods, 4 rows per unit-period (repeated cross-section, so a
* cluster-robust SE has something to see), onsets {4, 7, ., ., ., .}
*   M = 6!/(1! 1! 4!) = 30 distinct assignments
clear
set obs 6
gen long id = _n
gen double ons = .
replace ons = 4 in 1
replace ons = 7 in 2
expand 10
bysort id: gen int t = _n
expand 4
gen double y = 0.6*id + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p01_data, replace

*---------------------------------------------------------------- dumper
mata:
void wdmp(string scalar mname, string scalar fname)
{
    real matrix X
    real scalar fh, r
    if (st_matrix(mname) == J(0, 0, .)) {
        fh = fopen(fname, "w")
        fput(fh, "ABSENT")
        fclose(fh)
        return
    }
    X = st_matrix(mname)
    fh = fopen(fname, "w")
    fput(fh, invtokens(st_matrixcolstripe(mname)[., 2]', ","))
    for (r = 1; r <= rows(X); r++) fput(fh, invtokens(strofreal(X[r, .], "%26.17e"), ","))
    fclose(fh)
}
end
capture program drop dumpall
program define dumpall
  args tag
  foreach m in permdist_b permdist_se perm_comparator_b perm_comparator_se perm_p perm_n ///
               perm_p_lo perm_p_hi perm_p_se perm_p_ci perm_sd perm_ci perm_se_shift     ///
               perm_nullband b {
    cap mata: wdmp("e(`m')", "p01_`tag'_`m'.csv")
  }
  dis "TAG=`tag' NPOSS=" e(perm_npossible) " REPS=" e(perm_reps) " DONE=" e(perm_completed) ///
      " FAIL=" e(perm_failed) " EXACT=" e(perm_exact) " INCOMPLETE=<`e(perm_incomplete)'>" ///
      " ADDONE=<`e(perm_addone)'>" " OFFCENTRE=<`e(perm_offcentre)'>"
end

*---------------------------------------------------------------- runs
use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permsave(p01_asg) replace
dis "RC_default=" _rc
dumpall "default"
mat list e(perm_p)
mat list e(perm_n)

use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permabsp
dis "RC_absp=" _rc
dumpall "absp"

use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permrawp
dis "RC_rawp=" _rc
dumpall "rawp"

use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permside(left)
dis "RC_left=" _rc
dumpall "left"

use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permside(right)
dis "RC_right=" _rc
dumpall "right"

use p01_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permtstat
dis "RC_tstat=" _rc
dumpall "tstat"

*---------------------------------------------------------------- assignment file
use p01_asg, clear
desc, fullnames
export delimited using "p01_asg.csv", replace
dis "ALLDONE"
