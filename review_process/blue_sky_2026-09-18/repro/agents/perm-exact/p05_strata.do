* p05_strata.do -- EXACT ORACLE on a CROSSING-STRATA permexact design.
* NFPERM-02 (ado:7505-7517) records a defect that was fixed by keying each map-frame row
* to its generation-time matrix row: "a positional _n read at application time attaches
* assignments to the WRONG clusters whenever ascending cluster order crosses stratum order
* - silently wrong stratified inference (probed: crossing-strata exact designs reported
* p=0 with the observed assignment absent from their own enumeration)".
* This is the regression test for that fix, against a brute-force enumeration.
*
* Design: 8 clusters.  strata: clusters 1-4 are stratum 2, clusters 5-8 are stratum 1, so
* the (strata, cluster) generation order is 5 6 7 8 1 2 3 4 while frlink re-sorts to
* 1 2 3 4 5 6 7 8.  Onsets: cluster 1 -> 5 (stratum 2), cluster 5 -> 7 (stratum 1), rest
* never-treated.  Each stratum admits 4 placements, so M = 16.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

clear
set obs 8
gen long id = _n
gen byte str = cond(id <= 4, 2, 1)
gen double ons = .
replace ons = 5 if id == 1
replace ons = 7 if id == 5
expand 10
bysort id: gen int t = _n
expand 3
gen double y = 0.6*id + 0.15*t + 0.9*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p05_data, replace

*---------------------------------------------------------------- A: permexact
use p05_data, clear
cap noi wooldid y id t ons, cluster(id) permute(shuffle) permexact permstrata(str) permsave(p05_asg) replace
dis "RC_A=" _rc
dis "NPOSS=" e(perm_npossible) " REPS=" e(perm_reps) " DONE=" e(perm_completed) " FAIL=" e(perm_failed)
dis "INCOMPLETE=<`e(perm_incomplete)'>"
dis "ADDONE=<`e(perm_addone)'>"
mat list e(perm_p)
mat list e(perm_n)
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
foreach m in permdist_b perm_comparator_b perm_p perm_n perm_p_lo perm_p_hi b {
  cap mata: wdmp("e(`m')", "p05_A_`m'.csv")
}

* the saved assignments, to check each column respects its OWN stratum
use p05_asg, clear
export delimited using "p05_asg.csv", replace

*---------------------------------------------------------------- B: brute force
tempname fh
file open `fh' using "p05_B_draws.csv", write replace
file write `fh' "rep,a,b,name,value" _n
local rep 0
forvalues a = 1/4 {
  forvalues b = 5/8 {
    local ++rep
    use p05_data, clear
    gen double ons_r = .
    replace ons_r = 5 if id == `a'
    replace ons_r = 7 if id == `b'
    cap qui wooldid y id t ons_r, cluster(id) fastnose
    local rc = _rc
    if `rc' == 0 {
      tempname bb
      matrix `bb' = e(b)
      local nn : colfullnames `bb'
      local kk = colsof(`bb')
      forvalues j = 1/`kk' {
        local nmj : word `j' of `nn'
        mata: st_local("vj", strofreal(st_matrix("`bb'")[1, `j'], "%26.17e"))
        file write `fh' "`rep',`a',`b',`nmj',`vj'" _n
      }
    }
    dis "BDRAW rep=`rep' a=`a' b=`b' rc=`rc'"
  }
}
file close `fh'
dis "ALLDONE"
