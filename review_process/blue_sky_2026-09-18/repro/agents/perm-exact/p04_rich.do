* p04_rich.do -- EXACT ORACLE on a design with MANY permuted elements.
* 8 clusters, onsets {5, 7, ., ., ., ., ., .}  =>  M = 8!/(1!1!6!) = 56 assignments.
* Event-study horizons + 4 estimands + subgroups + jointtests, so the permutation kernel
* has dozens of columns, several of which are NOT estimable under every assignment
* (the "incomplete element" branch, where the exact share is taken over n < M).
* Part A runs permexact.  Part B enumerates all 56 assignments independently and refits
* on the engine's own replication option path (`fastnose'), dumping the full named e(b)
* of every assignment.  p04_oracle.py reconciles.
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
gen byte sg = mod(id, 2)
gen double ons = .
replace ons = 5 in 1
replace ons = 7 in 2
expand 10
bysort id: gen int t = _n
expand 3
gen double y = 0.6*id + 0.15*t + 0.9*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p04_data, replace

local COMMON esprelength(2) espostlength(3) att att_i att_it att_itime subgroups(sg) jointtests

*---------------------------------------------------------------- A: permexact
use p04_data, clear
cap noi wooldid y id t ons, cluster(id) `COMMON' permute(shuffle) permexact
dis "RC_A=" _rc
dis "NPOSS=" e(perm_npossible) " REPS=" e(perm_reps) " DONE=" e(perm_completed) " FAIL=" e(perm_failed)
dis "INCOMPLETE=<`e(perm_incomplete)'>"
dis "ADDONE=<`e(perm_addone)'>"
dis "OFFCENTRE=<`e(perm_offcentre)'>"

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
foreach m in permdist_b perm_comparator_b perm_p perm_n perm_p_lo perm_p_hi perm_p_se perm_sd b {
  cap mata: wdmp("e(`m')", "p04_A_`m'.csv")
}
mat list e(perm_p)
mat list e(perm_n)

*---------------------------------------------------------------- B: brute force
* the 56 assignments: cluster a takes onset 5, cluster b != a takes onset 7
tempname fh
file open `fh' using "p04_B_draws.csv", write replace
file write `fh' "rep,a,b,name,value" _n
local rep 0
forvalues a = 1/8 {
  forvalues b = 1/8 {
    if `a' != `b' {
      local ++rep
      use p04_data, clear
      gen double ons_r = .
      replace ons_r = 5 if id == `a'
      replace ons_r = 7 if id == `b'
      cap qui wooldid y id t ons_r, cluster(id) `COMMON' fastnose
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
}
file close `fh'
dis "ALLDONE"
