* p06_onc.do -- EXACT ORACLE on permexact + permonsetcohorts (T1, shipped 2026-09-16).
* i IS the onset grouping (cohort), cluster() is strictly FINER than i (the unit), so each
* replication rebuilds i from the onsets it drew (ado:7568-7638).  6 units carry the onset
* multiset {5, 7, ., ., ., .}, so M = 6!/(1!1!4!) = 30 assignments.
* Part B reproduces every one of those 30 replications from first principles: build the
* drawn onset, rebuild i as "the units sharing the drawn onset" using the OBSERVED cohort
* codes (1 for onset 5, 2 for onset 7, 3 for never-treated), refit on the engine's own
* replication option path, and compare bitwise.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

clear
set obs 6
gen long unit = _n
gen double ons = .
replace ons = 5 in 1
replace ons = 7 in 2
gen byte coh = cond(ons == 5, 1, cond(ons == 7, 2, 3))
expand 10
bysort unit: gen int t = _n
expand 3
gen double y = 0.6*unit + 0.15*t + 0.9*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p06_data, replace

*---------------------------------------------------------------- A
use p06_data, clear
cap noi wooldid y coh t ons, cluster(unit) permute(shuffle) permexact permonsetcohorts
dis "RC_A=" _rc
dis "NPOSS=" e(perm_npossible) " REPS=" e(perm_reps) " DONE=" e(perm_completed) " FAIL=" e(perm_failed)
dis "IREBUILT=" e(perm_irebuilt)
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
  cap mata: wdmp("e(`m')", "p06_A_`m'.csv")
}

*---------------------------------------------------------------- B: brute force
tempname fh
file open `fh' using "p06_B_draws.csv", write replace
file write `fh' "rep,a,b,name,value" _n
local rep 0
forvalues a = 1/6 {
  forvalues b = 1/6 {
    if `a' != `b' {
      local ++rep
      use p06_data, clear
      gen double ons_r = .
      replace ons_r = 5 if unit == `a'
      replace ons_r = 7 if unit == `b'
      * the T1 rebuild, done by hand with the OBSERVED cohort codes
      gen byte coh_r = cond(ons_r == 5, 1, cond(ons_r == 7, 2, 3))
      cap qui wooldid y coh_r t ons_r, cluster(unit) fastnose
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
