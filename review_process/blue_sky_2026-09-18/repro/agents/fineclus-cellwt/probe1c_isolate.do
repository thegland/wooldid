* probe1c_isolate.do -- the same ONE-CELL regime flip as probe1_straddle.do, on a
* design that ISOLATES the component keep-i loses.
*
* Idiosyncratic noise is turned down to sd 0.02 and the treatment effect is made a
* pure function of the CLUSTER (tau = 0.5*ind).  Then essentially all of att_i's
* bootstrap variance is cluster-composition variance, and the two regimes differ by
* exactly the thing the design note names: whether a twice-drawn cluster's units are
* counted twice (nest, `group(__wdbns, i)`) or once (keep-i, `__wdbnewi = i`).
*
* Prediction from the mechanism: the correct cluster bootstrap has variance
* Var_g(tau)/G; the presence bootstrap keeps only the 0-vs->=1 part, which for large G
* is (1-1/e)^-1 * (1/e) = 0.582 of it, i.e. an SE ratio of about 0.76.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

prog def mkdat
  args seed
  clear
  qui set obs 10
  gen long ind = _n
  qui expand 3
  bysort ind: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 4 if uw == 1
  qui replace ttre = 6 if uw == 2
  gen double tau = 0.5*ind if !missing(ttre)
  qui replace tau = 0 if missing(tau)
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  gen double eps = rnormal(0,0.02)
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.04*t + eps
end

tempname fh
file open `fh' using "p1c.csv", write replace
file write `fh' "arm,dseed,rseed,nest,rc,b_atti,se_atti,b_att,se_att" _n
file close `fh'

foreach ds in 201 202 203 {
  foreach arm in A B {
    qui mkdat `ds'
    if "`arm'" == "B" qui replace ind = 2 if unit == 1 & t == 6
    tempvar cg cmn cmx
    qui gegen double `cg' = group(ind)
    qui gegen double `cmn' = min(`cg'), by(unit)
    qui gegen double `cmx' = max(`cg'), by(unit)
    cap assert `cmn' == `cmx'
    local nst = (_rc == 0)
    qui drop `cg' `cmn' `cmx'
    foreach rs in 6001 6002 {
      dis "=== ARM `arm' ds `ds' rs `rs' nest=`nst' ==="
      cap noi wooldid y unit t ttre, att att_i cluster(ind) bootstrap(standard) reps(199) rseed(`rs') dots(0)
      local rc = _rc
      dis "RC_`arm'_`ds'_`rs'=`rc'"
      local bi = .
      local si = .
      local ba = .
      local sa = .
      if `rc' == 0 {
        mat b = e(b)
        mat V = e(V)
        local bi = b[1,colnumb(b,"main_att_i")]
        local si = sqrt(V[colnumb(b,"main_att_i"),colnumb(b,"main_att_i")])
        local ba = b[1,colnumb(b,"main_att")]
        local sa = sqrt(V[colnumb(b,"main_att"),colnumb(b,"main_att")])
      }
      file open `fh' using "p1c.csv", write append
      file write `fh' "`arm',`ds',`rs',`nst',`rc',`bi',`si',`ba',`sa'" _n
      file close `fh'
    }
  }
}
import delimited "p1c.csv", clear varnames(1) case(preserve)
list, sep(0) noobs
collapse (mean) se_atti se_att, by(arm)
list, sep(0) noobs
