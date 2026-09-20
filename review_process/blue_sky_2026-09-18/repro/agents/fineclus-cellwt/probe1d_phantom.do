* probe1d_phantom.do -- the PHANTOM STRADDLER.
*
* wooldid_bootrun decides the regime on `smp' (ado:8910-8916, 8991-8993), which is
* ONLY the user's `if' (==1) plus a non-missing cluster value.  It does NOT screen
* rows the estimator will drop -- a missing outcome, a zero weight, a period outside
* the usable window.  So a unit that carries a second cluster code on a row that NEVER
* ENTERS THE ESTIMATION SAMPLE still flips the whole run from the multiplicity
* bootstrap to the presence bootstrap.
*
* Arms A and D below have BIT-IDENTICAL estimation samples (the straddling row has a
* missing outcome in both) and bit-identical point estimates.  They differ only in the
* industry code recorded on that dropped row.
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
  * unit 1's period-3 outcome is MISSING in EVERY arm: the row never enters estimation
  qui replace y = . if unit == 1 & t == 3
end

tempname fh
file open `fh' using "p1d.csv", write replace
file write `fh' "arm,dseed,rseed,nest,rc,N,b_atti,se_atti,b_att,se_att" _n
file close `fh'

foreach ds in 201 202 {
  foreach arm in A D {
    qui mkdat `ds'
    * arm D: the SAME already-dropped row carries a different industry code
    if "`arm'" == "D" qui replace ind = 2 if unit == 1 & t == 3
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
      local nn = .
      if `rc' == 0 {
        mat b = e(b)
        mat V = e(V)
        local nn = e(N)
        local bi = b[1,colnumb(b,"main_att_i")]
        local si = sqrt(V[colnumb(b,"main_att_i"),colnumb(b,"main_att_i")])
        local ba = b[1,colnumb(b,"main_att")]
        local sa = sqrt(V[colnumb(b,"main_att"),colnumb(b,"main_att")])
      }
      file open `fh' using "p1d.csv", write append
      file write `fh' "`arm',`ds',`rs',`nst',`rc',`nn',`bi',`si',`ba',`sa'" _n
      file close `fh'
    }
  }
}
import delimited "p1d.csv", clear varnames(1) case(preserve)
list, sep(0) noobs
