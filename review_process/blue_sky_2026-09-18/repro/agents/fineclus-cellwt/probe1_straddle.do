* probe1_straddle.do -- GATE ARCHAEOLOGY, T2/T3.
* The removed guard (old ado:6364, error 459) refused bootstrap(standard) + a
* cell-weighted estimand when cluster() did not nest i.  Its replacement is a SILENT,
* data-dependent, ALL-OR-NOTHING regime switch in wooldid_bootrun (ado:8987-8998):
*   every i wholly inside one cluster -> wdb_nest=1, __wdbnewi = group(__wdbns, i)
*      (drawn clusters re-identified: the textbook cluster bootstrap of att_i)
*   otherwise                         -> wdb_nest=0, __wdbnewi = i
*      (i kept: a twice-drawn cluster arrives as duplicate ROWS, and the equal-cell
*       weighting counts that unit ONCE, so only 0-vs->=1 PRESENCE survives)
*
* This probe changes ONE CELL of ONE VARIABLE in a 180-row panel -- unit 1's last
* period is recorded under industry 2 instead of industry 1, which is what an
* industry reclassification mid-panel looks like -- and measures the bootstrap
* standard error of att_i on both sides of it.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

*------------------------------------------------------------------ DGP
* 10 industries x 3 units x T=6.  Dependence is at the UNIT level (unit random effect
* + AR(1)); cluster(ind) is a valid, coarser choice.  Effects are heterogeneous
* across units and correlated with industry, so att_i (equal weight per unit) has
* real cohort-composition variance and a weighting error has somewhere to show.
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
  gen double tau = 0.15*ind + 0.30*uw if !missing(ttre)
  qui replace tau = 0 if missing(tau)
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  qui bysort unit (t): gen double mu = rnormal(0,0.6) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.8)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.04*t + mu + eps
end

* the nesting predicate exactly as wooldid_bootrun computes it (ado:8987-8998)
prog def nestpred
  args ivar cvar
  tempvar cg cmn cmx
  qui gegen double `cg' = group(`cvar')
  qui gegen double `cmn' = min(`cg'), by(`ivar')
  qui gegen double `cmx' = max(`cg'), by(`ivar')
  cap assert `cmn' == `cmx'
  local nst = (_rc == 0)
  dis "NESTPRED i=`ivar' cluster=`cvar' wdb_nest=`nst'"
  c_local nst `nst'
end

tempname fh
file open `fh' using "p1.csv", write replace
file write `fh' "arm,dseed,rseed,nest,rc,b_atti,se_atti,b_att,se_att" _n
file close `fh'

foreach ds in 101 102 103 {
  foreach arm in A B C {
    qui mkdat `ds'
    if "`arm'" == "B" qui replace ind = 2 if unit == 1 & t == 6
    if "`arm'" == "C" qui drop if unit == 1
    nestpred unit ind
    local nst `nst'
    foreach rs in 5001 5002 {
      dis "=== ARM `arm' dseed `ds' rseed `rs' nest=`nst' ==="
      cap noi wooldid y unit t ttre, att att_i cluster(ind) bootstrap(standard) reps(99) rseed(`rs') dots(0)
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
      file open `fh' using "p1.csv", write append
      file write `fh' "`arm',`ds',`rs',`nst',`rc',`bi',`si',`ba',`sa'" _n
      file close `fh'
    }
  }
}

import delimited "p1.csv", clear varnames(1) case(preserve)
list, sep(0) noobs
dis "---- mean bootstrap SE by arm ----"
collapse (mean) se_atti se_att b_atti (first) nest, by(arm)
list, sep(0) noobs
