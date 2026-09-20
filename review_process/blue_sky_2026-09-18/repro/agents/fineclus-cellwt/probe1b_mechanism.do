* probe1b_mechanism.do -- the MECHANISM of the keep-i regime, deterministically,
* with no bootstrap and no Monte Carlo.
*
* Under wdb_nest=0 (cluster finer than, or cross-cutting, i) wooldid_bootrun sets
* __wdbnewi = i (ado:9349).  A cluster drawn TWICE therefore arrives as duplicate
* ROWS inside its units' existing cells.  Under wdb_nest=1 it sets
* __wdbnewi = group(__wdbns, i) (ado:9348), so the second copy becomes a fresh set
* of units.
*
* This file builds the two replicate frames BY HAND -- append one industry's rows to
* the panel, once with the original unit ids (keep-i) and once with fresh unit ids
* (re-id) -- and reports att_i on each.  If att_i is unchanged by the keep-i
* duplication, the equal-unit weighting is BLIND to cluster multiplicity, and the
* only thing left for the keep-i bootstrap to vary is which clusters are PRESENT.
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
  qui expand 4
  bysort ind: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 4 if uw == 1
  qui replace ttre = 6 if uw == 2
  gen double tau = 0.15*ind + 0.30*uw if !missing(ttre)
  qui replace tau = 0 if missing(tau)
  qui expand 8
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

prog def rep1
  args lab
  cap noi qui wooldid y unit t ttre, att att_i att_it
  local rc = _rc
  if `rc' != 0 {
    dis "`lab' RC=`rc'"
    exit
  }
  mat b = e(b)
  dis "`lab'  att=" %14.10f b[1,colnumb(b,"main_att")] "  att_i=" %14.10f b[1,colnumb(b,"main_att_i")] "  att_it=" %14.10f b[1,colnumb(b,"main_att_it")]
end

qui mkdat 101
tempfile D
qui save `D', replace

dis "===== BASE ====="
rep1 "BASE       "

dis "===== KEEP-I duplication: industry 1's rows appended with the SAME unit ids ====="
qui use `D', clear
qui expand 2 if ind == 1
rep1 "KEEPI_DUP1 "

dis "===== RE-ID duplication: industry 1's rows appended with FRESH unit ids ====="
qui use `D', clear
qui preserve
qui keep if ind == 1
qui replace unit = unit + 100000
tempfile C
qui save `C', replace
qui restore
qui append using `C'
rep1 "REID_DUP1  "

dis "===== KEEP-I duplication of FIVE industries (1-5), same unit ids ====="
qui use `D', clear
qui expand 2 if ind <= 5
rep1 "KEEPI_DUP5 "

dis "===== RE-ID duplication of FIVE industries (1-5), fresh unit ids ====="
qui use `D', clear
qui preserve
qui keep if ind <= 5
qui replace unit = unit + 100000
qui save `C', replace
qui restore
qui append using `C'
rep1 "REID_DUP5  "

dis "===== PRESENCE: industries 6-10 dropped (what keep-i CAN vary) ====="
qui use `D', clear
qui drop if ind > 5
rep1 "DROP_6_10  "
