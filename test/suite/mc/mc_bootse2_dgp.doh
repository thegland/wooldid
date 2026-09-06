* mc_bootse2_dgp.doh - the two designs of the F058 bootstrap calibration extension.
* Included by mc_bootse2.do and by smoke2/verify2.do so both see one definition.
*
* dgpA seed tau : mc_bootse.do's DESIGN A. 24 clusters of 2 units, T=6, clusters
*   strictly above the unit, cohorts ttre in {3,5,never} by mod(cl,3). Constant
*   effect tau on every treated post cell, so the true att is exactly tau under
*   any event-study window.
* dgpB seed mult : the UNEQUAL-cluster design. 12 states, state s holds s units
*   (78 units), T=8, odd states treated at t=5, tau_s = mult*0.1*n_s. Every treated
*   unit has the same 4 post periods, so the att, the att_it and any windowed
*   version of them all equal the unit-weighted mean of tau over treated units:
*   0.1*sum_{s odd} s^2 / sum_{s odd} s = 0.1*286/36 = 0.7944444...
* Both carry the same noise: unit effect, within-unit AR(1) (stationary start) and
* a cluster-by-period shock - the intra-cluster dependence that survives unit FEs
* and the reason a clustered variance is needed at all.

cap prog drop dgpA
prog def dgpA
  args seed tau
  clear
  qui set obs 48
  gen long unit = _n
  gen long cl   = ceil(unit/2)
  gen byte uw   = unit - (cl-1)*2
  gen int ttre = .
  qui replace ttre = 3 if mod(cl,3)==1
  qui replace ttre = 5 if mod(cl,3)==2
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  qui by unit: gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.7)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  tempvar clt
  qui egen long `clt' = group(cl t)
  sort `clt'
  gen double s = 0
  qui by `clt': replace s = rnormal(0,0.5) if _n==1
  qui by `clt': replace s = s[1]
  sort unit t
  gen double tau = `tau'
  gen double y = tau*(t>=ttre & ttre<.) + 0.05*t + mu + eps + s
end

cap prog drop dgpB
prog def dgpB
  args seed mult
  clear
  qui set obs 12
  gen long state = _n
  qui expand state
  bysort state: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 5 if mod(state,2)==1
  gen double tau = `mult'*0.1*state if mod(state,2)==1
  qui replace tau = 0 if missing(tau)
  qui expand 8
  bysort unit: gen int t = _n
  set seed `seed'
  qui by unit: gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.7)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  tempvar clt
  qui egen long `clt' = group(state t)
  sort `clt'
  gen double s = 0
  qui by `clt': replace s = rnormal(0,0.5) if _n==1
  qui by `clt': replace s = s[1]
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.05*t + mu + eps + s
end
