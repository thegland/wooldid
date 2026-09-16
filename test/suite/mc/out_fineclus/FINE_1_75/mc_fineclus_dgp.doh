* mc_fineclus_dgp.doh - the two designs of the finer-than-i feature readouts (T2/T3,
* prereg PREREG_FEATURE_FINECLUS_2026-09-16.md). Included by mc_fineclus.do (the F3
* coverage driver) and mc_fineclus_ctl.do (the F4 negative control) so both see one
* definition, the shape mc_bootse2_dgp.doh established.
*
* In both, i = state is a COHORT and cluster = unit is a level inside it, which is the
* pairing the feature is about: i says which cohorts the effects average over, cluster()
* says where dependence and assignment sit, and they are separate claims.
*
* dgpE - the F3 coverage design. i = state (40 cohorts), cluster = unit (160 units)
* is STRICTLY FINER than i; cluster = state IS i and is the nesting (regression) leg.
*
* The cohort effects are genuinely heterogeneous AND correlated with cohort size, so
* att_i (equal weight per state) and att (treated-observation weighted) differ by a
* wide margin - a DGP where they coincide cannot detect a weighting error:
*   treated states (s odd): mod(s,4)==1 -> 2 units, onset 3, tau = 0.4
*                           mod(s,4)==3 -> 6 units, onset 4, tau = 1.2
*   control states (s even): 4 units, never treated
* T = 6, 960 rows. 40 states -> 40 nesting clusters; 160 units -> 160 finer clusters;
* 20 of 40 states are treated, so neither leg is a few-treated-cluster design.
* Dependence is at the UNIT level only - a unit random effect plus a within-unit
* AR(1) - so clustering at unit is the correct level and state is valid but coarser.
*   dgpE seed mult [noise]     noise 0 turns the errors off (estimand readout)
cap prog drop dgpE
prog def dgpE
  args seed mult noise
  if "`noise'" == "" local noise 1
  clear
  qui set obs 40
  gen long state = _n
  gen int nu = 4
  qui replace nu = 2 if mod(state,4)==1
  qui replace nu = 6 if mod(state,4)==3
  qui expand nu
  bysort state: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 3 if mod(state,4)==1
  qui replace ttre = 4 if mod(state,4)==3
  gen double tau = `mult'*0.2*nu if mod(state,2)==1
  qui replace tau = 0 if missing(tau)
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  qui by unit: gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.7)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  qui replace mu = 0 if `noise' == 0
  qui replace eps = 0 if `noise' == 0
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.05*t + mu + eps
end
* dgpD - dgpC with larger states (n_s = s + 3, so 114 units over 12 states, T=8,
* 912 rows). i = state; cluster = unit is STRICTLY FINER than i, cluster = state IS
* i. Cohort effects are heterogeneous (tau_s = mult*0.1*s on the odd, treated
* states), so att_i (equal weight per state) differs from att: on the zero-noise
* draw att = 0.70761905 and att_i = 0.6 exactly. Dependence is at the UNIT level
* only - a unit random effect plus a within-unit AR(1) - so clustering at unit is
* the correct level and clustering at state is valid but coarser.
*   dgpD seed mult [noise]     noise 0 turns the errors off (estimand readout)
cap prog drop dgpD
prog def dgpD
  args seed mult noise
  if "`noise'" == "" local noise 1
  clear
  qui set obs 12
  gen long state = _n
  qui gen int nu = state + 3
  qui expand nu
  bysort state: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 4 if mod(state,4)==1
  qui replace ttre = 6 if mod(state,4)==3
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
  qui replace mu = 0 if `noise' == 0
  qui replace eps = 0 if `noise' == 0
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.05*t + mu + eps
end
