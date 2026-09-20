*! r03b — REFUTE the ATTRIBUTION in FCW-03: it says the cohort-composition term is
*! "the nesting re-id leg injecting" it.  wooldid only takes the re-id leg when a
*! cell-weighted estimand (or boottstat+cluster) is asked for (ado:378).  A run asking
*! for att ALONE at cluster()==i takes the COUNT path instead.  If the composition term
*! is a property of the re-id leg, se(att) must collapse on the count path.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkE
prog def mkE
  syntax , dseed(int) het(int) [esd(real .02)]
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs 40
  gen long state = _n
  gen int nu = cond(mod(state,4)==1, 2, cond(mod(state,4)==3, 6, 4))
  gen onset = cond(mod(state,2)==1, 4, .)
  gen double tau = cond(`het'==1, cond(mod(state,4)==1, .4, 1.2), .8)
  replace tau = 0 if missing(onset)
  expand nu
  sort state
  gen long unit = _n
  gen double uf = rnormal(0,.5)
  expand 6
  bysort unit (state): gen int t = _n
  gen byte D = (onset < . & t >= onset)
  gen double y = 2 + uf + .15*t + tau*D + rnormal(0,`esd')
  keep y state unit t onset
end

foreach h in 1 0 {
  dis "===== REID leg (att att_i): het=`h'"
  mkE , dseed(9001) het(`h')
  cap noi wooldid y state t onset, att att_i cluster(state) bootstrap(standard) reps(199) rseed(51) dots(0)
  dis "REID_RC=" _rc " het=`h'"
  mat V = e(V)
  dis "REID het=`h' se_att=" %12.0g sqrt(V[1,1]) " se_atti=" %12.0g sqrt(V[2,2])

  dis "===== COUNT path (att alone): het=`h'"
  mkE , dseed(9001) het(`h')
  cap noi wooldid y state t onset, att cluster(state) bootstrap(standard) reps(199) rseed(51) dots(0)
  dis "COUNT_RC=" _rc " het=`h'"
  mat V2 = e(V)
  dis "COUNT het=`h' se_att=" %12.0g sqrt(V2[1,1])
}
dis "predicted composition SE sqrt(Var_s/20) = " sqrt(.16/20)
dis "===== DONE"
