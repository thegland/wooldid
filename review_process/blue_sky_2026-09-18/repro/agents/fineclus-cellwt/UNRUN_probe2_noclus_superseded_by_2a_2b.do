* probe2_noclus.do -- GATE ARCHAEOLOGY, T2.
* The removed guard (old ado:6353, error 198: "observation-level resampling cannot
* re-identify duplicated treatment cells") refused bootstrap(standard) with a
* cell-weighted estimand and NO cluster().  It is gone.  Three questions:
*   Q1 does the run now proceed, and is the reps() resolution caution still shown?
*      (RF-2, 813c414, moved that caution INSIDE `if "`clus'" != ""' -- so an
*       unclustered bootstrap should now print nothing.)
*   Q2 is wooldid's unclustered bootstrap SE equal to an observation bootstrap I
*      write myself (faithful implementation)?
*   Q3 is that SE the RIGHT one -- how does it compare with the Monte Carlo sampling
*      SD of the same estimand on the same DGP?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

* i = unit, ONE row per (i,t): every treatment cell holds exactly one observation.
prog def mkd
  args seed
  clear
  qui set obs 40
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 3 if mod(unit,4)==1
  qui replace ttre = 5 if mod(unit,4)==2
  gen double tau = 0.3 + 0.05*unit if !missing(ttre)
  qui replace tau = 0 if missing(tau)
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  qui bysort unit (t): gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.8)
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.03*t + mu + eps
end

*=================================================== Q1: rc and the reps caution
qui mkd 901
dis "---- Q1a: unclustered cell-weighted bootstrap, reps(10) ----"
cap noi wooldid y unit t ttre, att att_i att_it bootstrap(standard) reps(10) rseed(11) dots(0)
dis "Q1A_RC=" _rc
dis "---- Q1b: SAME run with cluster(unit) added, reps(10) ----"
cap noi wooldid y unit t ttre, att att_i att_it cluster(unit) bootstrap(standard) reps(10) rseed(11) dots(0)
dis "Q1B_RC=" _rc

*=================================================== Q2: wooldid vs a hand bootstrap
qui mkd 901
tempfile D
qui save `D', replace

dis "---- Q2a: wooldid's own unclustered bootstrap, reps(399) ----"
cap noi wooldid y unit t ttre, att att_i att_it bootstrap(standard) reps(399) rseed(4242) dots(0)
dis "Q2A_RC=" _rc
mat bW = e(b)
mat VW = e(V)
foreach k in att att_i att_it {
  local sW_`k' = sqrt(VW[colnumb(bW,"`k'"),colnumb(bW,"`k'")])
  local bW_`k' = bW[1,colnumb(bW,"`k'")]
  dis "WOOLDID `k' b=" %12.8f `bW_`k'' "  se=" %12.8f `sW_`k''
}

dis "---- Q2b: my own observation bootstrap, 399 draws ----"
tempname H
postfile `H' long r double hatt double hatti double hattit using p2hand, replace
local nfail 0
forval r = 1/399 {
  qui use `D', clear
  set seed `=770000 + `r''
  bsample
  cap qui wooldid y unit t ttre, att att_i att_it
  if _rc != 0 {
    local nfail = `nfail' + 1
  }
  else {
    mat bh = e(b)
    post `H' (`r') (bh[1,colnumb(bh,"main_att")]) (bh[1,colnumb(bh,"main_att_i")]) (bh[1,colnumb(bh,"main_att_it")])
  }
}
postclose `H'
dis "HANDBOOT_FAILS=`nfail'"
use p2hand, clear
qui sum hatt
dis "HAND att    sd=" %12.8f r(sd) "  n=" r(N)
qui sum hatti
dis "HAND att_i  sd=" %12.8f r(sd) "  n=" r(N)
qui sum hattit
dis "HAND att_it sd=" %12.8f r(sd) "  n=" r(N)

*=================================================== Q3: Monte Carlo sampling SD
dis "---- Q3: MC sampling SD over 300 fresh datasets ----"
tempname M
postfile `M' long m double matt double matti double mattit using p2mc, replace
local mfail 0
forval m = 1/300 {
  qui mkd `=600000 + `m''
  cap qui wooldid y unit t ttre, att att_i att_it
  if _rc != 0 local mfail = `mfail' + 1
  else {
    mat bm = e(b)
    post `M' (`m') (bm[1,colnumb(bm,"main_att")]) (bm[1,colnumb(bm,"main_att_i")]) (bm[1,colnumb(bm,"main_att_it")])
  }
}
postclose `M'
dis "MC_FAILS=`mfail'"
use p2mc, clear
foreach k in att atti attit {
  qui sum m`k'
  dis "MC `k'  mean=" %12.8f r(mean) "  sd=" %12.8f r(sd) "  n=" r(N)
}
