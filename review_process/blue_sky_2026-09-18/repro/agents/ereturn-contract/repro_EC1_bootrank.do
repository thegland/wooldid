* EC-1: the bootstrap e(V) is rank-deficient at EVERY reps(), and at low reps `test`
* silently drops constraints and returns a wildly anti-conservative p.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
replace ttre = 11 if grp==4
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double dyn = 0
replace dyn = 0.5*(t-ttre+1) if post==1
gen double y = ufe + 0.25*t + 1.2*post + dyn + rnormal()

foreach R in 5 8 12 20 40 {
  cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod ///
      esprelength(3) espostlength(3) cluster(clu) bootstrap(standard) reps(`R') rseed(777)
  local rc = _rc
  mata: st_local("rk", strofreal(rank(st_matrix("e(V)"))))
  dis "reps=`R' rc=`rc' completed=" e(boot_completed) " failed=" e(boot_failed) ///
      " K=" colsof(e(b)) " rank(e(V))=`rk' properties=`e(properties)'"
  cap noi test main_att main_att_i ES_att_pre4 ES_att_pre3 ES_att_pre2 ES_att_post0 ES_att_post1
  dis "   test of 7 coefficients -> df=" r(df) "  p=" r(p)
  cap noi test ES_att_post0 ES_att_post1 ES_att_post2
  dis "   test of 3 coefficients -> df=" r(df) "  F=" r(F) "  p=" r(p)
}
* the control: the same panel with att ALONE is full rank at reps(60)
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    cluster(clu) bootstrap(standard) reps(60) rseed(11)
mata: st_local("rk", strofreal(rank(st_matrix("e(V)"))))
dis "att alone, reps=60: K=" colsof(e(b)) " rank=`rk'"
