clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
set obs 60
gen long id = _n
gen byte clu = mod(_n,20)+1
gen ttre = .
replace ttre = 5 if id<=20
replace ttre = 8 if id>20 & id<=40
expand 12
bys id: gen t = _n
gen double y = 0.3*mod(id,5) + 0.2*t + 1.2*(t>=ttre & ttre<.) + rnormal()
foreach r in 8 40 200 {
  dis "########## reps(`r') ##########"
  cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) cluster(clu) bootstrap(standard) reps(`r') rseed(99)
  dis "  RC=" _rc
  cap mat V = e(V)
  if _rc==0 {
    local K = colsof(V)
    mata: st_local("rk", strofreal(rank(st_matrix("V"))))
    mata: st_local("mev", strofreal(min(symeigenvalues(st_matrix("V")))))
    dis "  K=`K'  rank(e(V))=`rk'  min_eig=`mev'  boot_N_V=" e(boot_N_V)
    cap noi test ES_att_post0 ES_att_post1 ES_att_post2
    dis "  TESTp=" r(p) "  df=" r(df)
  }
}
dis "RANK_DONE"
