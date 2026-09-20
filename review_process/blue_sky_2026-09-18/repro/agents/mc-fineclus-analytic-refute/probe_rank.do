* probe_rank.do -- is the pre block of e(V_es_att) rank-deficient / badly conditioned at
* K=3 clusters per treated cell with cluster() STRICTLY FINER than i?  That is the
* candidate mechanism behind the 62% joint pre-trend rejection rate in mc1_ladder.do.
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
cap program drop gen3
program define gen3
  args K
  clear
  qui set obs 40
  gen long state = _n
  gen double s_g = rnormal()*0.8
  gen double ttre = cond(state<=`K', 5, cond(state<=2*`K', 7, .))
  qui expand 4
  bysort state: gen long unit = (state-1)*4 + _n
  gen double a_u = rnormal()*0.5
  qui expand 12
  bysort unit: gen int t = _n
  tempvar fst vv
  bysort state t: gen byte `fst' = (_n==1)
  gen double `vv' = rnormal()*0.4 if `fst'
  bysort state t (`fst'): replace `vv' = `vv'[_N]
  gen double v_st = `vv'
  gen double e = rnormal()*0.4
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_st + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  keep y state unit t ttre coh
end
foreach K in 3 10 {
  foreach iv in coh state {
    set seed 310041
    set sortseed 310041
    qui gen3 `K'
    cap qui wooldid y `iv' t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
    di "=== K=`K'  i=`iv'  rc=" _rc "  joint pre-trend p = " %8.6f e(joint_horiz)
    cap noi matrix VV = e(V_es_att)
    cap noi matrix PB = VV[1..4, 1..4]
    cap noi matrix list PB, format(%12.3e)
    cap noi mata: PBm = st_matrix("PB")
    cap noi mata: printf("   rank = %f   cond = %12.4e   det = %12.4e\n", rank(PBm), cond(PBm), det(PBm))
    cap noi mata: printf("   eigenvalues: "); mata: st_local("junk","")
    cap noi mata: e = symeigenvalues(makesymmetric(PBm)); e
  }
}
di "PROBE_RANK_DONE"
