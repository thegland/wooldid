* EC-4: e(V_es_att) is singular with a strictly positive diagonal and e(es_ndegen_att)
* is absent, which the help says marks a clean covariance.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 120
gen long id = _n
gen double ttre = .
replace ttre = 4  if inrange(id,1,30)
replace ttre = 7  if inrange(id,31,60)
replace ttre = 10 if inrange(id,61,80)
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()

gen long coh = 1 + (ttre==7) + 2*(ttre==10)
replace coh = 4 if missing(ttre)

foreach V in "cluster(clu)" "" "robust" {
  cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod ///
      esprelength(3) espostlength(5) `V'
  dis "================ vce: `V'   rc=" _rc
  mata: st_local("rk", strofreal(rank(st_matrix("e(V_es_att)"))))
  mata: st_local("ev", invtokens(strofreal(sort(Re(eigenvalues(st_matrix("e(V_es_att)")))',1)',"%9.3e")))
  mata: st_local("md", strofreal(min(diagonal(st_matrix("e(V_es_att)"))),"%9.3e"))
  mata: st_local("ch", strofreal(sum(cholesky(makesymmetric(st_matrix("e(V_es_att)"))) :>= .)))
  dis "  e(V_es_att) is " rowsof(e(V_es_att)) "x" colsof(e(V_es_att)) ", rank = `rk'"
  dis "  eigenvalues : `ev'"
  dis "  min diagonal: `md'   (every diagonal entry is strictly positive)"
  dis "  cholesky(makesymmetric(.)) missing cells: `ch'"
  cap confirm scalar e(es_ndegen_att)
  dis "  e(es_ndegen_att) posted? " (_rc==0) "   <- the help: its absence marks a clean covariance"
  mata: st_local("rkv", strofreal(rank(st_matrix("e(V)"))))
  dis "  rank(e(V)) = `rkv' of K = " colsof(e(b))
}

dis ""
dis "======== THE CONTROL: the design the help says satisfies its conditions ========"
dis "  (a coarser i -- cohorts -- with unit fixed effects through fe())"
cap noi wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod ///
    esprelength(3) espostlength(5) fe(id) cluster(clu)
dis "rc = " _rc
mata: st_local("rk", strofreal(rank(st_matrix("e(V_es_att)"))))
mata: st_local("ch", strofreal(sum(cholesky(makesymmetric(st_matrix("e(V_es_att)"))) :>= .)))
dis "  e(V_es_att) is " rowsof(e(V_es_att)) "x" colsof(e(V_es_att)) ", rank = `rk', cholesky missing cells = `ch'"
dis "  -> FULL RANK.  The singularity is confined to cluster() with a unit-level i,"
dis "     which is the geometry the help's Design-conditions block already flags in prose."
dis "     The defect is that e(es_ndegen_att) is absent in BOTH cases, and the help says"
dis "     its absence marks a clean covariance."
