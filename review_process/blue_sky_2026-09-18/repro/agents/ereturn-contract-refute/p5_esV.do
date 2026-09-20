* SKEPTIC probe of EC-4: is e(V_es_att) really singular, and is the singularity
* an honest property of the cluster-robust estimator or a computation fault?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606060
set sortseed 606060
which wooldid

cap prog drop esdiag
prog def esdiag
  args tag
  cap confirm matrix e(V_es_att)
  if _rc != 0 {
    dis "`tag': V_es_att ABSENT (rc=" _rc ")"
    exit
  }
  mat VE = e(V_es_att)
  mat BE = e(b_es_att)
  local nd = .
  cap local nd = e(es_ndegen_att)
  mata: V = st_matrix("VE")
  mata: st_local("rk", strofreal(rank(V)))
  mata: st_local("hh", strofreal(cols(V)))
  mata: st_local("mind", strofreal(min(diagonal(V))))
  mata: symeigensystem(V, X=., L=.)
  mata: st_local("ev", invtokens(strofreal(L, "%12.4e")))
  dis "`tag': H=`hh' rank=`rk' mindiag=`mind' ndegen=`nd'"
  dis "`tag':   eigenvalues: `ev'"
  dis "`tag':   colnames: `: colnames VE'"
  * PSD check: any eigenvalue materially negative relative to the largest?
  mata: st_local("negrel", strofreal(min(L)/max(L)))
  dis "`tag':   min/max eigenvalue ratio = `negrel'"
  * cholesky behaviour
  cap mat CH = cholesky(VE)
  dis "`tag':   cholesky rc=" _rc
  * rank of the whole e(V)
  mat VV = e(V)
  mata: st_local("rkv", strofreal(rank(st_matrix("VV"))))
  mata: st_local("kv", strofreal(cols(st_matrix("VV"))))
  dis "`tag':   rank(e(V))=`rkv' of `kv'; Nclust=" e(N_clust)
end

*========= my own panel: 3 treated cohorts + never-treated, unit-level i =========
clear
set obs 144
gen long id = _n
gen double ui = rnormal()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2007 if mod(id,4)==1
replace ttre = 2010 if mod(id,4)==2
expand 12
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 1.2 + 0.35*(t-ttre) if t>=ttre & ttre<.
gen byte c2  = ceil(id/2)
gen byte c6  = ceil(id/6)
gen byte c12 = ceil(id/12)
gen byte c36 = ceil(id/36)
gen long cid = id
save p5base, replace

local W "esfixedbaseperiod esprelength(3) espostlength(5)"

* --- does cluster granularity move the rank? ---
foreach C in cid c2 c6 c12 c36 {
  cap noisily wooldid y id t ttre, att `W' cluster(`C')
  dis "P5 CLU_`C' rc=" _rc
  if _rc==0 esdiag "CLU_`C'"
}
cap noisily wooldid y id t ttre, att `W' robust
dis "P5 ROBUST rc=" _rc
if _rc==0 esdiag "ROBUST"
cap noisily wooldid y id t ttre, att `W'
dis "P5 IID rc=" _rc
if _rc==0 esdiag "IID"

* --- does window width move it? ---
foreach P in 1 2 3 4 {
  cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(`P') espostlength(`P') cluster(c6)
  dis "P5 WIN_`P' rc=" _rc
  if _rc==0 esdiag "WIN_`P'"
}
foreach Q in 2 3 4 5 {
  cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(`Q') cluster(c6)
  dis "P5 POST_`Q' rc=" _rc
  if _rc==0 esdiag "POST_`Q'"
}

* --- is it the reghdfe path?  same VCE through the plain-reg route ---
cap noisily wooldid y id t ttre, att `W' customvce(cluster c6)
dis "P5 CUSTOMVCE rc=" _rc
if _rc==0 esdiag "CUSTOMVCE"
cap noisily wooldid y id t ttre, att `W' unconditionalse cluster(c6)
dis "P5 UNCOND rc=" _rc
if _rc==0 esdiag "UNCOND"

* --- the design the help endorses: i = cohort, unit FE in fe() ---
use p5base, clear
egen long coh = group(ttre), missing
cap noisily wooldid y coh t ttre, att `W' fe(id) cluster(c6)
dis "P5 COHORT_I rc=" _rc
if _rc==0 esdiag "COHORT_I"

* --- what is in the null space?  print the eigenvectors of the ~0 eigenvalues ---
use p5base, clear
cap noisily wooldid y id t ttre, att `W' cluster(c6)
if _rc==0 {
  mat VE = e(V_es_att)
  mata:
    V  = st_matrix("VE")
    nm = st_matrixcolstripe("VE")
    symeigensystem(V, X=., L=.)
    mx = max(L)
    for (j=1; j<=cols(L); j++) {
      if (L[1,j]/mx < 1e-10) {
        printf("NULLVEC eigen=%12.4e  loadings:\n", L[1,j])
        for (a=1; a<=rows(X); a++) printf("   %s %12.6f\n", nm[a,2], X[a,j])
      }
    }
  end
}

* --- does a downstream consumer notice?  test / lincom / pretrends handoff ---
if _rc==0 {
  cap noisily test ES_att_pre3 ES_att_pre2
  dis "P5 test-pre rc=" _rc " F=" r(F) " df=" r(df) " p=" r(p)
  cap noisily lincom ES_att_post5 - ES_att_post4
  dis "P5 lincom rc=" _rc " se=" r(se)
  dis "P5 es_ndegen_att=" e(es_ndegen_att)
  dis "P5 wdid_vcvstate=" e(wdid_vcvstate)
}
dis "ALLDONE_P5"
