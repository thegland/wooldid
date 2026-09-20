* r05: how far does the near-zero gteffects SE reach? cluster count and cluster level.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop sed2
prog def sed2
  syntax , tag(string)
  tempname G
  cap matrix `G' = e(gteffects)
  if _rc != 0 {
    dis "SED|`tag'|NOMAT"
    exit
  }
  local cn : colnames `G'
  local js : list posof "se" in cn
  local jp : list posof "p" in cn
  mata: G = st_matrix("`G'"); se = G[.,`js']; pp = G[.,`jp']; ///
        st_local("a", strofreal(rows(se))); ///
        st_local("b", strofreal(sum(se :>= .))); ///
        st_local("c", strofreal(sum(se :== 0))); ///
        st_local("d", strofreal(sum((se :< 1e-8) :& (se :> 0)))); ///
        st_local("e2", strofreal(min(select(pp, pp :< .)), "%10.3e"))
  dis "SED|`tag'|rows=`a'|se_missing=`b'|se_exact0=`c'|se_tiny=`d'|min_p=`e2'"
end

* ---- K onset cohorts + never-treated, 30 units per cohort, 24 periods ----
foreach K in 3 6 12 20 {
  clear
  local NN = 30*(`K'+1)
  set obs `NN'
  gen long unit = _n
  gen byte coh = 1 + mod(_n-1, `K'+1)
  gen double ttre = .
  forval c = 1/`K' {
    replace ttre = 2003 + `c' if coh == `c'
  }
  gen double ufe = rnormal()
  expand 24
  bysort unit: gen int t = 1999 + _n
  gen double tfe = 0.05*(t-1999)
  gen byte trt = (ttre<.) & (t >= ttre)
  gen double rel = t - ttre
  gen double y = 5 + ufe + tfe + cond(trt==1, 1.5 + 0.10*rel, 0) + rnormal()
  tempfile d`K'
  save `d`K'', replace

  use `d`K'', clear
  cap noi qui wooldid y coh t ttre, att esprelength(3) espostlength(3) gteffects cluster(coh)
  dis "K`K'_ati_RC=" _rc
  sed2, tag(K`K'_clusterAT_i)

  use `d`K'', clear
  cap noi qui wooldid y coh t ttre, att esprelength(3) espostlength(3) gteffects cluster(unit)
  dis "K`K'_below_RC=" _rc
  sed2, tag(K`K'_clusterBELOW_i)

  use `d`K'', clear
  gen long reg3 = 1 + mod(unit-1, 3)
  cap noi qui wooldid y coh t ttre, att esprelength(3) espostlength(3) gteffects cluster(reg3)
  dis "K`K'_cross3_RC=" _rc
  sed2, tag(K`K'_cluster3crosscut)
}

dis "===== R05 DONE ====="
