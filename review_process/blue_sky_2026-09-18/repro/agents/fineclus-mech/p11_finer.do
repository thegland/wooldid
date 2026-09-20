clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P11 -- the finer-than-i geometries, on a design where the ONSET IS
*        CONSTANT WITHIN i (p02's i=cgrp tests were invalid: onset varied
*        inside the group and wooldid refused them with r(459)).
*
*   96 firms in 12 cohort groups of 8; onset assigned BY GROUP.
*   i = cgrp (12)                     <- the estimand definition
*   cluster(firm)   96  STRICTLY FINER than i  -> keep-i regime
*   cluster(pair)   48  FINER                  -> keep-i regime
*   cluster(cgrp)   12  == i                   -> nesting regime
*   cluster(rowid) 768  one per observation    -> keep-i regime
*   no cluster()                               -> keep-i, plain bsample
*----------------------------------------------------------------------
prog def mkp
  syntax , [singleton(int 0) deadclus(int 0)]
  clear
  set seed 1234
  set sortseed 1234
  set obs 96
  gen long firm = _n
  gen long cgrp = ceil(_n/8)
  gen long pair = ceil(_n/2)
  gen int onset = .
  replace onset = 2004 if cgrp <= 3
  replace onset = 2006 if cgrp >= 4 & cgrp <= 6
  replace onset = 2008 if cgrp >= 7 & cgrp <= 9
  gen double gshk = .
  gen double gtau = .
  bysort cgrp (firm): replace gshk = rnormal()          if _n==1
  bysort cgrp (firm): replace gtau = 1 + rnormal()*0.6  if _n==1
  bysort cgrp (firm): replace gshk = gshk[1]
  bysort cgrp (firm): replace gtau = gtau[1]
  gen double fshk = rnormal()*0.8
  gen double ftau = rnormal()*0.5
  expand 8
  bysort firm (cgrp): gen int year = 2001 + _n - 1
  gen double e = rnormal()*0.4
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = gshk + fshk + 0.05*(year-2001) + (gtau + ftau)*D + e
  gen long rowid = _n
  gen byte one = 1
  if `singleton' == 1 {
    bysort firm (year): gen byte __k = _n
    drop if inlist(firm,7,19,31) & __k > 1
    drop __k
  }
  if `deadclus' == 1 replace y = . if firm <= 8
  sort firm year
end

prog def shw
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if inlist("`this'","main_att","main_att_i") {
      dis "P11 `tag' `this' b=" %18.14f `b'[1,`q'] " se=" %18.14f sqrt(`V'[`q',`q'])
    }
  }
end

*--- sanity: the design estimates at all -------------------------------
mkp
cap noisily wooldid y cgrp year onset, att att_i cluster(firm)
dis "P11 S0_ANALYTIC_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
shw , tag(S0_ANALYTIC)

*--- A. finer cluster, cell-weighted estimand, keep-i regime -----------
foreach CV in firm pair cgrp rowid {
  mkp
  cap noisily wooldid y cgrp year onset, att att_i cluster(`CV') bootstrap(standard) reps(100) rseed(77) dots(0)
  dis "P11 A_`CV'_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust) " DFR=" e(df_r)
  shw , tag(A_`CV')
}

*--- B. no cluster() at all (T2 2026-09-16 lifted the refusal) ---------
mkp
cap noisily wooldid y cgrp year onset, att att_i bootstrap(standard) reps(100) rseed(77) dots(0)
dis "P11 B_NOCLUS_RC=" _rc " N=" e(N)
shw , tag(B_NOCLUS)

*--- C. boottstat on the finer-cluster path (the P7 assert) ------------
foreach CV in firm pair cgrp {
  mkp
  cap noisily wooldid y cgrp year onset, att att_i cluster(`CV') bootstrap(standard) boottstat reps(50) rseed(77) dots(0)
  dis "P11 C_BTS_`CV'_RC=" _rc " NCLUST=" e(N_clust)
  shw , tag(C_BTS_`CV')
}

*--- D. the P7 assert's false-positive geometry: recommended pairing
*       (coarse i + fe(unit)) + finer cluster + boottstat + singletons
mkp , singleton(1)
cap noisily wooldid y cgrp year onset, att att_i fe(firm) cluster(firm) bootstrap(standard) boottstat reps(40) rseed(77) dots(0)
dis "P11 D_SINGLETON_BTS_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
shw , tag(D_SINGLETON_BTS)

mkp , singleton(1)
cap noisily wooldid y cgrp year onset, att att_i fe(firm) cluster(firm) bootstrap(standard) reps(40) rseed(77) dots(0)
dis "P11 D2_SINGLETON_NOBTS_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
shw , tag(D2_SINGLETON_NOBTS)

*--- E. the negative control the design names (WOOLDID_BOOTNOSUB) ------
global WOOLDID_BOOTNOSUB 1
mkp
cap noisily wooldid y cgrp year onset, att att_i cluster(firm) bootstrap(standard) boottstat reps(50) rseed(77) dots(0)
dis "P11 E_NOSUB_RC=" _rc
shw , tag(E_NOSUB)
global WOOLDID_BOOTNOSUB

*--- F. an ENTIRE cluster out of e(sample): e(df_r) vs e(N_clust)-1 ----
* wdb_C is counted on smp, which does not markout missing y, so the dead
* cluster is still a resampling unit and still sets bootdf(clusters).
mkp , deadclus(1)
cap noisily wooldid y cgrp year onset, att att_i cluster(firm)
dis "P11 F_PLAIN_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust)
mkp , deadclus(1)
cap noisily wooldid y cgrp year onset, att att_i cluster(firm) bootstrap(standard) reps(100) rseed(77) dots(0) bootdf(clusters)
dis "P11 F_BOOT_RC=" _rc " N=" e(N) " NCLUST=" e(N_clust) " DFR=" e(df_r)
dis "P11 F  e(N_clust)-1 = " e(N_clust)-1 "   vs   e(df_r) = " e(df_r)
shw , tag(F_DEADCLUS)

dis "ALLDONE"
