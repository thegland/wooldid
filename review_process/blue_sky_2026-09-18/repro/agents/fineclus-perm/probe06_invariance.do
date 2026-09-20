* probe06: is the permutation EXACTLY invariant, at fixed rseed, to row order,
* to relabelling i, and to relabelling the cluster ids?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* ---------- A. permexact: sorted reference distribution must be identical ----------
foreach v in ORIG ROWSHUF IRELAB CRELAB IREV {
  use "dataB.dta", clear
  if "`v'" == "ROWSHUF" {
    set seed 777
    gen double u = runiform()
    sort u
    drop u
  }
  if "`v'" == "IRELAB" replace coh = 100 - coh
  if "`v'" == "IREV"   replace coh = cond(coh==1,3,cond(coh==3,1,2))
  if "`v'" == "CRELAB" replace state = 60 - state*7
  cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts permsave(I_`v') replace
  dis "EXACT `v' rc=" _rc " M=" e(perm_npossible) " p=" e(perm_p_att) " sd=" %21x e(perm_sd_att) " b=" %21x _b[main_att]
}
foreach v in ROWSHUF IRELAB CRELAB IREV {
  use "I_ORIG.dta", clear
  keep _rep main_att
  drop if _rep == 0
  sort main_att
  gen long k = _n
  keep k main_att
  rename main_att a
  tempfile o
  save `o'
  use "I_`v'.dta", clear
  keep _rep main_att
  drop if _rep == 0
  sort main_att
  gen long k = _n
  keep k main_att
  merge 1:1 k using `o', nogen
  gen double d = abs(main_att - a)
  qui sum d
  dis "SORTEDDIFF `v': maxabs=" %12.4e r(max)
}

* ---------- B. sampled permute at a fixed rseed ----------
foreach v in ORIG ROWSHUF IRELAB CRELAB {
  use "dataB.dta", clear
  if "`v'" == "ROWSHUF" {
    set seed 777
    gen double u = runiform()
    sort u
    drop u
  }
  if "`v'" == "IRELAB" replace coh = 100 - coh
  if "`v'" == "CRELAB" replace state = 60 - state*7
  cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(40) rseed(12345) permonsetcohorts permsave(S_`v') replace
  dis "SAMP `v' rc=" _rc " p=" e(perm_p_att) " sd=" %21x e(perm_sd_att)
}
foreach v in ROWSHUF IRELAB CRELAB {
  use "S_ORIG.dta", clear
  keep _rep main_att
  rename main_att a
  tempfile o
  save `o'
  use "S_`v'.dta", clear
  keep _rep main_att
  merge 1:1 _rep using `o', nogen
  gen double d = abs(main_att - a)
  qui sum d
  dis "SEQDIFF `v': maxabs (rep-by-rep, same rseed)=" %12.4e r(max)
}
