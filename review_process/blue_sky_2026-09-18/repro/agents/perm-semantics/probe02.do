* probe02 -- permtstat / permtci / level(90): are the t-scale intervals and the
* inverted-|t| CI the objects the help defines?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
include "dgp.do"
tempfile base
save `base'
local OPTS cluster(unit) att esprelength(3) espostlength(3) jointtests permute(shuffle) permreps(99) rseed(7)

* ---------- A: permtstat + permtci, default equal-tailed ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permtstat permtci permsave(t_def) replace
dis "RC_TDEF=" _rc
local htag TDEF
local habsp 0
local hexact 0
include "harvest.do"
scalar TD_tlb = e(perm_tlb_att)
scalar TD_tub = e(perm_tub_att)
scalar TD_b   = _b[main_att]
scalar TD_se  = _se[main_att]
matrix TD_PD = e(permdist_b)
matrix TD_PS = e(permdist_se)
matrix TD_CM = e(perm_comparator_b)
matrix TD_CS = e(perm_comparator_se)
matrix TD_CI = e(perm_ci)
matrix TD_JW = e(perm_jointW)
scalar TD_jp = e(joint_horiz)
matrix TD_JWO = e(perm_jointWobs)

* ---------- B: permtstat + permabsp + permtci ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permtstat permtci permabsp permsave(t_abs) replace
dis "RC_TABS=" _rc
local htag TABS
local habsp 1
local hexact 0
include "harvest.do"
scalar TA_tlb = e(perm_tlb_att)
scalar TA_tub = e(perm_tub_att)
matrix TA_CI = e(perm_ci)
matrix TA_PD = e(permdist_b)
matrix TA_PS = e(permdist_se)

* permtci must coincide with the table interval for the pooled estimand under permabsp
mata:
  ci = st_matrix("TA_CI")
  b  = st_matrix("_OB")
  nm = st_matrixcolstripe("_OB")
  kk = 0
  for (j=1;j<=rows(nm);j++) if (nm[j,2]=="main_att") kk = j
  printf("TCI-vs-table col=%g  table [%21x , %21x]\n", kk, ci[1,kk], ci[2,kk])
  printf("TCI-vs-table        tci   [%21x , %21x]\n", st_numscalar("TA_tlb"), st_numscalar("TA_tub"))
  printf("TCI-vs-table reldif lb=%g ub=%g\n", reldif(ci[1,kk], st_numscalar("TA_tlb")), reldif(ci[2,kk], st_numscalar("TA_tub")))
  // independent recomputation of permtci: b +/- q_{lvl}(|t_r|) * se_obs
  D = st_matrix("TA_PD"); S = st_matrix("TA_PS")
  tt = J(0,1,.)
  for (r=1;r<=rows(D);r++) if (D[r,kk]<. & S[r,kk]<. & S[r,kk]>0) tt = tt \ abs(D[r,kk]/S[r,kk])
  q = rq(sort(tt,1), 0.95)
  printf("TCI recomputed from dump: q95|t|=%21x  n=%g\n", q, rows(tt))
  st_numscalar("TA_q", q)
end
dis "TCI own recompute: lb=" %21x (TD_b - TA_q*TD_se) " ub=" %21x (TD_b + TA_q*TD_se)
dis "TCI reported     : lb=" %21x TA_tlb " ub=" %21x TA_tub

* ---------- C: level(90), permabsp, no tstat ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permabsp level(90) permsave(l90_abs) replace
dis "RC_L90=" _rc
local htag L90ABS
local habsp 1
local hexact 0
include "harvest.do"

* ---------- D: level(90), default ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' level(90) permsave(l90_def) replace
dis "RC_L90D=" _rc
local htag L90DEF
local habsp 0
local hexact 0
include "harvest.do"

* ---------- E: same seed, are the ASSIGNMENTS identical with and without permtstat? ----
mata:
  st_local("cwd", pwd())
end
preserve
use t_def.dta, clear
qui ds onset_*
local ov `r(varlist)'
mata: TASG = st_data(., "`ov'")
use d_def.dta, clear
qui ds onset_*
local ov2 `r(varlist)'
mata: DASG = st_data(., "`ov2'")
mata: printf("ASSIGN identical across permtstat: %g (dims %g x %g vs %g x %g)\n", (TASG==DASG), rows(TASG), cols(TASG), rows(DASG), cols(DASG))
restore

* ---------- F: recompute the JOINT pre-trend W and its p from the dump ----------
mata:
  D  = st_matrix("TD_PD")
  cm = st_matrix("TD_CM")
  nm = st_matrixcolstripe("TD_PD")
  for (j=1;j<=rows(nm);j++) printf("col %g = %s\n", j, nm[j,2])
end
dis "PROBE02_DONE"
dis "PROBE02_DONE"
