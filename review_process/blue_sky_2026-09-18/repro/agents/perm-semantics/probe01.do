* probe01 -- SEMANTIC CONFORMANCE of every permute() p-value and interval, five
* conventions on one design and one seed.  Everything is recomputed from the saved
* null distribution using the HELP FILE's definitions.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
which wooldid
include "dgp.do"
tempfile base
save `base'

local OPTS cluster(unit) att esprelength(3) espostlength(3) jointtests permute(shuffle) permreps(99) rseed(7)

* ---------- run 1: default equal-tailed ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permsave(d_def) replace
dis "RC_DEF=" _rc
local htag DEF
local habsp 0
local hexact 0
include "harvest.do"
matrix DEF_PLO = e(perm_p_lo)
matrix DEF_PHI = e(perm_p_hi)
matrix DEF_PD  = e(permdist_b)
matrix DEF_P   = e(perm_p)
matrix DEF_CI  = e(perm_ci)
* permsave file must reproduce e(permdist_b)
preserve
use d_def.dta, clear
dis "DEF nobs in dump = " _N
mata: X = st_data(., "*"); printf("DEF dump cols=%g rows=%g\n", cols(X), rows(X))
restore

* ---------- run 2: permabsp ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permabsp permsave(d_abs) replace
dis "RC_ABS=" _rc
local htag ABS
local habsp 1
local hexact 0
include "harvest.do"
matrix ABS_PLO = e(perm_p_lo)
matrix ABS_PHI = e(perm_p_hi)
matrix ABS_PD  = e(permdist_b)

* ---------- run 3: permrawp ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permrawp permsave(d_raw) replace
dis "RC_RAW=" _rc
local htag RAW
local habsp 2
local hexact 0
include "harvest.do"
matrix RAW_PLO = e(perm_p_lo)
matrix RAW_PHI = e(perm_p_hi)
matrix RAW_PD  = e(permdist_b)

* ---------- run 4: permside(left) ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permside(left) permsave(d_left) replace
dis "RC_LEFT=" _rc
local htag LEFT
local habsp 3
local hexact 0
include "harvest.do"
matrix LEFT_P = e(perm_p)
matrix LEFT_PLO = e(perm_p_lo)
matrix LEFT_PHI = e(perm_p_hi)

* ---------- run 5: permside(right) ----------
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permside(right) permsave(d_right) replace
dis "RC_RIGHT=" _rc
local htag RIGHT
local habsp 4
local hexact 0
include "harvest.do"
matrix RIGHT_P = e(perm_p)
matrix RIGHT_PLO = e(perm_p_lo)
matrix RIGHT_PHI = e(perm_p_hi)

* ---------- cross-run invariants ----------
mata:
  a = st_matrix("DEF_PLO"); b = st_matrix("ABS_PLO"); c = st_matrix("RAW_PLO")
  d = st_matrix("LEFT_PLO"); e = st_matrix("RIGHT_PLO")
  printf("CROSS p_lo all equal: %g %g %g %g\n", a==b, a==c, a==d, a==e)
  a2= st_matrix("DEF_PHI"); b2= st_matrix("ABS_PHI"); c2= st_matrix("RAW_PHI")
  d2= st_matrix("LEFT_PHI"); e2= st_matrix("RIGHT_PHI")
  printf("CROSS p_hi all equal: %g %g %g %g\n", a2==b2, a2==c2, a2==d2, a2==e2)
  printf("CROSS draws identical across conventions: %g %g\n",
         st_matrix("DEF_PD")==st_matrix("ABS_PD"), st_matrix("DEF_PD")==st_matrix("RAW_PD"))
  printf("CROSS LEFT p == p_lo : %g   RIGHT p == p_hi : %g\n",
         st_matrix("LEFT_P")==d, st_matrix("RIGHT_P")==e2)
  // default p must be min(1, 2*min(p_lo,p_hi))
  pd = st_matrix("DEF_P"); lo = a; hi = a2
  mn = J(1, cols(pd), .)
  for (k=1;k<=cols(pd);k++) mn[1,k] = min((1, 2*min((lo[1,k], hi[1,k]))))
  printf("CROSS DEF p == min(1,2min(plo,phi)) : %g\n", (pd==mn))
  st_matrix("DEF_MN", mn)
end
matlist DEF_P, format(%12.8f)
matlist DEF_MN, format(%12.8f)
matlist LEFT_P, format(%12.8f)
matlist RIGHT_P, format(%12.8f)
dis "PROBE01_DONE"
