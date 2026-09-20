* probe03 -- permexact: are the shares plain, is the floor 2/M, does p ever hit 0,
* is the Monte-Carlo precision withheld, and does permanalyticse break membership?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
include "recheck.do"
* 8 clusters, 2 treated at one onset -> M = C(8,2) = 28 assignments
local dN 8
local dT 8
local dons 2004
local dnever 6
local dtau 2.5
local dsd 0.4
local dseed 991
include "dgp.do"
tab ttre, m
tempfile base
save `base'
local OPTS cluster(unit) att esprelength(2) espostlength(2) jointtests permute(shuffle) permexact

* ---- A: plain permexact ----
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permsave(x_def) replace
dis "RC_XDEF=" _rc
dis "npossible = " e(perm_npossible) "  completed=" e(perm_completed) " failed=" e(perm_failed)
dis "perm_addone = [" e(perm_addone) "]"
dis "perm_incomplete = [" e(perm_incomplete) "]"
dis "perm_exact = " e(perm_exact)
local htag XDEF
local habsp 0
local hexact 1
include "harvest.do"
matrix list e(perm_p), format(%12.8f)
matrix list e(perm_p_se)
matrix list e(perm_n)

* ---- B: permexact + permabsp ----
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permabsp permsave(x_abs) replace
dis "RC_XABS=" _rc
dis "perm_addone = [" e(perm_addone) "]"
local htag XABS
local habsp 1
local hexact 1
include "harvest.do"
matrix list e(perm_p), format(%12.8f)

* ---- C: permexact + permanalyticse (comparator is a fast-mode refit) ----
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permanalyticse permsave(x_ase) replace
dis "RC_XASE=" _rc
dis "perm_addone = [" e(perm_addone) "]"
dis "perm_incomplete = [" e(perm_incomplete) "]"
local htag XASE
local habsp 0
local hexact 1
include "harvest.do"
matrix list e(perm_p), format(%12.8f)
* is the comparator bit-identical to e(b)?
mata:
  b = st_matrix("_OB"); c = st_matrix("_CM")
  printf("XASE comparator==e(b): %g  maxreldif=%g\n", (b==c), max(reldif(b,c)))
  nm = st_matrixcolstripe("_OB")
  for (j=1;j<=cols(b);j++) printf("  col %g %s  b=%21x  cm=%21x\n", j, nm[j,2], b[1,j], c[1,j])
end

* ---- D: permexact + permtstat ----
use `base', clear
cap noi wooldid y unit t ttre, `OPTS' permtstat permsave(x_ts) replace
dis "RC_XTS=" _rc
dis "perm_addone = [" e(perm_addone) "]"
local htag XTS
local habsp 0
local hexact 1
include "harvest.do"
matrix list e(perm_p), format(%12.8f)
dis "PROBE03_DONE"
