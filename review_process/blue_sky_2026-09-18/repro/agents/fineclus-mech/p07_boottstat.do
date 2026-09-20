clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

*----------------------------------------------------------------------
* P07 -- boottstat on the physical re-id leg under BOTH regimes, plus the
*        P7 run-time assert (ado:9386-9414) and its documented negative
*        control (the internal global WOOLDID_BOOTNOSUB).
*  i = cgrp (12)        cluster(ind)  = 24  -> FINER than i  -> keep-i
*  i = cgrp (12)        cluster(firm) = 96  -> FINER than i  -> keep-i
*  i = firm (96)        cluster(ind)  = 24  -> COARSER       -> nesting
*----------------------------------------------------------------------
prog def dumpse
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local names : colnames `b'
  local K = colsof(`b')
  forvalues q = 1/`K' {
    local this : word `q' of `names'
    if inlist("`this'","main_att","main_att_i") {
      dis "P07 `tag' `this' b=" %18.14f `b'[1,`q'] " se=" %18.14f sqrt(`V'[`q',`q'])
    }
  }
end

use "p02_base.dta", clear
count
dis "P07 SETUP: firms="
qui levelsof firm, local(L)
dis "P07 nfirm=" `: word count `L''
qui levelsof ind, local(L2)
dis "P07 nind=" `: word count `L2''
qui levelsof cgrp, local(L3)
dis "P07 ncgrp=" `: word count `L3''

*--- A: keep-i regime, boottstat, cluster(ind) finer than i=cgrp ---------
use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att att_i cluster(ind) bootstrap(standard) boottstat reps(60) rseed(55) dots(0)
dis "P07 A_KEEPI_BTS_RC=" _rc " NCLUST=" e(N_clust)
dumpse , tag(A_KEEPI)

*--- B: keep-i regime, boottstat, cluster(firm) much finer ---------------
use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att att_i cluster(firm) bootstrap(standard) boottstat reps(60) rseed(55) dots(0)
dis "P07 B_KEEPIFINE_BTS_RC=" _rc " NCLUST=" e(N_clust)
dumpse , tag(B_KEEPIFINE)

*--- C: nesting regime, boottstat, i=firm cluster(ind) -------------------
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att att_i cluster(ind) bootstrap(standard) boottstat reps(60) rseed(55) dots(0)
dis "P07 C_NEST_BTS_RC=" _rc " NCLUST=" e(N_clust)
dumpse , tag(C_NEST)

*--- D: the documented negative control: suppress the subinstr -----------
global WOOLDID_BOOTNOSUB 1
use "p02_base.dta", clear
cap noisily wooldid y cgrp year onset, att att_i cluster(ind) bootstrap(standard) boottstat reps(60) rseed(55) dots(0)
dis "P07 D_NOSUB_KEEPI_RC=" _rc
dumpse , tag(D_NOSUB_KEEPI)
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att att_i cluster(ind) bootstrap(standard) boottstat reps(60) rseed(55) dots(0)
dis "P07 E_NOSUB_NEST_RC=" _rc
dumpse , tag(E_NOSUB_NEST)
global WOOLDID_BOOTNOSUB

*--- F: boottstat + cluster == i -----------------------------------------
use "p02_base.dta", clear
cap noisily wooldid y firm year onset, att att_i cluster(firm) bootstrap(standard) boottstat reps(60) rseed(55) dots(0)
dis "P07 F_CLUSEQI_BTS_RC=" _rc " NCLUST=" e(N_clust)
dumpse , tag(F_CLUSEQI)

*--- G: the P7 assert's false-positive geometry -------------------------
* the RECOMMENDED pairing (sthlp:255: coarse i + fe(unit)) plus the newly
* enabled finer-than-i cluster under boottstat, on an unbalanced panel
* where three firms have exactly one row (reghdfe singletons).
use "p02_base.dta", clear
bysort firm (year): gen byte __k = _n
drop if inlist(firm,7,19,31) & __k > 1
drop __k
save "p07_singleton.dta", replace
cap noisily wooldid y cgrp year onset, att att_i fe(firm) cluster(firm) bootstrap(standard) boottstat reps(40) rseed(55) dots(0)
dis "P07 G_SINGLETON_RC=" _rc " NCLUST=" e(N_clust)
dumpse , tag(G_SINGLETON)

*--- H: same geometry without boottstat (no subinstr, no assert) ---------
use "p07_singleton.dta", clear
cap noisily wooldid y cgrp year onset, att att_i fe(firm) cluster(firm) bootstrap(standard) reps(40) rseed(55) dots(0)
dis "P07 H_SINGLETON_NOBTS_RC=" _rc " NCLUST=" e(N_clust)
dumpse , tag(H_SINGLETON_NOBTS)

dis "P07REALLYDONE"
