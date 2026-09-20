clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

*==============================================================================
* base panel: 7 clusters, onsets {5, 5, 8, ., ., ., .}
*   M = 7!/(2!1!4!) = 5040/48 = 105 enumerated assignments
*==============================================================================
clear
set obs 7
gen long i = _n
gen ttre = .
replace ttre = 5 in 1
replace ttre = 5 in 2
replace ttre = 8 in 3
expand 12
bysort i: gen int t = _n
set seed 4242
gen double ui = .
bysort i (t): replace ui = rnormal() if _n == 1
bysort i (t): replace ui = ui[1]
gen double y = 2 + ui + 0.30*t + 1.20*(t >= ttre & !missing(ttre)) + rnormal()*0.5
gen int ycount = max(0, round(exp(0.20 + 0.08*t + 0.30*(t >= ttre & !missing(ttre)) + rnormal()*0.2)))
gen double w = 1 + 0.5*mod(i, 3)
gen double dose = 0.5 + 0.25*mod(i + t, 4)
tempfile base
save `base'

*==============================================================================
* PART 1 -- INTERLEAVED determinism.  This mimics the actual permutation loop:
* the rep-0 comparator fit, then many fits on REASSIGNED onsets, then the
* identity fit again.  If the second identity fit differs from the first in any
* bit, memb can be 0 under a complete enumeration.  This is the exact condition
* PERMEX-PX-2's reachability turns on.
*==============================================================================
prog def ilv
  syntax , tag(string) yv(string) opts(string asis)
  capture noisily {
    tempname b1 b2
    qui wooldid `yv' i t ttre , `opts'
    matrix `b1' = e(b)
    * 24 intervening fits on shuffled onsets, as the rep loop does
    forval r = 1/24 {
      preserve
      qui gen double __sh = .
      qui bysort i (t): replace __sh = runiform() if _n == 1
      qui bysort i (t): replace __sh = __sh[1]
      sort __sh i t
      qui by __sh: gen long __g = 1
      tempvar newons
      qui gegen double `newons' = mean(ttre), by(i)
      qui replace `newons' = .
      * assign the multiset {5,5,8,.,.,.,.} to units in a shuffled order
      qui egen long __u = group(__sh)
      qui replace `newons' = 5 if __u == 1
      qui replace `newons' = 5 if __u == 2
      qui replace `newons' = 8 if __u == 3
      cap qui wooldid `yv' i t `newons' , `opts'
      restore
    }
    qui wooldid `yv' i t ttre , `opts'
    matrix `b2' = e(b)
    local nd 0
    forval c = 1/`=colsof(`b1')' {
      if `b1'[1,`c'] != `b2'[1,`c'] local ++nd
    }
    dis "ILV_`tag'_K=" colsof(`b1') " ILV_`tag'_NDIFF=`nd'"
  }
  dis "ILV_`tag'_RC=" _rc
end

use `base', clear
ilv , tag(clusi)  yv(y)      opts(att cluster(i))
use `base', clear
ilv , tag(allest) yv(y)      opts(att att_i att_it att_itime cluster(i) esprelength(3) espostlength(3))
use `base', clear
ilv , tag(pois)   yv(ycount) opts(att cluster(i) poisson)
use `base', clear
ilv , tag(cont)   yv(y)      opts(att cluster(i) contreat(dose) contreatcontrols(basic))
use `base', clear
ilv , tag(fastn)  yv(y)      opts(att cluster(i) fastnose)

*==============================================================================
* PART 2 -- runtime detector.  element in e(perm_addone) but NOT in
* e(perm_incomplete)  <=>  exk=1 & memb=0  <=>  the state PERMEX-PX-2 claims.
*==============================================================================
prog def pxrep
  syntax , tag(string)
  local ao `"`e(perm_addone)'"'
  local in `"`e(perm_incomplete)'"'
  dis "PX_`tag'_M=" e(perm_reps) " DONE=" e(perm_completed) " FAIL=" e(perm_failed) " BLANK=" e(perm_blank)
  dis "PX_`tag'_ADDONE=[" `"`ao'"' "]"
  dis "PX_`tag'_INCOMP=[" `"`in'"' "]"
  local bad
  foreach el of local ao {
    local hit 0
    foreach e2 of local in {
      if "`el'" == "`e2'" local hit 1
    }
    if `hit' == 0 local bad `bad' `el'
  }
  dis "PX_`tag'_FALSECAUTION=[`bad']"
  tempname R
  capture matrix `R' = e(wdidmainresults)
  if _rc == 0 {
    local nr = rowsof(`R')
    local zz 0
    forval r = 1/`nr' {
      if `R'[`r',4] == 0 local ++zz
    }
    dis "PX_`tag'_NROW=`nr' PZERO=`zz'"
  }
end

use `base', clear
capture noisily {
  wooldid y i t ttre , att cluster(i) permute(shuffle) permexact
  pxrep , tag(A_base)
}
dis "PX_A_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att att_i att_it att_itime cluster(i) esprelength(3) espostlength(3) jointtests permute(shuffle) permexact
  pxrep , tag(B_wide)
}
dis "PX_B_RC=" _rc

use `base', clear
capture noisily {
  wooldid ycount i t ttre , att cluster(i) poisson permute(shuffle) permexact
  pxrep , tag(C_pois)
}
dis "PX_C_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre [aw=w] , att cluster(i) permute(shuffle) permexact
  pxrep , tag(D_awgt)
}
dis "PX_D_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att cluster(i) contreat(dose) contreatcontrols(basic) permute(shuffle) permexact
  pxrep , tag(E_cont)
}
dis "PX_E_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att cluster(i) permute(shuffle) permexact permtstat esprelength(3) espostlength(3)
  pxrep , tag(F_tstat)
}
dis "PX_F_RC=" _rc

*-- unbalanced: two clusters lose their late periods, so some reassigned onsets
*   leave some horizons unestimable -> exk=0 for those elements
use `base', clear
drop if inlist(i, 6, 7) & t >= 9
capture noisily {
  wooldid y i t ttre , att cluster(i) esprelength(3) espostlength(4) permute(shuffle) permexact
  pxrep , tag(G_unbal)
}
dis "PX_G_RC=" _rc

*-- strata
use `base', clear
gen byte stratum = 1 + (i > 3)
capture noisily {
  wooldid y i t ttre , att cluster(i) permute(shuffle) permexact permstrata(stratum)
  pxrep , tag(H_strat)
}
dis "PX_H_RC=" _rc

*-- a cluster with rows the design sample cannot reach: missing y on a block
use `base', clear
replace y = . if i == 4 & t <= 3
capture noisily {
  wooldid y i t ttre , att cluster(i) esprelength(3) espostlength(3) permute(shuffle) permexact
  pxrep , tag(I_holes)
}
dis "PX_I_RC=" _rc

dis "ALLDONE"
