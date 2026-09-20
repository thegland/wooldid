clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

*==============================================================================
* PART 0 -- a small staggered panel.  6 clusters, onsets {5, 8, ., ., ., .}
*           => M = 6!/(1!1!4!) = 30 enumerated assignments under permexact.
*==============================================================================
clear
set obs 6
gen long i = _n
gen ttre = .
replace ttre = 5 in 1
replace ttre = 8 in 2
expand 10
bysort i: gen int t = _n
set seed 111
gen double ui = .
bysort i (t): replace ui = rnormal() if _n == 1
bysort i (t): replace ui = ui[1]
gen double y = 2 + ui + 0.30*t + 1.20*(t >= ttre & !missing(ttre)) + rnormal()*0.5
gen int ycount = max(0, round(exp(0.20 + 0.10*t + 0.30*(t >= ttre & !missing(ttre)) + rnormal()*0.2)))
gen double w = 1 + 0.5*mod(i, 3)
gen double dose = 0.5 + 0.25*mod(i + t, 4)
tempfile base
save `base'

*==============================================================================
* PART 1 -- DETERMINISM OF THE FIT.
* memb=0 with exk=1 requires the identity replication's statistic to differ
* bitwise from the rep-0 comparator's.  Both are separate wooldid fits of the
* SAME data.  So: is a wooldid fit a bitwise-deterministic function of its
* input, path by path?  Any path where two identical calls differ in the last
* ulp is a path where memb can be 0 under a complete enumeration.
*==============================================================================
prog def detchk
  syntax , tag(string) opts(string asis)
  tempname b1 b2
  capture noisily {
    qui wooldid y i t ttre , `opts'
    matrix `b1' = e(b)
    qui wooldid y i t ttre , `opts'
    matrix `b2' = e(b)
    local K = colsof(`b1')
    local nd 0
    forval c = 1/`K' {
      if `b1'[1,`c'] != `b2'[1,`c'] local ++nd
    }
    dis "DET_`tag'_K=`K' DET_`tag'_NDIFF=`nd'"
  }
  dis "DET_`tag'_RC=" _rc
end

use `base', clear
detchk , tag(plain)   opts(att)
detchk , tag(fastn)   opts(att fastnose)
detchk , tag(clusi)   opts(att cluster(i))
detchk , tag(rob)     opts(att robust)
detchk , tag(aw)      opts(att)
detchk , tag(allest)  opts(att att_i att_it att_itime)

* weighted + poisson + customvce need their own calls
capture noisily {
  tempname c1 c2
  qui wooldid y i t ttre [aw=w] , att
  matrix `c1' = e(b)
  qui wooldid y i t ttre [aw=w] , att
  matrix `c2' = e(b)
  local nd 0
  forval c = 1/`=colsof(`c1')' {
    if `c1'[1,`c'] != `c2'[1,`c'] local ++nd
  }
  dis "DET_awgt_NDIFF=`nd'"
}
dis "DET_awgt_RC=" _rc

capture noisily {
  tempname d1 d2
  qui wooldid ycount i t ttre , att poisson
  matrix `d1' = e(b)
  qui wooldid ycount i t ttre , att poisson
  matrix `d2' = e(b)
  local nd 0
  forval c = 1/`=colsof(`d1')' {
    if `d1'[1,`c'] != `d2'[1,`c'] local ++nd
  }
  dis "DET_pois_NDIFF=`nd'"
}
dis "DET_pois_RC=" _rc

capture noisily {
  tempname e1 e2
  qui wooldid y i t ttre , att customvce(cluster i)
  matrix `e1' = e(b)
  qui wooldid y i t ttre , att customvce(cluster i)
  matrix `e2' = e(b)
  local nd 0
  forval c = 1/`=colsof(`e1')' {
    if `e1'[1,`c'] != `e2'[1,`c'] local ++nd
  }
  dis "DET_cvce_NDIFF=`nd'"
}
dis "DET_cvce_RC=" _rc

capture noisily {
  tempname f1 f2
  qui wooldid y i t ttre , att contreat(dose) cluster(i)
  matrix `f1' = e(b)
  qui wooldid y i t ttre , att contreat(dose) cluster(i)
  matrix `f2' = e(b)
  local nd 0
  forval c = 1/`=colsof(`f1')' {
    if `f1'[1,`c'] != `f2'[1,`c'] local ++nd
  }
  dis "DET_cont_NDIFF=`nd'"
}
dis "DET_cont_RC=" _rc

*==============================================================================
* PART 2 -- RUNTIME DETECTOR for the claimed state.
* An element sits in e(perm_addone) iff (M<. & !memb & sobs<. & n>=2).
* An element sits in e(perm_incomplete) iff (M<. & sum(s<.) != M), i.e. exk=0.
* So: an element in perm_addone but NOT in perm_incomplete has exk=1 & memb=0 --
* exactly the state PERMEX-PX-2 claims.  Print both lists for every design.
*==============================================================================
prog def pxrep
  syntax , tag(string)
  local ao `"`e(perm_addone)'"'
  local in `"`e(perm_incomplete)'"'
  dis "PX_`tag'_EXACT=" e(perm_exact)
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
    forval r = 1/`nr' {
      dis "PX_`tag'_P`r'=" %21x `R'[`r',4]
    }
  }
end

use `base', clear
capture noisily {
  wooldid y i t ttre , att permute(shuffle) permexact
  pxrep , tag(A1)
}
dis "PX_A1_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att permute(shuffle) permexact
  pxrep , tag(A2)
}
dis "PX_A2_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att att_i att_it att_itime cluster(i) permute(shuffle) permexact permtstat
  pxrep , tag(B_tstat)
}
dis "PX_B_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att cluster(i) permute(shuffle) permexact customvce(cluster i)
  pxrep , tag(C_cvce)
}
dis "PX_C_RC=" _rc

use `base', clear
capture noisily {
  wooldid ycount i t ttre , att poisson permute(shuffle) permexact
  pxrep , tag(D_pois)
}
dis "PX_D_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre [aw=w] , att permute(shuffle) permexact
  pxrep , tag(E_awgt)
}
dis "PX_E_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att cluster(i) contreat(dose) permute(shuffle) permexact
  pxrep , tag(F_cont)
}
dis "PX_F_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att esprelength(2) espostlength(2) jointtests permute(shuffle) permexact
  pxrep , tag(G_es)
}
dis "PX_G_RC=" _rc

use `base', clear
capture noisily {
  wooldid y i t ttre , att permute(shuffle) permexact permanalyticse cluster(i)
  pxrep , tag(H_ase)
}
dis "PX_H_RC=" _rc

*==============================================================================
* PART 3 -- an INCOMPLETE enumeration: does exk=0 & memb=1 behave as the entry
* says (plain share over the reduced set, incomplete caution, no addone)?
* Unbalanced panel: drop the late periods for two clusters so some reassigned
* onsets are not estimable.
*==============================================================================
use `base', clear
drop if inlist(i, 5, 6) & t >= 8
capture noisily {
  wooldid y i t ttre , att esprelength(2) espostlength(3) permute(shuffle) permexact
  pxrep , tag(I_unbal)
}
dis "PX_I_RC=" _rc

dis "ALLDONE"
