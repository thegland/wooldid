* Targeted suite for the issue-file permutation batch (NFPERM-05/06/09/15):
* - 05: permprobs() strict onset-level tokens (garbage/duplicate/extended-missing loud)
* - 06: probability rows within the 1e-6 tolerance are rescaled to sum exactly 1
* - 09: designs admitting exactly one assignment are refused (459)
* - 15: permute() x contreatcontrols(untreatedZero) refused (198)
* Run from test/: "<stata>" -e do test_permute2.do -> PERMUTE2_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop mkperm2
program define mkperm2
  clear
  set seed 90210
  set sortseed 90210
  qui set obs 36
  qui gen i = _n
  qui gen c = ceil(_n/3)
  qui gen ttre = cond(c <= 4, 4, cond(c <= 8, 6, .))
  qui expand 8
  qui bysort i: gen t = _n
  qui gen y = 2*(t >= ttre & ttre < .) + 0.05*i + rnormal()
end

* ---- P2a: duplicate listed onset level -> r(198), was accepted rc 0 with level 6
* silently absent from the randomization design ----
mkperm2
local pr "(ttre 4 4 .)"
forvalues cc = 1/12 {
  local pr `pr' (`cc' .3 .4 .3)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 198 {
  dis as error "P2a FAIL: duplicate level rc=" _rc " (expected 198)"
  local ++fails
}

* ---- P2b: garbage token (aliases the never-treated level) -> r(198) ----
mkperm2
local pr "(ttre 4 garbage .)"
forvalues cc = 1/12 {
  local pr `pr' (`cc' .3 .4 .3)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 198 {
  dis as error "P2b FAIL: garbage token rc=" _rc " (expected 198)"
  local ++fails
}

* ---- P2c: garbage token WITHOUT '.' listed (the caution-free silent class: garbage's
* probability column silently applied to never-treated) -> r(198) ----
mkperm2
local pr "(ttre 4 6 junk)"
forvalues cc = 1/12 {
  local pr `pr' (`cc' .3 .4 .3)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 198 {
  dis as error "P2c FAIL: junk-for-dot token rc=" _rc " (expected 198)"
  local ++fails
}

* ---- P2d: extended-missing token .a -> r(198) via the strict token rule ----
mkperm2
local pr "(ttre 4 6 .a)"
forvalues cc = 1/12 {
  local pr `pr' (`cc' .3 .4 .3)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 198 {
  dis as error "P2d FAIL: .a token rc=" _rc " (expected 198)"
  local ++fails
}

* ---- P2e: valid spec still runs ----
mkperm2
local pr "(ttre 4 6 .)"
forvalues cc = 1/12 {
  local pr `pr' (`cc' .3 .4 .3)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 0 {
  dis as error "P2e FAIL: valid spec rc=" _rc
  local ++fails
}

* ---- P3: NFPERM-06 engine-mode normalization - deficit row accepted and rescaled ----
mata: wooldid_permp("(ttre 4 .) (1 .5 .4999992) (2 .3 .7)", "ttre", "4 .", "1 2", "npP", "npL")
if "`permperr'" != "" {
  dis as error "P3 FAIL: deficit row rejected: `permperr'"
  local ++fails
}
else {
  mata: st_local("npbad", strofreal(max(abs(rowsum(st_matrix("npP")) :- 1)) > 1e-12))
  if "`npbad'" == "1" {
    dis as error "P3 FAIL: stored probability rows do not sum to 1 after normalization"
    local ++fails
  }
}
cap matrix drop npP
cap matrix drop npL

* ---- P4a: NFPERM-09 all-singleton strata (support 1) -> r(459), was rc 0 with p=1 ----
mkperm2
qui gen st = c
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permexact permstrata(st) esprelength(1) espostlength(1)
if _rc != 459 {
  dis as error "P4a FAIL: singleton strata rc=" _rc " (expected 459)"
  local ++fails
}

* ---- P4b: onset-homogeneous multi-cluster strata -> r(459) ----
mkperm2
qui gen st = cond(c <= 4, 1, cond(c <= 8, 2, 3))
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permexact permstrata(st) esprelength(1) espostlength(1)
if _rc != 459 {
  dis as error "P4b FAIL: homogeneous strata rc=" _rc " (expected 459)"
  local ++fails
}

* ---- P4c: all-deterministic probability rows targeting the observed assignment -> r(459) ----
mkperm2
local pr "(ttre 4 6 .)"
forvalues cc = 1/12 {
  if `cc' <= 4 local pr `pr' (`cc' 1 0 0)
  else if `cc' <= 8 local pr `pr' (`cc' 0 1 0)
  else local pr `pr' (`cc' 0 0 1)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 459 {
  dis as error "P4c FAIL: deterministic-at-observed rc=" _rc " (expected 459)"
  local ++fails
}

* ---- P4d: deterministic rows targeting a NON-observed assignment -> r(459) (formerly
* the worst subcase: ran with p = 1/(R+1)) ----
mkperm2
local pr "(ttre 4 6 .)"
forvalues cc = 1/12 {
  if `cc' <= 4 local pr `pr' (`cc' 0 1 0)
  else if `cc' <= 8 local pr `pr' (`cc' 1 0 0)
  else local pr `pr' (`cc' 0 0 1)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 459 {
  dis as error "P4d FAIL: deterministic-at-non-observed rc=" _rc " (expected 459)"
  local ++fails
}

* ---- P4e: acceptance - one mixed stratum among singletons still randomizes ----
mkperm2
qui gen st = c
qui replace st = 1 if c == 5
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permexact permstrata(st) esprelength(1) espostlength(1)
if _rc != 0 | e(perm_npossible) != 2 {
  dis as error "P4e FAIL: mixed stratum rc=" _rc " npos=" e(perm_npossible) " (expected rc 0, npos 2)"
  local ++fails
}

* ---- P4f: acceptance - one stochastic probability row still randomizes ----
mkperm2
local pr "(ttre 4 6 .)"
forvalues cc = 1/12 {
  if `cc' == 1 local pr `pr' (`cc' .5 .5 0)
  else if `cc' <= 4 local pr `pr' (`cc' 1 0 0)
  else if `cc' <= 8 local pr `pr' (`cc' 0 1 0)
  else local pr `pr' (`cc' 0 0 1)
}
set seed 1234
cap wooldid y i t ttre, att cluster(c) permute(probabilities) permprobs(`pr') permreps(25) esprelength(1) espostlength(1)
if _rc != 0 {
  dis as error "P4f FAIL: one-stochastic-row design rc=" _rc
  local ++fails
}

* ---- P5: NFPERM-15 untreatedZero x permute refused at parse (198), before any RNG ----
use "data/dta_cont", clear
set seed 1234
cap wooldid y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) permute(shuffle) permreps(20) esprelength(1) espostlength(1)
if _rc != 198 {
  dis as error "P5 FAIL: untreatedZero x permute rc=" _rc " (expected 198)"
  local ++fails
}
use "data/dta_cont", clear
set seed 1234
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) permute(shuffle) permreps(20) esprelength(1) espostlength(1)
if _rc != 0 {
  dis as error "P5b FAIL: contreat basic x permute rc=" _rc " (contreat in general stays allowed)"
  local ++fails
}

* ---- P6: NFPERM-02 keyed-row application - stratified assignments stay attached to
* their clusters when ascending cluster order crosses stratum order (frlink re-sorts
* the map frame; formerly the positional read scrambled assignments: exact designs
* reported p=0 with the observed assignment absent from their own enumeration) ----
cap program drop mkcross
program define mkcross
  args mode
  clear
  set obs 12
  qui gen i = _n
  qui gen c = ceil(_n/3)
  if "`mode'" == "cross" {
    qui gen st = cond(c == 2 | c == 3, 1, 2)
    qui gen ttre = cond(c == 2, 4, cond(c == 1, 6, .))
  }
  else {
    qui gen st = cond(c <= 2, 1, 2)
    qui gen ttre = cond(c == 1, 4, cond(c == 3, 6, .))
  }
  qui expand 10
  qui bysort i: gen t = _n
  qui gen y = 2*(t >= ttre & ttre < .) + 0.05*i + 0.3*sin(3*i + t)
end
mkcross cross
set seed 777
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permstrata(st) permexact esprelength(1) espostlength(1)
if _rc != 0 | e(perm_npossible) != 4 | !(e(perm_p_att) >= 1/4) {
  dis as error "P6a FAIL: crossing-strata exact rc=" _rc " npos=" e(perm_npossible) " p=" e(perm_p_att) " (observed assignment must be a member: p >= 1/4)"
  local ++fails
}
mkcross aligned
set seed 777
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permstrata(st) permexact esprelength(1) espostlength(1)
if _rc != 0 | !(e(perm_p_att) >= 1/4) {
  dis as error "P6b FAIL: aligned-strata exact rc=" _rc " p=" e(perm_p_att)
  local ++fails
}
* descending-contiguous variant: stratum 2 = low cluster IDs, stratum 1 = high
clear
set obs 12
qui gen i = _n
qui gen c = ceil(_n/3)
qui gen st = cond(c <= 2, 2, 1)
qui gen ttre = cond(c == 1, 4, cond(c == 3, 6, .))
qui expand 10
qui bysort i: gen t = _n
qui gen y = 2*(t >= ttre & ttre < .) + 0.05*i + 0.3*sin(3*i + t)
set seed 777
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permstrata(st) permexact esprelength(1) espostlength(1)
if _rc != 0 | !(e(perm_p_att) >= 1/4) {
  dis as error "P6c FAIL: descending-strata exact rc=" _rc " p=" e(perm_p_att)
  local ++fails
}

* ---- P7: NFPERM-03 exact-one mask - bare-flag if values 0/1/2/-1/. must match if==1
* exactly; the excluded rows are DESIGN-RELEVANT (they add a cluster + an onset level) ----
cap program drop mkflag
program define mkflag
  clear
  set obs 12
  qui gen i = _n
  qui gen c = ceil(_n/3)
  qui gen ttre = cond(c <= 2, 4, .)
  qui expand 8
  qui bysort i: gen t = _n
  qui gen y = 2*(t >= ttre & ttre < .) + 0.05*i + 0.3*sin(3*i + t)
  qui gen flag = 1
  * a fifth cluster carrying a NEW onset level, reachable only via flag values != 1
  qui set obs `=_N+4'
  forvalues k = 1/4 {
    qui replace i = 100 + `k' in `=_N-4+`k''
    qui replace c = 5 in `=_N-4+`k''
    qui replace ttre = 6 in `=_N-4+`k''
    qui replace t = `k' + 2 in `=_N-4+`k''
    qui replace y = 0.1*`k' in `=_N-4+`k''
  }
  qui replace flag = 2 if c == 5 & mod(i,2) == 0
  qui replace flag = . if c == 5 & mod(i,2) == 1
  qui replace flag = -1 if c == 5 & i == 101
end
mkflag
set seed 555
cap wooldid y i t ttre if flag, att cluster(c) permute(shuffle) permexact esprelength(1) espostlength(1)
local rcA = _rc
local pA = e(perm_p_att)
local nposA = e(perm_npossible)
mkflag
set seed 555
cap wooldid y i t ttre if flag == 1, att cluster(c) permute(shuffle) permexact esprelength(1) espostlength(1)
if `rcA' != _rc | `nposA' != e(perm_npossible) | `pA' != e(perm_p_att) {
  dis as error "P7 FAIL: bare-flag vs ==1 diverge: rc " `rcA' "/" _rc " npos " `nposA' "/" e(perm_npossible) " p " `pA' "/" e(perm_p_att)
  local ++fails
}

* ---- P7b: a cluster existing solely via flag!=1 rows no longer enters the design;
* with only one onset group left among if==1 clusters the freedom gate refuses r459 ----
clear
set obs 12
qui gen i = _n
qui gen c = ceil(_n/3)
qui gen ttre = cond(c <= 4, 4, .)
qui expand 8
qui bysort i: gen t = _n
qui gen y = 2*(t >= ttre & ttre < .) + 0.05*i + 0.3*sin(3*i + t)
qui gen flag = 1
qui replace flag = 2 if c == 4
qui replace ttre = 6 if c == 4
set seed 555
cap wooldid y i t ttre if flag, att cluster(c) permute(shuffle) permexact esprelength(1) espostlength(1)
local rc7b = _rc
set seed 555
cap wooldid y i t ttre if flag == 1, att cluster(c) permute(shuffle) permexact esprelength(1) espostlength(1)
if `rc7b' != _rc {
  dis as error "P7b FAIL: bare-flag rc=" `rc7b' " vs ==1 rc=" _rc " (must match; both are the post-fix design)"
  local ++fails
}

* ---- P8: NFPERM-10 comparator-validity invariant - a permutation run may either refuse
* loudly (missing/failed rep-0 comparator -> r459) or complete, but a completed run must
* NEVER carry a missing comparator (the silent-p=0/tail-corruption class). the noiseless
* fixture probes the nastiest corner (se ~ roundoff under permtstat) ----
clear
set obs 12
qui gen i = _n
qui gen c = ceil(_n/3)
qui gen ttre = cond(c <= 2, 4, .)
qui expand 6
qui bysort i: gen t = _n
qui gen y = 2*(t >= ttre & ttre < .) + 0.05*i
set seed 555
cap wooldid y i t ttre, att cluster(c) permute(shuffle) permexact permtstat esprelength(1) espostlength(1)
local rc8 = _rc
if `rc8' != 0 & `rc8' != 459 {
  dis as error "P8 FAIL: noiseless permtstat rc=" `rc8' " (expected 0 or the loud 459 refusal)"
  local ++fails
}
if `rc8' == 0 {
  cap confirm matrix e(perm_comparator)
  if _rc != 0 {
    dis as error "P8 FAIL: completed run without e(perm_comparator)"
    local ++fails
  }
  else if missing(e(perm_comparator)[1,1]) {
    dis as error "P8 FAIL: completed run with a MISSING comparator (the guard must refuse instead)"
    local ++fails
  }
}

* ---- P9: NFPERM-12 comparator export + exact reconstruction; NFPERM-11 population SD;
* NFPERM-01 no-false-fire (zero-failure exact run must NOT draw the failure caution) ----
mkperm2
qui keep if c <= 4 | c >= 9
set seed 999
qui log using perm2_p9.txt, text replace name(p2)
cap noi wooldid y i t ttre, att cluster(c) permute(shuffle) permexact esprelength(1) espostlength(1)
local rc9 = _rc
qui log close p2
if `rc9' != 0 {
  dis as error "P9 FAIL: rc=" `rc9'
  local ++fails
}
else {
  cap confirm matrix e(perm_comparator)
  if _rc != 0 {
    dis as error "P9 FAIL: e(perm_comparator) not returned"
    local ++fails
  }
  else {
    local cn : colnames e(perm_comparator)
    if "`cn'" != "att" {
      dis as error "P9 FAIL: comparator colnames [`cn'] != estimands [att]"
      local ++fails
    }
    tempname PD
    matrix `PD' = e(permdist)
    local recount 0
    forvalues r = 1/`=rowsof(`PD')' {
      if abs(`PD'[`r',1]) >= abs(e(perm_comparator)[1,1]) local ++recount
    }
    if reldif(`recount'/e(perm_completed), e(perm_p_att)) > 1e-12 {
      dis as error "P9 FAIL: reconstruction " `recount' "/" e(perm_completed) " != e(perm_p_att) " e(perm_p_att)
      local ++fails
    }
    * NFPERM-11: population SD of the enumerated distribution
    qui clear
    qui svmat double `PD', names(pd)
    qui sum pd1
    local psd = sqrt(r(Var)*(r(N)-1)/r(N))
    if reldif(`psd', e(perm_sd_att)) > 1e-12 {
      dis as error "P9 FAIL: exact-mode SD " e(perm_sd_att) " != population SD " `psd'
      local ++fails
    }
  }
  if e(perm_failed) == 0 & strpos(fileread("perm2_p9.txt"), "enumeration is incomplete") != 0 {
    dis as error "P9 FAIL: exact-failure caution false-fired on a zero-failure run"
    local ++fails
  }
}

if `fails' == 0 dis as result "PERMUTE2_ALL_PASS"
else dis as error "PERMUTE2_FAILS=`fails'"
