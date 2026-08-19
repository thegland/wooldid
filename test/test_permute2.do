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

if `fails' == 0 dis as result "PERMUTE2_ALL_PASS"
else dis as error "PERMUTE2_FAILS=`fails'"
