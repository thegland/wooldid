* test_adopath.do - wooldid loaded the way an installed user loads it: resolved from the
* adopath, never with `do'. Every other cert suite loads the ado with `do <path>', under which
* cleanmatrices' `clear mata' drops wooldid's own Mata kernels for the rest of the session (see
* test_stored S9 and the SI-17 note in test_permute2). Resolved from the adopath, Stata re-reads
* the file when a Mata function is missing, so the same runs recover; nothing pinned that until
* now (change audit 2026-09-22, cleanmatrices-cert-gap).
* Run from test/: stata-mp -e do test_adopath.do -> expect ADOPATH_ALL_PASS
*
* A0  the session resolves wooldid from ../wooldidinstall, not from an installed copy
* A1  plain run (event study, subgroups, gteffects), then the same run with cleanmatrices:
*     every compared e() matrix is identical and e(wdid_vcvstate) is the same analytic-blocks
*     (under do-loading the cleanmatrices run falls back to analytic-diagonal)
* A2  the plain run again after the cleanmatrices run: identical to the first
* A3  permute() + cleanmatrices completes and matches permute() without it, same seed
clear all
set more off
set graphics off
local fails 0
discard
adopath ++ "../wooldidinstall"

* ---- A0: resolution
qui findfile wooldid.ado
local a0fn `"`r(fn)'"'
dis `"A0: wooldid.ado resolves to `a0fn'"'
if strpos(`"`a0fn'"', "wooldidinstall") == 0 {
  dis as error "A0 FAIL: wooldid.ado resolves outside ../wooldidinstall; the legs below would test another build"
  dis as error "ADOPATH_FAILS=1"
  exit
}

clear
set obs 48
gen unit = _n
gen s = cond(mod(unit,2)==0, 1, 2)
gen ttre = cond(mod(unit,4)==0, 5, cond(mod(unit,4)==1, 7, .))
expand 12
bysort unit: gen t = _n
set seed 31
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + 0.3*s*(t>=ttre & ttre<.) + rnormal()
tempfile panel
save `panel'

* matrices compared across runs. cleanmatrices runs `matrix drop _all' and `clear mata', so the
* runs are kept with estimates store and read back afterwards
local a1mats b V wdidmainresults gteffects wdidsgresults_1 wdidsgresults_2 b_es_att V_es_att b_es_att_it V_es_att_it
cap program drop cmpest
program define cmpest, rclass
  args e1 e2 mats
  local worst 0
  local bad
  foreach m of local mats {
    qui estimates restore `e1'
    cap confirm matrix e(`m')
    if _rc != 0 {
      local bad `bad' `m'(absent-in-`e1')
      continue
    }
    tempname M1 M2
    matrix `M1' = e(`m')
    qui estimates restore `e2'
    cap confirm matrix e(`m')
    if _rc != 0 {
      local bad `bad' `m'(absent-in-`e2')
      continue
    }
    matrix `M2' = e(`m')
    if rowsof(`M1') != rowsof(`M2') | colsof(`M1') != colsof(`M2') {
      local bad `bad' `m'(shape)
      continue
    }
    local d = mreldif(`M1', `M2')
    if `d' > `worst' local worst = `d'
  }
  return scalar worst = `worst'
  return local bad `bad'
end

* ---- A1: plain vs cleanmatrices
use `panel', clear
cap noi wooldid y unit t ttre, att att_it cluster(unit) esprelength(2) espostlength(2) subgroups(s) gteffects
local rc1a = _rc
local vs1a "`e(wdid_vcvstate)'"
qui estimates store a_plain
use `panel', clear
cap noi wooldid y unit t ttre, att att_it cluster(unit) esprelength(2) espostlength(2) subgroups(s) gteffects cleanmatrices
local rc1b = _rc
local vs1b "`e(wdid_vcvstate)'"
qui estimates store a_clean
cmpest a_plain a_clean "`a1mats'"
local w1 = r(worst)
local b1 `r(bad)'
dis "A1: plain rc=`rc1a' state=`vs1a' | cleanmatrices rc=`rc1b' state=`vs1b' | worst mreldif=`w1' problems=[`b1']"
if `rc1a' != 0 | `rc1b' != 0 | "`vs1a'" != "analytic-blocks" | "`vs1b'" != "`vs1a'" | `w1' != 0 | "`b1'" != "" {
  dis as error "A1 FAIL: cleanmatrices under adopath loading must reproduce the plain run exactly"
  local ++fails
}

* ---- A2: the session is not left degraded
use `panel', clear
cap noi wooldid y unit t ttre, att att_it cluster(unit) esprelength(2) espostlength(2) subgroups(s) gteffects
local rc2 = _rc
local vs2 "`e(wdid_vcvstate)'"
qui estimates store a_after
cmpest a_plain a_after "`a1mats'"
local w2 = r(worst)
local b2 `r(bad)'
dis "A2: plain after cleanmatrices rc=`rc2' state=`vs2' | worst mreldif=`w2' problems=[`b2']"
if `rc2' != 0 | "`vs2'" != "`vs1a'" | `w2' != 0 | "`b2'" != "" {
  dis as error "A2 FAIL: a plain run after a cleanmatrices run must reproduce the first plain run"
  local ++fails
}

* ---- A3: permute() + cleanmatrices (the SI-17 combination) completes and matches
use `panel', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(4242) dots(0)
local rc3a = _rc
qui estimates store a_perm
use `panel', clear
cap noi wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(4242) dots(0) cleanmatrices
local rc3b = _rc
qui estimates store a_permc
cmpest a_perm a_permc "b wdidmainresults"
local w3 = r(worst)
local b3 `r(bad)'
dis "A3: permute rc=`rc3a' | permute + cleanmatrices rc=`rc3b' | worst mreldif=`w3' problems=[`b3']"
if `rc3a' != 0 | `rc3b' != 0 | `w3' != 0 | "`b3'" != "" {
  dis as error "A3 FAIL: permute() with cleanmatrices must complete under adopath loading and match the run without it"
  local ++fails
}

estimates clear
if `fails' == 0 dis as result "ADOPATH_ALL_PASS"
else dis as error "ADOPATH_FAILS=`fails'"
