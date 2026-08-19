* Targeted suite for the second issue-file low batch:
* - VAL-G13A: cfxplottypes() exact-token parse (garbage keyword-adjacent tokens and
*   fractional/missing numeric entries now error loudly; + strip and window-drop
*   semantics preserved; 2.0-spelling status quo preserved per the c48bdaa scoping)
* - VAL-G09C: F13 onset-constancy message uses group-defined-by-i terminology
* - VAL-G07 (docs-only resolution): pins the empirical facts the docs now state -
*   commensurate fractional t is bit-identical on the default path; the native-factor
*   paths (unconditionalse) reject it loudly
* Run from test/: "<stata>" -e do test_cfxtok.do -> CFXTOK_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 13579
set sortseed 13579
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

* ---- C1: fractional cfxplottypes entry -> loud r(198), was a silent no-op ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxplottypes(1.5)
if _rc != 198 {
  dis as error "C1 FAIL: cfxplottypes(1.5) rc=" _rc " (expected 198)"
  local ++fails
}

* ---- C2: postpost was substring-matched as post -> now unrecognized r(198) ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxplottypes(postpost)
if _rc != 198 {
  dis as error "C2 FAIL: cfxplottypes(postpost) rc=" _rc " (expected 198)"
  local ++fails
}

* ---- C3: valid tokens unchanged ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxplottypes(post pre 1 -1)
if _rc != 0 {
  dis as error "C3 FAIL: valid cfxplottypes rc=" _rc
  local ++fails
}

* ---- C4: leading + still stripped (c48bdaa regression guard) ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxplottypes(+2)
if _rc != 0 {
  dis as error "C4 FAIL: cfxplottypes(+2) rc=" _rc " (the + strip must survive the parser rewrite)"
  local ++fails
}

* ---- C5: pre without esf/pre-periods still silently dropped (documented precondition) ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) espostlength(2) makecfxplots cfxplottypes(pre post)
if _rc != 0 {
  dis as error "C5 FAIL: cfxplottypes(pre post) without espre rc=" _rc " (pre should drop silently)"
  local ++fails
}

* ---- C6: 2.0 spelling status quo (passes parse per the c48bdaa author scoping, then
* fails loudly downstream at name construction) ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxplottypes(2.0)
if _rc != 198 {
  dis as error "C6 FAIL: cfxplottypes(2.0) rc=" _rc " (expected the loud downstream 198, status quo)"
  local ++fails
}

* ---- C7: F13 message terminology (group defined by i) ----
clear
qui set obs 30
qui gen i = _n
qui gen ttre = cond(_n <= 10, 5, .)
qui expand 10
qui bysort i: gen t = _n
qui replace ttre = 6 if i == 3 & t >= 6
qui gen y = 2*(t >= ttre & ttre < .) + 0.1*i + rnormal()
qui log using cfxtok_c7.txt, text replace name(cfxtok)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
local rc7 = _rc
qui log close cfxtok
if `rc7' != 459 {
  dis as error "C7 FAIL: split-onset rc=" `rc7' " (expected 459)"
  local ++fails
}
else if strpos(fileread("cfxtok_c7.txt"), "group defined by") == 0 {
  dis as error "C7 FAIL: F13 message lacks group-defined-by-i terminology"
  local ++fails
}

* ---- C8: commensurate fractional t is bit-identical on the default path (the fact
* that made VAL-G07's proposed hard screen an over-rejection; docs-only resolution) ----
use "`DATA'/dta_binary", clear
qui replace t = t + 0.5
qui replace ttre = ttre + 0.5
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
tempname BF
matrix `BF' = e(b)
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
if mreldif(`BF', e(b)) != 0 {
  dis as error "C8 FAIL: fractional-t twin not bit-identical to integer twin, mreldif=" mreldif(`BF', e(b))
  local ++fails
}

* ---- C9: native-factor path rejects fractional t loudly (documented constraint) ----
use "`DATA'/dta_binary", clear
qui replace t = t + 0.5
qui replace ttre = ttre + 0.5
cap wooldid y i t ttre, att unconditionalse robust esprelength(1) espostlength(1)
if _rc != 452 {
  dis as error "C9 FAIL: uncond+fractional-t rc=" _rc " (expected loud 452)"
  local ++fails
}

if `fails' == 0 dis as result "CFXTOK_ALL_PASS"
else dis as error "CFXTOK_FAILS=`fails'"
