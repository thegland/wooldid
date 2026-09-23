* Targeted suite for the subgroup-split gate (author directive 2026-07-28): a unit
* split across subgroups over time is refused under default safety, and proceeds with
* share weighting plus a caution under safety(off) on the analytic path - while under
* permute() the same split is refused at BOTH safety settings (G5), as is a unit labelled on
* some rows and missing on others, never-treated units included (G6). Also pins the numeric-only
* cluster() syntax. Run from test/: "<stata>" -e do test_sgsplit.do -> SGSPLIT_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 24680
set sortseed 24680
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

* panel with one unit whose subgroup changes over time
cap program drop mksplit
program define mksplit
  clear
  set seed 24680
  set sortseed 24680
  qui set obs 30
  qui gen i = _n
  qui gen ttre = cond(_n <= 10, 5, .)
  qui expand 10
  qui bysort i: gen t = _n
  qui gen y = 2*(t >= ttre & ttre < .) + 0.1*i + rnormal()
  qui gen sg = mod(i, 2)
  qui replace sg = 1 - sg if i == 3 & t >= 6
end

* ---- G1: split unit + default safety -> informative r(459) ----
mksplit
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 459 {
  dis as error "G1 FAIL: rc=" _rc " (expected 459)"
  local ++fails
}

* ---- G2: split unit + safety(off) -> runs with the share-weighting caution ----
mksplit
qui log using sgsp_g2.txt, text replace name(sgsp)
cap noi wooldid y i t ttre, att att_i cluster(i) esprelength(2) espostlength(2) subgroups(sg) safety(off)
local rc2 = _rc
qui log close sgsp
if `rc2' != 0 {
  dis as error "G2 FAIL: rc=" `rc2'
  local ++fails
}
else if strpos(fileread("sgsp_g2.txt"), "splits some units") == 0 | strpos(fileread("sgsp_g2.txt"), "share weighting") == 0 {
  dis as error "G2 FAIL: share-weighting caution missing under safety(off)"
  local ++fails
}

* ---- G3: unit-constant subgroups unaffected, both safety settings ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 {
  dis as error "G3 FAIL: clean sg run rc=" _rc
  local ++fails
}
tempname b3
matrix `b3' = e(b)
use "`DATA'/dta_binary", clear
qui log using sgsp_g3.txt, text replace name(sgsp)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) safety(off)
qui log close sgsp
if _rc != 0 | mreldif(`b3', e(b)) != 0 {
  dis as error "G3 FAIL: safety(off) clean sg run rc=" _rc " mreldif=" mreldif(`b3', e(b))
  local ++fails
}
if strpos(fileread("sgsp_g3.txt"), "splits some units") != 0 {
  dis as error "G3 FAIL: split caution fired on unit-constant subgroups"
  local ++fails
}

* ---- G4: string cluster variables are refused at parse (numeric-only syntax) ----
use "`DATA'/dta_binary", clear
qui gen str8 sclus = "c" + string(i)
cap wooldid y i t ttre, att cluster(sclus)
if _rc != 109 {
  dis as error "G4 FAIL: rc=" _rc " (expected 109 from varlist numeric)"
  local ++fails
}
cap wooldid y i t ttre, att cluster(i)
if _rc != 0 {
  dis as error "G4 FAIL: numeric cluster run rc=" _rc
  local ++fails
}

* ---- G5: under permute() the split-unit refusal is UNLIFTABLE by safety(off) ----
* Ruled by the author 2026-09-21. The worker's own screen (G2) keeps its escape; the
* permutation screen does not, because each replication assigns one onset per cluster and a
* unit in two subgroups has no defined subgroup in a replication that treats it. The refusal
* is raised in wooldid_permctl during validation, before any replication runs, so this leg is
* cheap. Both safety settings must give 459.
mksplit
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(20)
if _rc != 459 {
  dis as error "G5 FAIL: permute + split sg under default safety rc=" _rc " (expected 459)"
  local ++fails
}
mksplit
qui log using sgsp_g5.txt, text replace name(sgsp5)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(20) safety(off)
local rc5 = _rc
qui log close sgsp5
if `rc5' != 459 {
  dis as error "G5 FAIL: safety(off) lifted the permute split-sg refusal, rc=" `rc5' " (expected 459)"
  local ++fails
}
else if strpos(fileread("sgsp_g5.txt"), "multiple subgroup variable levels") == 0 {
  dis as error "G5 FAIL: the permute split-sg refusal did not name itself under safety(off)"
  local ++fails
}
cap erase sgsp_g5.txt

* ---- G6: under permute() a unit labelled on some rows and missing on others is refused over
* the whole design sample, never-treated units included, and safety(off) does not lift it. The
* observed run's own screen reads treated units only, so a never-treated unit's partial label
* is legal there; a replication that treated the unit would enter it into its subgroup by its
* labelled share of sample weight. Refused in wooldid_permctl before any replication runs.
cap program drop mkpartial
program define mkpartial
  clear
  set seed 24680
  set sortseed 24680
  qui set obs 30
  qui gen i = _n
  qui gen ttre = cond(_n <= 10, 5, .)
  qui expand 10
  qui bysort i: gen t = _n
  qui gen y = 2*(t >= ttre & ttre < .) + 0.1*i + rnormal()
  qui gen sg = mod(i, 2)
  qui replace sg = . if i == 25 & t >= 6
end
mkpartial
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg)
if _rc != 0 {
  dis as error "G6 FAIL: a never-treated unit's partial label refused without permute(), rc=" _rc " (expected 0)"
  local ++fails
}
mkpartial
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(20)
if _rc != 459 {
  dis as error "G6 FAIL: permute + partially labelled never-treated unit rc=" _rc " (expected 459)"
  local ++fails
}
mkpartial
qui log using sgsp_g6.txt, text replace name(sgsp6)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(20) safety(off)
local rc6 = _rc
qui log close sgsp6
if `rc6' != 459 {
  dis as error "G6 FAIL: safety(off) lifted the permute partial-label refusal, rc=" `rc6' " (expected 459)"
  local ++fails
}
else if strpos(fileread("sgsp_g6.txt"), "Error: the subgroup variable sg is missing for some in-sample row(s) of units of i that carry a subgroup value on their other rows.") == 0 {
  dis as error "G6 FAIL: the permute partial-label refusal did not name itself"
  local ++fails
}
cap erase sgsp_g6.txt
* the same unit labelled on every row runs
mkpartial
qui replace sg = mod(i, 2)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(20)
if _rc != 0 {
  dis as error "G6 FAIL: permute with the unit fully labelled rc=" _rc " (expected 0)"
  local ++fails
}

if `fails' == 0 dis "SGSPLIT_ALL_PASS"
else dis as error "SGSPLIT_FAILS=`fails'"
exit, clear
