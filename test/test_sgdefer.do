* test_sgdefer.do - the subgroup level-margins deferral (EFF3) and its boundary.
* Run from test/: stata-mp -e do test_sgdefer.do
*
* At the BINARY subgroup sites wooldid reads the per-subgroup level margins off the
* pwcompare(effects) pass instead of issuing a second identical margins call. At the
* CONTREAT sites it deliberately does NOT, and D1/D2 below are why.
*
* The pwcompare pass does not emit margins' `o.' omitted markers. At the binary sites
* only the second half of the vector (targetflag==1) is consumed and the markers live in
* the discarded first half, so nothing downstream depends on them. At the contreat sites
* the WHOLE vector is consumed, and the truerow/sghaspre filter uses exactly those
* markers to drop a subgroup's pre-treatment row when the per-sg continuous
* pre-treatment margins come back all-omitted. Deferring there let an unmarked row
* survive the filter: e(b) gained an all-zero pretreatment_SG<k>_ame column, every later
* coefficient shifted one position into the wrong name, and e(wdidsgresults_*) gained a
* row of raw missings. That is the same shifted/duplicated-colname symptom the
* cont_sg_tidose battery config was added to pin.
*
* D1 the all-omitted contreat shape produces NO pretreatment_SG*_ame column
* D2 e(b) column names are unique (no shifted/duplicated names) on that shape
* D3 the deferral is live on the binary side: multi-subgroup runs still report per-sg
*    levels, and a single-subgroup run (where the pwcompare pass never runs, so the
*    dedicated call must be retained) reports them too
* D4 a subgroup whose pre-period pwcompare fails still gets its level margins (fallback)
* D5 subgroups whose treated margin is degenerate do not gain spurious rows
* Prints SGDEFER_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

* ---- the cont_sg_tidose shape: time-invariant dose makes the per-sg continuous
* pre-treatment margins come back all-omitted.
clear
set seed 447102
set obs 50
gen i = _n
gen ttre = cond(_n<=13, 5, cond(_n<=26, 8, cond(_n<=38, 10, .)))
expand 13
bysort i: gen t = _n
gen dose = 0.5 + mod(i,6) + 0.25*mod(i,7)/7
gen post = (t >= ttre & ttre < .)
gen y = 0.4*dose*post + 0.02*i + 0.08*t + sin(3*i + 2*t)
gen sg = mod(i,3)
drop post
cap noi wooldid y i t ttre, contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg)
local rcA = _rc
matrix bA = e(b)
local cnA : colnames bA

* ---- D1: no per-sg pretreatment ame column on the all-omitted shape
local d1bad 0
foreach nm in pretreatment_SG0_ame pretreatment_SG1_ame pretreatment_SG2_ame {
  if strpos("`cnA'", "`nm'") > 0 local ++d1bad
}
di "D1_spurious_pretreat_ame_cols=" `d1bad'
if `rcA' != 0 | `d1bad' > 0 {
  di as error "D1 FAIL: rc `rcA'; the all-omitted contreat shape produced `d1bad' spurious pretreatment_SG*_ame column(s) - the contreat deferral must stay off"
  local ++fails
}

* ---- D2: column names unique (the shifted/duplicated-colname symptom)
local nA : word count `cnA'
local uniqA 0
local seen
foreach c of local cnA {
  if strpos(" `seen' ", " `c' ") == 0 {
    local seen `seen' `c'
    local ++uniqA
  }
}
di "D2_bcols=" `nA' " unique=" `uniqA'
if `nA' != `uniqA' {
  di as error "D2 FAIL: e(b) has `nA' columns but only `uniqA' unique names - names have shifted"
  local ++fails
}

* ---- D3: the binary deferral is live and still reports per-sg levels, with and
* without the pwcompare pass (>1 vs ==1 realized subgroup)
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg3 = mod(i,3)+1
cap qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg3) safety(off)
local rc3 = _rc
local mm3 : e(matrices)
local n3 = 0
foreach m of local mm3 {
  if strpos("`m'", "wdidsgresults_") == 1 local ++n3
}
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg1 = 7
cap qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg1) safety(off)
local rc1 = _rc
local mm1 : e(matrices)
local n1 = 0
foreach m of local mm1 {
  if strpos("`m'", "wdidsgresults_") == 1 local ++n1
}
di "D3_multi_sg_matrices=" `n3' " single_sg_matrices=" `n1'
if `rc3' != 0 | `n3' != 3 | `rc1' != 0 | `n1' != 1 {
  di as error "D3 FAIL: multi-sg rc `rc3' n `n3' (want 3); single-sg rc `rc1' n `n1' (want 1)"
  local ++fails
}

* ---- D4: pre-period pwcompare failure still yields level margins via the fallback
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sg2 = 1
qui replace sg2 = 2 if mod(i,7) == 0
qui drop if sg2 == 2 & t < ttre - 1 & ttre != .
cap qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) esfixedbaseperiod esrelativeto(-1) subgroups(sg2) safety(off)
local rc4 = _rc
local mm4 : e(matrices)
local n4 = 0
foreach m of local mm4 {
  if strpos("`m'", "wdidsgresults_") == 1 local ++n4
}
di "D4_prefail_rc=" `rc4' " sg_matrices=" `n4'
if `rc4' != 0 | `n4' < 1 {
  di as error "D4 FAIL: the pre-pwcompare fallback lost the level margins (rc `rc4', n `n4')"
  local ++fails
}

* ---- D5: a never-treated subgroup must not gain spurious rows
use "`REPO'/test/data/dta_binary.dta", clear
qui gen sgx = 1
qui replace sgx = 2 if mod(i,5) == 0
qui replace ttre = . if sgx == 2
cap qui wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sgx) safety(off)
local rc5 = _rc
matrix b5 = e(b)
local cn5 : colnames b5
local n5 : word count `cn5'
local u5 = 0
local seen5
foreach c of local cn5 {
  if strpos(" `seen5' ", " `c' ") == 0 {
    local seen5 `seen5' `c'
    local ++u5
  }
}
di "D5_rc=" `rc5' " bcols=" `n5' " unique=" `u5'
if `rc5' != 0 | `n5' != `u5' {
  di as error "D5 FAIL: never-treated subgroup produced duplicated e(b) names (rc `rc5', `n5' cols, `u5' unique)"
  local ++fails
}

if `fails' == 0 di as result "SGDEFER_ALL_PASS"
else di as error "SGDEFER_FAILS=`fails'"
