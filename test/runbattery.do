* wooldid regression battery. Run from a fresh subdirectory of test/ (see test/README.md):
*   mkdir test/run1 && cd test/run1 && "<stata>" -e do ../runbattery.do
* Loads the CURRENT repo wooldid.ado, runs 34 configs, dumps all e() content
* (exact %21x hex) + rc + post-run varlist to <cfg>.csv in cwd.
clear all
set more off
set maxvar 20000
set seed 98765
set sortseed 98765

local REPO "../.."
local DATA "../data"

do "`REPO'/wooldidinstall/wooldid.ado"

cap program drop dumpres
program define dumpres
  args cfg rc
  tempname fh
  file open `fh' using "`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    if _rc != 0 continue
    local R = rowsof(`M')
    local C = colsof(`M')
    local rn : rowfullnames `M'
    local cn : colfullnames `M'
    file write `fh' `"matnames|`m'|`R'x`C'|`rn'|`cn'"' _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  file close `fh'
end

* helper to reset state between configs
cap program drop freshstart
program define freshstart
  args dtafile
  ereturn clear
  use "`dtafile'", clear
end

* ---------------- configs ----------------

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats jointtests
local rc = _rc
dumpres basic `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres allest_aw `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2) adversarial
local rc = _rc
dumpres esf_adv `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att_it cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests
local rc = _rc
dumpres sg_attit `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att_itime cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests
local rc = _rc
dumpres sg_attitime `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) impvalues
local rc = _rc
dumpres impv `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid ycount i t ttre, att cluster(i) poisson
local rc = _rc
dumpres pois `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [pw=w], att cluster(i) esprelength(2) espostlength(2) summarystats
local rc = _rc
dumpres pw_ss `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [aw=w], att robust unconditionalse esprelength(2) espostlength(2)
local rc = _rc
dumpres aw_uncond `rc'

freshstart "`DATA'/dta_tiny"
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
local rc = _rc
dumpres tiny `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(-999)
local rc = _rc
dumpres esptrip `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) esfixedbaseperiod esprelength(3) espostlength(3) cluster(i)
local rc = _rc
dumpres cont_basic `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres cont_uz `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2) summarystats
local rc = _rc
dumpres cont_uz_ss `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre [aw=w], contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres cont_lat `rc'

freshstart "`DATA'/dta_cont"
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) summarystats
local rc = _rc
dumpres cont_ss `rc'

freshstart "`DATA'/dta_simple1"
cap noi wooldid y i t ttre, att cluster(i) subgroup(ttre)
local rc = _rc
dumpres simple1 `rc'

freshstart "`DATA'/dta_simple1"
cap noi wooldid y i t ttre, att cluster(i) subgroup(ttre) esprelength(3) espostlength(3) esfixedbaseperiod
local rc = _rc
dumpres simple_esf `rc'

* E1 coverage (added 2026-07-26, blessed pre-E1): the contreat SG-comparison block had
* zero battery coverage; sg_missif pins missing-sg + bare-flag-if through the SG block;
* sg_degpre pins the degenerate no-pre-period state (absent pre joints, post-only pwsgs)
freshstart "`DATA'/dta_cont"
qui gen sg = mod(i,2)+1
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests
local rc = _rc
dumpres cont_sg `rc'

freshstart "`DATA'/dta_binary"
qui gen sgm = mod(i,3)
qui replace sgm = . if mod(i,11) == 0
cap noi wooldid y i t ttre if w > 0, att att_it cluster(i) esprelength(2) espostlength(2) subgroups(sgm) jointtests
local rc = _rc
dumpres sg_missif `rc'

freshstart "`DATA'/dta_binary"
qui sum t, meanonly
local t2 = r(min) + 1
qui gen sgd2 = 1
qui replace sgd2 = 2 if mod(i,7) == 0
qui replace ttre = `t2' if mod(i,7) == 0
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sgd2) jointtests
local rc = _rc
dumpres sg_degpre `rc'

* B16 coverage (added 2026-07-26, blessed pre-fix): the esf adversarial reference-row
* planting and the altmethod_adversarial pre-trend cut had no battery coverage; these
* four pin the affected cells (adv_altb: below-window reference regime with altmethod;
* adv_r2001: degenerate 1-row fit, rc 2001 at bless time; adv_gapa/adv_gapb: gap-year
* data where the reference-gap hunt aborts (rc 198) or plants at the wrong horizon)
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-3) esprelength(2) espostlength(2) adversarial altmethod_adversarial
local rc = _rc
dumpres adv_altb `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(1) espostlength(2) adversarial altmethod_adversarial
local rc = _rc
dumpres adv_r2001 `rc'

freshstart "`DATA'/dta_path_relyrgap"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-1) esprelength(4) espostlength(4) adversarial
local rc = _rc
dumpres adv_gapa `rc'

freshstart "`DATA'/dta_path_relyrgap"
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-4) esprelength(3) espostlength(3) adversarial
local rc = _rc
dumpres adv_gapb `rc'

* dupcol coverage (added 2026-07-26, blessed pre-fix): continuous + subgroups where the
* time-invariant dose makes the per-sg continuous pre-treatment margins come back
* all-omitted, tripping the sg rowname filter (silently shifted/duplicated colnames in
* e(b)/e(V)/e(wdidsgresults_*) on the pre-fix build: 45 columns, 42 unique names)
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
ereturn clear
cap noi wooldid y i t ttre, contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg)
local rc = _rc
dumpres cont_sg_tidose `rc'

freshstart "`DATA'/dta_path_relyrgap"
gen dose2 = 0.5 + mod(i,5) + 0.2*mod(i,3)
gen sgr = mod(i,3)
cap noi wooldid y i t ttre, contreat(dose2) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sgr)
local rc = _rc
dumpres cont_sg_gap `rc'

* F3 coverage (added 2026-07-27, blessed on the feature build): pins e(gteffects) content
* and the off-state e() surface of a weighted two-estimand run
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) esprelength(2) espostlength(2) gteffects
local rc = _rc
dumpres gt_basic `rc'

* F6 coverage (added 2026-07-27, blessed on the feature build): pins the permutation
* e() surface; the seed is set INSIDE the config so the dump is independent of
* battery ordering (draws consume the session RNG stream)
freshstart "`DATA'/dta_binary"
set seed 246813
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20)
local rc = _rc
dumpres perm_basic `rc'

* New-feature benchmark pins (added 2026-08-24, blessed in-commit; PREREG_NEWFEAT
* chunk): close the battery's coverage gaps for fastnose, pretrendpower,
* anticipation, permtstat/permtci, and permute(probabilities). Appended AFTER all
* prior configs so earlier dumps are bit-unaffected; permute configs set their own
* seed (perm_basic precedent).
freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fastnose gteffects
local rc = _rc
dumpres fn_gt `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod pretrendpower(50 80)
local rc = _rc
dumpres ptp_two `rc'

freshstart "`DATA'/dta_binary"
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) anticipation(1)
local rc = _rc
dumpres antic_es `rc'

freshstart "`DATA'/dta_binary"
set seed 135791
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permtstat permtci
local rc = _rc
dumpres perm_tstat `rc'

* 10-cluster subset spanning all three onset blocks (2005: i 1-25, 2008: 26-50,
* never: 51-80); permprobs needs one row per in-sample cluster value
freshstart "`DATA'/dta_binary"
qui keep if inlist(i, 1, 2, 3, 26, 27, 28, 51, 52, 53, 54)
set seed 975310
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(30) permprobs((ttre 2005 2008 .) (1 .3 .3 .4) (2 .3 .3 .4) (3 .3 .3 .4) (26 .3 .3 .4) (27 .3 .3 .4) (28 .3 .3 .4) (51 .3 .3 .4) (52 .3 .3 .4) (53 .3 .3 .4) (54 .3 .3 .4))
local rc = _rc
dumpres perm_probs `rc'

di "BATTERY_DONE"
exit, clear
