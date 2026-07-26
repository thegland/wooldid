* wooldid regression battery. Run from a fresh subdirectory of test/ (see test/README.md):
*   mkdir test/run1 && cd test/run1 && "<stata>" -e do ../runbattery.do
* Loads the CURRENT repo wooldid.ado, runs 18 configs, dumps all e() content
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

di "BATTERY_DONE"
exit, clear
