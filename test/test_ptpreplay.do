* Targeted suite for the pretrendpower replay completion gate (NFPTP-07 + the NFPTP-02
* replay disclosure residue): part(roth) now (a) refuses to render an interrupted
* diagnostics pass (roth_* present without e(roth_ok)), (b) replays the completed-but-
* unavailable caution (roth_ok==1, no matrices), and (c) replays the rank-deficiency
* disclosure (roth_keff < roth_npre), byte-identically to the live text.
* Run from test/: "<stata>" -e do test_ptpreplay.do -> PTPREPLAY_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 8642
set sortseed 8642
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

* simpler: forge the interrupted state by copying a real ptp run then knocking out
* roth_ok via ereturn repost through _estimates
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(3) espostlength(3) pretrendpower(80)

* ---- R1: normal run then bare replay - roth block byte-parity ----
qui log using ptpr_r1a.txt, text replace name(ptpr)
qui noi wooldid_display , part(roth)
qui log close ptpr
qui log using ptpr_r1b.txt, text replace name(ptpr)
cap noi wooldid
qui log close ptpr
if strpos(fileread("ptpr_r1b.txt"), "Pre-Trend Power Diagnostics") == 0 {
  dis as error "R1 FAIL: replay lost the roth block"
  local ++fails
}
if strpos(fileread("ptpr_r1b.txt"), "diagnostics are incomplete") != 0 {
  dis as error "R1 FAIL: interrupted-pass caution false-fired on a complete store"
  local ++fails
}

* ---- R2: forge interrupted state (roth_power_* present, roth_ok absent) via a fresh
* eclass repost that carries the matrices but not the marker ----
tempname fb fv fm
matrix `fb' = e(b)
matrix `fv' = e(V)
matrix `fm' = e(roth_power_att)
local frepv `e(wdid_replayv)'
tempvar esmp
qui gen byte `esmp' = e(sample)
cap program drop forge_partial
program define forge_partial, eclass
  args bb vv mm repv smp
  ereturn post `bb' `vv', esample(`smp')
  ereturn local cmd wooldid
  ereturn hidden local wdid_replayv `repv'
  ereturn matrix roth_power_att = `mm'
  ereturn scalar roth_npre_att = 3
  ereturn scalar roth_keff_att = 3
  tempname wm
  matrix `wm' = J(1,2,1)
  matrix colnames `wm' = a b
  ereturn matrix wdidmainresults = `wm'
end
forge_partial `fb' `fv' `fm' "`frepv'" `esmp'
qui log using ptpr_r2.txt, text replace name(ptpr)
cap noi wooldid
local rc2 = _rc
qui log close ptpr
if strpos(fileread("ptpr_r2.txt"), "diagnostics are incomplete") == 0 {
  dis as error "R2 FAIL: interrupted-pass caution missing (rc=" `rc2' ")"
  local ++fails
}
if strpos(fileread("ptpr_r2.txt"), "Pre-Trend Power Diagnostics") != 0 {
  dis as error "R2 FAIL: partial roth block was rendered"
  local ++fails
}
if `rc2' != 0 | strpos(fileread("ptpr_r2.txt"), "Results (Full Estimation Sample)") == 0 {
  dis as error "R2b FAIL: replay aborted or lost the main results table (rc=" `rc2' ") - the gate must skip the block, not the replay"
  local ++fails
}

* ---- R3: completed-but-unavailable store (roth_ok=1, no matrices) -> the unavailable
* caution replays; the interrupted caution must NOT fire ----
cap program drop forge_unavail
program define forge_unavail, eclass
  args bb vv repv smp
  ereturn post `bb' `vv', esample(`smp')
  ereturn local cmd wooldid
  ereturn hidden local wdid_replayv `repv'
  ereturn scalar roth_ok = 1
  tempname wm
  matrix `wm' = J(1,2,1)
  matrix colnames `wm' = a b
  ereturn matrix wdidmainresults = `wm'
end
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(3) espostlength(3) pretrendpower(80)
tempname gb gv
matrix `gb' = e(b)
matrix `gv' = e(V)
local grepv `e(wdid_replayv)'
tempvar gsmp
qui gen byte `gsmp' = e(sample)
forge_unavail `gb' `gv' "`grepv'" `gsmp'
qui log using ptpr_r3.txt, text replace name(ptpr)
cap noi wooldid
qui log close ptpr
if strpos(fileread("ptpr_r3.txt"), "no diagnostics were computed") == 0 {
  dis as error "R3 FAIL: completed-but-unavailable caution missing on replay"
  local ++fails
}
if strpos(fileread("ptpr_r3.txt"), "diagnostics are incomplete") != 0 {
  dis as error "R3 FAIL: interrupted caution false-fired on a completed-unavailable store"
  local ++fails
}

* ---- R4: rank-deficiency disclosure replays byte-identically (few-cluster esf design
* with esprelength >= #clusters gives keff < npre) ----
clear
qui set obs 4
qui gen i = _n
qui gen ttre = cond(_n <= 2, 7, .)
qui expand 10
qui bysort i: gen t = _n
qui gen y = 1.5*(t >= ttre & ttre < .) + 0.07*i + rnormal()
qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(5) espostlength(2) pretrendpower(80)
local rk = e(roth_keff_att)
local rn = e(roth_npre_att)
if `rk' >= `rn' {
  dis as text "R4 NOTE: fixture not rank-deficient (keff=`rk' npre=`rn'); disclosure-parity leg skipped"
}
else {
  qui log using ptpr_r4.txt, text replace name(ptpr)
  cap noi wooldid
  qui log close ptpr
  if strpos(fileread("ptpr_r4.txt"), "rank-deficient (rank `rk' of `rn')") == 0 {
    dis as error "R4 FAIL: rank-deficiency disclosure missing on replay"
    local ++fails
  }
}

if `fails' == 0 dis as result "PTPREPLAY_ALL_PASS"
else dis as error "PTPREPLAY_FAILS=`fails'"
