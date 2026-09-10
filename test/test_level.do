* test_level.do - level() on the analytic path (F054 + F044, 2026-09-06).
* Run from test/: stata-mp -e do test_level.do -> expect LEVEL_ALL_PASS
*
* L1 analytic level(90): lb_90ci/ub_90ci columns, e(level) = 90, estimates/se identical to
*    the default run, the att interval = est -/+ invttail(df_r, .05) se, c(level) restored
* L2 the default run keeps lb_95ci/ub_95ci and posts e(level) = 95
* L3 a FAILING call with level(90) still restores c(level)
* L4 gteffects at level(90): lb_90ci columns and the same critical value
* L5 pretrendpower at level(90): the slope line prints "90% CI", and so does the replay
* L6 level(92.5) labels lb_92p5ci
* L7 subgroups: the sg tables and the pairwise table carry lb_90ci
* L8 bootstrap level(90) unchanged in kind: lb_90ci, e(boot_level) = 90, e(level) = 90,
*    and c(level) restored after the replication loop
* L9 the stray-companion screen no longer refuses level() without an engine, and the
*    range check still refuses level(5)
* L11 impvalues at level(90): the Ibragimov-Muller interval follows level() (bug sweep BS-5);
*    b/se/t/p bit-identical to the default run, the 95% interval pinned, columns lb_90ci
* L12 pretrendpower: the NIS pre-test alpha is pinned at 5% and does NOT follow level()
*    (bug sweep BS-11): e(roth_nissize_att) and every e(roth_power_att) cell identical at
*    level(90), the default-level numbers pinned, the slope CI still level-dependent, and
*    mdt_nis strictly positive at level(50) where it used to print an exact 0
clear all
set more off
do "../wooldidinstall/wooldid.ado"
local fails 0
set level 95
clear
set obs 12
gen unit = _n
gen s = ceil(unit/4)
gen ttre = cond(unit<=4, 4, cond(unit<=8, 6, .))
expand 8
bysort unit: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
tempfile panel
save `panel'

* ---- L1 / L2
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2)
matrix M95 = e(wdidmainresults)
local c95 : colnames M95
local lv95 = e(level)
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) level(90)
matrix M90 = e(wdidmainresults)
local c90 : colnames M90
local lv90 = e(level)
local dfr = e(df_r)
local clv = c(level)
local r1 = rowsof(M90)
local est = M90[`r1', 1]
local se = M90[`r1', 2]
local lbhand = `est' - invttail(`dfr', .05)*`se'
local ubhand = `est' + invttail(`dfr', .05)*`se'
mata: st_local("dbse", strofreal(max(abs(st_matrix("M90")[., 1..2] - st_matrix("M95")[., 1..2]))))
mata: st_local("dint", strofreal(max(abs(st_matrix("M90")[., 5..6] - st_matrix("M95")[., 5..6]))))
di "L1 cols90=[`c90'] e(level)=`lv90' c(level)=`clv' df_r=`dfr' lb=" %12.8f M90[`r1',5] " hand=" %12.8f `lbhand' " ub=" %12.8f M90[`r1',6] " hand=" %12.8f `ubhand' " |b,se diff|=`dbse' |int diff|=`dint'"
if !`: list posof "lb_90ci" in c90' | !`: list posof "ub_90ci" in c90' | `lv90' != 90 | `clv' != 95 ///
   | abs(M90[`r1',5] - `lbhand') > 1e-8 | abs(M90[`r1',6] - `ubhand') > 1e-8 | `dbse' != 0 | `dint' == 0 {
  di as error "L1 FAIL: level(90) is not honoured on the analytic path"
  local ++fails
}
di "L2 cols95=[`c95'] e(level)=`lv95'"
if !`: list posof "lb_95ci" in c95' | `lv95' != 95 {
  di as error "L2 FAIL: the default run does not carry lb_95ci / e(level) = 95"
  local ++fails
}

* ---- L3 a failing call restores c(level)
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) level(90) bootdf(clusters)
local rc3 = _rc
local clv3 = c(level)
di "L3 rc=`rc3' c(level)=`clv3'"
if `rc3' == 0 | `clv3' != 95 {
  di as error "L3 FAIL: a failing call with level(90) left c(level) = `clv3' (rc `rc3')"
  local ++fails
}

* ---- L4 gteffects
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) gteffects level(90)
local cg : colnames e(gteffects)
matrix G = e(gteffects)
local dfr4 = e(df_r)
local g4 = G[1, 5] - invttail(`dfr4', .05)*G[1, 6]
di "L4 gteffects cols=[`cg'] lb=" %12.8f G[1,8] " hand=" %12.8f `g4'
if !`: list posof "lb_90ci" in cg' | abs(G[1,8] - `g4') > 1e-8 {
  di as error "L4 FAIL: gteffects cells do not carry the 90% interval"
  local ++fails
}

* ---- L5 pretrendpower slope line at 90%, and its replay
use `panel', clear
qui log using level_l5.txt, text replace name(lv)
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) esfixedbaseperiod pretrendpower(50) level(90)
local rc5 = _rc
cap noi wooldid
local rc5b = _rc
qui log close lv
local lg = fileread("level_l5.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local n90 = 0
local pos = strpos(`"`lg'"', "; 90% CI ")
while `pos' > 0 {
  local ++n90
  local lg = substr(`"`lg'"', `pos' + 9, .)
  local pos = strpos(`"`lg'"', "; 90% CI ")
}
di "L5 rcs=`rc5'/`rc5b' slope lines at 90%: `n90'"
if `rc5' != 0 | `rc5b' != 0 | `n90' < 2 {
  di as error "L5 FAIL: the Roth slope line does not carry level() (fit + replay)"
  local ++fails
}
cap erase level_l5.txt

* ---- L6 non-integer level
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) level(92.5)
local c6 : colnames e(wdidmainresults)
di "L6 cols=[`c6'] e(level)=" e(level)
if !`: list posof "lb_92p5ci" in c6' | e(level) != 92.5 {
  di as error "L6 FAIL: level(92.5) is not labelled lb_92p5ci"
  local ++fails
}

* ---- L7 subgroups
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) subgroups(s) level(90)
local c7s : colnames e(wdidsgresults_1)
local c7p : colnames e(pwsgs)
di "L7 sg=[`c7s'] pw=[`c7p']"
if !`: list posof "lb_90ci" in c7s' | !`: list posof "lb_90ci" in c7p' {
  di as error "L7 FAIL: subgroup tables do not carry lb_90ci"
  local ++fails
}

* ---- L8 bootstrap
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) level(90)
local c8 : colnames e(wdidmainresults)
di "L8 cols=[`c8'] e(boot_level)=" e(boot_level) " e(level)=" e(level) " c(level)=" c(level)
if !`: list posof "lb_90ci" in c8' | e(boot_level) != 90 | e(level) != 90 | c(level) != 95 {
  di as error "L8 FAIL: bootstrap level(90) surface changed or c(level) not restored"
  local ++fails
}

* ---- L9 stray screen and range
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) level(90)
local rc9a = _rc
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) level(5)
local rc9b = _rc
di "L9 rcs=`rc9a' (want 0) / `rc9b' (want 198) c(level)=" c(level)
if `rc9a' != 0 | `rc9b' != 198 | c(level) != 95 {
  di as error "L9 FAIL: level() acceptance/range screen"
  local ++fails
}

* ---- L10 (CS-2, consistency sweep 2026-09-06): an UNSPECIFIED level() means c(level) on
* every path. Under `set level 90' with no level() option a bootstrap run and a permute run
* must equal the explicit level(90) runs bit for bit (columns, intervals, e(boot_level) /
* e(perm_level)), instead of labelling lb_90ci columns that hold 95% intervals
set level 90
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7)
matrix L10B = e(wdidmainresults)
local l10bc : colnames L10B
local l10bl = e(boot_level)
matrix L10BP = e(boot_ci_percentile)
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(7)
matrix L10P = e(wdidmainresults)
local l10pc : colnames L10P
local l10pl = e(perm_level)
set level 95
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(10) rseed(7) level(90)
matrix L10Bx = e(wdidmainresults)
matrix L10BPx = e(boot_ci_percentile)
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(7) level(90)
matrix L10Px = e(wdidmainresults)
di "L10 set-level-90: boot cols=[`l10bc'] e(boot_level)=`l10bl' mreldif vs level(90)=" mreldif(L10B, L10Bx) "/" mreldif(L10BP, L10BPx) " | perm cols=[`l10pc'] e(perm_level)=`l10pl' mreldif=" mreldif(L10P, L10Px) " c(level)=" c(level)
if !`: list posof "lb_90ci" in l10bc' | `l10bl' != 90 | mreldif(L10B, L10Bx) != 0 | mreldif(L10BP, L10BPx) != 0 ///
   | !`: list posof "lb_90ci" in l10pc' | `l10pl' != 90 | mreldif(L10P, L10Px) != 0 | c(level) != 95 {
  di as error "L10 FAIL: an unspecified level() under set level 90 does not match level(90) on the engines"
  local ++fails
}

* ---- L11 impvalues: the hand-built IM interval follows level() (BS-5)
use data/dta_binary.dta, clear
qui wooldid y i t ttre, att cluster(i) impvalues esprelength(2) espostlength(2)
matrix L11A = e(wdidmainresults)
use data/dta_binary.dta, clear
qui wooldid y i t ttre, att cluster(i) impvalues esprelength(2) espostlength(2) level(90)
matrix L11B = e(wdidmainresults)
local l11bc : colnames L11B
local l11lv = e(level)
local l11rn : rownames L11A
local k : list posof "att" in l11rn
local hw95 = (L11A[`k',6] - L11A[`k',5])/(2*L11A[`k',2])
local hw90 = (L11B[`k',6] - L11B[`k',5])/(2*L11B[`k',2])
mata: st_local("d14", strofreal(max(abs(st_matrix("L11B")[., 1..4] - st_matrix("L11A")[., 1..4]))))
di "L11 impvalues att b=" %10.8f L11A[`k',1] " se=" %10.8f L11A[`k',2] " 95 int=[" %10.8f L11A[`k',5] "," %10.8f L11A[`k',6] "] hw/se=" %9.6f `hw95' " vs invttail(49,.025)=" %9.6f invttail(49,.025) " | level(90) int=[" %10.8f L11B[`k',5] "," %10.8f L11B[`k',6] "] hw/se=" %9.6f `hw90' " vs invttail(49,.05)=" %9.6f invttail(49,.05) " cols=[`l11bc'] e(level)=`l11lv' |b..p diff|=`d14'"
if abs(L11A[`k',5] - .61519121) > 1e-7 | abs(L11A[`k',6] - .87494362) > 1e-7 | abs(`hw95' - invttail(49,.025)) > 1e-6 ///
   | abs(`hw90' - invttail(49,.05)) > 1e-6 | abs(L11B[`k',5] - .63671413) > 1e-7 | abs(L11B[`k',6] - .85342071) > 1e-7 ///
   | `d14' != 0 | !`: list posof "lb_90ci" in l11bc' | `l11lv' != 90 {
  di as error "L11 FAIL: the impvalues interval does not follow level()"
  local ++fails
}

* ---- L12 the Roth NIS pre-test alpha is pinned at 5% and ignores level() (BS-11)
set level 95
use data/dta_binary.dta, clear
qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) pretrendpower(80)
matrix L12A = e(roth_power_att)
local n12a = e(roth_nissize_att)
local s12a = e(roth_slope_ll_att)
use data/dta_binary.dta, clear
qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) pretrendpower(80) level(90)
matrix L12B = e(roth_power_att)
local n12b = e(roth_nissize_att)
local s12b = e(roth_slope_ll_att)
use data/dta_binary.dta, clear
qui wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(2) espostlength(2) pretrendpower(50) level(50)
matrix L12C = e(roth_power_att)
local m12 = L12C[1,3]
mata: st_local("d12", strofreal(max(abs(st_matrix("L12B") - st_matrix("L12A")))))
di "L12 nissize " %12.8f `n12a' " vs " %12.8f `n12b' " | mdt_nis " %12.8f L12A[1,3] " vs " %12.8f L12B[1,3] " | mdt_wald " %12.8f L12A[1,2] " | max|roth_power diff|=`d12' | slope_ll " %12.8f `s12a' " vs " %12.8f `s12b' " | level(50) mdt_nis=" %12.8f `m12' " c(level)=" c(level)
if abs(`n12a' - `n12b') > 1e-12 | `d12' > 1e-12 | abs(`n12a' - .07961397) > 1e-7 ///
   | abs(L12A[1,3] - .20181734) > 1e-7 | abs(L12A[1,5] - .69718718) > 1e-7 | abs(L12A[1,2] - .22199476) > 1e-7 ///
   | abs(`s12a' - `s12b') < 1e-8 | `m12' <= 0 | missing(`m12') {
  di as error "L12 FAIL: the Roth NIS pre-test alpha is not pinned at 5% (or the default-level numbers moved)"
  local ++fails
}
set level 95

if `fails' == 0 di as result "LEVEL_ALL_PASS"
else di as error "LEVEL_FAILS=`fails'"
