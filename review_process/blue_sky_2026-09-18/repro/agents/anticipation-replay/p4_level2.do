*==============================================================================
* p4_level2.do -- second wave on level().
*  (i)  EXACT EQUIVALENCE: `level(L)' with the session at 95  ==  `set level L'
*       with no level() option.  Both are documented to mean the same thing
*       (sthlp:540, "the default is c(level)"), so every posted quantity must
*       agree bit for bit.
*  (ii) level() must move ONLY intervals: with a fixed rseed, every non-interval
*       quantity must be bit-identical across level(80) and level(95).
*  (iii) interval identity on the z/poisson, contreat, simultaneousbands,
*       pretrendpower, permute and bootstrap surfaces.
*  (iv)  non-integer level labels on the bootstrap pairwise-subgroup table.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
include "`W'/dumpe.do"
which wooldid

local dN 90
local dT 12
local donsets "2004 2006 2008"
include "`W'/mkdata.do"
gen double ypos = exp(y/3)
save "`W'/d4.dta", replace

* ==========================================================================
* (i) level(L) == set level L
* ==========================================================================
cap prog drop lvlpair
prog def lvlpair
  syntax , tag(string) lev(string) opts(string) w(string)
  dis as result "==== LVLPAIR `tag' at `lev' : `opts' ===="
  qui use "`w'/d4.dta", clear
  set level 95
  cap noi wooldid `opts' level(`lev')
  dis "OPT_RC_`tag'=" _rc
  if _rc == 0 dumpe , f("`w'/lv_`tag'_opt.txt")
  qui use "`w'/d4.dta", clear
  set level `lev'
  cap noi wooldid `opts'
  dis "SET_RC_`tag'=" _rc
  if _rc == 0 dumpe , f("`w'/lv_`tag'_set.txt")
  set level 95
end

lvlpair , tag(L1) lev(80) w("`W'") opts(y i t ttre, att att_it att_i att_itime robust espre(3) espost(3) jointtests subgroups(sg) gteffects)
lvlpair , tag(L2) lev(80) w("`W'") opts(y i t ttre, att cluster(i) espre(2) espost(2) bootstrap(standard) reps(25) rseed(4242) dots(0) jointtests)
lvlpair , tag(L3) lev(80) w("`W'") opts(y i t ttre, att cluster(i) espre(2) espost(2) permute(shuffle) permreps(25) rseed(4242) dots(0) permtstat permtci)
lvlpair , tag(L4) lev(90) w("`W'") opts(y i t ttre, att robust espre(3) espost(3) esfixedbaseperiod pretrendpower(80))
lvlpair , tag(L5) lev(90) w("`W'") opts(ypos i t ttre, att poisson robust espre(2) espost(2))
lvlpair , tag(L6) lev(90) w("`W'") opts(y i t ttre, att robust contreat(dose) contreatcontrols(basic) espre(2) espost(2) makecfxplots)

* ==========================================================================
* (ii) level() must move ONLY intervals (fixed rseed)
* ==========================================================================
cap prog drop lvlinv
prog def lvlinv
  syntax , tag(string) opts(string) w(string)
  dis as result "==== LVLINV `tag' : `opts' ===="
  qui use "`w'/d4.dta", clear
  cap noi wooldid `opts' level(80)
  dis "INV80_RC_`tag'=" _rc
  if _rc == 0 dumpe , f("`w'/iv_`tag'_80.txt")
  qui use "`w'/d4.dta", clear
  cap noi wooldid `opts' level(95)
  dis "INV95_RC_`tag'=" _rc
  if _rc == 0 dumpe , f("`w'/iv_`tag'_95.txt")
end

lvlinv , tag(I1) w("`W'") opts(y i t ttre, att att_it robust espre(3) espost(3) jointtests adversarial)
lvlinv , tag(I2) w("`W'") opts(y i t ttre, att cluster(i) espre(2) espost(2) bootstrap(standard) reps(25) rseed(4242) dots(0))
lvlinv , tag(I3) w("`W'") opts(y i t ttre, att cluster(i) espre(2) espost(2) permute(shuffle) permreps(25) rseed(4242) dots(0))

* ==========================================================================
* (iii)/(iv) label + interval identity on the replication surfaces
* ==========================================================================
qui use "`W'/d4.dta", clear
cap noi wooldid y i t ttre, att cluster(i) espre(2) espost(2) subgroups(sg) ///
    bootstrap(standard) reps(30) rseed(99) dots(0) level(92.5)
dis "BT925_RC=" _rc
dis "BT925_level=" e(level)
foreach M in wdidmainresults pwsgs wdidsgresults_1 {
  cap confirm matrix e(`M')
  if _rc == 0 {
    tempname X
    mat `X' = e(`M')
    local cn : colnames `X'
    dis "BT925 `M' colnames: `cn'"
  }
}
cap noi ereturn list

* permute at a non-integer level
qui use "`W'/d4.dta", clear
cap noi wooldid y i t ttre, att cluster(i) espre(2) espost(2) subgroups(sg) ///
    permute(shuffle) permreps(25) rseed(99) dots(0) level(92.5)
dis "PM925_RC=" _rc
dis "PM925_level=" e(level)
foreach M in wdidmainresults pwsgs perm_ci {
  cap confirm matrix e(`M')
  if _rc == 0 {
    tempname X
    mat `X' = e(`M')
    local cn : colnames `X'
    dis "PM925 `M' colnames: `cn'"
  }
}

* simultaneousbands, analytic
qui use "`W'/d4.dta", clear
cap noi wooldid y i t ttre, att robust espre(3) espost(3) simultaneousbands level(80)
dis "SB_RC=" _rc
cap noi ereturn list

dis "P4_DONE"

* ==========================================================================
* (v) LABEL ROUND-TRIP: the worker labels its tables from `wdlvl' directly,
*     but the gteffects Mata kernel rebuilds the label from the round-trip
*     100 - 200*((100-L)/200).  Any level where the two disagree gives a
*     gteffects column name that does not match the main table's.
* ==========================================================================
dis as result "==== LABEL ROUND-TRIP SWEEP ===="
mata:
  real scalar L, gtal
  string scalar wk, gt
  for (L = 10; L <= 99.99; L = L + 0.01) {
    L = round(L*100)/100
    wk = (mod(L,1)==0 ? strofreal(L) : subinstr(strofreal(L,"%9.0g"), ".", "p", .))
    gtal = (100 - L)/200
    gt = subinstr(strofreal(100 - 200*gtal, "%9.0g"), ".", "p", .)
    if (wk != gt) printf("LABELMISMATCH L=%9.4f worker=%s gteffects=%s\n", L, wk, gt)
  }
  printf("LABELSWEEP_DONE\n")
end

* ==========================================================================
* (vi) gteffects cell CI identity at a non-default level
* ==========================================================================
qui use "`W'/d4.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2) gteffects level(80)
dis "GT80_RC=" _rc
dis "GT80_df_r=" e(df_r)
cap noi {
  tempname G
  mat `G' = e(gteffects)
  local cn : colnames `G'
  dis "GT80 colnames: `cn'"
  local crit = invttail(e(df_r), .10)
  local worst = 0
  forvalues r = 1/`=rowsof(`G')' {
    local b = `G'[`r',5]
    local s = `G'[`r',6]
    local lo = `G'[`r',8]
    local hi = `G'[`r',9]
    if `b' < . & `s' < . & `lo' < . {
      local e1 = `b' - `crit'*`s'
      local d = abs(`lo'-`e1')/max(abs(`e1'),1e-12)
      if `d' > `worst' local worst = `d'
    }
  }
  dis "GT80_worst_rel_dev=" %12.4e `worst'
}
dis "GTCHK_RC=" _rc

dis "P4_EXTRA_DONE"
