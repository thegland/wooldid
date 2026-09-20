* probe11b.do -- cell-level (gteffects) and subgroup surfaces under the two bootstrap legs,
* on designs with enough clusters.  Under cluster() the two legs draw the SAME resamples
* (probe1/cmp1), so every cell number must agree.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 1212
local gN 60
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui gen long cl12 = mod(unit,12) + 1
save base11b.dta, replace

cap prog drop mcmp
prog def mcmp
  syntax , m1(string) m2(string) tag(string) [maxcol(integer 0) maxrow(integer 0)]
  cap confirm matrix `m1'
  local r1 = _rc
  cap confirm matrix `m2'
  local r2 = _rc
  if `r1' | `r2' {
    dis "MCMP|`tag'|MISSING m1rc=`r1' m2rc=`r2'"
    exit
  }
  local nr = min(rowsof(`m1'), rowsof(`m2'))
  local nc = min(colsof(`m1'), colsof(`m2'))
  if `maxcol' > 0 & `maxcol' < `nc' local nc `maxcol'
  if `maxrow' > 0 & `maxrow' < `nr' local nr `maxrow'
  dis "MCMP|`tag'|dims m1=" rowsof(`m1') "x" colsof(`m1') " m2=" rowsof(`m2') "x" colsof(`m2') " compared " `nr' "x" `nc'
  local worst = 0
  local wr 0
  local wc 0
  local nmiss 0
  forvalues r = 1/`nr' {
    forvalues c = 1/`nc' {
      local a = `m1'[`r',`c']
      local b = `m2'[`r',`c']
      if missing(`a') & missing(`b') continue
      if missing(`a') | missing(`b') {
        local ++nmiss
        continue
      }
      local d = reldif(`a', `b')
      if `d' > `worst' {
        local worst = `d'
        local wr `r'
        local wc `c'
      }
    }
  }
  local cn : colfullnames `m1'
  local wcn : word `wc' of `cn'
  dis "MCMP|`tag'|MAXRELDIF=" %12.4e `worst' "|at row " `wr' " col " `wc' " (" "`wcn'" ")|one-sided-missing cells=" `nmiss'
end

local R 40
local SD 7070

dis "=========== NESTING regime: i=unit, cluster(cl12) (12 clusters nest units) ==========="
use base11b.dta, clear
cap noi wooldid ync unit t ttre, att gteffects cluster(cl12) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "N1_RC=" _rc
mat GN1 = e(gteffects)
mat BN1 = e(b)
mat VN1 = e(V)
dis "N1|b_att=" %21x BN1[1,colnumb(BN1,"main_att")] "|se_att=" %21x sqrt(VN1[colnumb(BN1,"main_att"),colnumb(BN1,"main_att")])
dis "N1|gt dims=" rowsof(GN1) "x" colsof(GN1)

use base11b.dta, clear
cap noi wooldid ync unit t ttre, att att_i gteffects cluster(cl12) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "N2_RC=" _rc
mat GN2 = e(gteffects)
mat BN2 = e(b)
mat VN2 = e(V)
dis "N2|b_att=" %21x BN2[1,colnumb(BN2,"main_att")] "|se_att=" %21x sqrt(VN2[colnumb(BN2,"main_att"),colnumb(BN2,"main_att")])
dis "N2|gt dims=" rowsof(GN2) "x" colsof(GN2)
cap noi mcmp , m1(GN1) m2(GN2) tag("NESTING gteffects COUNT vs REID") maxcol(11)
dis "MC1_RC=" _rc

dis "=========== KEEP-i regime: i=cohort, cluster(unit) (cluster finer than i) ==========="
use base11b.dta, clear
cap noi wooldid ync cohort t ttre, att gteffects cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "F1_RC=" _rc
mat GF1 = e(gteffects)
mat BF1 = e(b)
mat VF1 = e(V)
dis "F1|b_att=" %21x BF1[1,colnumb(BF1,"main_att")] "|se_att=" %21x sqrt(VF1[colnumb(BF1,"main_att"),colnumb(BF1,"main_att")])

use base11b.dta, clear
cap noi wooldid ync cohort t ttre, att att_i gteffects cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "F2_RC=" _rc
mat GF2 = e(gteffects)
mat BF2 = e(b)
mat VF2 = e(V)
dis "F2|b_att=" %21x BF2[1,colnumb(BF2,"main_att")] "|se_att=" %21x sqrt(VF2[colnumb(BF2,"main_att"),colnumb(BF2,"main_att")])
cap noi mcmp , m1(GF1) m2(GF2) tag("KEEP-i gteffects COUNT vs REID") maxcol(11)
dis "MC2_RC=" _rc

dis "=========== SUBGROUPS under the two legs, cluster(cl12) ==========="
use base11b.dta, clear
cap noi wooldid ync unit t ttre, att subgroups(sg2) cluster(cl12) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "S1_RC=" _rc
mat BS1 = e(b)
mat VS1 = e(V)
cap mat PS1 = e(pwsgs)
dis "S1|colnames=" "`: colfullnames BS1'"
local K = colsof(BS1)
forvalues k = 1/`K' {
  local nm : word `k' of `: colfullnames BS1'
  dis "S1|`nm'|b=" %21x BS1[1,`k'] "|se=" %21x sqrt(VS1[`k',`k'])
}

use base11b.dta, clear
cap noi wooldid ync unit t ttre, att att_i subgroups(sg2) cluster(cl12) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "S2_RC=" _rc
mat BS2 = e(b)
mat VS2 = e(V)
cap mat PS2 = e(pwsgs)
local K = colsof(BS2)
forvalues k = 1/`K' {
  local nm : word `k' of `: colfullnames BS2'
  dis "S2|`nm'|b=" %21x BS2[1,`k'] "|se=" %21x sqrt(VS2[`k',`k'])
}
cap noi mcmp , m1(PS1) m2(PS2) tag("subgroup pairwise contrasts COUNT vs REID") maxrow(1)
dis "MC3_RC=" _rc

dis "PROBE11C_DONE"
