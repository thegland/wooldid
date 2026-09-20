*==============================================================================
* p3_replay.do -- EXACT-IDENTITY test of `replay' (bare wooldid after an
*                 estimation) against the live display.
* Every wooldid call is issued from inside a PROGRAM so Stata does not echo it
* into the captured log; the live and replay outputs then differ only in
* wooldid's own output.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

local dN 120
local dT 12
local donsets "2004 2006 2008"
include "`W'/mkdata.do"
gen double ypos = exp(y/3)
save "`W'/d3.dta", replace

cap prog drop RUNIT
prog def RUNIT
  args nothing
  wooldid $CMDBODY
end
cap prog drop REPIT
prog def REPIT
  wooldid
end

cap prog drop dopair
prog def dopair
  syntax , tag(string) body(string) w(string)
  global CMDBODY `body'
  dis as result "==== REPLAY CONFIG `tag' : wooldid `body' ===="
  qui log using "`w'/rep_`tag'_live.txt", replace text name(LL)
  cap noi RUNIT
  local rc = _rc
  qui log close LL
  qui log using "`w'/rep_`tag'_replay.txt", replace text name(RR)
  cap noi REPIT
  local rc2 = _rc
  qui log close RR
  dis "REPLAY_`tag'_RCS=`rc',`rc2'"
end

dopair , tag(r01) w("`W'") body(y i t ttre, att att_it att_i att_itime cluster(i) espre(3) espost(3) jointtests adversarial)
dopair , tag(r02) w("`W'") body(y i t ttre, att subgroups(sg) robust espre(3) espost(3) jointtests adversarial)
dopair , tag(r03) w("`W'") body(y i t ttre, att robust contreat(dose) espre(2) espost(2))
dopair , tag(r04) w("`W'") body(ypos i t ttre, att poisson robust espre(2) espost(2))
dopair , tag(r05) w("`W'") body(y i t ttre, att att_it robust espre(3) espost(3) pretrendpower(50 80))
dopair , tag(r06) w("`W'") body(y i t ttre, att robust espre(3) espost(3) gteffects summarystats)
dopair , tag(r07) w("`W'") body(y i t ttre, att_it cluster(i) espre(2) espost(2) impvalues)
dopair , tag(r08) w("`W'") body(y i t ttre, att cluster(i) espre(3) espost(3) permute(shuffle) permreps(25) rseed(77) dots(0) jointtests)
dopair , tag(r09) w("`W'") body(y i t ttre, att cluster(i) espre(3) espost(3) bootstrap(standard) reps(25) rseed(77) dots(0) jointtests pretrendpower(80))
dopair , tag(r10) w("`W'") body(y i t ttre, att cluster(i) espre(3) espost(3) bootstrap(standard) reps(25) rseed(77) dots(0) simultaneousbands)
dopair , tag(r11) w("`W'") body(y i t ttre, att robust espre(3) espost(3) level(80) jointtests)
dopair , tag(r12) w("`W'") body(y i t ttre, att att_it att_i att_itime robust espre(3) espost(3) adversarial jointtests subgroups(sg) gteffects)

* ---------------------------------------------------------------------------
* STALENESS: what does the user see after a FAILED wooldid?
* ---------------------------------------------------------------------------
dis as result "==== STALENESS ===="
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
dis "S_goodrc=" _rc
dis "S_ecmd=[" e(cmd) "]  e(N)=" e(N)
cap noi wooldid y i t ttre, att robust controls(nosuchvariable)
dis "S_failrc=" _rc
dis "S_LASTRC=[$WOOLDID_LASTRC]"
dis "S_ecmd_after_fail=[" e(cmd) "]"
qui log using "`W'/rep_stale_refused.txt", replace text name(SS)
cap noi REPIT
dis "S_replay_rc=" _rc
qui log close SS

* defeat the guard with macro drop, as recon reports
macro drop WOOLDID_LASTRC
qui log using "`W'/rep_stale_defeated.txt", replace text name(TT)
cap noi REPIT
dis "S_replay_defeated_rc=" _rc
qui log close TT

* what is in e() after a failure that gets FAR enough to leave reghdfe's model?
cap noi wooldid y i t ttre if t < 2002, att robust espre(3) espost(3)
dis "S2_failrc=" _rc
dis "S2_ecmd=[" e(cmd) "]"
dis "S2_LASTRC=[$WOOLDID_LASTRC]"
macro drop WOOLDID_LASTRC
qui log using "`W'/rep_stale_reghdfe.txt", replace text name(UU)
cap noi REPIT
dis "S2_replay_rc=" _rc
qui log close UU

dis "P3_DONE"
