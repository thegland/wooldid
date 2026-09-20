*==============================================================================
* p5_replay2.do -- adversarial replay probes that nobody suggested.
*  (a) histogramcohorteffects + bootstrap()/permute(): the shell prints a caution
*      (ado:299, ado:511) that wooldid_display never prints, so replay should
*      lose it.  EXACT test: diff live vs replay.
*  (b) replay after the DATA has been replaced / cleared.
*  (c) replay after estimates store / restore.
*  (d) replay after another estimation command has been run in between.
*  (e) replay twice in a row (idempotence).
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

local dN 100
local dT 12
local donsets "2004 2006 2008"
include "`W'/mkdata.do"
save "`W'/d5.dta", replace

cap prog drop RUNIT
prog def RUNIT
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
  dis as result "==== P5 CONFIG `tag' : wooldid `body' ===="
  qui use "`w'/d5.dta", clear
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

* (a) the shell-only cautions
dopair , tag(h01) w("`W'") body(y i t ttre, att cluster(i) espre(2) espost(2) histogramcohorteffects bootstrap(standard) reps(25) rseed(5) dots(0))
dopair , tag(h02) w("`W'") body(y i t ttre, att cluster(i) espre(2) espost(2) histogramcohorteffects permute(shuffle) permreps(25) rseed(5) dots(0))
dopair , tag(h03) w("`W'") body(y i t ttre, att cluster(i) espre(2) espost(2) histogramcohorteffects)

* (b) replay after the data is gone / replaced
qui use "`W'/d5.dta", clear
global CMDBODY y i t ttre, att att_it robust espre(3) espost(3) jointtests subgroups(sg)
qui log using "`W'/rep_b01_live.txt", replace text name(LL)
cap noi RUNIT
qui log close LL
clear
qui log using "`W'/rep_b01_replay.txt", replace text name(RR)
cap noi REPIT
dis "B01_replay_rc=" _rc
qui log close RR

* (c) estimates store / restore, then replay
qui use "`W'/d5.dta", clear
global CMDBODY y i t ttre, att att_it robust espre(3) espost(3) jointtests
qui log using "`W'/rep_c01_live.txt", replace text name(LL)
cap noi RUNIT
qui log close LL
cap noi estimates store WD1
dis "C01_store_rc=" _rc
* an unrelated estimation in between
qui reg y x1 x2
dis "C01_ecmd_after_reg=[" e(cmd) "]"
* restore FIRST, then replay -- the clean path
cap noi estimates restore WD1
dis "C01_ecmd_after_restore=[" e(cmd) "]"
dis "C01_LASTRC_before_replay=[$WOOLDID_LASTRC]"
qui log using "`W'/rep_c01_replay.txt", replace text name(RR)
cap noi REPIT
dis "C01_replay_rc=" _rc
qui log close RR

* (c2) THE TRAP: type a bare `wooldid' while some OTHER command owns e().
*      What does the user see, and what does it do to the replay marker?
qui use "`W'/d5.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
dis "C02_setup_rc=" _rc
cap noi estimates store WD2
qui reg y x1 x2
dis "C02_ecmd=[" e(cmd) "]"
dis "C02_LASTRC_before=[$WOOLDID_LASTRC]"
cap noi wooldid
dis "C02_bare_after_reg_rc=" _rc
dis "C02_LASTRC_after=[$WOOLDID_LASTRC]"
cap noi estimates restore WD2
dis "C02_ecmd_after_restore=[" e(cmd) "]"
cap noi wooldid
dis "C02_replay_after_restore_rc=" _rc

* (c3) two stored wooldid models: does restore+replay show the RIGHT one?
qui use "`W'/d5.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
cap noi estimates store M1
qui log using "`W'/rep_c03_live.txt", replace text name(LL)
cap noi REPIT
qui log close LL
cap noi wooldid y i t ttre, att_it cluster(i) espre(3) espost(3)
dis "C03_second_rc=" _rc
cap noi estimates restore M1
qui log using "`W'/rep_c03_replay.txt", replace text name(RR)
cap noi REPIT
dis "C03_replay_rc=" _rc
qui log close RR

* (e) replay idempotence: replay twice, diff the two replays
qui use "`W'/d5.dta", clear
global CMDBODY y i t ttre, att robust espre(2) espost(2) jointtests
cap noi RUNIT
qui log using "`W'/rep_e01_live.txt", replace text name(LL)
cap noi REPIT
qui log close LL
qui log using "`W'/rep_e01_replay.txt", replace text name(RR)
cap noi REPIT
qui log close RR

* (f) what a FAILED run leaves behind, in detail
qui use "`W'/d5.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
dis "F_good_ecmd=[" e(cmd) "] N=" e(N)
* a failure that gets past parse and dies inside estimation
cap noi wooldid y i t ttre, att robust espre(2) espost(2) fe(nosuchvar)
dis "F_fail_rc=" _rc
dis "F_LASTRC=[$WOOLDID_LASTRC]"
dis "F_ecmd_now=[" e(cmd) "]  e(N)=" e(N)
cap noi REPIT
dis "F_replay_refused_rc=" _rc
macro drop WOOLDID_LASTRC
qui log using "`W'/rep_f01_defeated.txt", replace text name(RR)
cap noi REPIT
dis "F_replay_defeated_rc=" _rc
qui log close RR

* a failure DEEP in estimation: singleton/collinear blow-up
qui use "`W'/d5.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2) controls(x1 x1)
dis "F2_rc=" _rc
dis "F2_ecmd=[" e(cmd) "] e(cmd2)=[" e(cmd2) "]"
dis "F2_LASTRC=[$WOOLDID_LASTRC]"

dis "P5_DONE"
