*==============================================================================
* p7_fixups.do -- re-runs of the configs whose option syntax was wrong in
*                 p2/p3 (contreat needs contreatcontrols(); pretrendpower()
*                 needs esfixedbaseperiod), plus the extra replay probes.
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
local dT 14
local donsets "2005 2008 2010"
local dantic 2
include "`W'/mkdata.do"
save "`W'/d7.dta", replace
include "`W'/dumpe.do"

cap prog drop runpair
prog def runpair
  syntax , tag(string) k(integer) opts(string) w(string) [dta(string)]
  if "`dta'" == "" local dta d7
  dis as result "==== CONFIG `tag'  k=`k'  opts: `opts' ===="
  qui use "`w'/`dta'.dta", clear
  cap noi wooldid `opts' anticipation(`k')
  local rcA = _rc
  dis "A_RC_`tag'=`rcA'"
  if `rcA' == 0 dumpe , f("`w'/out_`tag'_A.txt")
  qui use "`w'/`dta'.dta", clear
  qui replace ttre = ttre - `k'
  cap noi wooldid `opts'
  local rcB = _rc
  dis "B_RC_`tag'=`rcB'"
  if `rcB' == 0 dumpe , f("`w'/out_`tag'_B.txt")
  dis "PAIR_`tag'_RCS=`rcA',`rcB'"
end

runpair , tag(f06) k(2) w("`W'") opts(y i t ttre, att robust contreat(dose) contreatcontrols(basic) espre(2) espost(2))
runpair , tag(f06b) k(2) w("`W'") opts(y i t ttre, att robust contreat(dose) contreatcontrols(byTreated) lattice(all 3) espre(2) espost(2) makecfxplots)
runpair , tag(f14) k(2) w("`W'") opts(y i t ttre, att robust espre(3) espost(3) esfixedbaseperiod pretrendpower(80 50))
runpair , tag(f14b) k(2) w("`W'") opts(y i t ttre, att att_it robust espre(3) espost(3) esfixedbaseperiod esrelativeto(-2) pretrendpower(80))

* ---------------------------------------------------------------- replay ----
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
  dis as result "==== REPLAY CONFIG `tag' : wooldid `body' ===="
  qui use "`w'/d7.dta", clear
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

dopair , tag(r03) w("`W'") body(y i t ttre, att robust contreat(dose) contreatcontrols(basic) espre(2) espost(2))
dopair , tag(r05) w("`W'") body(y i t ttre, att att_it robust espre(3) espost(3) esfixedbaseperiod pretrendpower(50 80))
dopair , tag(r05b) w("`W'") body(y i t ttre, att att_it robust espre(3) espost(3) esfixedbaseperiod pretrendpower(50 80) level(92.5))
dopair , tag(r09) w("`W'") body(y i t ttre, att cluster(i) espre(3) espost(3) esfixedbaseperiod bootstrap(standard) reps(25) rseed(77) dots(0) jointtests pretrendpower(80))
dopair , tag(r13) w("`W'") body(y i t ttre, att robust espre(3) espost(3) esfixedbaseperiod pretrendpower(80) level(80))
dopair , tag(r14) w("`W'") body(y i t ttre, att robust contreat(dose) contreatcontrols(byTreated) lattice(all 3) espre(2) espost(2) makecfxplots jointtests)

dis "P7_DONE"

* ---------------------------------------------------- impvalues at level() --
qui use "`W'/d7.dta", clear
cap noi wooldid y i t ttre, att_it cluster(i) espre(2) espost(2) impvalues level(80)
dis "IMP80_RC=" _rc
dis "IMP80_level=" e(level) "  df_r=" e(df_r) "  Nclus=" e(Nclus)
cap noi {
  tempname M
  mat `M' = e(wdidmainresults)
  local cn : colnames `M'
  dis "IMP80 colnames: `cn'"
  local crit = invttail(e(df_r), .10)
  dis "IMP80 crit(df=" e(df_r) ")=" %12.10f `crit'
  local worst = 0
  forvalues r = 1/`=rowsof(`M')' {
    local b = `M'[`r',1]
    local s = `M'[`r',2]
    local lo = `M'[`r',5]
    if `b' < . & `s' < . & `lo' < . {
      local e1 = `b' - `crit'*`s'
      local d = abs(`lo'-`e1')/max(abs(`e1'),1e-12)
      if `d' > `worst' local worst = `d'
    }
  }
  dis "IMP80_worst_rel_dev=" %12.4e `worst'
  mat list `M', format(%12.8f)
}
dis "IMPCHK_RC=" _rc
dis "P7_ALLDONE"
