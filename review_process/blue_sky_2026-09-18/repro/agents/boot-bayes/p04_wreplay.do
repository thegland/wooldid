* p04: does the Dirichlet multiplier MULTIPLY the user weight (not replace it)?
*      Exact replay: rebuild manw = pop * w_c from bootsave() and re-run the
*      replication's own command by hand; compare to e(bootdist) bitwise.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a571-42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

use p01_panel, clear
gen double bigpop = pop * cond(st==3, 500, 1)
save p04_panel, replace

foreach am in equal size {
  use p04_panel, clear
  cap noi wooldid y st year onset [aw=bigpop], att cluster(st) bootstrap(bayes) bootalpha(`am') reps(8) rseed(999) bootsave(p04_`am') replace dots(0)
  dis "RUN_`am'_RC=" _rc
  matrix BD = e(bootdist)
  matrix OB = e(b)
  dis "OBS_b_`am' = " %21x OB[1,1]

  use p04_`am', clear
  qui ds w_*
  local wv `r(varlist)'
  qui keep _rep `wv'
  qui reshape long w_, i(_rep) j(st)
  rename w_ wdraw
  tempfile draws
  qui save `draws'

  forval r = 1/3 {
    preserve
      use `draws', clear
      qui keep if _rep == `r'
      qui drop _rep
      tempfile dr`r'
      qui save `dr`r''
    restore
    use p04_panel, clear
    qui merge m:1 st using `dr`r'', keep(match) nogen
    qui gen double manw = bigpop * wdraw
    cap noi wooldid y st year onset [aw=manw], att cluster(st) fastnose esrelativeto(-1)
    local mrc = _rc
    matrix MB = e(b)
    dis "`am' rep`r' MANUAL_RC=`mrc' manual_b=" %21x MB[1,1] " boot_b=" %21x BD[`r',1] " DIFF=" %21x (MB[1,1]-BD[`r',1])
  }
  * negative control: multiplier REPLACING the user weight
  use p04_panel, clear
  qui merge m:1 st using `dr1', keep(match) nogen
  cap noi wooldid y st year onset [aw=wdraw], att cluster(st) fastnose esrelativeto(-1)
  matrix RB = e(b)
  dis "`am' NEGCTRL replace-weight b=" %21x RB[1,1] " vs boot " %21x BD[1,1]
}
dis "PROBE_DONE"
