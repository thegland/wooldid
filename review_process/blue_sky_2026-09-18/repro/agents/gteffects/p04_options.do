* p04_options.do -- aggregation identity under option combinations a real user would use.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime
local ESW esprelength(3) espostlength(3)

cap prog drop frun
prog def frun
  syntax anything(name=tg) , CMD(string asis)
  ereturn clear
  cap noi `cmd'
  dis "`tg'_RC=" _rc
  if _rc != 0 exit
  gtcheck , tag(`tg') pre es
  gtdefcheck , tag(`tg')
end

mkpanel , n(80) t(10) onsets("2004 2006 2008") ysd(0.4) seed(4242)
gen byte st = mod(unit, 8) + 1
gen byte sg = mod(unit, 2)
bys unit (year): gen double wt = 0.3 + 2*runiform() if _n==1
bys unit (year): replace wt = wt[1]
* customweightmultiplier that VARIES WITHIN a cohort-by-period cell
gen double cwmv = 0.5 + 3*runiform()
* customweightmultiplier that is CONSTANT within a cohort-by-period cell (cohort-year level)
bys cohp year (unit): gen double cwmc = 0.5 + 3*runiform() if _n==1
bys cohp year (unit): replace cwmc = cwmc[1]
tempfile base
qui save `base'

frun F1_subgroups  , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' subgroups(sg))
use `base', clear
frun F2_fe         , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' fe(unit))
use `base', clear
frun F3_controls   , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' controls(x1))
use `base', clear
frun F4_intercon   , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' interactivecontrols(x1))
use `base', clear
frun F5_poisson    , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' poisson)
use `base', clear
frun F6_poisexp    , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' poisson poisexpresults)
use `base', clear
frun F7_semielas   , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' semielasticity)
use `base', clear
frun F8_impv       , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' impvalues cluster(st))
use `base', clear
frun F9_uncond     , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' unconditionalse)
use `base', clear
frun F10_clus      , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' cluster(st))
use `base', clear
frun F11_pw        , cmd(wooldid y cohp year ttre [pw=wt], `EST' gteffects `ESW')
use `base', clear
frun F12_esrelto   , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' esrelativeto(-2) esfixedbaseperiod)
use `base', clear
frun F13_oos       , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' oostreatedcontrols(2))
use `base', clear
frun F14_shortwin  , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(1) espostlength(1))
use `base', clear
frun F15_cwm_within, cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' customweightmultiplier(cwmv))
use `base', clear
frun F16_cwm_cell  , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' customweightmultiplier(cwmc))
use `base', clear
frun F17_aw_cwm    , cmd(wooldid y cohp year ttre [aw=wt], `EST' gteffects `ESW' customweightmultiplier(cwmv))
use `base', clear
frun F18_espost0   , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(0))
use `base', clear
frun F19_ccc       , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' coarsecohortcontrols(x1))
use `base', clear
frun F20_jt        , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' jointtests adversarial esfixedbaseperiod)

dis "DONE_P04"
