* p03_oracle.do -- INDEPENDENT oracle for the e(gteffects) estimate/se/p/ci columns.
* The stage-1 regression is stored by `verbose' as ___wooldidfullmodel; we restore it and
* recompute every cell's b/se/p/CI from its own e(b)/e(V)/e(df_r).  Any divergence is a
* wrong number in the printed cell table, not just an internal inconsistency.
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

cap prog drop orun
prog def orun
  syntax anything(name=tg) , CMD(string asis) [LVL(real 95)]
  ereturn clear
  cap noi `cmd'
  dis "`tg'_RC=" _rc
  if _rc != 0 exit
  cap confirm matrix e(gteffects)
  if _rc != 0 {
    dis as result "OR|`tg'|NOTABLE"
    exit
  }
  tempname GTS
  matrix `GTS' = e(gteffects)
  gtcheck , tag(`tg') pre es
  qui est restore ___wooldidfullmodel
  dis as result "OR|`tg'|REGCMD|`e(cmd)'|VCE|`e(vce)'|NCLUST|`e(N_clust)'"
  cap noi mata: gtoracle("`tg'", "`GTS'", `lvl')
  dis "ORACLE_RC=" _rc
end

mkpanel , n(80) t(10) onsets("2004 2006 2008") ysd(0.4) seed(4242)
gen byte st = mod(unit, 8) + 1
gen double ypos = exp(y/4)
tempfile base
qui save `base'

orun O1_iid       , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose)
use `base', clear
orun O2_robust    , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose robust)
use `base', clear
orun O3_clusunit  , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose cluster(unit))
use `base', clear
orun O4_clusstate , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose cluster(st))
use `base', clear
orun O5_level90   , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose level(90)) lvl(90)
use `base', clear
orun O6_fastnose  , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose fastnose)
use `base', clear
orun O7_uncond    , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose unconditionalse)
use `base', clear
orun O8_impv      , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose impvalues cluster(st))
use `base', clear
orun O9_poisson   , cmd(wooldid ypos cohp year ttre, `EST' gteffects `ESW' verbose poisson)
use `base', clear
orun O10_iunit    , cmd(wooldid y unit year ttre, `EST' gteffects `ESW' verbose cluster(st))
use `base', clear
orun O11_level975 , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' verbose level(97.5)) lvl(97.5)

dis "DONE_P03"
