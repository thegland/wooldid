* p09_options2.do -- the option arms p04 could not reach (my earlier specification
* errors), plus the "interactivecontrols with a degenerate slope" case the help says
* reduces exactly to controls(i.t#c.x) and where the cells should therefore be EXACT.
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
gen double ypos = exp(y/4)
* a strictly unit-level control: with i = unit it cannot vary within a cell
bys unit (year): gen double xu = rnormal() if _n==1
bys unit (year): replace xu = xu[1]
tempfile base
qui save `base'

frun H1_sg_iunit   , cmd(wooldid y unit year ttre, `EST' gteffects `ESW' subgroups(sg))
use `base', clear
frun H2_poisson    , cmd(wooldid ypos cohp year ttre, `EST' gteffects `ESW' poisson)
use `base', clear
frun H3_poisexp    , cmd(wooldid ypos cohp year ttre, `EST' gteffects `ESW' poisson poisexpresults)
use `base', clear
frun H4_uncond     , cmd(wooldid y cohp year ttre, att gteffects `ESW' unconditionalse)
use `base', clear
frun H5_customvce  , cmd(wooldid y cohp year ttre, att gteffects `ESW' customvce(cluster st))
use `base', clear
frun H6_ic_unitlev , cmd(wooldid y unit year ttre, `EST' gteffects `ESW' interactivecontrols(xu))
use `base', clear
frun H7_ccc        , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' coarsecohortcontrols(c.x1))
use `base', clear
frun H8_ant_oos    , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' anticipation(2) oostreatedcontrols(3))
use `base', clear
frun H9_tt_clus    , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' timetrends cluster(st))
use `base', clear
frun H10_attonly   , cmd(wooldid y cohp year ttre, att gteffects `ESW')
use `base', clear
frun H11_icell_coh , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' interactivecontrols(x1))
use `base', clear
frun H12_ctrl_cell , cmd(wooldid y cohp year ttre, `EST' gteffects `ESW' controls(i.year##c.x1))

dis "DONE_P09"
