* p11_window.do -- does the inwindow column name exactly the horizons the run reports?
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

cap prog drop wrun
prog def wrun
  syntax anything(name=tg) , CMD(string asis)
  ereturn clear
  cap noi `cmd'
  dis "`tg'_RC=" _rc
  if _rc != 0 exit
  gtcheck , tag(`tg') pre es
  gtwincheck , tag(`tg')
end

mkpanel , n(80) t(12) onsets("2003 2005") ysd(0.4) seed(9090)
tempfile base
qui save `base'

wrun W1_default  , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2))
use `base', clear
wrun W2_esf      , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) esfixedbaseperiod)
use `base', clear
wrun W3_esf_rel2 , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2))
use `base', clear
wrun W4_post0    , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(0))
use `base', clear
wrun W5_pre0     , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(0) espostlength(3))
use `base', clear
wrun W6_wide     , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(6) espostlength(8))
use `base', clear
wrun W7_antic    , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) anticipation(2))
use `base', clear
wrun W8_esfantic , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(2) espostlength(2) anticipation(2) esfixedbaseperiod)

mkpanel , n(80) t(9) onsets("2004 2008") ysd(0.4) tstep(2) t0(2000) seed(9091)
wrun W9_biennial , cmd(wooldid y cohp year ttre, `EST' gteffects esprelength(4) espostlength(4))

dis "DONE_P11"
