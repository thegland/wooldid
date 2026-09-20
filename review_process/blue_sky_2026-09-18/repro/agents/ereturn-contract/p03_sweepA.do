clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioA.txt"
cap erase "vioA.txt"
include "mkdata.doi"
include "ckc.doi"
include "ckc2.doi"
mkD1
mkD2
mkD3
mkD4
dis "DATA_OK"

* record the full set of e() names each config posts, for the documented-vs-present diff
cap prog drop nameslog
prog def nameslog
  args tag
  tempname fh
  file open `fh' using "names.txt", write text append
  foreach s in `: e(scalars)' {
    file write `fh' "N|`tag'|SCA|`s'" _n
  }
  foreach s in `: e(macros)' {
    file write `fh' "N|`tag'|MAC|`s'" _n
  }
  foreach s in `: e(matrices)' {
    file write `fh' "N|`tag'|MAT|`s'" _n
  }
  file close `fh'
end

cap prog drop RUN
prog def RUN
  args tag rc
  dis "==================== `tag' (rc=`rc') ===================="
  if "`rc'" != "0" {
    ckcw INF "`tag'" SKIPPED "run failed rc=`rc'; e() is STALE, assertions skipped"
    dis "SKIPPED|`tag'|rc=`rc'"
    exit
  }
  cap noi ckc "`tag'"
  dis "CKC_RC=" _rc
  cap noi ckc2 "`tag'"
  dis "CKC2_RC=" _rc
  nameslog "`tag'"
end

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"

* =================== D1: unit panel ===================
use D1.dta, clear
cap noi wooldid y id t ttre, att
local lrc = _rc
dis "A01_RC=`lrc'"
RUN A01_att `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime
local lrc = _rc
dis "A02_RC=`lrc'"
RUN A02_all4 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime esprelength(3) espostlength(3)
local lrc = _rc
dis "A03_RC=`lrc'"
RUN A03_es `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES'
local lrc = _rc
dis "A04_RC=`lrc'"
RUN A04_esfix `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu)
local lrc = _rc
dis "A05_RC=`lrc'"
RUN A05_clu `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(id)
local lrc = _rc
dis "A06_RC=`lrc'"
RUN A06_cluid `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' robust
local lrc = _rc
dis "A07_RC=`lrc'"
RUN A07_robust `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' controls(xctl)
local lrc = _rc
dis "A08_RC=`lrc'"
RUN A08_ctl `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' [aw=awt]
local lrc = _rc
dis "A09_RC=`lrc'"
RUN A09_aw `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' [pw=pwt]
local lrc = _rc
dis "A10_RC=`lrc'"
RUN A10_pw `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' level(90)
local lrc = _rc
dis "A11_RC=`lrc'"
RUN A11_lvl90 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' level(97.5)
local lrc = _rc
dis "A12_RC=`lrc'"
RUN A12_lvl975 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' level(99.9)
local lrc = _rc
dis "A13_RC=`lrc'"
RUN A13_lvl999 `lrc'

cap noi wooldid ypos id t ttre, att att_i att_it att_itime `ES' poisson
local lrc = _rc
dis "A14_RC=`lrc'"
RUN A14_poisson `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' unconditionalse
local lrc = _rc
dis "A15_RC=`lrc'"
RUN A15_uncond `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' customvce(cluster clu)
local lrc = _rc
dis "A16_RC=`lrc'"
RUN A16_cvceclu `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' customvce(robust)
local lrc = _rc
dis "A17_RC=`lrc'"
RUN A17_cvcerob `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' fastnose
local lrc = _rc
dis "A18_RC=`lrc'"
RUN A18_fastnose `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' impvalues cluster(clu)
local lrc = _rc
dis "A19_RC=`lrc'"
RUN A19_impvalues `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' gteffects
local lrc = _rc
dis "A20_RC=`lrc'"
RUN A20_gteffects `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' jointtests adversarial
local lrc = _rc
dis "A21_RC=`lrc'"
RUN A21_jointadv `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sg) cluster(clu)
local lrc = _rc
dis "A22_RC=`lrc'"
RUN A22_sg2 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' subgroups(sg3) cluster(clu)
local lrc = _rc
dis "A23_RC=`lrc'"
RUN A23_sg3 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' anticipation(1)
local lrc = _rc
dis "A24_RC=`lrc'"
RUN A24_antic `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) esrelativeto(-2)
local lrc = _rc
dis "A25_RC=`lrc'"
RUN A25_relto2 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' timetrends
local lrc = _rc
dis "A26_RC=`lrc'"
RUN A26_timetrends `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' simultaneousbands
local lrc = _rc
dis "A27_RC=`lrc'"
RUN A27_simbands `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' pretrendpower(0.5)
local lrc = _rc
dis "A28_RC=`lrc'"
RUN A28_ptp `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' oldsyntax
local lrc = _rc
dis "A29_RC=`lrc'"
RUN A29_oldsyntax `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' interactivecontrols(xctl)
local lrc = _rc
dis "A30_RC=`lrc'"
RUN A30_interctl `lrc'

cap noi wooldid y id t ttre, att_itime
local lrc = _rc
dis "A31_RC=`lrc'"
RUN A31_itimeonly `lrc'

cap noi wooldid y id t ttre, att_i att_itime `ES'
local lrc = _rc
dis "A32_RC=`lrc'"
RUN A32_subsetest `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' summarystats histogramcohorteffects
local lrc = _rc
dis "A33_RC=`lrc'"
RUN A33_summ `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) subgroups(sg) jointtests adversarial
local lrc = _rc
dis "A34_RC=`lrc'"
RUN A34_sgjoint `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' nofastcollapse
local lrc = _rc
dis "A35_RC=`lrc'"
RUN A35_nofastcol `lrc'

dis "SWEEPA_DONE"
