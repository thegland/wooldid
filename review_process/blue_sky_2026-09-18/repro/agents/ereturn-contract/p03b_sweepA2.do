clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioA2.txt"
cap erase "vioA2.txt"
include "mkdata.doi"
include "ckc.doi"
include "ckc2.doi"
dis "DATA_REUSED"

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

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"
use D1.dta, clear
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

dis "SWEEPA2_DONE"
