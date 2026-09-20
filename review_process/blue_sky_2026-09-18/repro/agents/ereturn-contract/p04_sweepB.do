clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioB.txt"
cap erase "vioB.txt"
include "ckc.doi"
include "ckc2.doi"
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
cap prog drop PERMX
prog def PERMX
  args tag
  if "`e(perm_method)'" == "" exit
  tempname M P
  cap confirm matrix e(wdidmainresults)
  if _rc exit
  matrix `M' = e(wdidmainresults)
  local rn : rowfullnames `M'
  local cn : colnames `M'
  local pc : list posof "p" in cn
  cap confirm matrix e(perm_p)
  local hasP = (_rc==0)
  if `hasP' matrix `P' = e(perm_p)
  local pn ""
  if `hasP' local pn : colnames `P'
  foreach est in `e(wdid_estimands)' {
    local sc = .
    cap local sc = e(perm_p_`est')
    local trow = 0
    local i = 0
    foreach r of local rn {
      local ++i
      if `"`r'"' == "main:`est'" local trow = `i'
    }
    local tval = .
    if `trow' > 0 & `pc' > 0 local tval = `M'[`trow',`pc']
    local vval = .
    if `hasP' {
      local q : list posof "main_`est'" in pn
      if `q' > 0 local vval = `P'[1,`q']
    }
    dis "PERMX|`tag'|`est'|scalar=`sc'|table=`tval'|vector=`vval'"
    if !missing(`sc') & !missing(`tval') {
      if abs(`sc'-`tval') > 1e-12 ckcw VIO "`tag'" PERMPSCA "e(perm_p_`est')=`sc' but e(wdidmainresults)[main:`est',p]=`tval'"
    }
    if !missing(`vval') & !missing(`tval') {
      if abs(`vval'-`tval') > 1e-12 ckcw VIO "`tag'" PERMPVEC "e(perm_p)[main_`est']=`vval' but table p=`tval'"
    }
  }
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
  cap noi PERMX "`tag'"
  dis "PERMX_RC=" _rc
  nameslog "`tag'"
end

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"

* =================== D2: state-year administrative panel ===================
use D2.dta, clear
cap noi wooldid y st t ttre, att att_i att_it att_itime `ES' cluster(st)
local lrc = _rc
dis "B01_RC=`lrc'"
RUN B01_d2 `lrc'

cap noi wooldid y st t ttre, att att_i att_it att_itime `ES' cluster(st) subgroups(reg)
local lrc = _rc
dis "B02_RC=`lrc'"
RUN B02_d2sg `lrc'

cap noi wooldid y st t ttre, att att_i att_it att_itime `ES' cluster(clu) subgroups(reg) jointtests adversarial
local lrc = _rc
dis "B03_RC=`lrc'"
RUN B03_d2sgjoint `lrc'

* =================== D3: repeated cross-section, genuine continuous dose ===================
use D3.dta, clear
cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic)
local lrc = _rc
dis "B04_RC=`lrc'"
RUN B04_ct_basic `lrc'

cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(untreatedZero)
local lrc = _rc
dis "B05_RC=`lrc'"
RUN B05_ct_uz `lrc'

cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) cluster(clu)
local lrc = _rc
dis "B06_RC=`lrc'"
RUN B06_ct_clu `lrc'

cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) subgroups(sg)
local lrc = _rc
dis "B07_RC=`lrc'"
RUN B07_ct_sg `lrc'

cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort)
local lrc = _rc
dis "B08_RC=`lrc'"
RUN B08_ct_bycoh `lrc'

cap noi wooldid y cty t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) lattice(all 3)
local lrc = _rc
dis "B09_RC=`lrc'"
RUN B09_ct_lattice `lrc'

cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) semielasticity
local lrc = _rc
dis "B10_RC=`lrc'"
RUN B10_ct_semi `lrc'

* =================== bootstrap ===================
use D1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(40) rseed(777)
local lrc = _rc
dis "B11_RC=`lrc'"
RUN B11_boot40 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(8) rseed(777)
local lrc = _rc
dis "B12_RC=`lrc'"
RUN B12_boot8 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) subgroups(sg) bootstrap(standard) reps(40) rseed(777)
local lrc = _rc
dis "B13_RC=`lrc'"
RUN B13_bootsg40 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(40) rseed(777) boottstat
local lrc = _rc
dis "B14_RC=`lrc'"
RUN B14_boottstat `lrc'

cap noi wooldid y id t ttre, att cluster(clu) bootstrap(bayes) reps(40) rseed(777)
local lrc = _rc
dis "B15_RC=`lrc'"
RUN B15_bootbayes `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) bootstrap(standard) reps(40) rseed(777) simultaneousbands
local lrc = _rc
dis "B16_RC=`lrc'"
RUN B16_bootbands `lrc'

* =================== permute ===================
cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) permute(shuffle) permreps(40) rseed(777)
local lrc = _rc
dis "B17_RC=`lrc'"
RUN B17_perm40 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) subgroups(sg) permute(shuffle) permreps(40) rseed(777)
local lrc = _rc
dis "B18_RC=`lrc'"
RUN B18_permsg `lrc'

cap noi wooldid y id t ttre, att_itime cluster(clu) permute(shuffle) permreps(40) rseed(777)
local lrc = _rc
dis "B19_RC=`lrc'"
RUN B19_permitime `lrc'

cap noi wooldid y id t ttre, att_i att_itime cluster(clu) permute(shuffle) permreps(40) rseed(777)
local lrc = _rc
dis "B20_RC=`lrc'"
RUN B20_permsubset `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime `ES' cluster(clu) permute(shuffle) permreps(40) rseed(777) permanalyticse permtstat
local lrc = _rc
dis "B21_RC=`lrc'"
RUN B21_permanase `lrc'

dis "SWEEPB_DONE"
