*==============================================================================
* p15: known_bydesign section 7, item 11 - "Check whether fastnose actually posts
* e(b_es_att)/e(numpre_att)".  sthlp:1325 still opens the "the covariance is not
* returned but the coefficient vector still is" list with {cmd:fastnose}.
* Run WITHOUT permute(), so the permute-specific withholding is not in play.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear

cap prog drop ichk
prog def ichk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) `opts'
  dis "I11|`tag'|RC=" _rc
  if _rc != 0 exit
  foreach m in b_es_att V_es_att numpre_att b_es_att_i V_es_att_i {
    cap confirm matrix e(`m')
    local hm = (_rc==0)
    if `hm' == 0 {
      cap confirm scalar e(`m')
      local hm = 2*(_rc==0)
    }
    dis "I11|`tag'|`m'=" `hm' "   (0=absent, 1=matrix, 2=scalar)"
  }
  dis "I11|`tag'|numpre_att value=" e(numpre_att) "|vcvstate=`e(wdid_vcvstate)'"
end
ichk plain
ichk fastnose  , opts(fastnose)
ichk impvalues , opts(impvalues)
dis "P15_DONE"
