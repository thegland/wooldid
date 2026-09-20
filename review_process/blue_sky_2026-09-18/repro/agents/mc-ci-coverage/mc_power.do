clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* harness power check: analytic SEs with cluster() AT i (=cg, 4 clusters) are documented
* to UNDER-cover.  If my harness cannot see that, a clean result elsewhere means nothing.
local R 250
local FN 40
local FT 8
local FSEED 11
include "$WD/mkfix.doi"
replace y = base
qui wooldid y cg t ttre, att att_i cluster(cg) esprelength(2) espostlength(2) level(95)
mat TB = e(b)
local tn : colnames TB
tempname fh
file open `fh' using "$WD/mc_power.txt", write replace text
file write `fh' "draw|arm|param|b|se" _n
forvalues c = 1/`=colsof(TB)' {
  local nm : word `c' of `tn'
  file write `fh' "0|truth|`nm'|" (TB[1,`c']) "|." _n
}
forvalues r = 1/`R' {
  local FSEED = 300000 + `r'
  local FN 40
  local FT 8
  include "$WD/mkfix.doi"
  cap qui wooldid y cg t ttre, att att_i cluster(cg) esprelength(2) espostlength(2) level(95)
  if _rc==0 {
    mat BB = e(b)
    mat VV = e(V)
    local bn : colnames BB
    forvalues c = 1/`=colsof(BB)' {
      local nm : word `c' of `bn'
      file write `fh' "`r'|clusI|`nm'|" (BB[1,`c']) "|" (sqrt(VV[`c',`c'])) _n
    }
    if `r'==1 dis "POWARM|G=" e(N_clust) "|df_r=" e(df_r)
  }
}
file close `fh'
dis "MC_POWER_DONE"
