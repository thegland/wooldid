clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* ===== simultaneous (sup-t) band: does it cover EVERY horizon at once at level()? =====
local R 400
local FN 40
local FT 10
local FSEED 11
include "$WD/mkfix.doi"
replace y = base
qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) level(95) simultaneousbands
mat TB = e(b_es_att)
local tn : colnames TB
dis "BANDTRUTH|`tn'"
mat list TB, format(%14.9f)

tempname fh
file open `fh' using "$WD/mc_band.txt", write replace text
file write `fh' "draw|lvl|horizon|lb|ub|truth|crit" _n
forvalues c = 1/`=colsof(TB)' {
  local nm : word `c' of `tn'
  file write `fh' "0|.|`nm'|.|.|" (TB[1,`c']) "|." _n
}
forvalues r = 1/`R' {
  local FSEED = 700000 + `r'
  local FN 40
  local FT 10
  include "$WD/mkfix.doi"
  foreach L in 90 95 {
    cap qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) level(`L') simultaneousbands
    if _rc==0 {
      cap confirm matrix e(esband_att)
      if _rc==0 {
        mat BB = e(esband_att)
        local bc : colnames BB
        local cr = e(esband_crit_att)
        forvalues h = 1/`=colsof(BB)' {
          local hn : word `h' of `bc'
          file write `fh' "`r'|`L'|`hn'|" (BB[1,`h']) "|" (BB[2,`h']) "|.|" (`cr') _n
        }
      }
    }
  }
  if mod(`r',100)==0 dis "BANDPROG|`r'"
}
file close `fh'
dis "MC_BAND_DONE"
