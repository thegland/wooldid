clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* ---- MC arm 1: analytic intervals, i=cg, cluster(unit) finer than i ----
* writes long-format rows:  draw|arm|param|b|se|truth
local R   500
local FN  40
local FT  8
local OPT "att att_i att_it att_itime"

* ---- truth on the noiseless outcome ----
local FSEED 11
include "$WD/mkfix.doi"
replace y = base
qui wooldid y cg t ttre, `OPT' cluster(unit) esprelength(2) espostlength(2) level(95)
mat TB = e(b)
local tn : colnames TB
local te : coleq TB
dis "TRUTHNAMES|`te'"
dis "TRUTHNAMES2|`tn'"
mat list TB, format(%14.9f)

tempname fh
file open `fh' using "$WD/mc_analytic.txt", write replace text
file write `fh' "draw|arm|param|b|se" _n
forvalues c = 1/`=colsof(TB)' {
  local nm : word `c' of `tn'
  local eq : word `c' of `te'
  file write `fh' "0|truth|`eq':`nm'|" (TB[1,`c']) "|." _n
}

forvalues r = 1/`R' {
  * regenerate fixture with a per-draw noise seed
  local FSEED = 100000 + `r'
  local FN 40
  local FT 8
  include "$WD/mkfix.doi"
  foreach arm in clus robu {
    if "`arm'" == "clus" local vceopt cluster(unit)
    if "`arm'" == "robu" local vceopt robust
    cap qui wooldid y cg t ttre, `OPT' `vceopt' esprelength(2) espostlength(2) level(95)
    if _rc == 0 {
      mat BB = e(b)
      mat VV = e(V)
      local bn : colnames BB
      local be : coleq BB
      forvalues c = 1/`=colsof(BB)' {
        local nm : word `c' of `bn'
        local eq : word `c' of `be'
        local sv = sqrt(VV[`c',`c'])
        file write `fh' "`r'|`arm'|`eq':`nm'|" (BB[1,`c']) "|" (`sv') _n
      }
      if `r' == 1 {
        dis "ARM|`arm'|G=" e(N_clust) "|df_r=" e(df_r) "|K=" colsof(BB)
      }
    }
    else {
      file write `fh' "`r'|`arm'|RC|" (_rc) "|." _n
    }
  }
  if mod(`r',50)==0 dis "MCPROG|`r'"
}
file close `fh'
dis "MC_ANALYTIC_DONE"
