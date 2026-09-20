clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* usage: -e do mc_tstat.do <chunk> <ndraws> <reps>
local CH `1'
local R  `2'
local REPS `3'
if "`CH'"=="" local CH 1
if "`R'"==""  local R 5
if "`REPS'"=="" local REPS 99

local FN 24
local FT 8
local FSEED 11
include "$WD/mkfix.doi"
replace y = base
qui wooldid y cg t ttre, att cluster(unit) esprelength(1) espostlength(1) level(95)
mat TB = e(b)
local tn : colnames TB

tempname fh
file open `fh' using "$WD/mctstat_`CH'.txt", write replace text
file write `fh' "draw|param|arm|b|se|lb|ub|truth|bin" _n
forvalues c = 1/`=colsof(TB)' {
  local nm : word `c' of `tn'
  file write `fh' "0|`nm'|truth|" (TB[1,`c']) "|.|.|.|" (TB[1,`c']) "|." _n
}
forvalues r = 1/`R' {
  local dr = (`CH'-1)*1000 + `r'
  local FSEED = 800000 + `dr'
  local FN 24
  local FT 8
  include "$WD/mkfix.doi"
  cap qui wooldid y cg t ttre, att cluster(unit) esprelength(1) espostlength(1) ///
       level(95) bootstrap(standard) reps(`REPS') rseed(`=900000+`dr'') boottstat dots(0)
  if _rc != 0 {
    file write `fh' "`dr'|RC|rc|" (_rc) "|.|.|.|.|." _n
    continue
  }
  mat Q = e(wdidmainresults)
  local rn : rownames Q
  local re : roweq Q
  forvalues i = 1/`=rowsof(Q)' {
    local rr : word `i' of `rn'
    local ee : word `i' of `re'
    local bb = Q[`i',1]
    local lo = Q[`i',5]
    local hi = Q[`i',6]
    local inb = cond(`lo'<. & `hi'<. & `bb'<., cond(`bb'>=`lo' & `bb'<=`hi',1,0), .)
    file write `fh' "`dr'|`ee':`rr'|tstud|" (`bb') "|" (Q[`i',2]) "|" (`lo') "|" (`hi') "|.|" (`inb') _n
  }
  dis "MCT|`dr'|done"
}
file close `fh'
dis "MCTSTAT_DONE_`CH'"
