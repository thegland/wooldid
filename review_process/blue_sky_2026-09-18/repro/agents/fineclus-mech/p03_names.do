clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use "p01_A.dta", clear
qui wooldid y firm year onset, att_i cluster(ind)
mat b = e(b)
local nm : colnames b
local eq : coleq b
dis "COLNAMES: `nm'"
dis "COLEQ: `eq'"
mat list b
ereturn list
dis "ALLDONE"
