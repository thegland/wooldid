* p11 -- does the DISPLAYED e(pwsgs) row key resolve against the STORED one?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid
local gN 240
local gT 9
local gonsets "2003 2004 2006"
local gysd 0.60
include "`GEN'/gen_core.do"
egen long grp = group(cohort region)
gen int sgq = 1
qui replace sgq = 3  if region == 2
qui replace sgq = 7  if region == 3
qui replace sgq = 99 if region == 4
qui wooldid y grp t ttre , att cluster(unit) fe(i.unit) subgroups(sgq) esprelength(2) espostlength(2)
matrix PWX = e(pwsgs)
dis "=== displayed by matrix list ==="
mat list PWX
dis "=== stored stripe ==="
local rn : rownames PWX
local re : roweq PWX
dis "rownames: `rn'"
dis "roweq   : `re'"
dis "=== rownumb() with the DISPLAYED key ==="
local k1 att:3vs1.sgq
local k2 att:3vs1bn.sgq
local k3 att:7vs1.sgq
dis "rownumb displayed-key `k1' = " rownumb(PWX,"`k1'")
dis "rownumb stored-key    `k2' = " rownumb(PWX,"`k2'")
dis "rownumb plain-key     `k3' = " rownumb(PWX,"`k3'")
cap noi dis "value by displayed key = " PWX["`k1'","estimate"]
dis "SUB_RC=" _rc
cap noi dis "value by stored key    = " PWX["`k2'","estimate"]
dis "SUB2_RC=" _rc
cap noi dis "value by plain key     = " PWX["`k3'","estimate"]
dis "SUB3_RC=" _rc
