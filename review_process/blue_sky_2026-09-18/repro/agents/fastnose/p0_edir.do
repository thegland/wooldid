clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 1
sysuse auto, clear
qui reg price mpg weight
dis "A1:[`e(scalars)']"
local z "`e(scalars)'"
dis "A2:[`z'] LEN=" length("`z'")
local z2 : e(scalars)
dis "A3:[`z2']"
cap mata: st_local("m1", invtokens(st_dir("e()", "numscalar", "*")'))
dis "A4_rc=" _rc " [`m1']"
cap mata: st_local("m2", invtokens(st_dir("e()", "numeric scalar", "*")'))
dis "A5_rc=" _rc " [`m2']"
cap mata: st_local("m3", invtokens(st_dir("e()", "macro", "*")'))
dis "A6_rc=" _rc " [`m3']"
cap mata: st_local("m4", invtokens(st_dir("e()", "matrix", "*")'))
dis "A7_rc=" _rc " [`m4']"
cap mata: st_local("m5", invtokens(st_dir("e()", "scalar", "*")'))
dis "A8_rc=" _rc " [`m5']"
cap mata: st_local("m6", invtokens(st_dir("e()", "string", "*")'))
dis "A9_rc=" _rc " [`m6']"
dis "DONE_P0"
