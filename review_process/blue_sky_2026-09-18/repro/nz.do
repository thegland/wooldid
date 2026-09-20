clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set graphics off
set seed 909090
set sortseed 909090
set obs 60
gen i = _n
gen ttre = .
replace ttre = 6 if i<=20
replace ttre = 9 if i>20 & i<=40
expand 12
bys i: gen t = _n
gen y = 0.2*i/10 + 0.3*t + rnormal() + 1.2*(t>=ttre)
dis "=== POOLED (default) reference, espre(3) espost(3) ==="
wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3)
dis "COLNAMES_POOLED:"
mat b = e(b)
local nm : colnames b
dis "`nm'"
dis "NUMPRE=" e(numpre_att)
mat list b, format(%12.6f)
dis "--- b_es_att ---"
cap mat list e(b_es_att), format(%12.6f)
dis "=== FIXED reference for contrast ==="
wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-1)
mat b2 = e(b)
local nm2 : colnames b2
dis "COLNAMES_FIXED: `nm2'"
mat list b2, format(%12.6f)
dis "NZ_DONE"
