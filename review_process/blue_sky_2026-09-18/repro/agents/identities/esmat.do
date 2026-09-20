clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 8
replace ttre = 6 if i > 8 & i <= 16
expand 8
bys i: gen t = _n
gen ui = rnormal()
gen ut = rnormal()/3
bys i (t): replace ui = ui[1]
bys t (i): replace ut = ut[1]
gen y = ui + ut + cond(t>=ttre & ttre<., 0.5+0.2*(t-ttre), 0) + rnormal()
drop if i == 3 & t == 7
dis "SENT1"
qui wooldid y i t ttre , att att_it esprelength(2) espostlength(6) cluster(i)
ereturn list
dis "SENT2"
dis "b_es_att exists? " _rc
cap confirm matrix e(b_es_att)
dis "confirm b_es_att rc=" _rc
cap confirm matrix e(V_es_att)
dis "confirm V_es_att rc=" _rc
matrix list e(b_es_att)
dis "SENT3"
