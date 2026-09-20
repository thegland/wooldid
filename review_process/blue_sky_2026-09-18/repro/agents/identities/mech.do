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
gen byte one = 1
drop if i == 3 & t == 7
drop if i == 11 & t == 8
drop if i == 12 & t >= 7
drop if i == 20 & t == 2
save base, replace

* ---------- 1: does e(b_es_att) post when gteffects is on? ----------
dis "SENT_1"
qui wooldid y i t ttre , att att_it gteffects esprelength(2) espostlength(6) cluster(i)
cap confirm matrix e(b_es_att)
dis "gteffects-on: confirm e(b_es_att) rc=" _rc
cap confirm matrix e(V_es_att)
dis "gteffects-on: confirm e(V_es_att) rc=" _rc
dis "SENT_1_END"

* ---------- 2: one-level subgroups vs no subgroups ----------
dis "SENT_2"
qui wooldid y i t ttre , att att_it esprelength(2) espostlength(6) cluster(i)
mat A = e(b)
qui est store noSG
cap noi wooldid y i t ttre , att att_it esprelength(2) espostlength(6) cluster(i) subgroups(one)
dis "one-level subgroups rc=" _rc
cap mat B = e(b)
cap mat list A, format(%21x)
cap mat list B, format(%21x)
dis "SENT_2_END"

* ---------- 3: fastnose points vs full points ----------
dis "SENT_3"
use base, clear
qui wooldid y i t ttre , att att_it att_i att_itime esprelength(2) espostlength(6) cluster(i)
mat F0 = e(b)
qui wooldid y i t ttre , att att_it att_i att_itime esprelength(2) espostlength(6) cluster(i) fastnose
mat F1 = e(b)
mata: a=st_matrix("F0"); b=st_matrix("F1"); printf("fastnose maxreldiff = %g ; cols %g vs %g\n", max(abs(a-b):/(abs(a):+1e-30)), cols(a), cols(b))
dis "SENT_3_END"

* ---------- 4: single cohort, balanced, unweighted ----------
dis "SENT_4"
clear
set obs 30
gen i = _n
gen ttre = cond(i<=15, 5, .)
expand 8
bys i: gen t = _n
gen ui = rnormal()
bys i (t): replace ui = ui[1]
gen y = ui + cond(t>=ttre & ttre<., 0.7, 0) + rnormal()
qui wooldid y i t ttre , att att_it att_i att_itime esprelength(2) espostlength(4) cluster(i)
mat S = e(b)
mat list S, format(%21x)
dis "SENT_4_END"

* ---------- 5: single cohort, UNBALANCED (drop rows for some units) ----------
dis "SENT_5"
drop if i == 2 & t == 8
drop if i == 5 & t >= 7
qui wooldid y i t ttre , att att_it att_i att_itime esprelength(2) espostlength(4) cluster(i)
mat S2 = e(b)
mat list S2, format(%21x)
dis "SENT_5_END"
