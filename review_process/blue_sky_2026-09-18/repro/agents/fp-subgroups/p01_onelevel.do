clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

* ---------------- data ----------------
local gN 80
local gT 8
local gt0 2000
local gonsets "2003 2005"
local gysd 0.30
include "`GEN'/gen_core.do"

codebook sgall sg2, compact
tab sg2, m

* ---------------- run A: no subgroups ----------------
qui wooldid y cohort t ttre , att att_i att_it att_itime esprelength(2) espostlength(2) cluster(cohort)
tempname A_b A_V A_main A_bes A_Ves
matrix `A_b' = e(b)
matrix `A_V' = e(V)
matrix `A_main' = e(wdidmainresults)
matrix `A_bes' = e(b_es_att)
matrix `A_Ves' = e(V_es_att)
local A_N = e(N)
local A_df = e(df_r)

* ---------------- run B: one-level subgroups ----------------
cap noi wooldid y cohort t ttre , att att_i att_it att_itime esprelength(2) espostlength(2) cluster(cohort) subgroups(sgall)
dis "B_RC=" _rc
tempname B_b B_V B_main B_bes B_Ves
matrix `B_b' = e(b)
matrix `B_V' = e(V)
matrix `B_main' = e(wdidmainresults)
matrix `B_bes' = e(b_es_att)
matrix `B_Ves' = e(V_es_att)
local B_N = e(N)
local B_df = e(df_r)

dis "A_N=`A_N' B_N=`B_N'  A_df=`A_df' B_df=`B_df'"

dis "=== A colnames e(b) ==="
local an : colnames `A_b'
dis "`an'"
dis "=== B colnames e(b) ==="
local bn : colnames `B_b'
dis "`bn'"

dis "=== A wdidmainresults ==="
mat list `A_main', format(%21x)
dis "=== B wdidmainresults ==="
mat list `B_main', format(%21x)

dis "=== A b_es_att ==="
mat list `A_bes', format(%21x)
dis "=== B b_es_att ==="
mat list `B_bes', format(%21x)

dis "=== A V_es_att ==="
mat list `A_Ves', format(%21x)
dis "=== B V_es_att ==="
mat list `B_Ves', format(%21x)

dis "=== B ereturn list (names only) ==="
qui wooldid y cohort t ttre , att att_i att_it att_itime esprelength(2) espostlength(2) cluster(cohort) subgroups(sgall)
ereturn list, all
