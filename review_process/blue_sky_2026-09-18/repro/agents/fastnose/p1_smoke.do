clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

dis "WHICH_START"
which wooldid
dis "WHICH_END"

local gN 80
local gT 10
local gonsets "2004 2006 2008"
local gysd 0.30
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"

dis "NOBS=" _N
tab cohort, m

* ---------- analytic ----------
dis "=== ANALYTIC RUN ==="
ereturn clear
cap noi wooldid y unit t ttre , esprelength(3) espostlength(4)
dis "AN_RC=" _rc
dis "AN_SCALARS: `e(scalars)'"
dis "AN_MACROS: `e(macros)'"
dis "AN_MATRICES: `e(matrices)'"
dis "AN_numpre_att=" e(numpre_att)
cap confirm matrix e(b_es_att)
dis "AN_HAS_b_es_att_rc=" _rc
cap confirm matrix e(V_es_att)
dis "AN_HAS_V_es_att_rc=" _rc

* ---------- fastnose ----------
dis "=== FASTNOSE RUN ==="
ereturn clear
cap noi wooldid y unit t ttre , esprelength(3) espostlength(4) fastnose
dis "FN_RC=" _rc
dis "FN_SCALARS: `e(scalars)'"
dis "FN_MACROS: `e(macros)'"
dis "FN_MATRICES: `e(matrices)'"
dis "FN_numpre_att=" e(numpre_att)
cap confirm matrix e(b_es_att)
dis "FN_HAS_b_es_att_rc=" _rc
cap confirm matrix e(V_es_att)
dis "FN_HAS_V_es_att_rc=" _rc
dis "FN_wdid_fastnose=" e(wdid_fastnose)
dis "DONE_P1"
