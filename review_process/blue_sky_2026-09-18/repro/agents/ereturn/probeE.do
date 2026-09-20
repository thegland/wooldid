clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn"
local ES esprelength(3) espostlength(3)

use wd_base.dta, clear
qui wooldid y coh t ttre , att att_i att_it att_itime esfixedbaseperiod `ES' contreat(dose) contreatcontrols(untreatedZero)
tempname BB VV
matrix `BB' = e(b)
matrix `VV' = e(V)
local cn : colnames `BB'
local k 0
dis "ROW|NAME|b|Vdiag"
foreach c of local cn {
  local ++k
  dis "ROW|`c'|" %12.0g `BB'[1,`k'] "|" %12.4e `VV'[`k',`k']
}
dis "--- eigenvalues of e(V) ---"
mata: ev = Re(eigenvalues(st_matrix("e(V)"))); printf("mineig=%g maxeig=%g nneg=%g\n", min(ev), max(ev), sum(ev :< -1e-8))
dis "--- symmetry ---"
mata: V=st_matrix("e(V)"); printf("maxasym=%g\n", max(abs(V-V')))
dis "--- V_es_ame diagonal ---"
matrix list e(V_es_ame), format(%12.4e)
dis "--- b_es_ame ---"
matrix list e(b_es_ame), format(%12.4g)
dis "--- es_ndegen_ame = " e(es_ndegen_ame)
dis "--- numpre_ame = " e(numpre_ame)

dis "=== TEST 1: all-good ame post horizons ==="
cap noi test ES_ame_post0 ES_ame_post1 ES_ame_post2
dis "T1_p=" r(p) " df=" r(df)
dis "=== TEST 2: mixed (sentinel + good) ==="
cap noi test ES_ame_pre4 ES_ame_post0
dis "T2_p=" r(p) " df=" r(df)
dis "=== TEST 3: all sentinel ==="
cap noi test ES_ame_pre4 ES_ame_pre3 ES_ame_pre2
dis "T3_p=" r(p) " df=" r(df)
dis "=== TEST 4: cross-family (zero off-diagonal) ==="
cap noi test main_att = main_att_i
dis "T4_p=" r(p)
dis "V(main_att,main_att_i) = " `VV'[1,3]
dis "=== TEST 5: cross-family ES vs pooled ==="
cap noi lincom main_att - ES_att_post0
dis "T5_se=" r(se)
dis "=== estat / predict consumers ==="
cap noi estat vce
dis "ESTATVCE_RC=" _rc
cap noi matrix li e(V)
dis "=== does cholesky of e(V) work (a generic consumer) ==="
mata: rc = 0; V=st_matrix("e(V)"); C = cholesky(makesymmetric(V)); printf("cholesky_hasmissing=%g\n", missing(C))

dis "=== wdidmainresults row names vs e(b) colnames ==="
local rn : rowfullnames e(wdidmainresults)
dis "WDIDROWS=`rn'"
local bn : colnames e(b)
dis "BCOLS=`bn'"

dis "=== update option ==="
cap noi wooldid , update
dis "UPDATE_RC=" _rc

dis "PROBEE_DONE"
