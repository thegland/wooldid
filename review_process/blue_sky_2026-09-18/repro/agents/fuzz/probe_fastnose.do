clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
use "`D'/ds_canon.dta", clear
ereturn clear
wooldid y unit t ttre, att att_i att_itime cluster(region) fastnose esprelength(2) espostlength(2)
dis "FN_RC=" _rc
ereturn list
dis "== dim b: " colsof(e(b))
capture matrix VV = e(V)
dis "VRC=" _rc
if _rc==0 dis "== dim V: " rowsof(VV) " x " colsof(VV)
dis "PROPS=`e(properties)'"
dis "--- downstream consumers ---"
capture noisily test _b[att]=0
dis "TEST_RC=" _rc
capture noisily lincom _b[att]
dis "LINCOM_RC=" _rc
capture noisily estimates store fnrun
dis "ESTSTORE_RC=" _rc
capture noisily estimates replay fnrun
dis "REPLAY_RC=" _rc
capture noisily nlcom _b[att]*2
dis "NLCOM_RC=" _rc
dis "--- does e(b_es_att)/e(numpre_att) exist (help item 11)? ---"
capture confirm matrix e(b_es_att)
dis "BESATT_RC=" _rc
capture confirm scalar e(numpre_att)
dis "NUMPRE_RC=" _rc
dis "MATRICES: `e(matrices)'"
dis "SCALARS: `e(scalars)'"
