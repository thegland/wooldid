* ============================================================================
* P4 -- CONSEQUENCE SURFACE of the contreat() sentinel / zero-diagonal.
* Self-contained.  Asks: what does a real user's postestimation actually do
* with (a) the -1e100 in e(V) and (b) the exact 0 in e(V_es_ame)?
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
dis "=== WHICH ==="
which wooldid

* ---- DGP: canonical "exposure x post" continuous DiD ------------------------
clear
set obs 100
gen long unit = _n
gen byte coh = ceil(unit/20)
gen int ttre = .
replace ttre = 4  if coh==1
replace ttre = 6  if coh==2
replace ttre = 8  if coh==3
replace ttre = 10 if coh==4
gen double dose0 = 0.5 + 1.5*runiform()
gen double ufe   = rnormal()
expand 12
bys unit: gen int t = _n
gen double dose = dose0*(ttre<. & t>=ttre)
gen double y = ufe + 0.2*t + 0.8*dose + rnormal()

dis "=== ESTIMATE (jointtests + adversarial + pretrendpower) ==="
cap noisily wooldid y coh t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero) jointtests adversarial pretrendpower(80)
dis "EST_RC=" _rc

* ---------------------------------------------------------------- marker scan
dis "=== MARKER SCAN: is there ANY e() scalar recording the e(V) sentinel? ==="
dis "M_es_ndegen_ame=" e(es_ndegen_ame)
dis "M_es_ndegen_ame_it=" e(es_ndegen_ame_it)
dis "M_es_ndegen_att=" e(es_ndegen_att)
dis "M_wdid_vcvstate=" e(wdid_vcvstate)
dis "M_wdid_vblocks=" e(wdid_vblocks)
dis "M_numpre_ame=" e(numpre_ame)
dis "--- full hidden ereturn list follows (grep for anything sentinel-shaped) ---"
ereturn list, all

* ---------------------------------------------------------------- native table
dis "=== e(wdidmainresults): is pretreatment:ame 0 or missing? ==="
matrix WM = e(wdidmainresults)
local wr : rowfullnames WM
local jj 0
foreach rn of local wr {
  local ++jj
  dis "WMR | " %-22s "`rn'" " | est=" %12.0g WM[`jj',1] " | se=" %12.0g WM[`jj',2]
}

* ---------------------------------------------------------------- esttab
dis "=== esttab (what a user pastes into a paper) ==="
cap noisily esttab ., se b(%9.4f) nostar wide
dis "ESTTAB_RC=" _rc
cap noisily esttab ., ci b(%9.4f) nostar
dis "ESTTABCI_RC=" _rc

* ---------------------------------------------------------------- test/lincom
dis "=== test / lincom against the -1e100 rows ==="
cap noisily test ES_ame_pre4 ES_ame_pre3 ES_ame_pre2
dis "T1_RC=" _rc " p=" r(p) " df=" r(df)
cap noisily test ES_ame_pre4 ES_ame_post0
dis "T2_RC=" _rc " p=" r(p) " df=" r(df)
cap noisily test ES_att_pre4 ES_att_pre3 ES_att_pre2
dis "T3_RC=" _rc " p=" r(p) " df=" r(df)
cap noisily lincom ES_ame_pre4
dis "L1_RC=" _rc " est=" r(estimate) " se=" r(se)
cap noisily lincom ES_ame_post0 - ES_ame_pre4
dis "L2_RC=" _rc " est=" r(estimate) " se=" r(se)
cap noisily lincom main_ame + 0.0000001*ES_ame_pre4
dis "L3_RC=" _rc " est=" r(estimate) " se=" r(se)

* --------------------------------------------------- matrix-algebra consumers
dis "=== matrix consumers on e(V) ==="
matrix VV = e(V)
mata: V = st_matrix("VV"); st_local("mineig", strofreal(min(Re(eigenvalues(makesymmetric(V))))))
dis "MINEIG=" `mineig'
mata: st_local("nchmiss", strofreal(sum(cholesky(makesymmetric(st_matrix("VV"))) :>= .)))
dis "CHOL_MISSING_CELLS=" `nchmiss'
cap noisily estat vce
dis "ESTATVCE_RC=" _rc

* --------------------------------------------- the e(V_es_ame) copy, reposted
dis "=== the OTHER copy: ereturn post e(b_es_ame), e(V_es_ame) [help says valid for lincom/test] ==="
matrix BE = e(b_es_ame)
matrix VE = e(V_es_ame)
matrix list BE
matrix list VE, format(%12.0g)
local npre = e(numpre_ame)
dis "NUMPRE_AME=" `npre'
tempname hold
estimates store wdrun
cap noisily ereturn post BE VE
dis "REPOST_RC=" _rc
cap noisily ereturn display
dis "EDISP_RC=" _rc
cap noisily lincom pre4
dis "R_L1_RC=" _rc " est=" r(estimate) " se=" r(se) " lb=" r(lb) " ub=" r(ub)
cap noisily lincom pre4 + pre3 + pre2
dis "R_L2_RC=" _rc " est=" r(estimate) " se=" r(se) " lb=" r(lb) " ub=" r(ub)
cap noisily test pre4 pre3 pre2
dis "R_T1_RC=" _rc " p=" r(p) " df=" r(df)

* --------------------------------------------------------------- honestdid
dis "=== honestdid handoff on the degenerate covariance ==="
cap noisily honestdid, b(BE) vcov(VE) numpre(`npre') mvec(0(0.5)1) parallel(0)
dis "HONEST_RC=" _rc
cap noisily matrix list r(HonestEventStudy)

* --------------------------------------------------- estimates save / restore
dis "=== does the sentinel survive estimates save/use, and what marks it? ==="
estimates restore wdrun
cap noisily estimates save "p4_saved", replace
dis "SAVE_RC=" _rc
clear all
cap noisily estimates use "p4_saved"
dis "USE_RC=" _rc
matrix VB = e(V)
matrix BB = e(b)
local bn : colfullnames BB
local j 0
local nsent 0
foreach nm of local bn {
  local ++j
  if VB[`j',`j'] < 0 & !missing(VB[`j',`j']) local ++nsent
}
dis "AFTERUSE_NSENT=" `nsent'
dis "AFTERUSE_es_ndegen_ame=" e(es_ndegen_ame)
dis "AFTERUSE_cmd=" e(cmd)
cap noisily test ES_ame_pre4 ES_ame_pre3 ES_ame_pre2
dis "AFTERUSE_T_RC=" _rc " p=" r(p) " df=" r(df)

dis "=== ALLDONE ==="
