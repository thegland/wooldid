* =============================================================================
* R2 -- SKEPTIC: adversarial rank claim + sentinel-marker claim, my own data.
*
* Kernel: I transcribe the alt=0 / tt=0 branch of wooldid's own wooldid_advW
* (ado:14105-14155) verbatim into skeptic_advW(), because the package's mata
* functions are not addressable from a caller's session.  The transcription is
* VALIDATED by reproducing wooldid's own printed adversarial p-values for the
* att family (a non-degenerate covariance) before it is used on the ame family.
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r2_adv.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
which wooldid

mata:
real scalar skeptic_advW(real colvector b, real matrix Vin, real colvector rel,
                         real scalar esrel, real scalar qout)
{
    real matrix V, X, U, T, M
    real rowvector lam, sv
    real colvector z, g, r, sel
    real scalar k1, scale, tol, rankM, q, W
    qout = .
    if (hasmissing(Vin) | hasmissing(b) | hasmissing(rel)) return(.)
    V = (Vin + Vin')/2
    U = .
    lam = .
    symeigensystem(V, U, lam)
    scale = max(abs(lam))
    if (scale <= 0 | scale >= .) return(.)
    tol = scale*1e-12
    sel = selectindex(lam' :> tol)
    k1 = rows(sel)
    if (k1 == 0) return(.)
    T = diag(1:/sqrt(lam[sel]')) * U[., sel]'
    z = T*b
    X = rel :- esrel
    M = T*X
    sv = svdsv(M)
    rankM = sum(sv :> max(sv)*1e-12)
    q = k1 - rankM
    if (q <= 0) return(.)
    g = svsolve(M, z)
    r = z - M*g
    W = r'r
    if (W >= .) return(.)
    qout = q
    return(W)
}
end

cap prog drop advprobe
prog def advprobe
  syntax , tag(string) fam(string)
  cap confirm matrix e(V_es_`fam')
  if _rc {
    dis "ADV_`tag'_`fam'_NOV=1"
    exit
  }
  tempname B V
  matrix `B' = e(b_es_`fam')
  matrix `V' = e(V_es_`fam')
  scalar SDFR = e(df_r)
  scalar SPRT = e(adversarial_pval_`fam')
  local k = colsof(`B')
  local npre = e(numpre_`fam')
  matrix REL_ = J(`k',1,.)
  forval j = 1/`k' {
    if `j' <= `npre'  matrix REL_[`j',1] = -(`npre'+2-`j')
    if `j' >  `npre'  matrix REL_[`j',1] = `j'-`npre'-1
  }
  dis "RELVEC_`tag'_`fam':"
  matrix list REL_, noheader
  * full
  cap noisily mata: qf=.; Wf=skeptic_advW(st_matrix("`B'")', st_matrix("`V'"), st_matrix("REL_"), -1, qf); st_numscalar("Wf_",Wf); st_numscalar("qf_",qf)
  local rcf = _rc
  * post-only
  local p1 = `npre'+1
  matrix BP_ = `B'[1,`p1'..`k']
  matrix VP_ = `V'[`p1'..`k', `p1'..`k']
  matrix RP_ = REL_[`p1'..`k', 1]
  cap noisily mata: qp=.; Wp=skeptic_advW(st_matrix("BP_")', st_matrix("VP_"), st_matrix("RP_"), -1, qp); st_numscalar("Wp_",Wp); st_numscalar("qp_",qp)
  local rcp = _rc
  cap dis "ADVW_`tag'_`fam'_FULL=" %21x Wf_ " q=" qf_
  cap dis "ADVW_`tag'_`fam'_POST=" %21x Wp_ " q=" qp_
  cap dis "ADVW_`tag'_`fam'_BITEQ=" (Wf_==Wp_ & qf_==qp_)
  cap dis "ADVP_`tag'_`fam'_MINE_FULL=" %14.12f Ftail(qf_, SDFR, Wf_/qf_)
  cap dis "ADVP_`tag'_`fam'_MINE_POST=" %14.12f Ftail(qp_, SDFR, Wp_/qp_)
  cap dis "ADVP_`tag'_`fam'_PRINTED="   %14.12f SPRT
  cap dis "ADVP_`tag'_`fam'_MATCH=" (reldif(Ftail(qf_, SDFR, Wf_/qf_), SPRT) < 1e-9)
end

foreach CFG in IID CLHOSP {
  use "r1_panel.dta", clear
  local vopt ""
  if "`CFG'" == "CLHOSP" local vopt cluster(hid)
  dis "############### R2 CFG=`CFG' ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(untreatedZero) adversarial jointtests `vopt'
  dis "R2RC_`CFG'=" _rc
  if _rc continue
  dis "DFR_`CFG'=" e(df_r)
  advprobe , tag(`CFG') fam(att)
  advprobe , fam(ame) tag(`CFG')

  * ---- PS-62 splice integrity: e(V) ame post block vs e(V_es_ame) post block
  tempname FV EV
  matrix `FV' = e(V)
  matrix `EV' = e(V_es_ame)
  local fn : colfullnames `FV'
  local ix = 0
  local i2 = 0
  forval j = 1/`=colsof(`FV')' {
    local nm : word `j' of `fn'
    if "`nm'" == "ES_ame_post0" local ix = `j'
    if "`nm'" == "ES_ame_post2" local i2 = `j'
  }
  dis "SPLICE_`CFG'_idx=" `ix' " .. " `i2'
  if `ix' > 0 & `i2' >= `ix' {
    matrix SUBFV = `FV'[`ix'..`i2', `ix'..`i2']
    matrix SUBEV = `EV'[4..6, 4..6]
    matrix DIFFS = SUBFV - SUBEV
    dis "SPLICE_`CFG'_maxabsdiff=" %21x mreldif(SUBFV, SUBEV)
    matrix list SUBFV, format(%21x)
    matrix list SUBEV, format(%21x)
  }

  * ---- hv-noVsentinel-marker
  if "`CFG'" == "CLHOSP" {
    dis "=== EREturn LIST (CLHOSP) ==="
    ereturn list
    estimates save "r2_clhosp", replace
  }
  if "`CFG'" == "IID" estimates save "r2_iid", replace
}

* ======== does the sentinel count survive estimates use? =================
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cap noisily estimates use "r2_clhosp"
dis "RESTORE_RC=" _rc
dis "RESTORED_ndegen_ame=" e(es_ndegen_ame)
cap confirm matrix e(V)
if !_rc {
  tempname RV
  matrix `RV' = e(V)
  local ns 0
  forval j = 1/`=colsof(`RV')' {
    if `RV'[`j',`j'] < 0 & !missing(`RV'[`j',`j']) local ++ns
  }
  dis "RESTORED_NSENT=" `ns'
}
dis "=== restored ereturn list ==="
ereturn list

dis "R2_DONE=1"
