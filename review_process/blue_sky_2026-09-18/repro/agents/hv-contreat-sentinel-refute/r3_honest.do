* =============================================================================
* R3 -- SKEPTIC: (a) adversarial rank, computed with a transcription of
*                    wooldid's own kernel, printed from inside mata;
*                (b) the honestdid handoff at M = 0 / 0.5 / 1, AME vs ATT;
*                (c) THE DECISIVE CONTROL: the same M-invariance from a
*                    hand-built b/V with a zero pre block and NO wooldid at all
*                    -- does the invariance belong to wooldid or to honestdid?
*                (d) the Roth `pretrends' handoff on the same inputs.
*
* Data: r1_panel.dta, built by r1_own.do (hospital-year admin panel, i = state,
*       dose = hospital exposure x post).  r1_own.do is self-contained.
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r3_honest.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
which wooldid

mata:
// verbatim transcription of the alt==0 / tt==0 branch of wooldid_advW
// (build wooldid.ado:14105-14155)
real scalar sk_advW(real colvector b, real matrix Vin, real colvector rel,
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

void sk_advreport(string scalar tag, string scalar fam, string scalar bn,
                  string scalar vn, real scalar npre, real scalar dfr,
                  real scalar printed)
{
    real colvector b, rel, relp, bp
    real matrix V, Vp
    real scalar k, i, qf, qp, Wf, Wp, pf, pp
    b = st_matrix(bn)'
    V = st_matrix(vn)
    k = rows(b)
    rel = J(k, 1, .)
    for (i = 1; i <= k; i++) {
        if (i <= npre) rel[i] = -(npre + 2 - i)
        else           rel[i] = i - npre - 1
    }
    qf = .
    Wf = sk_advW(b, V, rel, -1, qf)
    bp   = b[(npre+1)::k]
    Vp   = V[(npre+1)::k, (npre+1)::k]
    relp = rel[(npre+1)::k]
    qp = .
    Wp = sk_advW(bp, Vp, relp, -1, qp)
    printf("{txt}ADV_%s_%s FULL W=%21x q=%g\n", tag, fam, Wf, qf)
    printf("{txt}ADV_%s_%s POST W=%21x q=%g\n", tag, fam, Wp, qp)
    printf("{txt}ADV_%s_%s BITEQ=%g\n", tag, fam, (Wf == Wp) & (qf == qp))
    pf = Ftail(qf, dfr, Wf/qf)
    pp = Ftail(qp, dfr, Wp/qp)
    printf("{txt}ADV_%s_%s p_full=%18.15f p_post=%18.15f printed=%18.15f\n",
           tag, fam, pf, pp, printed)
    printf("{txt}ADV_%s_%s MATCHPRINTED=%g\n", tag, fam,
           (abs(pf - printed) < 1e-12))
    // exactness of the pre block
    printf("{txt}ADV_%s_%s PREBLOCK_ALLEXACTZERO=%g\n", tag, fam,
           all(vec(V[1::npre, .]) :== 0))
}
end

* ======================= estimate once per VCE ==============================
foreach CFG in IID CLHOSP {
  use "r1_panel.dta", clear
  local vopt ""
  if "`CFG'" == "CLHOSP" local vopt cluster(hid)
  dis "############### R3 CFG=`CFG' ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(untreatedZero) adversarial jointtests pretrendpower(80) `vopt'
  dis "R3RC_`CFG'=" _rc
  if _rc continue

  matrix BA = e(b_es_ame)
  matrix VA = e(V_es_ame)
  matrix BT = e(b_es_att)
  matrix VT = e(V_es_att)
  local npA = e(numpre_ame)
  local npT = e(numpre_att)
  local dfr = e(df_r)
  local pA  = e(adversarial_pval_ame)
  local pT  = e(adversarial_pval_att)
  dis "R3_`CFG'_ndegen_ame=" e(es_ndegen_ame) " ndegen_att=" e(es_ndegen_att)
  dis "R3_`CFG'_dfr=" `dfr'

  mata: sk_advreport("`CFG'", "ame", "BA", "VA", `npA', `dfr', `pA')
  mata: sk_advreport("`CFG'", "att", "BT", "VT", `npT', `dfr', `pT')

  * ------------------------- honestdid handoff -----------------------------
  foreach FAM in ame att {
    if "`FAM'" == "ame" {
      matrix HB = BA
      matrix HV = VA
      local NP `npA'
    }
    else {
      matrix HB = BT
      matrix HV = VT
      local NP `npT'
    }
    foreach D in rm sd {
      dis "--- honestdid `FAM' `CFG' delta(`D') ---"
      cap noisily honestdid, b(HB) vcov(HV) numpre(`NP') mvec(0(0.5)1) delta(`D') parallel(0) matasave(HR_`FAM')
      dis "HON_`CFG'_`FAM'_`D'_RC=" _rc
      if _rc == 0 {
        cap noisily mata: printf("{txt}HONCI_`CFG'_`FAM'_`D'\n"); HR_`FAM'.CI
        cap noisily mata: st_matrix("HCI", HR_`FAM'.CI)
        cap matrix list HCI, format(%21x)
        cap mata: printf("{txt}HONINV_`CFG'_`FAM'_`D'_rows234identical=%g\n", (HR_`FAM'.CI[2,2]==HR_`FAM'.CI[3,2]) & (HR_`FAM'.CI[3,2]==HR_`FAM'.CI[4,2]) & (HR_`FAM'.CI[2,3]==HR_`FAM'.CI[3,3]) & (HR_`FAM'.CI[3,3]==HR_`FAM'.CI[4,3]))
      }
    }
  }

  * ------------------------- pretrends handoff -----------------------------
  foreach FAM in ame att {
    if "`FAM'" == "ame" {
      matrix PB = BA
      matrix PV = VA
      local NP `npA'
    }
    else {
      matrix PB = BT
      matrix PV = VT
      local NP `npT'
    }
    dis "--- pretrends power 0.5 `FAM' `CFG' ---"
    cap noisily pretrends power 0.5, b(PB) vcov(PV) numpre(`NP')
    dis "PRET_`CFG'_`FAM'_RC=" _rc
    cap dis "PRET_`CFG'_`FAM'_slope=" r(slope)
    dis "--- pretrends slope(0.05) `FAM' `CFG' ---"
    cap noisily pretrends, slope(0.05) b(PB) vcov(PV) numpre(`NP') nocoefplot
    dis "PRETS_`CFG'_`FAM'_RC=" _rc
  }
}

* ===========================================================================
* (c) THE CONTROL: no wooldid at all.  A hand-built 6-coefficient event study
*     whose 3 pre-period rows are exactly 0 with an exactly-0 covariance block,
*     and whose 3 post rows are an ordinary well-conditioned block.  If
*     honestdid is M-invariant HERE too, the invariance is a property of
*     honestdid given a degenerate pre block, not of anything wooldid computed.
* ===========================================================================
clear
matrix SB = (0, 0, 0, 1.2, 1.0, 0.8)
matrix colnames SB = pre3 pre2 pre1 post0 post1 post2
matrix SV = J(6,6,0)
matrix SV[4,4] = 0.09
matrix SV[5,5] = 0.10
matrix SV[6,6] = 0.12
matrix SV[4,5] = 0.03
matrix SV[5,4] = 0.03
matrix SV[4,6] = 0.02
matrix SV[6,4] = 0.02
matrix SV[5,6] = 0.04
matrix SV[6,5] = 0.04
matrix colnames SV = pre3 pre2 pre1 post0 post1 post2
matrix rownames SV = pre3 pre2 pre1 post0 post1 post2
dis "=== SYNTHETIC CONTROL (no wooldid) ==="
matrix list SB
matrix list SV
cap noisily honestdid, b(SB) vcov(SV) numpre(3) mvec(0(0.5)1) delta(rm) parallel(0) matasave(SR)
dis "SYN_RM_RC=" _rc
cap noisily mata: SR.CI
cap noisily mata: printf("{txt}SYN_RM_INVARIANT=%g\n", (SR.CI[2,2]==SR.CI[3,2]) & (SR.CI[3,2]==SR.CI[4,2]) & (SR.CI[2,3]==SR.CI[3,3]) & (SR.CI[3,3]==SR.CI[4,3]))
cap noisily honestdid, b(SB) vcov(SV) numpre(3) mvec(0(0.5)1) delta(sd) parallel(0) matasave(SR2)
dis "SYN_SD_RC=" _rc

* ... and the SAME synthetic with a non-degenerate pre block, as a sanity check
matrix SV2 = SV
matrix SV2[1,1] = 0.08
matrix SV2[2,2] = 0.07
matrix SV2[3,3] = 0.06
cap noisily honestdid, b(SB) vcov(SV2) numpre(3) mvec(0(0.5)1) delta(rm) parallel(0) matasave(SR3)
dis "SYN2_RM_RC=" _rc
cap noisily mata: SR3.CI
cap noisily mata: printf("{txt}SYN2_RM_INVARIANT=%g\n", (SR3.CI[2,2]==SR3.CI[3,2]) & (SR3.CI[3,2]==SR3.CI[4,2]))

dis "R3_DONE=1"
