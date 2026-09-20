* =============================================================================
* SELF-CONTAINED REPRO for the four seeded claims, on the SKEPTIC's own data.
* Builds its own panel; needs nothing but the pinned build, honestdid and pretrends.
*
* Covers, with a tag per claim:
*   F5     -- e(V) diag = -1e100 with b = 0 on the pretreatment ame rows while
*             e(V_es_ame) carries EXACTLY 0 on the same horizons
*   ADV    -- the adversarial ame statistic is bitwise the post-only statistic
*   HONEST -- honestdid delta(rm) returns an M-invariant set on e(b_es_ame)/e(V_es_ame)
*   CTRL   -- ... and returns the SAME M-invariance on a hand-built b/V with a zero
*             pre block and NO wooldid anywhere (so the invariance is honestdid's)
*   PRET   -- Roth's pretrends, the other advertised handoff, REFUSES loudly here
*   SENT   -- no e() scalar records the e(V) sentinel count, before or after
*             estimates save/use
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_seeded_claims.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
dis "=== build under test ==="
which wooldid

* ------------------------------------------------------------------- DGP ---
* hospital-year administrative panel, i = state (12), 30 hospitals per state,
* T = 2008..2019, staggered onsets 2012/2014/2016, 3 never-expanding states,
* dose = hospital exposure switched on at that state's own onset.
clear
set obs 12
gen byte state = _n
gen int expyr = .
replace expyr = 2012 if inlist(state,1,2,3,4)
replace expyr = 2014 if inlist(state,5,6,7)
replace expyr = 2016 if inlist(state,8,9)
expand 30
bys state: gen int hosp = _n
gen long hid = state*100 + hosp
gen double expo = 0.05 + 0.35*rbeta(2,3)
gen double hfe = rnormal(0,1.2)
expand 12
bys hid: gen int year = 2007 + _n
gen double dose = expo*(expyr<. & year>=expyr)
gen double y = 4 + hfe + 0.03*(year-2007) + 0.8*dose + rnormal(0,0.9)
compress
save "rs_panel.dta", replace

mata:
// verbatim transcription of the alt==0 / tt==0 branch of wooldid's own
// wooldid_advW (build wooldid.ado:14105-14155)
real scalar rs_advW(real colvector b, real matrix Vin, real colvector rel,
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
void rs_adv(string scalar fam, string scalar bn, string scalar vn,
            real scalar npre, real scalar dfr, real scalar printed)
{
    real colvector b, rel, relp, bp
    real matrix V, Vp
    real scalar k, i, qf, qp, Wf, Wp, pf, pp
    b = st_matrix(bn)'
    V = st_matrix(vn)
    k = rows(b)
    rel = J(k,1,.)
    for (i=1; i<=k; i++) rel[i] = (i <= npre ? -(npre+2-i) : i-npre-1)
    qf = .
    Wf = rs_advW(b, V, rel, -1, qf)
    bp = b[(npre+1)::k]
    Vp = V[(npre+1)::k, (npre+1)::k]
    relp = rel[(npre+1)::k]
    qp = .
    Wp = rs_advW(bp, Vp, relp, -1, qp)
    pf = Ftail(qf, dfr, Wf/qf)
    pp = Ftail(qp, dfr, Wp/qp)
    printf("{txt}ADV_%s: W_full=%21x q=%g | W_postonly=%21x q=%g | BITWISE_EQUAL=%g\n",
           fam, Wf, qf, Wp, qp, (Wf==Wp)&(qf==qp))
    printf("{txt}ADV_%s: my_p_full=%18.15f  my_p_postonly=%18.15f  wooldid_printed=%18.15f  MATCH=%g\n",
           fam, pf, pp, printed, (abs(pf-printed)<1e-12))
    printf("{txt}ADV_%s: whole pre BLOCK of V exactly zero (rows AND cross-covariances) = %g\n",
           fam, all(vec(V[1::npre, .]) :== 0))
}
end

* ============================ the run ======================================
use "rs_panel.dta", clear
dis "###### wooldid, contreat + untreatedZero, cluster(hid) ######"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(untreatedZero) adversarial jointtests pretrendpower(80) cluster(hid)
dis "RUN_RC=" _rc

* ---------------------------------------------------------------- F5 ------
tempname V B
matrix `B' = e(b)
matrix `V' = e(V)
local bn : colfullnames `B'
local ns 0
local j 0
foreach nm of local bn {
  local ++j
  if `V'[`j',`j'] < 0 & !missing(`V'[`j',`j']) {
    local ++ns
    dis "F5_SENTINEL|" %-22s "`nm'" "| b=" %14.10f `B'[1,`j'] "| Vdiag=" %12.0g `V'[`j',`j']
  }
}
dis "F5_nsentinels=" `ns' " of " colsof(`B')
matrix VA = e(V_es_ame)
matrix BA = e(b_es_ame)
local nz 0
forval j = 1/`=colsof(VA)' {
  if VA[`j',`j'] == 0 local ++nz
}
dis "F5_V_es_ame_exact_zero_diagonals=" `nz' " of " colsof(VA)
dis "F5_es_ndegen_ame=" e(es_ndegen_ame)
dis "F5_two_representations_disagree=" (`ns' > 0 & `nz' > 0)

* --------------------------------------------------------------- ADV ------
matrix VT = e(V_es_att)
matrix BT = e(b_es_att)
local dfr = e(df_r)
mata: rs_adv("ame", "BA", "VA", `=e(numpre_ame)', `dfr', `=e(adversarial_pval_ame)')
mata: rs_adv("att", "BT", "VT", `=e(numpre_att)', `dfr', `=e(adversarial_pval_att)')

* ------------------------------------------------------------ HONEST ------
local npA = e(numpre_ame)
local npT = e(numpre_att)
dis "###### honestdid on e(b_es_ame)/e(V_es_ame) ######"
cap noisily honestdid, b(BA) vcov(VA) numpre(`npA') mvec(0(0.5)1) delta(rm) parallel(0) matasave(HA)
dis "HONEST_ame_rc=" _rc
cap noisily mata: printf("{txt}HONEST_ame_M_INVARIANT=%g\n", (HA.CI[2,2]==HA.CI[3,2])&(HA.CI[3,2]==HA.CI[4,2])&(HA.CI[2,3]==HA.CI[3,3])&(HA.CI[3,3]==HA.CI[4,3]))
dis "###### honestdid on the att family, same run (comparator) ######"
cap noisily honestdid, b(BT) vcov(VT) numpre(`npT') mvec(0(0.5)1) delta(rm) parallel(0) matasave(HT)
dis "HONEST_att_rc=" _rc
cap noisily mata: printf("{txt}HONEST_att_M_INVARIANT=%g\n", (HT.CI[2,2]==HT.CI[3,2])&(HT.CI[3,2]==HT.CI[4,2]))

* --------------------------------------------------------------- PRET -----
dis "###### pretrends (the other advertised handoff) on the same two families ######"
cap noisily pretrends power 0.5, b(BA) vcov(VA) numpre(`npA')
dis "PRET_ame_rc=" _rc
cap noisily pretrends power 0.5, b(BT) vcov(VT) numpre(`npT')
dis "PRET_att_rc=" _rc

* --------------------------------------------------------------- SENT -----
dis "###### does any e() scalar record the sentinel count? ######"
ereturn list
estimates save "rs_saved", replace
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
estimates use "rs_saved"
tempname RV
matrix `RV' = e(V)
local ns2 0
forval j = 1/`=colsof(`RV')' {
  if `RV'[`j',`j'] < 0 & !missing(`RV'[`j',`j']) local ++ns2
}
dis "SENT_after_restore_nsentinels=" `ns2'
dis "SENT_after_restore_es_ndegen_ame=" e(es_ndegen_ame)
dis "###### restored ereturn list -- look for a scalar naming the sentinel count ######"
ereturn list

* ==================== CTRL: honestdid WITHOUT wooldid ======================
clear
matrix SB = (0, 0, 0, 1.2, 1.0, 0.8)
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
dis "###### CONTROL: hand-built b/V, zero pre block, NO wooldid ######"
cap noisily honestdid, b(SB) vcov(SV) numpre(3) mvec(0(0.5)1) delta(rm) parallel(0) matasave(SR)
dis "CTRL_rc=" _rc
cap noisily mata: printf("{txt}CTRL_M_INVARIANT_no_wooldid=%g\n", (SR.CI[2,2]==SR.CI[3,2])&(SR.CI[3,2]==SR.CI[4,2])&(SR.CI[2,3]==SR.CI[3,3])&(SR.CI[3,3]==SR.CI[4,3]))
matrix SV2 = SV
matrix SV2[1,1] = 0.08
matrix SV2[2,2] = 0.07
matrix SV2[3,3] = 0.06
cap noisily honestdid, b(SB) vcov(SV2) numpre(3) mvec(0(0.5)1) delta(rm) parallel(0) matasave(SR2)
cap noisily mata: printf("{txt}CTRL2_M_INVARIANT_nondegenerate_preblock=%g\n", (SR2.CI[2,2]==SR2.CI[3,2])&(SR2.CI[3,2]==SR2.CI[4,2]))

dis "RS_DONE=1"
