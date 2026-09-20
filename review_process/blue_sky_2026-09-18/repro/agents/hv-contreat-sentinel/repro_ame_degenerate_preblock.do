* =============================================================================
* SELF-CONTAINED REPRO
*
* Continuous-treatment DiD with the canonical "exposure x post" dose
* (dose_it = D_i * 1{t >= g_i}).  wooldid posts e(b_es_ame)/e(V_es_ame) in which
* the three PRE-treatment horizons are exact structural zeros with EXACTLY-ZERO
* variance (and -1e100 sentinels in e(V)).
*
* Two consequences measured here, neither of which is announced at the point of
* use:
*  (1) honestdid -- the handoff wooldid's own help advertises -- returns a
*      robust confidence set that is IDENTICAL for M = 0, 0.5 and 1.0, i.e. it
*      reports the ame effect as completely insensitive to violations of
*      parallel trends.  It is not: the pre-period covariance is zero, so the
*      relative-magnitudes bound M * max|pre violation| is 0 for every M.
*      The att family on the SAME run behaves normally (interval widens with M).
*  (2) the reported "adversarial" ame p-value is bitwise identical to the one
*      computed from the THREE POST-TREATMENT HORIZONS ALONE: a pre-trend
*      diagnostic that uses no pre-treatment information.  On the same run the
*      ame joint pre-trend test correctly prints "." and pretrendpower()
*      correctly refuses with a caution.
*
* Run:  stata-mp -e do repro_ame_degenerate_preblock.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
dis "=== WHICH (must be the pinned build) ==="
which wooldid

* ---- my own minimum-distance line statistic: independent reimplementation of
*      the help's description of the adversarial test (GLS under a generalized
*      inverse of e(V_es_*); numerator df = rank(V) - 1) ----------------------
mata:
real scalar mdline(real colvector b, real matrix V, real colvector rel,
                   real scalar esrel, real scalar qout)
{
  real matrix Vs, U, T, M
  real rowvector lam, sv
  real colvector sel, z, X, g, r
  real scalar scale, tol, k1, rankM, q
  qout = .
  Vs = (V + V')/2
  U = .
  lam = .
  symeigensystem(Vs, U, lam)
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
  qout = q
  return(r'r)
}
end

* --------------------------------------------------------------------- DGP
* 100 units in 5 cohorts (4 treated, onset 4/6/8/10; 1 never treated), T = 1..12.
* dose = unit-specific intensity, switched on at that unit's own onset.
clear
set obs 100
gen long unit = _n
gen byte coh  = ceil(unit/20)
gen int ttre  = .
replace ttre = 4  if coh==1
replace ttre = 6  if coh==2
replace ttre = 8  if coh==3
replace ttre = 10 if coh==4
gen double dose0 = 0.5 + 1.5*runiform()
gen double ufe   = rnormal()
expand 12
bys unit: gen int t = _n
gen double dose = dose0*(ttre<. & t>=ttre)
gen double y    = ufe + 0.2*t + 0.8*dose + rnormal()
save "repro_panel.dta", replace

foreach VCE in IID ROBUST CLUSTERUNIT {
  use "repro_panel.dta", clear
  local vopt ""
  if "`VCE'" == "ROBUST"      local vopt customvce(robust)
  if "`VCE'" == "CLUSTERUNIT" local vopt cluster(unit)

  dis "############################ VCE = `VCE' ############################"
  cap noisily wooldid y coh t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(untreatedZero) adversarial jointtests pretrendpower(80) `vopt'
  dis "EST_RC_`VCE'=" _rc
  if _rc != 0 continue

  matrix BA = e(b_es_ame)
  matrix VA = e(V_es_ame)
  matrix BT = e(b_es_att)
  matrix VT = e(V_es_att)
  scalar DFR = e(df_r)
  local npA = e(numpre_ame)
  local npT = e(numpre_att)
  dis "NDEGEN_ame_`VCE'=" e(es_ndegen_ame) "   NDEGEN_att_`VCE'=" e(es_ndegen_att)
  dis "e(b_es_ame) and e(V_es_ame):"
  matrix list BA
  matrix list VA, format(%12.0g)

  * ---- (2) adversarial: full 6 horizons vs the 3 post horizons alone --------
  matrix REL  = (-4\-3\-2\0\1\2)
  matrix RELP = (0\1\2)
  matrix BAP  = (BA[1,4], BA[1,5], BA[1,6])
  matrix VAP  = VA[4..6,4..6]
  mata: q1 = .; W1 = mdline(st_matrix("BA")', st_matrix("VA"), st_matrix("REL"), -1, q1); st_numscalar("W_full", W1); st_numscalar("q_full", q1)
  mata: q2 = .; W2 = mdline(st_matrix("BAP")', st_matrix("VAP"), st_matrix("RELP"), -1, q2); st_numscalar("W_post", W2); st_numscalar("q_post", q2)
  mata: q3 = .; W3 = mdline(st_matrix("BT")', st_matrix("VT"), st_matrix("REL"), -1, q3); st_numscalar("W_att", W3); st_numscalar("q_att", q3)
  dis "ADV_`VCE'_W_full="  %21x W_full " q=" q_full
  dis "ADV_`VCE'_W_post="  %21x W_post " q=" q_post
  dis "ADV_`VCE'_BITEQUAL=" (W_full == W_post & q_full == q_post)
  dis "ADV_`VCE'_p_full="     Ftail(q_full, DFR, W_full/q_full)
  dis "ADV_`VCE'_p_postonly=" Ftail(q_post, DFR, W_post/q_post)
  dis "ADV_`VCE'_p_att="      Ftail(q_att,  DFR, W_att/q_att)  " (q=" q_att ")"

  * ---- (1) the advertised honestdid handoff --------------------------------
  foreach D in rm sd {
    dis "--- honestdid AME  vce=`VCE' delta(`D') ---"
    matrix HB = BA
    matrix HV = VA
    cap noisily honestdid, b(HB) vcov(HV) numpre(`npA') mvec(0(0.5)1) delta(`D') parallel(0)
    dis "HONEST_AME_`VCE'_`D'_RC=" _rc
    dis "--- honestdid ATT  vce=`VCE' delta(`D') (comparator: behaves normally) ---"
    matrix HB2 = BT
    matrix HV2 = VT
    cap noisily honestdid, b(HB2) vcov(HV2) numpre(`npT') mvec(0(0.5)1) delta(`D') parallel(0)
    dis "HONEST_ATT_`VCE'_`D'_RC=" _rc
  }
}
dis "=== ALLDONE ==="
