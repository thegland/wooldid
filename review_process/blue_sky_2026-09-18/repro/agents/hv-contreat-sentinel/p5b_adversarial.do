* ============================================================================
* P5b -- INDEPENDENT reimplementation of wooldid's adversarial (minimum-distance
*        line) statistic, to determine how much of the reported ame p-value is
*        driven by the pre-treatment horizons when those horizons are
*        structurally degenerate (exact-zero variance).
* Self-contained.
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
dis "=== WHICH ==="
which wooldid

* ---- my own minimum-distance line statistic (mirrors the help's description:
*      GLS under a generalized inverse, numerator df = rank(V) - 1) -----------
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

cap noisily wooldid y coh t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero) adversarial
dis "EST_RC=" _rc

matrix BA = e(b_es_ame)
matrix VA = e(V_es_ame)
matrix BT = e(b_es_att)
matrix VT = e(V_es_att)
scalar DFR = e(df_r)
dis "DFR=" DFR
matrix list BA
matrix list VA, format(%12.0g)

matrix REL  = (-4\-3\-2\0\1\2)
matrix RELP = (0\1\2)
matrix BAP  = (BA[1,4], BA[1,5], BA[1,6])
matrix VAP  = VA[4..6,4..6]
matrix list BAP
matrix list VAP, format(%12.0g)

mata:
  q1 = .
  W1 = mdline(st_matrix("BA")', st_matrix("VA"), st_matrix("REL"), -1, q1)
  st_numscalar("W_full", W1); st_numscalar("q_full", q1)
  q2 = .
  W2 = mdline(st_matrix("BAP")', st_matrix("VAP"), st_matrix("RELP"), -1, q2)
  st_numscalar("W_post", W2); st_numscalar("q_post", q2)
  q3 = .
  W3 = mdline(st_matrix("BT")', st_matrix("VT"), st_matrix("REL"), -1, q3)
  st_numscalar("W_att", W3); st_numscalar("q_att", q3)
end

dis "W_FULL_AME=" %21x W_full "  q_full=" q_full
dis "W_POSTONLY_AME=" %21x W_post "  q_post=" q_post
dis "W_ATT=" %21x W_att "  q_att=" q_att
dis "P_FULL_AME="     Ftail(q_full, DFR, W_full/q_full)
dis "P_POSTONLY_AME=" Ftail(q_post, DFR, W_post/q_post)
dis "P_ATT="          Ftail(q_att,  DFR, W_att/q_att)
dis "BITEQUAL_W=" (W_full == W_post) " BITEQUAL_q=" (q_full == q_post)
dis "REPORTED_ame_p_IN_LOG_ABOVE (grep 'single line')"

* ---- and what a user would think the ame pre-trend test did -----------------
dis "=== honestdid handoff ==="
matrix HB = BA
matrix HV = VA
local np = e(numpre_ame)
dis "NUMPRE_AME=" `np'
cap noisily honestdid, b(HB) vcov(HV) numpre(`np') mvec(0(0.5)1) parallel(0)
dis "HONEST_AME_RC=" _rc
matrix HB2 = BT
matrix HV2 = VT
local np2 = e(numpre_att)
cap noisily honestdid, b(HB2) vcov(HV2) numpre(`np2') mvec(0(0.5)1) parallel(0)
dis "HONEST_ATT_RC=" _rc

dis "=== ALLDONE ==="
