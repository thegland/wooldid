* ============================================================================
* P5 -- the adversarial ame p-value on a covariance whose pre-treatment block
*       is structurally degenerate.  Does the reported "adversarial pre-trend"
*       p use any pre-treatment information at all?
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

cap noisily wooldid y coh t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) contreat(dose) contreatcontrols(untreatedZero) adversarial jointtests
dis "EST_RC=" _rc

dis "=== what e() records about the adversarial test ==="
cap noisily matrix list e(adversarial_line_att)
cap noisily matrix list e(adversarial_line_ame)
foreach s in adversarialp_att adversarialp_ame adv_q_att adv_q_ame {
  dis "E_`s'=" e(`s')
}
dis "--- ereturn list, all (adversarial-related names) ---"
ereturn list, all

dis "=== the raw inputs ==="
matrix BA = e(b_es_ame)
matrix VA = e(V_es_ame)
matrix BT = e(b_es_att)
matrix VT = e(V_es_att)
matrix list BA
matrix list VA, format(%12.0g)
local dfr = e(df_r)
dis "DFR=" `dfr'

* relyr vector for the 6 horizons: -4 -3 -2 0 1 2 ; esrelativeto = -1
matrix REL = (-4\-3\-2\0\1\2)
matrix RELP = (0\1\2)
matrix BAP = (BA[1,4], BA[1,5], BA[1,6])
matrix VAP = VA[4..6,4..6]

dis "=== call wooldid's OWN kernel: FULL 6x6 vs POST-ONLY 3x3 ==="
mata:
  qf = .
  Wf = wooldid_advW(st_matrix("BA")', st_matrix("VA"), st_matrix("REL"), -1, 0, 0, qf)
  st_local("Wfull", strofreal(Wf, "%21x"))
  st_local("qfull", strofreal(qf))
  qp = .
  Wp = wooldid_advW(st_matrix("BAP")', st_matrix("VAP"), st_matrix("RELP"), -1, 0, 0, qp)
  st_local("Wpost", strofreal(Wp, "%21x"))
  st_local("qpost", strofreal(qp))
  qt = .
  Wt = wooldid_advW(st_matrix("BT")', st_matrix("VT"), st_matrix("REL"), -1, 0, 0, qt)
  st_local("Watt", strofreal(Wt, "%21x"))
  st_local("qatt", strofreal(qt))
end
dis "W_FULL_AME=`Wfull'  q=`qfull'"
dis "W_POSTONLY_AME=`Wpost'  q=`qpost'"
dis "W_ATT=`Watt'  q=`qatt'"
dis "IDENTICAL_FULL_VS_POSTONLY=" ("`Wfull'" == "`Wpost'" & "`qfull'" == "`qpost'")
local pf = Ftail(`qfull', `dfr', (real("`Wfull'"))/`qfull')
dis "P_FROM_FULL=" `pf'
* reconstruct the reported p directly
mata:
  q1 = .
  W1 = wooldid_advW(st_matrix("BA")', st_matrix("VA"), st_matrix("REL"), -1, 0, 0, q1)
  st_numscalar("WW", W1)
  st_numscalar("QQ", q1)
end
dis "REPORTED_STYLE_P=" Ftail(QQ, `dfr', WW/QQ)

dis "=== honestdid handoff on the degenerate ame covariance ==="
matrix HB = e(b_es_ame)
matrix HV = e(V_es_ame)
local np = e(numpre_ame)
dis "NUMPRE=" `np'
cap noisily honestdid, b(HB) vcov(HV) numpre(`np') mvec(0(0.5)1) parallel(0)
dis "HONEST_AME_RC=" _rc
matrix HB2 = e(b_es_att)
matrix HV2 = e(V_es_att)
local np2 = e(numpre_att)
cap noisily honestdid, b(HB2) vcov(HV2) numpre(`np2') mvec(0(0.5)1) parallel(0)
dis "HONEST_ATT_RC=" _rc

dis "=== ALLDONE ==="
