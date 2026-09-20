* ===========================================================================
* p05 -- DEGENERATE-LIMIT IDENTITY: a BLOCK-SEPARABLE triple-diff design.
*
* Construction: the subgroup sg partitions the units; each subgroup has its OWN
* onset cohorts (sg0 -> 2003,2005 ; sg1 -> 2004,2006) and its OWN never-treated
* units, and the run adds sg-specific time effects fe(i.sg##i.t).  The unit FEs
* are nested in sg.  Every column of the design matrix therefore belongs to
* exactly one subgroup, so the OLS problem is EXACTLY block diagonal in sg.
*
* Three things MUST then hold, and none of them is an approximation:
*  (S1) each subgroup's estimate == the estimate from running wooldid on that
*       subgroup's SUBSAMPLE alone (same FEs, same everything).
*  (S2) cluster(unit) with units nested in sg makes the cluster-robust
*       covariance block diagonal too, so cov(b_sg0, b_sg1) == 0 and the
*       contrast SE MUST equal sqrt(se0^2 + se1^2) exactly.
*  (S3) the "triple diff" the help promises == b_sg1 - b_sg0 == e(pwsgs).
*
* Also: (S4) suppresssgprimaryeffects must leave e(pwsgs) numerically unchanged.
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------------- purpose-built separable panel ----------------
local NU 200
local T0 2000
local TT 10
qui set obs `NU'
gen long unit = _n
gen byte sg = mod(unit,2)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if sg==0 & u < 0.30
qui replace ttre = 2005 if sg==0 & u >= 0.30 & u < 0.60
qui replace ttre = 2004 if sg==1 & u < 0.30
qui replace ttre = 2006 if sg==1 & u >= 0.30 & u < 0.60
gen double a_i = rnormal()
expand `TT'
bysort unit: gen int t = `T0' + _n - 1
* sg-SPECIFIC time effects (so fe(i.sg##i.t) is the true model)
gen double d_sgt = 0.4*sin(t/2) + 0.9*sg*cos(t/3)
gen byte post = (ttre < . & t >= ttre)
* subgroup-specific treatment effects, plus horizon variation
gen double truetau = 0
qui replace truetau = 1.0 + 0.6*sg + 0.15*(t-ttre) if post==1
gen double y = a_i + d_sgt + truetau + 0.5*rnormal()

tab sg ttre, m

* ---------------- run 1: pooled, with subgroups ----------------
cap noi wooldid y unit t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.sg##i.t) subgroups(sg) esprelength(2) espostlength(2)
dis "POOLED_RC=" _rc
matrix P_SG0 = e(wdidsgresults_0)
matrix P_SG1 = e(wdidsgresults_1)
matrix P_PW  = e(pwsgs)
matrix P_MAIN = e(wdidmainresults)
scalar P_N = e(N)
scalar P_DFR = e(df_r)
scalar P_JATT = e(joint_att_allsgssame)
dis "POOLED N=" P_N " dfr=" P_DFR

* ---------------- run 2: suppresssgprimaryeffects ----------------
cap noi wooldid y unit t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.sg##i.t) subgroups(sg) esprelength(2) espostlength(2) suppresssgprimaryeffects
dis "SUPPRESS_RC=" _rc
matrix S_PW = e(pwsgs)
scalar S_JATT = e(joint_att_allsgssame)

* ---------------- run 3/4: the two SUBSAMPLE runs ----------------
cap noi wooldid y unit t ttre if sg==0 , att att_it att_i att_itime cluster(unit) ///
    fe(i.t) esprelength(2) espostlength(2)
dis "SUB0_RC=" _rc
matrix A_SG0 = e(wdidmainresults)
scalar A0_N = e(N)
scalar A0_DFR = e(df_r)

cap noi wooldid y unit t ttre if sg==1 , att att_it att_i att_itime cluster(unit) ///
    fe(i.t) esprelength(2) espostlength(2)
dis "SUB1_RC=" _rc
matrix A_SG1 = e(wdidmainresults)
scalar A1_N = e(N)
scalar A1_DFR = e(df_r)

dis "N: pooled=" P_N "  sub0=" A0_N "  sub1=" A1_N "  sum=" A0_N+A1_N
dis "dfr: pooled=" P_DFR "  sub0=" A0_DFR "  sub1=" A1_DFR

* ---------------------------------------------------------------------------
mata:
  P0 = st_matrix("P_SG0")
  P1 = st_matrix("P_SG1")
  A0 = st_matrix("A_SG0")
  A1 = st_matrix("A_SG1")
  PW = st_matrix("P_PW")
  SW = st_matrix("S_PW")
  pws = st_matrixrowstripe("P_PW")
  sws = st_matrixrowstripe("S_PW")
  ests = ("att","att_it","att_i","att_itime")

  printf("\n===== S1: subgroup estimate vs SUBSAMPLE run =====\n")
  printf("%-12s %22s %22s %11s | %22s %22s %11s\n","estimand","sg0 in-pooled","sg0 subsample","reldiff","sg1 in-pooled","sg1 subsample","reldiff")
  for (e=1; e<=8; e++) {
    lbl = (e<=4 ? ests[e] : "pre_"+ests[e-4])
    d0 = (P0[e,1]-A0[e,1])/abs(A0[e,1])
    d1 = (P1[e,1]-A1[e,1])/abs(A1[e,1])
    printf("%-12s %22.17f %22.17f %11.2e | %22.17f %22.17f %11.2e\n", lbl, P0[e,1], A0[e,1], d0, P1[e,1], A1[e,1], d1)
  }
  printf("\n--- and the SEs ---\n")
  for (e=1; e<=8; e++) {
    lbl = (e<=4 ? ests[e] : "pre_"+ests[e-4])
    d0 = (P0[e,2]-A0[e,2])/abs(A0[e,2])
    d1 = (P1[e,2]-A1[e,2])/abs(A1[e,2])
    printf("%-12s se0 %18.14f vs %18.14f  rd=%11.2e | se1 %18.14f vs %18.14f  rd=%11.2e\n", lbl, P0[e,2], A0[e,2], d0, P1[e,2], A1[e,2], d1)
  }

  printf("\n===== S2/S3: contrast, and cov==0 under block separability =====\n")
  for (r=1; r<=rows(PW); r++) {
    eq = pws[r,1]
    e = 0
    for (k=1; k<=4; k++) if (ests[k]==eq) e = k
    if (e == 0) continue
    s0 = P0[e,2]
    s1 = P1[e,2]
    sc = PW[r,2]
    indep = sqrt(s0^2 + s1^2)
    cov = (s0^2 + s1^2 - sc^2)/2
    rho = cov/(s0*s1)
    printf("%-10s %-14s est=%20.16f  b1-b0=%20.16f  d=%10.2e\n", eq, pws[r,2], PW[r,1], P1[e,1]-P0[e,1], PW[r,1]-(P1[e,1]-P0[e,1]))
    printf("            se_c=%18.14f  sqrt(s0^2+s1^2)=%18.14f  reldiff=%10.2e  implied rho=%12.8f\n", sc, indep, (sc-indep)/indep, rho)
  }

  printf("\n===== S4: suppresssgprimaryeffects leaves e(pwsgs) alone =====\n")
  printf("rows: full=%f suppressed=%f\n", rows(PW), rows(SW))
  for (r=1; r<=min((rows(PW),rows(SW))); r++) {
    printf("%-10s %-14s | full est=%20.16f se=%18.14f p=%18.14f\n", pws[r,1], pws[r,2], PW[r,1], PW[r,2], PW[r,4])
    printf("%-10s %-14s | supp est=%20.16f se=%18.14f p=%18.14f   d_est=%10.2e d_se=%10.2e\n", sws[r,1], sws[r,2], SW[r,1], SW[r,2], SW[r,4], SW[r,1]-PW[r,1], SW[r,2]-PW[r,2])
  }
  printf("joint att: full=%21.17f suppressed=%21.17f\n", st_numscalar("P_JATT"), st_numscalar("S_JATT"))
end
