* ===========================================================================
* p04 -- INDEPENDENT REIMPLEMENTATION of the across-subgroup homogeneity test
*        with FOUR non-contiguous, value-labelled subgroup levels (1 3 7 99).
*
* METHOD.  wooldid reports, for each estimand:
*    - four level estimates b_v and their SEs s_v          (e(wdidsgresults_v))
*    - six pairwise contrasts and their SEs s_ab           (e(pwsgs))
*    - one "all SGs the same" p-value                      (e(joint_*_allsgssame))
* The six pairwise variances pin the six off-diagonal covariances EXACTLY:
*    C_ab = (s_a^2 + s_b^2 - s_ab^2)/2
* so the full 4x4 level covariance V is recoverable from the printed numbers.
* We then rebuild the Wald statistic of H0: b1=b3=b7=b99 ourselves and compare
* the p-value.  Disagreement is a defect in the joint test, the pairwise table,
* or the level table -- all three are checked against each other.
*
* Also checks: pwsgs row LABEL vs. DIRECTION (does "AvsB" mean b_A - b_B?),
* and the e() naming contract for non-contiguous subgroup values.
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

local gN 240
local gT 9
local gt0 2000
local gonsets "2003 2004 2006"
local gysd 0.60
local gtsd 0.50
local gaisd 1
local gtaumode hetero
local gtauhet 0.7
include "`GEN'/gen_core.do"

egen long grp = group(cohort region)
gen int sgq = 1
qui replace sgq = 3  if region == 2
qui replace sgq = 7  if region == 3
qui replace sgq = 99 if region == 4
cap label drop SGQ
label define SGQ 1 "North" 3 "South" 7 "East" 99 "Outlying territories"
label values sgq SGQ
qui gen double __tauold = truetau
qui replace truetau = truetau * (1 + 0.5*(sgq==3) + 1.2*(sgq==7) - 0.4*(sgq==99))
qui replace y = y - __tauold + truetau
drop __tauold

tab sgq, m

cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) subgroups(sgq) esprelength(2) espostlength(2)
dis "RUN_RC=" _rc

scalar DFR = e(df_r)
dis "DFR=" DFR

matrix PWX  = e(pwsgs)
matrix SGX1 = e(wdidsgresults_1)
matrix SGX3 = e(wdidsgresults_3)
matrix SGX7 = e(wdidsgresults_7)
matrix SGX99 = e(wdidsgresults_99)
scalar JP_att       = e(joint_att_allsgssame)
scalar JP_att_it    = e(joint_att_it_allsgssame)
scalar JP_att_i     = e(joint_att_i_allsgssame)
scalar JP_att_itime = e(joint_att_itime_allsgssame)
scalar JP_pre_att   = e(joint_pre_att_allsgssame)
dis "JP_att = " %20.17f JP_att
dis "JP_att_it = " %20.17f JP_att_it
dis "JP_att_i = " %20.17f JP_att_i
dis "JP_att_itime = " %20.17f JP_att_itime
dis "JP_pre_att = " %20.17f JP_pre_att

dis "=== e(b) colnames (NAMING CONTRACT) ==="
matrix EBX = e(b)
local ebn : colnames EBX
dis "`ebn'"

ereturn list

dis ""
dis "=== e(pwsgs) ==="
mat list PWX

* ---------------------------------------------------------------------------
* REIMPLEMENTATION (Mata)
* ---------------------------------------------------------------------------
mata:
  PW   = st_matrix("PWX")
  pws  = st_matrixrowstripe("PWX")     // col1 = eq, col2 = name
  lv   = ("1","3","7","99")
  SG   = J(4,0,.)
  SGb  = J(4,2,.)
  dfr  = st_numscalar("DFR")

  ests = ("att","att_it","att_i","att_itime")
  for (e=1; e<=4; e++) {
    nm = ests[e]
    printf("\n############ estimand: %s ############\n", nm)
    b = J(4,1,.)
    s = J(4,1,.)
    for (k=1; k<=4; k++) {
      M = st_matrix("SGX" + lv[k])
      b[k] = M[e,1]
      s[k] = M[e,2]
      printf("  level %s: b=%21.17f  se=%21.17f\n", lv[k], b[k], s[k])
    }
    V = diag(s:^2)
    npair = 0
    badrho = 0
    for (r=1; r<=rows(PW); r++) {
      if (pws[r,1] != nm) continue
      rn   = pws[r,2]
      core = subinstr(rn, ".sgq", "")
      pos  = strpos(core, "vs")
      A    = ustrregexra(substr(core, 1, pos-1), "[^0-9\-]", "")
      B    = ustrregexra(substr(core, pos+2), "[^0-9\-]", "")
      ia = 0
      ib = 0
      for (k=1; k<=4; k++) {
        if (lv[k] == A) ia = k
        if (lv[k] == B) ib = k
      }
      printf("  [dbg] eq=<%s> rn=<%s> core=<%s> pos=%f A=<%s> B=<%s> ia=%f ib=%f\n", pws[r,1], rn, core, pos, A, B, ia, ib)
      if (ia == 0 | ib == 0) continue
      pb = PW[r,1]
      ps = PW[r,2]
      want = b[ia] - b[ib]
      printf("  pair %s A=%s(%f) B=%s(%f) est=%21.17f bA-bB=%21.17f d=%10.3e se=%18.14f\n", rn, A, ia, B, ib, pb, want, pb-want, ps)
      cab = (s[ia]^2 + s[ib]^2 - ps^2)/2
      V[ia,ib] = cab
      V[ib,ia] = cab
      rho = cab/(s[ia]*s[ib])
      printf("      implied rho = %14.10f\n", rho)
      if (abs(rho) > 1) {
        printf("      !!!! IMPLIED CORRELATION OUT OF RANGE\n")
        badrho = 1
      }
      npair++
    }
    printf("  pairs found: %f (expect 6)\n", npair)
    if (npair != 6) continue
    printf("  eigenvalues of rebuilt V: ")
    ev = Re(eigenvalues(V))
    for (k=1; k<=4; k++) printf("%12.6e ", ev[k])
    printf("\n")

    R  = (1,-1,0,0 \ 1,0,-1,0 \ 1,0,0,-1)
    Rb = R*b
    W  = (Rb' * invsym(R*V*R') * Rb)
    pF = Ftail(3, dfr, W/3)
    pX = chi2tail(3, W)
    jp = st_numscalar("JP_" + nm)
    printf("  rebuilt Wald  = %18.12f\n", W)
    printf("  p (F 3,%f)  = %21.17f\n", dfr, pF)
    printf("  p (chi2 3)    = %21.17f\n", pX)
    printf("  wooldid joint = %21.17f\n", jp)
    printf("  dF = %10.3e   dX = %10.3e\n", pF-jp, pX-jp)
  }
end
