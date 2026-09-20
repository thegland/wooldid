* ===========================================================================
* p06 -- SUBGROUP EVENT STUDY, PRE-TREND TESTS, AND THE e(V) SPLICE.
*
*  (E1) per-subgroup pre-trend joint test: rebuild W = b_pre' V_pre^-1 b_pre from
*       e(b_es_att_s<v>) / e(V_es_att_s<v>) and compare p = Ftail(npre,dfr,W/npre)
*       with e(joint_horiz_s<v>).  Same for the pooled e(joint_horiz).
*  (E2) e(V) splice: the ES block for subgroup v inside e(V) must equal
*       e(V_es_att_s<v>) cell for cell (bydesign PS-62 says the ES blocks are
*       spliced in; this checks the SUBGROUP blocks actually were).
*  (E3) ES aggregation: pooled ES_att at horizon h == the weight-weighted mean of
*       the subgroup ES_att at h, weights recomputed from the data.
*  (E4) the per-subgroup ES rows in e(wdidsgresults_<v>) == e(b_es_att_s<v>).
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
local gT 11
local gt0 2000
local gonsets "2005 2006 2008"
local gysd 0.60
local gtsd 0.50
local gaisd 1
local gtaumode hetero
local gtauhet 0.7
include "`GEN'/gen_core.do"

egen long grp = group(cohort region)
gen byte sgb = (region <= 2)

cap noi wooldid y grp t ttre , att att_it att_i att_itime cluster(unit) ///
    fe(i.unit) subgroups(sgb) esprelength(3) espostlength(2) jointtests
dis "RUN_RC=" _rc

qui gen byte insmp = e(sample)
scalar DFR = e(df_r)
matrix VV   = e(V)
matrix BB   = e(b)
matrix MAIN = e(wdidmainresults)
matrix SGR0 = e(wdidsgresults_0)
matrix SGR1 = e(wdidsgresults_1)
matrix BES   = e(b_es_att)
matrix VES   = e(V_es_att)
matrix BES0  = e(b_es_att_s0)
matrix VES0  = e(V_es_att_s0)
matrix BES1  = e(b_es_att_s1)
matrix VES1  = e(V_es_att_s1)
matrix BESI  = e(b_es_att_it)
matrix VESI  = e(V_es_att_it)
matrix BESI0 = e(b_es_att_it_s0)
matrix VESI0 = e(V_es_att_it_s0)
matrix BESI1 = e(b_es_att_it_s1)
matrix VESI1 = e(V_es_att_it_s1)
matrix ESC   = e(escomp)
matrix ESC0  = e(escomp_s0)
matrix ESC1  = e(escomp_s1)
scalar JH    = e(joint_horiz)
scalar JH0   = e(joint_horiz_s0)
scalar JH1   = e(joint_horiz_s1)
scalar JHI   = e(joint_horiz_i)
scalar JHI0  = e(joint_horiz_i_s0)
scalar JHI1  = e(joint_horiz_i_s1)
scalar NP    = e(numpre_att)
scalar NP0   = e(numpre_att_s0)
scalar NP1   = e(numpre_att_s1)

dis "DFR=" DFR "  numpre=" NP " / " NP0 " / " NP1
dis "joint_horiz    =" %20.17f JH
dis "joint_horiz_s0 =" %20.17f JH0
dis "joint_horiz_s1 =" %20.17f JH1
dis "joint_horiz_i    =" %20.17f JHI
dis "joint_horiz_i_s0 =" %20.17f JHI0
dis "joint_horiz_i_s1 =" %20.17f JHI1

dis "=== escomp colnames ==="
local ecn : colnames ESC
dis "`ecn'"
mat list ESC
mat list ESC0
mat list ESC1
mat list BES
mat list BES0
mat list BES1

* ---- independently recompute the per-horizon att weights from the data ----
cap drop __*
qui gen int __h = t - ttre if insmp & ttre < .
qui gen double __one = 1 if insmp & ttre < .
levelsof __h if insmp & ttre < . & __h >= -3 & __h <= 2, local(HLIST)
matrix HW = J(20,4,.)
local hr 0
foreach h of local HLIST {
  local hr = `hr' + 1
  qui summarize __one if insmp & __h == `h' & sgb==0, meanonly
  local w0 = r(sum)
  qui summarize __one if insmp & __h == `h' & sgb==1, meanonly
  local w1 = r(sum)
  matrix HW[`hr',1] = `h'
  matrix HW[`hr',2] = `w0'
  matrix HW[`hr',3] = `w1'
  matrix HW[`hr',4] = `w0'+`w1'
}
mat list HW

* ---------------------------------------------------------------------------
mata:
  dfr = st_numscalar("DFR")

  // ---------- E1 : pre-trend joint tests ----------
  printf("\n===== E1: pre-trend joint test rebuilt from e(b_es)/e(V_es) =====\n")
  fam = ("att",    "att_s0",  "att_s1",  "att_it",  "att_it_s0", "att_it_s1")
  sc  = ("JH",     "JH0",     "JH1",     "JHI",     "JHI0",      "JHI1")
  bm  = ("BES",    "BES0",    "BES1",    "BESI",    "BESI0",     "BESI1")
  vm  = ("VES",    "VES0",    "VES1",    "VESI",    "VESI0",     "VESI1")
  for (f=1; f<=6; f++) {
    b  = st_matrix(bm[f])
    V  = st_matrix(vm[f])
    cs = st_matrixcolstripe(bm[f])
    // pre columns only
    idx = J(1,0,.)
    for (c=1; c<=cols(b); c++) if (strpos(cs[c,2],"pre")==1) idx = (idx, c)
    bp = b[1,idx]'
    Vp = V[idx,idx]
    np = cols(idx)
    W  = bp' * invsym(Vp) * bp
    pF = Ftail(np, dfr, W/np)
    pX = chi2tail(np, W)
    jp = st_numscalar(sc[f])
    printf("%-12s npre=%f  W=%16.10f  pF=%21.17f  pX=%21.17f  posted=%21.17f  dF=%10.2e dX=%10.2e\n", fam[f], np, W, pF, pX, jp, pF-jp, pX-jp)
  }

  // ---------- E2 : e(V) splice for the SUBGROUP ES blocks ----------
  printf("\n===== E2: e(V) subgroup ES block vs e(V_es_att_s<v>) =====\n")
  Vfull = st_matrix("VV")
  cs    = st_matrixcolstripe("VV")
  nmv   = cs[.,2]
  for (g=0; g<=1; g++) {
    V = st_matrix("VES" + strofreal(g))
    cse = st_matrixcolstripe("VES" + strofreal(g))
    k = cols(V)
    maxd = 0
    nfound = 0
    for (a=1; a<=k; a++) {
      for (bq=1; bq<=k; bq++) {
        na = "ES_att_SG" + strofreal(g) + "_" + cse[a,2]
        nb = "ES_att_SG" + strofreal(g) + "_" + cse[bq,2]
        ia = 0
        ib = 0
        for (c=1; c<=rows(nmv); c++) {
          if (nmv[c]==na) ia = c
          if (nmv[c]==nb) ib = c
        }
        if (ia==0 | ib==0) continue
        nfound++
        d = abs(Vfull[ia,ib] - V[a,bq])
        if (d > maxd) maxd = d
      }
    }
    printf("  sg%f: cells matched=%f of %f, max |e(V) - e(V_es)| = %10.3e\n", g, nfound, k*k, maxd)
    printf("      e(V_es_att_s%f) colnames: ", g)
    for (a=1; a<=k; a++) printf("%s ", cse[a,2])
    printf("\n")
  }

  // ---------- E3 : ES aggregation across subgroups ----------
  printf("\n===== E3: pooled ES_att horizon == weighted mean of subgroup ES =====\n")
  B  = st_matrix("BES")
  B0 = st_matrix("BES0")
  B1 = st_matrix("BES1")
  csP = st_matrixcolstripe("BES")
  cs0 = st_matrixcolstripe("BES0")
  cs1 = st_matrixcolstripe("BES1")
  HW  = st_matrix("HW")
  ESC0 = st_matrix("ESC0")
  ESC1 = st_matrix("ESC1")
  // escomp: last column is relyr, col5 = sumw_att
  for (c=1; c<=cols(B); c++) {
    nmP = csP[c,2]
    // relative period from the name
    if (strpos(nmP,"pre")==1) h = -strtoreal(substr(nmP,4))
    else h = strtoreal(substr(nmP,5))
    i0 = 0
    i1 = 0
    for (a=1; a<=cols(B0); a++) if (cs0[a,2]==nmP) i0 = a
    for (a=1; a<=cols(B1); a++) if (cs1[a,2]==nmP) i1 = a
    if (i0==0 | i1==0) {
      printf("  %-8s : not present in both subgroups (i0=%f i1=%f)\n", nmP, i0, i1)
      continue
    }
    // weights from my own recomputation
    w0 = .
    w1 = .
    for (r=1; r<=rows(HW); r++) {
      if (HW[r,1]==h) {
        w0 = HW[r,2]
        w1 = HW[r,3]
      }
    }
    agg = (w0*B0[1,i0] + w1*B1[1,i1])/(w0+w1)
    // weights from escomp
    e0 = .
    e1 = .
    for (r=1; r<=rows(ESC0); r++) if (ESC0[r,cols(ESC0)]==h) e0 = ESC0[r,5]
    for (r=1; r<=rows(ESC1); r++) if (ESC1[r,cols(ESC1)]==h) e1 = ESC1[r,5]
    agg2 = (e0*B0[1,i0] + e1*B1[1,i1])/(e0+e1)
    printf("  %-8s h=%3.0f pooled=%20.16f  wmean(data)=%20.16f rd=%10.2e | wmean(escomp)=%20.16f rd=%10.2e   [w0=%g w1=%g e0=%g e1=%g]\n", nmP, h, B[1,c], agg, (agg-B[1,c])/abs(B[1,c]), agg2, (agg2-B[1,c])/abs(B[1,c]), w0, w1, e0, e1)
  }

  // ---------- E4 : subgroup results table vs e(b_es_att_s<v>) ----------
  printf("\n===== E4: e(wdidsgresults_v) ES rows vs e(b_es_att_s<v>) =====\n")
  for (g=0; g<=1; g++) {
    R  = st_matrix("SGR" + strofreal(g))
    rs = st_matrixrowstripe("SGR" + strofreal(g))
    Bg = st_matrix("BES" + strofreal(g))
    cg = st_matrixcolstripe("BES" + strofreal(g))
    for (a=1; a<=cols(Bg); a++) {
      nmg = cg[a,2]
      ir = 0
      for (r=1; r<=rows(R); r++) if (rs[r,1]=="ES_att" & rs[r,2]==nmg) ir = r
      if (ir==0) {
        printf("  sg%f %-8s : NOT FOUND in wdidsgresults\n", g, nmg)
        continue
      }
      printf("  sg%f %-8s table=%20.16f  b_es=%20.16f  d=%10.2e | table_se=%18.14f  sqrt(V_es_diag)=%18.14f  d=%10.2e\n", g, nmg, R[ir,1], Bg[1,a], R[ir,1]-Bg[1,a], R[ir,2], sqrt(st_matrix("VES"+strofreal(g))[a,a]), R[ir,2]-sqrt(st_matrix("VES"+strofreal(g))[a,a]))
    }
  }
end
