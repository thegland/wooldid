* ===========================================================================
* p08 -- RAGGED SUBGROUP EVENT-STUDY HORIZONS, with an exact LABEL oracle.
*
* The code carries a positional filter for the case the comment calls
* "not all sgs have all relyrs, necessarily" (wooldid_sglbl).  A positional
* filter that slips would mislabel a subgroup's event-study row -- a silent
* wrong number, not a crash.  This probe makes that detectable:
*
*   - NOISELESS DGP (no error term), so every estimate is exact.
*   - the true effect at relative period h is  base(sg) + slope(sg)*h,
*     the SAME for every cohort inside a subgroup.  So the event-study row
*     labelled post<h> for subgroup v MUST come back as base(v)+slope(v)*h
*     to machine precision.  Any label slip shows up as a whole-number-sized
*     error, not a rounding one.
*   - the two subgroups are given DISJOINT onset windows so that sg1 cannot
*     reach the longest post horizon that sg0 can: the ragged case.
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------------- purpose-built ragged panel ----------------
qui set obs 200
gen long unit = _n
gen byte sg = (unit > 100)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if sg==0 & u < 0.35
qui replace ttre = 2004 if sg==0 & u >= 0.35 & u < 0.70
qui replace ttre = 2007 if sg==1 & u < 0.35
qui replace ttre = 2008 if sg==1 & u >= 0.35 & u < 0.70
gen double a_i = rnormal()
expand 11
bysort unit: gen int t = 2000 + _n - 1
gen double g_t = 0.7*sin(t/2)
gen int h = t - ttre
gen byte post = (ttre < . & t >= ttre)
* the true effect depends ONLY on the relative period and the subgroup
gen double truetau = 0
qui replace truetau = 1.00 + 0.10*h if post==1 & sg==0
qui replace truetau = 2.00 - 0.25*h if post==1 & sg==1
gen double y = a_i + g_t + truetau

tab sg ttre, m
tab h sg if post==1

cap noi wooldid y unit t ttre , att att_it cluster(unit) subgroups(sg) ///
    esprelength(4) espostlength(4)
dis "RUN_RC=" _rc

matrix MAIN = e(wdidmainresults)
matrix SG0  = e(wdidsgresults_0)
matrix SG1  = e(wdidsgresults_1)
matrix BES  = e(b_es_att)
matrix BES0 = e(b_es_att_s0)
matrix BES1 = e(b_es_att_s1)
matrix ESC  = e(escomp)
matrix ESC0 = e(escomp_s0)
matrix ESC1 = e(escomp_s1)
scalar NP  = e(numpre_att)
scalar NP0 = e(numpre_att_s0)
scalar NP1 = e(numpre_att_s1)
dis "numpre: pooled=" NP " sg0=" NP0 " sg1=" NP1

dis "=== e(b) colnames ==="
matrix EBX = e(b)
local ebn : colnames EBX
dis "`ebn'"

dis "=== b_es pooled / sg0 / sg1 ==="
mat list BES
mat list BES0
mat list BES1
dis "=== escomp pooled / sg0 / sg1 (last col = relyr) ==="
mat list ESC
mat list ESC0
mat list ESC1
dis "=== the two subgroup results tables ==="
mat list SG0
mat list SG1

ereturn list

* ---- data-side weights per (horizon, subgroup) ----
qui gen byte insmp = 1
matrix HW = J(30,3,.)
local hr 0
forval hh = -6/4 {
  local hr = `hr' + 1
  qui count if post==1 & h == `hh' & sg==0
  local w0 = r(N)
  qui count if post==1 & h == `hh' & sg==1
  local w1 = r(N)
  qui count if h == `hh' & sg==0 & ttre < .
  local a0 = r(N)
  qui count if h == `hh' & sg==1 & ttre < .
  local a1 = r(N)
  matrix HW[`hr',1] = `hh'
  matrix HW[`hr',2] = `a0'
  matrix HW[`hr',3] = `a1'
}
mat list HW

mata:
  printf("\n===== LABEL ORACLE: subgroup ES row post<h> must equal base+slope*h =====\n")
  base = (1.00, 2.00)
  slope = (0.10, -0.25)
  for (g=0; g<=1; g++) {
    B  = st_matrix("BES" + strofreal(g))
    cs = st_matrixcolstripe("BES" + strofreal(g))
    printf("-- subgroup %f --\n", g)
    for (c=1; c<=cols(B); c++) {
      nm = cs[c,2]
      if (strpos(nm,"post")==1) {
        hh = strtoreal(substr(nm,5))
        want = base[g+1] + slope[g+1]*hh
        printf("   %-8s h=%3.0f   wooldid=%20.16f   truth=%20.16f   d=%10.3e\n", nm, hh, B[1,c], want, B[1,c]-want)
      }
      else {
        hh = -strtoreal(substr(nm,4))
        printf("   %-8s h=%3.0f   wooldid=%20.16f   truth=%20.16f   d=%10.3e\n", nm, hh, B[1,c], 0, B[1,c]-0)
      }
    }
  }

  printf("\n===== POOLED ES vs subgroup-weighted mean (ragged: post4 only in sg0) =====\n")
  B  = st_matrix("BES")
  cs = st_matrixcolstripe("BES")
  B0 = st_matrix("BES0")
  c0 = st_matrixcolstripe("BES0")
  B1 = st_matrix("BES1")
  c1 = st_matrixcolstripe("BES1")
  HW = st_matrix("HW")
  for (c=1; c<=cols(B); c++) {
    nm = cs[c,2]
    if (strpos(nm,"post")==1) hh = strtoreal(substr(nm,5))
    else hh = -strtoreal(substr(nm,4))
    i0 = 0
    i1 = 0
    for (a=1; a<=cols(B0); a++) if (c0[a,2]==nm) i0 = a
    for (a=1; a<=cols(B1); a++) if (c1[a,2]==nm) i1 = a
    w0 = 0
    w1 = 0
    for (r=1; r<=rows(HW); r++) {
      if (HW[r,1]==hh) {
        w0 = HW[r,2]
        w1 = HW[r,3]
      }
    }
    if (i0==0) w0 = 0
    if (i1==0) w1 = 0
    num = 0
    if (i0>0) num = num + w0*B0[1,i0]
    if (i1>0) num = num + w1*B1[1,i1]
    agg = num/(w0+w1)
    printf("  %-8s h=%3.0f  pooled=%20.16f  wmean=%20.16f  d=%10.3e   [in sg0=%f w0=%g | in sg1=%f w1=%g]\n", nm, hh, B[1,c], agg, agg-B[1,c], i0, w0, i1, w1)
  }
end
