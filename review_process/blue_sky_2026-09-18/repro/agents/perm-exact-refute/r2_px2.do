* r2_px2.do -- SKEPTIC probe of PERMEX-PX-2: can wooldid be driven into
*   exk = 1 (complete enumeration) AND memb = 0 (observed row not in the reference set)?
* If it can, the ado:8297 caution ("their p-values use the add-one form") fires while
* ado:15992 uses the PLAIN share, and a p of exactly 0 becomes attainable.
* Three triggers, plus a control that proves memb detection works from the command line.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- reporter: is the rep-0 comparator statistic bitwise present in the draws?
prog def membchk
  syntax , tag(string)
  cap confirm matrix e(permdist_b)
  if _rc != 0 {
    dis "`tag'_NOPERMDIST"
    exit
  }
  tempname PB CB
  matrix `PB' = e(permdist_b)
  matrix `CB' = e(perm_comparator_b)
  local cn : colnames `CB'
  dis "`tag'_comparator_cols: `cn'"
  matrix list `CB', format(%21x)
  mata: wdp_membchk("`PB'", "`CB'", "`tag'")
end
mata:
void wdp_membchk(string scalar pbn, string scalar cbn, string scalar tag)
{
  real matrix P
  real rowvector c
  real scalar r, k, nhit, kk, j
  string rowvector pn, cnames
  P = st_matrix(pbn)
  c = st_matrix(cbn)
  pn = st_matrixcolstripe(pbn)[., 2]'
  cnames = st_matrixcolstripe(cbn)[., 2]'
  for (k = 1; k <= cols(c); k++) {
    kk = 0
    for (j = 1; j <= cols(P); j++) {
      if (pn[j] == cnames[k]) kk = j
    }
    if (kk == 0) {
      printf("%s_NOCOL %s\n", tag, cnames[k])
    }
    else {
      nhit = 0
      for (r = 1; r <= rows(P); r++) {
        if (P[r, kk] == c[k]) nhit++
      }
      printf("%s_MEMB col=%s hits=%f nonmiss=%f rows=%f\n", tag, cnames[k], nhit, sum(P[., kk] :< .), rows(P))
      printf("%s_CMP col=%s cmp=%21x\n", tag, cnames[k], c[k])
    }
  }
}
end

* =====================================================================
* T1: permonsetcohorts + permexact, cluster STRICTLY finer than i
*     i = cohort (the onset grouping), cluster = unit.  6 units -> M = 6!/(2!1!3!) = 60
* =====================================================================
clear
set obs 6
gen long unit = _n
gen double onset = .
replace onset = 4 in 1
replace onset = 4 in 2
replace onset = 6 in 3
gen byte coh = 3
replace coh = 1 if onset == 4
replace coh = 2 if onset == 6
gen double ue = rnormal()*0.5
expand 8
bysort unit: gen int year = _n
gen double te = rnormal()*0.3
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.15*year + ue + te + 0.9*post
tempfile d1
save `d1', replace

dis "===== T1: permonsetcohorts + permexact ====="
capture noisily wooldid y coh year onset, att espostlength(2) cluster(unit) permute(shuffle) permexact permonsetcohorts permstepdown dots(0)
dis "T1_RC=" _rc
dis "T1_perm_exact=" e(perm_exact)
dis "T1_perm_irebuilt=" e(perm_irebuilt)
dis "T1_perm_reps=" e(perm_reps)
cap confirm matrix e(perm_n)
if _rc==0 matrix list e(perm_n)
cap confirm matrix e(perm_p)
if _rc==0 matrix list e(perm_p), format(%20.17f)
membchk , tag(T1)

* =====================================================================
* T2: permtstat + permexact on the same design (statistic is b/se; s_obs and the
*     identity draw's s come from two separate fits, so an ulp would show here)
* =====================================================================
use `d1', clear
dis "===== T2: permtstat + permexact ====="
capture noisily wooldid y coh year onset, att espostlength(2) cluster(unit) permute(shuffle) permexact permonsetcohorts permtstat dots(0)
dis "T2_RC=" _rc
dis "T2_perm_exact=" e(perm_exact)
cap confirm matrix e(perm_p)
if _rc==0 matrix list e(perm_p), format(%20.17f)
membchk , tag(T2)

* =====================================================================
* T3: control -- INCOMPLETE enumeration (exk = 0) with the observed row complete
*     (memb = 1).  Unbalanced panel: unit 5 ends at year 5, so assignments that give
*     it the early onset make the late horizons inestimable.
* =====================================================================
clear
set obs 5
gen long unit = _n
gen double onset = .
replace onset = 4 in 1
replace onset = 7 in 2
gen double ue = rnormal()*0.5
expand 8
bysort unit: gen int year = _n
drop if unit == 5 & year > 5
gen double te = rnormal()*0.3
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.15*year + ue + te + 0.9*post
dis "===== T3: incomplete enumeration control ====="
capture noisily wooldid y unit year onset, att espostlength(4) cluster(unit) permute(shuffle) permexact dots(0)
dis "T3_RC=" _rc
dis "T3_perm_exact=" e(perm_exact)
dis "T3_perm_reps=" e(perm_reps)
cap confirm matrix e(perm_n)
if _rc==0 matrix list e(perm_n)
cap confirm matrix e(perm_p)
if _rc==0 matrix list e(perm_p), format(%20.17f)
membchk , tag(T3)
dis "DONE_R2"
