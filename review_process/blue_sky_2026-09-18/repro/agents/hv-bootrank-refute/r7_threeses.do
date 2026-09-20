*=======================================================================
* r7_threeses.do -- how many different standard errors does ONE wooldid
* bootstrap run carry for the SAME coefficient?  Printed table, e(V)
* (listwise over every column), e(V_boot), and e(V_es_att) (listwise over
* that family's columns only, the honestdid entry point).  Plus a
* no-subgroups CONTROL on the same panel, same reps, same seed.
* SELF-CONTAINED: needs gen_skep.do beside it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/gen_skep.do"

mata:
void wd_three(string scalar tag)
{
  real matrix V, VB, VE
  real rowvector d, db
  string matrix cn, rn
  real scalar k
  V  = st_matrix("e(V)")
  cn = st_matrixcolstripe("e(V)")
  d  = sqrt(diagonal(V))'
  VB = st_matrix("e(V_boot)")
  db = sqrt(diagonal(VB))'
  printf("%s_NV=%f  K=%f\n", tag, st_numscalar("e(boot_N_V)"), cols(V))
  for (k=1; k<=cols(V); k++) {
     if (cn[k,2] == "main_att" | substr(cn[k,2],1,10) == "ES_att_pos" | substr(cn[k,2],1,10) == "ES_att_pre") {
        printf("%s_EV %-22s sqrt_diag_eV=%10.6f  sqrt_diag_eVboot=%10.6f\n", tag, cn[k,2], d[k], db[k])
     }
  }
  VE = st_matrix("e(V_es_att)")
  if (rows(VE) > 0) {
     rn = st_matrixrowstripe("e(V_es_att)")
     for (k=1; k<=rows(VE); k++) {
        printf("%s_VES %-22s sqrt_diag_eVes_att=%10.6f\n", tag, rn[k,2], sqrt(VE[k,k]))
     }
  }
}
end

* ---------------- A: with subgroups ----------------
cap noisily wooldid y state t ttre, att esfixedbaseperiod esprelength(5) ///
    espostlength(10) cluster(state) subgroups(reg) bootstrap(standard) ///
    reps(60) rseed(6201)
dis "A_RC=" _rc
dis "A_N_V=" e(boot_N_V)
dis "A_COMPLETED=" e(boot_completed)
mata: wd_three("A")
dis "=== printed table A ==="
cap noisily mat list e(wdidmainresults), format(%12.6f)

* ---------------- B: CONTROL, identical but no subgroups() ----------------
cap noisily wooldid y state t ttre, att esfixedbaseperiod esprelength(5) ///
    espostlength(10) cluster(state) bootstrap(standard) ///
    reps(60) rseed(6201)
dis "B_RC=" _rc
dis "B_N_V=" e(boot_N_V)
dis "B_COMPLETED=" e(boot_completed)
mata: wd_three("B")
dis "=== printed table B ==="
cap noisily mat list e(wdidmainresults), format(%12.6f)
dis "R7_DONE"
