*=======================================================================
* r2_pilot.do -- skeptic's pilot: reproduce (or not) the listwise/per-column
* divergence on MY OWN data, verify the arithmetic offline from
* e(bootdist), and time the run so the big version can be sized.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/gen_skep.do"

timer clear 1
timer on 1
cap noisily wooldid y state t ttre, att esfixedbaseperiod esprelength(5) ///
    espostlength(10) cluster(state) subgroups(reg) bootstrap(standard) ///
    reps(30) rseed(9001)
dis "P_RC=" _rc
timer off 1
timer list 1
dis "P_COMPLETED=" e(boot_completed)
dis "P_FAILED=" e(boot_failed)
dis "P_N_V=" e(boot_N_V)
dis "P_VPOSTED=" e(boot_v_posted)
dis "P_DFR=" e(df_r)
dis "P_VCVSTATE=" e(wdid_vcvstate)

cap confirm matrix e(bootdist)
dis "P_BOOTDIST_PRESENT=" (_rc==0)
tempname BD VV NAV BM
cap matrix `BD' = e(bootdist)
dis "P_BD_ROWS=" rowsof(`BD') " COLS=" colsof(`BD')
matrix `VV' = e(V)
matrix `NAV' = e(boot_navail)
dis "P_K=" colsof(`VV')

* ---- offline: redo BOTH moment passes from e(bootdist) and compare ----
mata:
  A   = st_matrix("`BD'")
  V   = st_matrix("`VV'")
  nav = st_matrix("`NAV'")
  cs  = st_matrixcolstripe("`BD'")
  K   = cols(A)
  B   = rows(A)
  printf("OFF_B=%f K=%f\n", B, K)
  // per-column sd over that column's own available rows
  se_pc = J(1,K,.)
  cnt   = J(1,K,0)
  for (k=1;k<=K;k++) {
     u = selectindex(A[.,k] :< .)
     cnt[k] = length(u)
     if (length(u) >= 2) se_pc[k] = sqrt(variance(A[u,k]))
  }
  printf("OFF_NAVAIL_MIN=%f MAX=%f\n", min(cnt), max(cnt))
  printf("OFF_NAVAIL_MATCHES_E=%f\n", all(cnt :== nav))
  // listwise complete rows over ALL columns (no exclusions in this run, checked below)
  ok  = rowmissing(A) :== 0
  nlw = sum(ok)
  printf("OFF_NLW=%f  (e(boot_N_V)=%s)\n", nlw, st_global("e(boot_N_V)"))
  CC  = select(A, ok)
  dev = CC :- mean(CC)
  Vlw = (dev'*dev)/(nlw-1)
  d_e = sqrt(diagonal(V))'
  d_o = sqrt(diagonal(Vlw))'
  fin = selectindex((d_e :< .) :* (d_e :> 0))
  printf("OFF_MAXRELDIFF_DIAG=%g over %f finite cols\n",
         max(abs(d_e[fin] :/ d_o[fin] :- 1)), cols(fin))
  // the ratio that matters: exported sd / printed sd, per column
  r = d_o :/ se_pc
  rs = sort(r[fin]', 1)
  nf = rows(rs)
  printf("OFF_RATIO n=%f  min=%6.4f  p25=%6.4f  med=%6.4f  p75=%6.4f  max=%6.4f\n",
         nf, rs[1], rs[max((1,ceil(.25*nf)))], rs[max((1,ceil(.5*nf)))],
         rs[max((1,ceil(.75*nf)))], rs[nf])
  printf("OFF_RATIO_BELOW1=%f of %f\n", sum(rs :< 1), nf)
  st_matrix("RAT", r)
  st_matrix("SEPC", se_pc)
  st_matrix("SELW", d_o)
  st_matrix("CNT", cnt)
  // print the worst 12 and the main att row
  idx = order(r[fin]', 1)
  for (j=1; j<=min((12, cols(fin))); j++) {
     kk = fin[idx[j]]
     printf("WORST %-28s navail=%4.0f  printed_se=%10.6f  exported_se=%10.6f  ratio=%6.4f\n",
            cs[kk,2], cnt[kk], se_pc[kk], d_o[kk], r[kk])
  }
  for (k=1;k<=K;k++) {
     if (cs[k,2] == "main_att" | cs[k,2] == "main_att_it") {
        printf("KEY %-28s navail=%4.0f  printed_se=%10.6f  exported_se=%10.6f  ratio=%6.4f\n",
               cs[k,2], cnt[k], se_pc[k], d_o[k], r[k])
     }
  }
end

dis "=== cautions? grep the log for 'Caution' ==="
dis "=== the footer line, verbatim, is above ==="
cap noisily lincom main_att
dis "P_LINCOM_se=" r(se)
dis "P_LINCOM_RC=" _rc
cap noisily mat list e(wdidmainresults), format(%12.6f)
dis "R2_PILOT_DONE"
