*! n03 (NEW PROBE, not seeded) — which sample counts the bootstrap's CLUSTERS?
*! wooldid_bootrun builds its cluster map on `smp' = (user if == 1) & nonmissing cluster
*! (ado:8910-8925) and hands wdb_C back as wdb_r_nclus (ado:9537), which becomes the
*! default bootdf(clusters) => e(df_r) = G - 1 (ado:423) and drives the "only N clusters"
*! and "only N treated clusters" cautions (ado:6502, 6520).  None of those counts is the
*! ESTIMATION sample.  A county whose outcome is never recorded therefore still buys a
*! degree of freedom and still suppresses the few-cluster caution.
*! Realistic shape: a 14-county panel where the outcome is collected in 7 counties only.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mk
prog def mk
  syntax , nblank(int)
  clear
  set seed 4242
  set sortseed 4242
  set obs 14
  gen long cty = _n
  gen onset = cond(cty <= 4, 4, cond(cty <= 8, 6, .))
  expand 3
  sort cty
  gen long unit = _n
  gen double uf = rnormal(0,1)
  expand 7
  bysort unit (cty): gen int t = _n
  gen byte D = (onset < . & t >= onset)
  gen double y = 2 + uf + .1*t + 1*D + rnormal(0,1)
  * the outcome is not collected in the last `nblank' counties (still onset-labelled,
  * still in the frame, still carrying a cluster code)
  * blank the outcome in `nblank' counties, keeping treated AND never-treated counties measured
  if `nblank' == 7  replace y = . if inlist(cty, 2, 3, 6, 7, 10, 11, 12)
  if `nblank' == 12 replace y = . if !inlist(cty, 1, 9)
end

dis "===== C0: all 14 counties measured (reference)"
mk , nblank(0)
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(99) rseed(13) dots(0)
dis "C0_RC=" _rc
dis "C0 N=" e(N) "  df_r=" e(df_r) "  boot_N_V=" e(boot_N_V)
qui levelsof cty if e(sample), local(L)
dis "C0 clusters in e(sample) = " `:word count `L''

dis "===== C7: outcome collected in 7 of 14 counties"
mk , nblank(7)
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(99) rseed(13) dots(0)
dis "C7_RC=" _rc
dis "C7 N=" e(N) "  df_r=" e(df_r) "  boot_N_V=" e(boot_N_V)
qui levelsof cty if e(sample), local(L)
dis "C7 clusters in e(sample) = " `:word count `L''

dis "===== C7A: the ANALYTIC run on the same data, for the df it uses"
mk , nblank(7)
cap noi wooldid y unit t onset, att att_i cluster(cty)
dis "C7A_RC=" _rc
dis "C7A N=" e(N) "  df_r=" e(df_r)

dis "===== C12: outcome collected in 2 of 14 counties (edge)"
mk , nblank(12)
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(99) rseed(13) dots(0)
dis "C12_RC=" _rc
dis "C12 N=" e(N) "  df_r=" e(df_r)
dis "===== DONE"
