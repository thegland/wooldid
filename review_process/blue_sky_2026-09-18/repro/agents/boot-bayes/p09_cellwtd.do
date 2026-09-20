* p09: bayes + a cell-weighted estimand under safety(off) - what actually comes out?
* The caution says the SE "can be too small". Quantify against bootstrap(standard),
* which resamples cells, and against the analytic clustered SE.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

* effect heterogeneous ACROSS clusters - the regime the caution names
use p01_panel, clear
set seed 4711
qui bysort st (year): gen double shk = rnormal() if _n==1
qui bysort st (year): replace shk = shk[1]
qui replace y = y + 1.2*shk*post
save p09_panel, replace

dis "=== analytic ==="
cap noi wooldid y st year onset, att att_i att_it att_itime cluster(st)
dis "AN_RC=" _rc
mat li e(wdidmainresults)
dis "=== bootstrap(standard) ==="
use p09_panel, clear
cap noi wooldid y st year onset, att att_i att_it att_itime cluster(st) bootstrap(standard) reps(300) rseed(99) dots(0)
dis "ST_RC=" _rc
mat li e(wdidmainresults)
matrix BS = e(b)
matrix VS = e(V)
dis "=== bootstrap(bayes) safety(off) ==="
use p09_panel, clear
cap noi wooldid y st year onset, att att_i att_it att_itime cluster(st) bootstrap(bayes) reps(300) rseed(99) dots(0) safety(off)
dis "BY_RC=" _rc
mat li e(wdidmainresults)
matrix BB = e(b)
matrix VB = e(V)
mata {
  vs = diagonal(st_matrix("VS")); vb = diagonal(st_matrix("VB"))
  nm = st_matrixcolstripe("VS")
  for (k=1; k<=rows(vs); k++) printf("%-14s  se_standard=%10.6f  se_bayes=%10.6f  ratio=%7.4f\n", nm[k,2], sqrt(vs[k]), sqrt(vb[k]), sqrt(vb[k]/vs[k]))
}
dis "PROBE_DONE"
