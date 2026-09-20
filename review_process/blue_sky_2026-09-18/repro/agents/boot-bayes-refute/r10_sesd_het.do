* r10 (own hunt): QUANTIFY the analytic-SE blindness.
* review_process/se_coverage_2026-09-14 measured SE/SD = 0.66-0.75 for i = unit id with
* cluster(i); its two DGPs carry NO unit-level treatment-effect heterogeneity (tau is
* constant within cohort).  Re-run the same shape with a UNIT-level effect shock of
* increasing size and watch the reported analytic SE stand still while the sampling sd
* grows.  Analytic runs only - fast.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid
tempname fh
file open `fh' using "r10_sesd.csv", write replace
file write `fh' "arm,draw,b,se,cellest,cellse" _n

foreach hsd in 0 0.5 1.5 {
  forvalues d = 1/80 {
    set seed `=700000 + `d'*29'
    set sortseed `=700000 + `d'*29'
    clear
    set obs 30
    gen int st = _n
    gen double a_i = rnormal()
    gen double shk = rnormal()
    gen int onset = .
    replace onset = 2004 if st <=  8
    replace onset = 2006 if st >  8 & st <= 15
    replace onset = 2008 if st > 15 & st <= 21
    expand 10
    bysort st: gen int year = 2000 + _n
    gen byte post = (year >= onset) & !missing(onset)
    gen double y = a_i + 0.03*(year-2000) + (1.0 + `hsd'*shk)*post + rnormal()*0.7
    cap noi qui wooldid y st year onset, att cluster(st) gteffects
    if _rc == 0 {
      matrix A = e(wdidmainresults)
      matrix G = e(gteffects)
      * e(gteffects) cols: 1 i, 2 onset, 3 caltime, 4 relyr, 5 estimate, 6 se
      * pick the cell i == 1, relyr == 0 (state 1, its first treated year)
      local ce = .
      local cs = .
      forvalues r = 1/`=rowsof(G)' {
        if G[`r',1] == 1 & G[`r',4] == 0 {
          local ce = G[`r',5]
          local cs = G[`r',6]
        }
      }
      file write `fh' "`hsd',`d'," (A[1,1]) "," (A[1,2]) "," (`ce') "," (`cs') _n
    }
    else dis "FAIL arm=`hsd' draw=`d' rc=" _rc
  }
  dis "ARM `hsd' DONE"
}
file close `fh'
dis "R10_DONE"
