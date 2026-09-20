clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

dis as result "==== LABEL ROUND-TRIP SWEEP ===="
mata:
void wdlabsweep() {
  real scalar L, k, gtal, nbad
  string scalar wk, gt
  nbad = 0
  for (k = 1000; k <= 9999; k++) {
    L = k/100
    wk = (mod(L,1)==0 ? strofreal(L) : subinstr(strofreal(L,"%9.0g"), ".", "p", .))
    gtal = (100 - L)/200
    gt = subinstr(strofreal(100 - 200*gtal, "%9.0g"), ".", "p", .)
    if (wk != gt) {
      nbad++
      if (nbad <= 20) printf("LABELMISMATCH L=%9.4f worker=%s gteffects=%s\n", L, wk, gt)
    }
  }
  printf("LABELSWEEP_DONE nbad=%f\n", nbad)
}
wdlabsweep()
end

local dN 60
local dT 10
local donsets "2004 2006"
include "`W'/mkdata.do"

foreach L in 80 92.5 99.99 10 33.33 66.67 {
  cap noi wooldid y i t ttre, att robust espre(2) espost(2) gteffects level(`L')
  dis "GT_RC_`L'=" _rc
  if _rc == 0 {
    tempname G M
    mat `G' = e(gteffects)
    mat `M' = e(wdidmainresults)
    local gc : colnames `G'
    local mc : colnames `M'
    dis "L=`L'  main cols: `mc'"
    dis "L=`L'  gteffects cols: `gc'"
    local crit = invttail(e(df_r), (1-`L'/100)/2)
    local worst = 0
    forvalues r = 1/`=rowsof(`G')' {
      local b = `G'[`r',5]
      local s = `G'[`r',6]
      local lo = `G'[`r',8]
      if `b' < . & `s' < . & `lo' < . {
        local e1 = `b' - `crit'*`s'
        local d = abs(`lo'-`e1')/max(abs(`e1'),1e-12)
        if `d' > `worst' local worst = `d'
      }
    }
    dis "L=`L'  gteffects_worst_rel_dev=" %12.4e `worst' "   crit=" %14.10f `crit'
  }
}
dis "P8_DONE"
