* r10 (UNSEEDED): which replication subset does the REPORTED bootstrap standard error
* use?  The help says e(V) is the "complete-replication listwise" bootstrap covariance.
* Recompute three candidates from e(bootdist) and see which one the printed column 2
* of e(wdidmainresults)/e(wdidsgresults_*) actually equals.
*   SDP = per-column SD (each column's own non-missing rows)
*   SDL = listwise SD over ALL columns of e(bootdist)
*   SDE = listwise SD over the columns whose OBSERVED estimate is non-missing
* Percent differences are reported as true relative differences, not Stata's reldif().
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid
set seed 242424
set sortseed 242424
clear
set obs 24
gen int state = _n
gen int onset = .
qui replace onset = 2002 if state <= 6
qui replace onset = 2003 if state >= 7  & state <= 12
qui replace onset = 2004 if state >= 13 & state <= 17
gen double a_i = rnormal()
gen int nun = 2 + mod(state,3)
expand nun
bys state: gen int unit = _n + 100*state
gen double u_i = rnormal()*0.6
expand 9
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t-onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
gen byte sg3 = mod(state,3)

cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   subgroups(sg3) bootstrap(standard) reps(399) dots(0)
dis "RC=" _rc
mat BD = e(bootdist)
mat BB = e(b)
local cfn : colfullnames BB
dis "bootdist rows=" rowsof(BD) " cols=" colsof(BD)
mat REP = J(1, colsof(BB), .)
mat WM = e(wdidmainresults)
local kk 0
forval r = 1/`=rowsof(WM)' {
  local ++kk
  mat REP[1,`kk'] = WM[`r',2]
}
foreach sg of local e_sglist {
}
local sgl `e(wdid_sglist)'
foreach sg of local sgl {
  cap mat WS = e(wdidsgresults_`sg')
  if _rc == 0 {
    forval r = 1/`=rowsof(WS)' {
      local ++kk
      mat REP[1,`kk'] = WS[`r',2]
    }
  }
}
mata:
  A  = st_matrix("BD")
  ob = st_matrix("BB")
  R = rows(A); K = cols(A)
  // per column
  sdp = J(1,K,.); np = J(1,K,.)
  for (k=1;k<=K;k++) {
    v = select(A[.,k], A[.,k] :< .)
    np[k] = rows(v)
    if (rows(v) > 1) sdp[k] = sqrt(variance(v))
  }
  // listwise over all columns
  ok = J(R,1,1)
  for (r=1;r<=R;r++) for (k=1;k<=K;k++) if (A[r,k] >= .) ok[r] = 0
  B = select(A, ok)
  sdl = J(1,K,.)
  for (k=1;k<=K;k++) sdl[k] = sqrt(variance(B[.,k]))
  st_matrix("SDP", sdp)
  st_matrix("SDL", sdl)
  st_matrix("NP", np)
  st_local("nl", strofreal(rows(B)))
end
dis "N_LISTWISE_ALLCOLS=" `nl'
local wp 0
local wl 0
forval k = 1/`=colsof(BB)' {
  local nm : word `k' of `cfn'
  local dp = abs(REP[1,`k']-SDP[1,`k'])/SDP[1,`k']
  local dl = abs(REP[1,`k']-SDL[1,`k'])/SDL[1,`k']
  if `dp' > `wp' local wp = `dp'
  if `dl' > `wl' local wl = `dl'
  dis "SD `nm'  reported=" %12.9f REP[1,`k'] "  percol=" %12.9f SDP[1,`k'] " (n=" NP[1,`k'] ")" ///
      "  listwise=" %12.9f SDL[1,`k'] "  pctdiff_percol=" %9.5f 100*`dp' "  pctdiff_listwise=" %9.5f 100*`dl'
}
dis "WORST_PCTDIFF_PERCOL=" %10.5f 100*`wp'
dis "WORST_PCTDIFF_LISTWISE=" %10.5f 100*`wl'
dis "R10_DONE"
