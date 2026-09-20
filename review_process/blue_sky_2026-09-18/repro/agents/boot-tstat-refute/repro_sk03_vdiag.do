* r11 (UNSEEDED follow-up to r10): the printed bootstrap SE is the PER-COLUMN SD of the
* draws (r10: exact to 0.00000 percent).  The run's footer says e(V) is the
* "complete-replication listwise" bootstrap covariance, and the bootdf() help says a
* test/lincom on e(V) reproduces the Wald family.  Does sqrt(diag(e(V))) equal the
* printed SE, or the listwise SD?  If it equals the listwise one, the printed row and
* e(V) disagree, and lincom after the command does not reproduce the printed interval.
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
   subgroups(sg3) bootstrap(standard) reps(399) bootdisplay(wald) dots(0)
dis "RC=" _rc
mat V  = e(V)
mat BD = e(bootdist)
mat BB = e(b)
local cfn : colfullnames BB
mat REP = J(1, colsof(BB), .)
mat WM = e(wdidmainresults)
local kk 0
forval r = 1/`=rowsof(WM)' {
  local ++kk
  mat REP[1,`kk'] = WM[`r',2]
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
  R = rows(A); K = cols(A)
  sdp = J(1,K,.)
  for (k=1;k<=K;k++) {
    v = select(A[.,k], A[.,k] :< .)
    if (rows(v) > 1) sdp[k] = sqrt(variance(v))
  }
  ok = J(R,1,1)
  for (r=1;r<=R;r++) for (k=1;k<=K;k++) if (A[r,k] >= .) ok[r] = 0
  B = select(A, ok)
  sdl = J(1,K,.)
  for (k=1;k<=K;k++) sdl[k] = sqrt(variance(B[.,k]))
  st_matrix("SDP", sdp)
  st_matrix("SDL", sdl)
end
local wV 0
forval k = 1/`=colsof(BB)' {
  local nm : word `k' of `cfn'
  local sv = sqrt(V[`k',`k'])
  local d1 = abs(`sv'-REP[1,`k'])/REP[1,`k']
  local d2 = abs(`sv'-SDP[1,`k'])/SDP[1,`k']
  local d3 = abs(`sv'-SDL[1,`k'])/SDL[1,`k']
  if `d1' > `wV' local wV = `d1'
  dis "V `nm'  printed_se=" %12.9f REP[1,`k'] "  sqrt(V_kk)=" %12.9f `sv' ///
      "  percolSD=" %12.9f SDP[1,`k'] "  listwiseSD=" %12.9f SDL[1,`k'] ///
      "  pct(V vs printed)=" %9.5f 100*`d1' "  pct(V vs percol)=" %9.5f 100*`d2' "  pct(V vs listwise)=" %9.5f 100*`d3'
}
dis "WORST_PCT_V_VS_PRINTED=" %10.5f 100*`wV'
* does lincom on e(V) reproduce the printed row?
dis "--- lincom check on the worst coefficient ---"
cap noi lincom _b[pretreatment_SG0_att]
dis "LINCOM_RC=" _rc
dis "R11_DONE"
