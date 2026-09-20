* r09 (UNSEEDED): two independent checks on the bootstrap rewiring surface.
*  A) wooldid_bootwire's own comment says only the FIRST results-matrix walk is
*     guarded, and that a later walk that mis-lands "shifts every subgroup standard
*     error one column with no error and no caution".  Test it by recomputing every
*     reported bootstrap SE -- main table AND subgroup tables -- from e(bootdist)
*     directly (listwise-complete SD, which is what e(V) is documented to hold) and
*     comparing name-for-name.
*  B) bootsave() under bootalpha(size): the help says the _rep = 0 row of the w_
*     columns holds 1/N, N the number of observations in the estimation sample.
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
tempfile base
save `base', replace

dis "### A: bootstrap + subgroups(sg3), recompute every reported SE from e(bootdist)"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   subgroups(sg3) bootstrap(standard) reps(399) dots(0)
dis "A_RC=" _rc
mat BD = e(bootdist)
mat BB = e(b)
local cfn : colfullnames BB
dis "A_bootdist rows=" rowsof(BD) " cols=" colsof(BD) " ; e(b) cols=" colsof(BB)
dis "A_boot_N_complete=" e(boot_N_complete)
* stack the reported SEs in e(b) column order: main table then sg blocks in sglist order
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
dis "A_stacked=" `kk' " of " colsof(BB)
* listwise-complete SD of each bootdist column
mata:
  A = st_matrix("BD")
  R = rows(A); K = cols(A)
  ok = J(R,1,1)
  for (r=1;r<=R;r++) for (k=1;k<=K;k++) if (A[r,k] >= .) ok[r] = 0
  B = select(A, ok)
  n = rows(B)
  sd = J(1,K,.)
  for (k=1;k<=K;k++) sd[k] = sqrt(variance(B[.,k]))
  st_matrix("MYSD", sd)
  st_local("ncomp", strofreal(n))
end
dis "A_my_complete=" `ncomp'
local worst 0
local worstnm ""
local shifted 0
forval k = 1/`=colsof(BB)' {
  local nm : word `k' of `cfn'
  local rd = reldif(REP[1,`k'], MYSD[1,`k'])
  if `rd' > `worst' {
    local worst = `rd'
    local worstnm "`nm'"
  }
  local tag ""
  if `rd' > 1e-8 {
    local tag " <-- MISMATCH"
    * is it the NEIGHBOUR's sd?
    if `k' > 1 {
      if reldif(REP[1,`k'], MYSD[1,`k'-1]) < 1e-10 local tag " <-- EQUALS PREVIOUS COLUMN SD"
    }
    if `k' < colsof(BB) {
      if reldif(REP[1,`k'], MYSD[1,`k'+1]) < 1e-10 local tag " <-- EQUALS NEXT COLUMN SD"
    }
    local ++shifted
  }
  dis "A_SE `nm'  reported=" %14.10f REP[1,`k'] "  recomputed=" %14.10f MYSD[1,`k'] "  reldif=" %10.3e `rd' "`tag'"
}
dis "A_WORST_RELDIF=" %12.4e `worst' "  on `worstnm'   A_NMISMATCH=" `shifted'

dis "### B: bootsave under bootalpha(size) -- the _rep=0 w_ row should be 1/N"
use `base', clear
cap noi wooldid y state t ttre , cluster(state) att esprelength(1) espostlength(1) ///
   bootstrap(bayes) bootalpha(size) reps(99) bootsave(bsz) replace dots(0)
dis "B_RC=" _rc
local NN = e(N)
dis "B_N=" `NN' "  1/N=" %18.14f 1/`NN'
preserve
use bsz, clear
qui ds w_*
local wl `r(varlist)'
local w1 : word 1 of `wl'
qui sum `w1' if _rep==0
dis "B_w_rep0=" %18.14f r(mean) "  nwcols=" wordcount("`wl'")
qui egen double _rs = rowtotal(`wl') if _rep>0
qui sum _rs
dis "B_rowsum_reps mean=" %18.14f r(mean) " min=" %18.14f r(min) " max=" %18.14f r(max)
restore
dis "### C: same but bootalpha(equal) -- _rep=0 w_ row should be 1/G"
use `base', clear
cap noi wooldid y state t ttre , cluster(state) att esprelength(1) espostlength(1) ///
   bootstrap(bayes) bootalpha(equal) reps(99) bootsave(bse) replace dots(0)
dis "C_RC=" _rc
dis "C_nclust=" e(N_clust) "  1/G=" %18.14f 1/e(N_clust)
preserve
use bse, clear
qui ds w_*
local wl `r(varlist)'
local w1 : word 1 of `wl'
qui sum `w1' if _rep==0
dis "C_w_rep0=" %18.14f r(mean) "  nwcols=" wordcount("`wl'")
restore
dis "R09_DONE"
