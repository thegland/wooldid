*==========================================================================
* REPRO (uncertain finding): with cluster() at the level of `i', every
*        subgroups() row's analytic standard error equals the pooled one to
*        15 significant digits -- on every event-study row, on both subgroups.
*
*  Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_sgse_clusteri.do
*
*  Three comparisons on the SAME data:
*    V1  cluster(state), state == i, subgroups(reg2)
*    V2  robust        , same subgroups
*    V3  two separate runs, one per subgroup, cluster(state)
*  V1's three standard errors are the same number.  V2's are not.  V3 shows
*  what the subgroup-specific standard errors actually look like.
*
*  Why this sits on the boottstat surface: e(boot_ase_obs) -- the percentile-t
*  studentization denominator, and the scale of the printed interval
*  b - q*se_obs -- is exactly this quantity, taken from column 2 of
*  e(wdidsgresults_*) before the rewiring (wooldid_bootwire, ado:9606-9625).
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

clear
set obs 24
gen int state = _n
gen int onset = .
forval s = 1/18 {
  local w = mod(`s'-1, 3) + 1
  local o : word `w' of 2002 2003 2004
  qui replace onset = `o' in `s'
}
gen double a_i = rnormal()*1.0
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand 7
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t - onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
gen byte reg2 = mod(state, 2)

dis _n(2) "########## V1: cluster(state), state == i ##########"
qui wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(reg2)
mat V1M = e(wdidmainresults)
mat V1A = e(wdidsgresults_0)
mat V1B = e(wdidsgresults_1)
dis "  main att  b=" %12.9f V1M[1,1] "  se=" %12.9f V1M[1,2] "  (" %21x V1M[1,2] ")"
dis "  sg0  att  b=" %12.9f V1A[1,1] "  se=" %12.9f V1A[1,2] "  (" %21x V1A[1,2] ")"
dis "  sg1  att  b=" %12.9f V1B[1,1] "  se=" %12.9f V1B[1,2] "  (" %21x V1B[1,2] ")"
dis "  se ratios to pooled: " %14.11f V1A[1,2]/V1M[1,2] "   " %14.11f V1B[1,2]/V1M[1,2]
mata:
  A = st_matrix("V1M"); B = st_matrix("V1A"); C = st_matrix("V1B")
  rs = st_matrixrowstripe("V1M")
  n = 0
  for (r = 1; r <= rows(A); r++) {
    if (A[r,2] >= .) continue
    same = (abs(B[r,2]/A[r,2] - 1) < 1e-12) & (abs(C[r,2]/A[r,2] - 1) < 1e-12)
    n = n + same
    printf("  row %-14s|%-10s  pooled=%12.9f  sg0=%12.9f  sg1=%12.9f  %s\n",
       rs[r,1], rs[r,2], A[r,2], B[r,2], C[r,2], (same ? "IDENTICAL" : ""))
  }
  printf("\nV1_ROWS_WITH_IDENTICAL_SE = %f of %f\n", n, rows(A))
end

dis _n(2) "########## V2: robust, same data, same subgroups ##########"
qui wooldid y state t ttre , robust att esprelength(2) espostlength(2) subgroups(reg2)
mat V2M = e(wdidmainresults)
mat V2A = e(wdidsgresults_0)
mat V2B = e(wdidsgresults_1)
dis "  main att se=" %12.9f V2M[1,2] "   sg0 se=" %12.9f V2A[1,2] "   sg1 se=" %12.9f V2B[1,2]
dis "  se ratios to pooled: " %14.11f V2A[1,2]/V2M[1,2] "   " %14.11f V2B[1,2]/V2M[1,2]

dis _n(2) "########## V3: separate runs, one per subgroup, cluster(state) ##########"
qui wooldid y state t ttre if reg2==0 , cluster(state) att esprelength(2) espostlength(2)
dis "  reg2==0 alone: b=" %12.9f e(wdidmainresults)[1,1] "  se=" %12.9f e(wdidmainresults)[1,2]
qui wooldid y state t ttre if reg2==1 , cluster(state) att esprelength(2) espostlength(2)
dis "  reg2==1 alone: b=" %12.9f e(wdidmainresults)[1,1] "  se=" %12.9f e(wdidmainresults)[1,2]
