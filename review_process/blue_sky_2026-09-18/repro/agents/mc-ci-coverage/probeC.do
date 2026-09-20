clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
set seed 20260916
set sortseed 20260916
include "$WD/shapelib.doi"

* independent re-implementation of the three interval families from the raw draws
mata:
real scalar qs(real colvector v, real scalar p)
{
  real scalar n, x
  n = rows(v)
  x = n*p
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1,n))])/2)
  return(v[max((min((ceil(x),n)),1))])
}
void recheck(string scalar tag, real scalar lvl, real scalar dfden)
{
  real matrix A, cin, cip, cib
  real rowvector obs
  real colvector v, use
  real scalar K, k, al, zc, nn, se, z0, mx1, mx2, mx3
  string matrix cs
  A   = st_matrix("e(bootdist)")
  obs = st_matrix("e(b)")
  cin = st_matrix("e(boot_ci_normal)")
  cip = st_matrix("e(boot_ci_percentile)")
  cib = st_matrix("e(boot_ci_bc)")
  cs  = st_matrixcolstripe("e(bootdist)")
  K = cols(A)
  al = (100 - lvl)/200
  zc = (dfden < . ? invt(dfden, 1-al) : invnormal(1-al))
  mx1 = 0; mx2 = 0; mx3 = 0
  for (k = 1; k <= K; k++) {
    use = selectindex(A[.,k] :< .)
    nn = length(use)
    if (nn < 2) continue
    v = sort(A[use,k],1)
    se = sqrt(variance(v))
    if (se <= 0 | se >= .) continue
    mx1 = max((mx1, abs(cin[1,k] - (obs[k]-zc*se)), abs(cin[2,k] - (obs[k]+zc*se))))
    mx2 = max((mx2, abs(cip[1,k] - qs(v,al)), abs(cip[2,k] - qs(v,1-al))))
    z0 = invnormal(sum(v :<= obs[k])/nn)
    if (z0 < .) {
      mx3 = max((mx3, abs(cib[1,k] - qs(v, normal(2*z0+invnormal(al)))),
                      abs(cib[2,k] - qs(v, normal(2*z0+invnormal(1-al))))))
    }
  }
  printf("RECHK|%s|lvl=%g|K=%g|dnormal=%12.4e|dpctl=%12.4e|dbc=%12.4e\n", tag, lvl, K, mx1, mx2, mx3)
}
end

* ------------- 40-unit fixture (cheap jackknife) -------------
local FN 40
local FT 8
local FSEED 11
include "$WD/mkfix.doi"
save "$WD/fixC.dta", replace

foreach L in 95 90 99 87.5 {
  local Lt = subinstr("`L'",".","p",.)
  cap noi wooldid y cg t ttre, att cluster(unit) reps(399) rseed(90210) ///
       bootstrap(standard) bootdisplay(bca) level(`L') esprelength(2) espostlength(2)
  dis "RC|C1_`Lt'|" _rc
  if _rc == 0 {
    dis "ELEVEL|C1_`Lt'|" e(level) "|completed=" e(boot_completed) "|failed=" e(boot_failed) "|jk=" e(boot_jk_completed)
    wdck C1_`Lt' e(wdidmainresults)
    mata: recheck("C1_`Lt'", `L', .)
    mat list e(boot_ci_normal), format(%14.9f)
    mat list e(boot_ci_percentile), format(%14.9f)
    mat list e(boot_ci_bc), format(%14.9f)
    cap confirm matrix e(boot_ci_bca)
    if _rc==0 {
      mat list e(boot_ci_bca), format(%14.9f)
      mat list e(boot_accel), format(%14.9f)
    }
    mat list e(boot_z0), format(%14.9f)
    * half-interval / lb>=ub scan over every family
    foreach fam in normal percentile bc bca {
      cap confirm matrix e(boot_ci_`fam')
      if _rc==0 {
        mat CC = e(boot_ci_`fam')
        local cn : colnames CC
        local nbad = 0
        local nhalf = 0
        local nout = 0
        mat BO = e(b)
        forvalues k = 1/`=colsof(CC)' {
          local lo = CC[1,`k']
          local hi = CC[2,`k']
          local bb = BO[1,`k']
          local nm : word `k' of `cn'
          if `lo'<. & `hi'<. {
            if !(`lo'<`hi') {
              local ++nbad
              dis "VIOL_CILBUB|C1_`Lt'|`fam'|`nm'|" `lo' "|" `hi'
            }
            if `bb'<. {
              if (`bb'<`lo') | (`bb'>`hi') {
                local ++nout
                dis "VIOL_BOUT|C1_`Lt'|`fam'|`nm'|b=" `bb' "|[" `lo' "," `hi' "]"
              }
            }
          }
          if (`lo'<. & `hi'>=.) | (`lo'>=. & `hi'<.) {
            local ++nhalf
            dis "VIOL_CIHALF|C1_`Lt'|`fam'|`nm'|lb=" `lo' "|ub=" `hi'
          }
        }
        dis "CIFAM|C1_`Lt'|`fam'|K=" colsof(CC) "|lbub=`nbad'|half=`nhalf'|bout=`nout'"
      }
    }
  }
}
dis "PROBEC_DONE"
