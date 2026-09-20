clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* usage: stata-mp -e do mcboot.do <chunk> <ndraws> <reps>
local CH  `1'
local R   `2'
local REPS `3'
if "`CH'"   == "" local CH 1
if "`R'"    == "" local R 5
if "`REPS'" == "" local REPS 199

mata:
real scalar qsx(real colvector v, real scalar p)
{
  real scalar n, x
  n = rows(v)
  x = n*p
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1,n))])/2)
  return(v[max((min((ceil(x),n)),1))])
}
// families from the raw draws; dfden for the normal leg
void bootfam(real scalar lvl, real scalar dfden, string scalar outname)
{
  real matrix A, out
  real rowvector obs, z0v, av
  real colvector v, use
  real scalar K, k, al, zc, nn, se, z0, a, w, pl, pu
  A   = st_matrix("e(bootdist)")
  obs = st_matrix("e(b)")
  z0v = st_matrix("e(boot_z0)")
  av  = J(1, cols(A), .)
  if (st_matrix("e(boot_accel)") != J(0,0,.)) av = st_matrix("e(boot_accel)")
  K = cols(A)
  al = (100 - lvl)/200
  zc = (dfden < . ? invt(dfden, 1-al) : invnormal(1-al))
  out = J(9, K, .)
  for (k = 1; k <= K; k++) {
    use = selectindex(A[.,k] :< .)
    nn = length(use)
    if (nn < 2) continue
    v = sort(A[use,k],1)
    se = sqrt(variance(v))
    if (se <= 0 | se >= .) continue
    out[9,k] = se
    out[1,k] = obs[k] - zc*se
    out[2,k] = obs[k] + zc*se
    out[3,k] = qsx(v, al)
    out[4,k] = qsx(v, 1-al)
    z0 = (cols(z0v) >= k ? z0v[k] : .)
    if (z0 >= .) continue
    out[5,k] = qsx(v, normal(2*z0+invnormal(al)))
    out[6,k] = qsx(v, normal(2*z0+invnormal(1-al)))
    a = (cols(av) >= k ? av[k] : .)
    if (a >= .) continue
    w = z0 + invnormal(al)
    if (1 - a*w <= 0) continue
    out[7,k] = qsx(v, normal(z0 + w/(1-a*w)))
    w = z0 + invnormal(1-al)
    if (1 - a*w <= 0) continue
    out[8,k] = qsx(v, normal(z0 + w/(1-a*w)))
  }
  st_matrix(outname, out)
}
end

local FN 40
local FT 8
local FSEED 11
include "$WD/mkfix.doi"
replace y = base
qui wooldid y cg t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) level(95)
mat TB = e(b)
local tn : colnames TB

tempname fh
file open `fh' using "$WD/mcboot_`CH'.txt", write replace text
file write `fh' "draw|param|lvl|fam|lb|ub|b|bsd" _n
forvalues c = 1/`=colsof(TB)' {
  local nm : word `c' of `tn'
  file write `fh' "0|`nm'|.|truth|.|.|" (TB[1,`c']) "|." _n
}

forvalues r = 1/`R' {
  local dr = (`CH'-1)*1000 + `r'
  local FSEED = 200000 + `dr'
  local FN 40
  local FT 8
  include "$WD/mkfix.doi"
  cap qui wooldid y cg t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) ///
        level(95) bootstrap(standard) reps(`REPS') rseed(`=500000+`dr'') bootdisplay(bca) dots(0)
  if _rc != 0 {
    file write `fh' "`dr'|RC|.|rc|.|.|" (_rc) "|." _n
    continue
  }
  mat BB = e(b)
  local bn : colnames BB
  * implied dfden from the native normal interval
  mat CN = e(boot_ci_normal)
  mat BD = e(bootdist)
  local zimp = .
  * native check at 95
  mata: bootfam(95, 39, "NAT")
  local dnat = 0
  foreach fam in normal percentile bc bca {
    cap confirm matrix e(boot_ci_`fam')
    if _rc == 0 {
      mat NC = e(boot_ci_`fam')
      local row1 = cond("`fam'"=="normal",1,cond("`fam'"=="percentile",3,cond("`fam'"=="bc",5,7)))
      forvalues k = 1/`=colsof(NC)' {
        if NC[1,`k'] < . & NAT[`row1',`k'] < . {
          local d = abs(NC[1,`k'] - NAT[`row1',`k']) + abs(NC[2,`k'] - NAT[`=`row1'+1',`k'])
          if `d' > `dnat' local dnat = `d'
        }
      }
    }
  }
  file write `fh' "`dr'|NATIVEDIFF|.|chk|.|.|" (`dnat') "|" (e(boot_completed)) _n
  foreach L in 80 90 95 99 {
    mata: bootfam(`L', 39, "FE")
    forvalues k = 1/`=colsof(BB)' {
      local nm : word `k' of `bn'
      local bb = BB[1,`k']
      local sd = FE[9,`k']
      local i = 0
      foreach fam in normal percentile bc bca {
        local r1 = 1 + 2*`i'
        local r2 = 2 + 2*`i'
        local ++i
        file write `fh' "`dr'|`nm'|`L'|`fam'|" (FE[`r1',`k']) "|" (FE[`r2',`k']) "|" (`bb') "|" (`sd') _n
      }
    }
  }
  dis "MCB|`dr'|done"
}
file close `fh'
dis "MCBOOT_DONE_`CH'"
