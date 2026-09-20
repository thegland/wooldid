* p04b_brute.do -- part B of the p04 oracle: enumerate all 56 assignments independently
* and refit on the engine's OWN replication option path.  The engine's repopts
* (ado:125-163) carry the estimand flags, subgroups(), esprelength()/espostlength(),
* cluster() and `fastnose' -- but NOT jointtests (fastnose computes no covariance).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid

local REPOPTS att att_it att_i att_itime subgroups(sg) esprelength(2) espostlength(3) fastnose

tempname fh
file open `fh' using "p04_B_draws.csv", write replace
file write `fh' "rep,a,b,name,value" _n
local rep 0
forvalues a = 1/8 {
  forvalues b = 1/8 {
    if `a' != `b' {
      local ++rep
      use p04_data, clear
      gen double ons_r = .
      replace ons_r = 5 if id == `a'
      replace ons_r = 7 if id == `b'
      cap qui wooldid y id t ons_r, cluster(id) `REPOPTS'
      local rc = _rc
      if `rc' == 0 {
        tempname bb
        matrix `bb' = e(b)
        local nn : colfullnames `bb'
        local kk = colsof(`bb')
        forvalues j = 1/`kk' {
          local nmj : word `j' of `nn'
          mata: st_local("vj", strofreal(st_matrix("`bb'")[1, `j'], "%26.17e"))
          file write `fh' "`rep',`a',`b',`nmj',`vj'" _n
        }
      }
      dis "BDRAW rep=`rep' a=`a' b=`b' rc=`rc' k=`kk'"
    }
  }
}
file close `fh'
dis "ALLDONE"
