* p03_oracle.do -- INDEPENDENT BRUTE-FORCE ENUMERATION.
* The design has M = 6!/(1!1!4!) = 30 distinct onset assignments: cluster a takes onset 4,
* cluster b != a takes onset 7, the rest are never-treated.  Enumerate all 30 here, refit
* wooldid once per assignment on the SAME rep-option path the engine uses (`fastnose' for
* the coefficient statistic, no `fastnose' for the studentized one) and dump the resulting
* reference distributions.  p03_oracle.py then forms every documented p-value from them and
* compares with what `permexact' reported in p01_conv.do.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"

which wooldid

*---------------------------------------------------------------- the same data as p01
clear
set obs 6
gen long id = _n
gen double ons = .
replace ons = 4 in 1
replace ons = 7 in 2
expand 10
bysort id: gen int t = _n
expand 4
gen double y = 0.6*id + 0.15*t + 1.2*(ons < . & t >= ons) + 0.5*rnormal()
compress
save p03_data, replace

*---------------------------------------------------------------- enumerate
tempname pf
postfile `pf' int rep int ca int cb double bfast double bfull double sefull ///
    using p03_draws, replace
local rep 0
forvalues a = 1/6 {
  forvalues b = 1/6 {
    if `a' != `b' {
      local ++rep
      use p03_data, clear
      gen double ons_r = .
      replace ons_r = 4 if id == `a'
      replace ons_r = 7 if id == `b'
      * scalars, never locals: a local's decimal round-trip re-parses one ulp off
      scalar sbf = .
      scalar sbl = .
      scalar ssl = .
      cap qui wooldid y id t ons_r, cluster(id) fastnose
      if _rc == 0 scalar sbf = _b[main_att]
      cap qui wooldid y id t ons_r, cluster(id)
      if _rc == 0 {
        scalar sbl = _b[main_att]
        scalar ssl = _se[main_att]
      }
      post `pf' (`rep') (`a') (`b') (sbf) (sbl) (ssl)
      dis "DRAW rep=`rep' a=`a' b=`b' bfast=" sbf " bfull=" sbl " sefull=" ssl
    }
  }
}
postclose `pf'

use p03_draws, clear
* the observed assignment is a=1 (onset 4), b=2 (onset 7)
list, noobs sepby(ca)
outsheet using "p03_draws.csv", comma replace

mata:
X = st_data(., ("rep", "ca", "cb", "bfast", "bfull", "sefull"))
fh = fopen("p03_draws_hi.csv", "w")
fput(fh, "rep,ca,cb,bfast,bfull,sefull")
for (r = 1; r <= rows(X); r++) fput(fh, invtokens(strofreal(X[r, .], "%26.17e"), ","))
fclose(fh)
end
dis "ALLDONE"
