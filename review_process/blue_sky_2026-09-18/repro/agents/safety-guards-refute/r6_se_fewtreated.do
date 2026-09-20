* R6 -- side measurement provoked by the panel shape the si8 false positive fires on:
* i = unit, cluster(unit), MANY never-treated controls and FEW treated units.
* Is the analytic ATT SE inside the project's documented understatement range
* (SE/SD 0.66-0.75, coverage 0.79-0.88 at G=40) or far outside it?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

local R = 60
mat RES = J(`R', 3, .)
forvalues r = 1/`R' {
  qui {
    clear
    set obs 300
    gen long unit = _n
    gen byte tg = (unit <= 20)
    expand 8
    bysort unit: gen int t = 1999 + _n
    gen double ttre = .
    replace ttre = 2002 + mod(unit,3) if tg == 1
    gen double ui = rnormal()
    bysort unit (t): replace ui = ui[1]
    gen double y = 2 + 0.15*(t-1999) + ui + rnormal() + 1.5*(t >= ttre & ttre < .)
  }
  cap qui wooldid y unit t ttre, att cluster(unit)
  if _rc == 0 {
    mat b = e(b)
    mat V = e(V)
    mat RES[`r',1] = b[1,1]
    mat RES[`r',2] = sqrt(V[1,1])
    mat RES[`r',3] = abs(b[1,1] - 1.5) <= invttail(e(df_r), 0.025)*sqrt(V[1,1])
  }
}
clear
svmat RES
rename RES1 bhat
rename RES2 sehat
rename RES3 cov95
qui sum bhat
local sdb = r(sd)
local nrep = r(N)
qui sum sehat
local msе = r(mean)
qui sum sehat
local mse = r(mean)
qui sum cov95
local cov = r(mean)
dis "R6_NREP=" `nrep'
dis "R6_SD_OF_BHAT=" %9.6f `sdb'
dis "R6_MEAN_ANALYTIC_SE=" %9.6f `mse'
dis "R6_SE_OVER_SD=" %9.6f `mse'/`sdb'
dis "R6_COVERAGE95=" %9.4f `cov'
dis "R6_DONE"
