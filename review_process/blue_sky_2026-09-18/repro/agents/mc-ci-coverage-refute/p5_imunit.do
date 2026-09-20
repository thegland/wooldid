*=========================================================================
* p5_imunit.do -- does the impvalues under-coverage extend to the ORDINARY
* geometry  i = unit, cluster(unit), where every cluster owns its own
* coefficients and Ibragimov-Muller is textbook-legitimate?
* Arms on the SAME draws: impvalues cluster(unit) vs analytic cluster(unit).
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p5_imunit.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local R 100
local OUT "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage-refute/p5_imunit.txt"

cap prog drop mkdat
prog def mkdat
  args nunit sd
  clear
  set seed `sd'
  set sortseed `sd'
  set obs `nunit'
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ttre = .
  replace ttre = 2004 if cg==1
  replace ttre = 2006 if cg==2
  replace ttre = 2008 if cg==3
  gen double ui   = rnormal(0,0.70)
  gen double hmul = 1 + rnormal(0,0.25)
  expand 12
  bys unit: gen int t = 2000 + _n - 1
  gen double e = .
  bys unit (t): replace e = rnormal(0,0.55) if _n==1
  bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
  gen double truetau = 0
  replace truetau = 1.00*hmul if cg==1 & t>=ttre
  replace truetau = 1.15*hmul if cg==2 & t>=ttre
  replace truetau = 0.85*hmul if cg==3 & t>=ttre
  gen double y = ui + 0.25*(t-2000) + truetau + e
end

file open fh using "`OUT'", write replace text
file write fh "rep,arm,b,se,lb,ub,dfr,truth_att,truth_ceq" _n
file close fh

forvalues r = 1/`R' {
  mkdat 60 `=70000+`r''
  qui su truetau if !missing(ttre) & t>=ttre
  local truth_att = r(mean)
  tempvar um
  qui egen double `um' = mean(cond(!missing(ttre) & t>=ttre, truetau, .)), by(unit)
  qui bys unit: gen byte __f = (_n==1)
  qui su `um' if !missing(ttre) & __f==1, meanonly
  local truth_ceq = r(mean)
  drop __f
  tempfile dd
  qui save "`dd'"
  foreach arm in I C {
    use "`dd'", clear
    if "`arm'"=="I" cap qui wooldid y unit t ttre, att cluster(unit) impvalues level(95)
    if "`arm'"=="C" cap qui wooldid y unit t ttre, att cluster(unit) level(95)
    if _rc == 0 {
      mat T = e(wdidmainresults)
      local dfr = e(df_r)
      file open fh using "`OUT'", write append text
      file write fh "`r',`arm'," (T[1,1]) "," (T[1,2]) "," (T[1,5]) "," (T[1,6]) "," (`dfr') "," (`truth_att') "," (`truth_ceq') _n
      file close fh
    }
    else {
      file open fh using "`OUT'", write append text
      file write fh "`r',`arm',.,.,.,.,.,`truth_att',`truth_ceq'" _n
      file close fh
    }
  }
  if mod(`r',20)==0 di "PROGRESS rep=`r'"
}
di "ALLDONE_P5"
