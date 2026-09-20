*=========================================================================
* p2_mcimp.do -- SKEPTIC Monte Carlo for mcci-01.
* Three arms on the SAME draws:
*   A impvalues cluster(unit)  (cluster FINER than i=cg) + ES post horizons
*   B impvalues cluster(cg)    (cluster == i : textbook Ibragimov-Muller)
*   C analytic  cluster(unit)
* Decisive statistic: mean(se) versus sd(b) ACROSS DRAWS for each arm, which
* needs no view about which estimand is targeted, plus coverage of each arm's
* interval around that arm's own Monte Carlo mean AND around the DGP truth.
* Own DGP, unlike the seeded prober's: AR(1) errors, within-cohort effect
* heterogeneity, unit random intercepts not absorbed by i=cg.
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p2_mcimp.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local R 140
local OUT "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage-refute/p2_mcimp.txt"

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

* ---- one diagnostic run: what rows does the table carry? ----------------
mkdat 60 10001
qui wooldid y cg t ttre, att cluster(unit) impvalues esprelength(0) espostlength(2) level(95)
mat DIAG = e(wdidmainresults)
di "DIAG_rows=" rowsof(DIAG) " cols=" colsof(DIAG)
local rn : rownames DIAG
di "DIAG_rownames: `rn'"
local cn : colnames DIAG
di "DIAG_colnames: `cn'"

* ---- the Monte Carlo ----------------------------------------------------
file open fh using "`OUT'", write replace text
file write fh "rep,arm,b,se,lb,ub,dfr,truth_att,truth_ceq" _n
file close fh

forvalues r = 1/`R' {
  mkdat 60 `=20000+`r''

  * --- DGP truths on THIS draw --------------------------------------------
  * (1) att  : cell-weighted mean of truetau over treated post-onset rows
  qui su truetau if !missing(ttre) & t>=ttre
  local truth_att = r(mean)
  * (2) cluster-equally-weighted mean : per treated unit, mean of its own
  *     post-onset truetau; then simple mean over treated units
  tempvar um
  qui egen double `um' = mean(cond(!missing(ttre) & t>=ttre, truetau, .)), by(unit)
  qui su `um' if !missing(ttre), meanonly
  qui bys unit: gen byte __f = (_n==1)
  qui su `um' if !missing(ttre) & __f==1, meanonly
  local truth_ceq = r(mean)
  drop __f
  tempfile dd
  qui save "`dd'"

  foreach arm in A B C {
    use "`dd'", clear
    local ok = 1
    if "`arm'"=="A" cap qui wooldid y cg t ttre, att cluster(unit) impvalues esprelength(0) espostlength(2) level(95)
    if "`arm'"=="B" cap qui wooldid y cg t ttre, att cluster(cg)   impvalues level(95)
    if "`arm'"=="C" cap qui wooldid y cg t ttre, att cluster(unit) level(95)
    if _rc != 0 local ok = 0
    if `ok' == 1 {
      mat T = e(wdidmainresults)
      local dfr = e(df_r)
      local rn2 : rownames T
      * att row is row 1 on every arm
      local bb = T[1,1]
      local ss = T[1,2]
      local ll = T[1,5]
      local uu = T[1,6]
      file open fh using "`OUT'", write append text
      file write fh "`r',`arm'," (`bb') "," (`ss') "," (`ll') "," (`uu') "," (`dfr') "," (`truth_att') "," (`truth_ceq') _n
      file close fh
      * arm A also writes the event-study horizon rows
      if "`arm'"=="A" {
        local nr = rowsof(T)
        forvalues j = 2/`nr' {
          local thisname : word `j' of `rn2'
          local b2 = T[`j',1]
          local s2 = T[`j',2]
          local l2 = T[`j',5]
          local u2 = T[`j',6]
          file open fh using "`OUT'", write append text
          file write fh "`r',ES_`thisname'," (`b2') "," (`s2') "," (`l2') "," (`u2') "," (`dfr') "," (`truth_att') "," (`truth_ceq') _n
          file close fh
        }
      }
    }
    else {
      file open fh using "`OUT'", write append text
      file write fh "`r',`arm',.,.,.,.,.,`truth_att',`truth_ceq'" _n
      file close fh
    }
  }
  if mod(`r',20)==0 di "PROGRESS rep=`r'"
}

di "ALLDONE_MC"
