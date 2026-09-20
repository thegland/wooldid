*=========================================================================
* p3_mcband.do -- SKEPTIC Monte Carlo for mcci-02 (simultaneousbands).
* DGP chosen so that EVERY weighting scheme gives the same horizon truth:
* the true dynamic effect depends on the RELATIVE period only, identically
* across cohorts and units, so pre-horizon truth is 0 and post-horizon
* truth is tau_k exactly.  That removes any estimand ambiguity from the
* coverage measurement and leaves only the band's calibration.
*
* Recorded per draw per horizon: b, se, band lb/ub, truth, plus the common
* critical value e(esband_crit_att) and e(df_r).  From those we compute
*   - measured SIMULTANEOUS coverage of the reported band
*   - the ORACLE critical value  = the level() quantile across draws of
*     max_h |b_h - truth_h| / se_h
*   - pointwise coverage per horizon at invttail(df_r, .05)
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p3_mcband.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local R 250
local OUT "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage-refute/p3_mcband.txt"

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
  gen double ui = rnormal(0,0.70)
  expand 12
  bys unit: gen int t = 2000 + _n - 1
  gen double e = .
  bys unit (t): replace e = rnormal(0,0.55) if _n==1
  bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
  * effect depends on RELATIVE period only, identical across cohorts/units
  gen double rel = t - ttre
  gen double truetau = 0
  replace truetau = 0.60 if rel==0
  replace truetau = 0.90 if rel==1
  replace truetau = 1.10 if rel>=2 & rel<.
  gen double y = ui + 0.25*(t-2000) + truetau + e
end

* ---- one diagnostic run --------------------------------------------------
mkdat 12 55001
qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(2) simultaneousbands level(90)
di "DIAG_engine=" e(esband_engine)
di "DIAG_crit="   e(esband_crit_att)
di "DIAG_level="  e(esband_level)
di "DIAG_ndraw="  e(esband_ndraw)
di "DIAG_dfr="    e(df_r)
mat BB = e(esband_att)
di "DIAG_bandcols=" colsof(BB)
local bcn : colnames BB
di "DIAG_bandcolnames: `bcn'"
mat list BB
mat T = e(wdidmainresults)
local rn : rownames T
di "DIAG_tablerows: `rn'"

* ---- the Monte Carlo ------------------------------------------------------
file open fh using "`OUT'", write replace text
file write fh "rep,horiz,b,se,blb,bub,truth,crit,dfr" _n
file close fh

forvalues r = 1/`R' {
  mkdat 12 `=55000+`r''
  cap qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(2) simultaneousbands level(90)
  if _rc == 0 {
    local crit = e(esband_crit_att)
    local dfr  = e(df_r)
    cap mat BB = e(esband_att)
    if _rc == 0 {
      mat T = e(wdidmainresults)
      local rn : rownames T
      local bcn : colnames BB
      local nb = colsof(BB)
      forvalues j = 1/`nb' {
        local hz : word `j' of `bcn'
        local relv = BB[3,`j']
        local blb  = BB[1,`j']
        local bub  = BB[2,`j']
        * find this horizon in the main table by matching relyr (col 7)
        local bb = .
        local ss = .
        forvalues k = 1/`=rowsof(T)' {
          if T[`k',7] == `relv' & T[`k',7] < . {
            local bb = T[`k',1]
            local ss = T[`k',2]
          }
        }
        * truth: 0 for pre horizons; 0.60 / 0.90 / 1.10 for rel 0 / 1 / >=2
        local tr = 0
        if `relv' == 0 local tr = 0.60
        if `relv' == 1 local tr = 0.90
        if `relv' >= 2 & `relv' < . local tr = 1.10
        file open fh using "`OUT'", write append text
        file write fh "`r'," (`relv') "," (`bb') "," (`ss') "," (`blb') "," (`bub') "," (`tr') "," (`crit') "," (`dfr') _n
        file close fh
      }
    }
  }
  if mod(`r',25)==0 di "PROGRESS rep=`r'"
}

di "ALLDONE_BAND"
