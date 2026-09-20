* r06_whyrank.do -- SELF-CONTAINED. Why is the ANALYTIC ptp pre block rank-deficient
* on design A while the bootstrap block is full rank?  If the deficiency disappears
* under robust / a finer cluster, it is the DEGENERATE CLUSTERED analytic covariance
* (wooldid's own documented cluster-at-or-above-i degeneracy), not a genuine linear
* dependency among the pre coefficients - and the bootstrap's full rank is correct.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 30
gen state = _n
gen onset = .
replace onset = 7 + mod(_n, 6) if mod(_n,3)==0
expand 2
bysort state: gen cty = _n
gen county = state*10 + cty
expand 16
bysort county: gen t = _n
bysort county (t): gen double ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
drop rel te
gen long cellid = county*100 + t

cap prog drop rk
prog def rk
  args tag
  tempname b V S
  matrix `b' = e(b_es_att)
  matrix `V' = e(V_es_att)
  local ln : colnames `b'
  local k 0
  local idx
  local j 0
  foreach l of local ln {
    local ++j
    if substr("`l'",1,3)=="pre" {
      local ++k
      local idx `idx' `j'
    }
  }
  matrix `S' = J(`k',`k',0)
  local r 0
  foreach i1 of local idx {
    local ++r
    local c 0
    foreach i2 of local idx {
      local ++c
      matrix `S'[`r',`c'] = `V'[`i1',`i2']
    }
  }
  matrix `S' = (`S'+`S'')/2
  mata: st_local("eigs", invtokens(strofreal(sort(symeigenvalues(st_matrix("`S'"))',1)',"%12.4e")))
  tempname Si
  matrix `Si' = invsym(`S')
  dis "WHY|`tag'|k=`k'|keff=" `k'-diag0cnt(`Si') "|keff_stored=" e(roth_keff_att) "|eig=`eigs'"
end

local OPT "esfixedbaseperiod esprelength(6) espostlength(4) pretrendpower(80)"
cap noi wooldid y county t onset, cluster(state) `OPT'
dis "RC_cl_state=" _rc
rk CLUSTER_STATE
cap noi wooldid y county t onset, cluster(county) `OPT'
dis "RC_cl_county=" _rc
rk CLUSTER_COUNTY
cap noi wooldid y county t onset, robust `OPT'
dis "RC_robust=" _rc
rk ROBUST
cap noi wooldid y county t onset, `OPT'
dis "RC_iid=" _rc
rk IID
dis "DONE"
