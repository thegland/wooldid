* r12b_corr.do -- SELF-CONTAINED. Pin the exact singularity found in r07: print the
* analytic pre-period covariance block and its correlation matrix under iid, robust
* and two cluster levels, on an ordinary panel.
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
tab onset, m


local OPT "esfixedbaseperiod esprelength(6) espostlength(4)"
qui wooldid y county t onset, `OPT'
matrix S_iid = e(V_es_att)[1..6,1..6]
qui wooldid y county t onset, robust `OPT'
matrix S_rob = e(V_es_att)[1..6,1..6]
qui wooldid y county t onset, cluster(state) `OPT'
matrix S_cls = e(V_es_att)[1..6,1..6]
matrix b_es = e(b_es_att)
matrix list b_es
mata:
  nmv = ("S_iid","S_rob","S_cls")
  for (q=1; q<=3; q++) {
    S = st_matrix(nmv[q])
    d = sqrt(diagonal(S))
    printf("\n==== %s ====\n", nmv[q])
    "variances (pre7 pre6 pre5 pre4 pre3 pre2):"
    diagonal(S)'
    "correlation matrix:"
    S :/ (d*d')
    "eigenvalues (ascending):"
    sort(symeigenvalues(makesymmetric(S))', 1)'
  }
end
dis "DONE"
