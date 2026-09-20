clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
* Homogeneous treatment effect => EVERY estimand equals tau exactly => truth unambiguous
local TAU 1
local R 200
postfile pf arm draw double b double se double lb double ub using res, replace
forvalues d = 1/`R' {
  quietly {
    clear
    set seed `=90000 + `d''
    set sortseed `=90000 + `d''
    set obs 40
    gen long st = _n
    gen double a_st = rnormal(0,2)            // sub-i permanent heterogeneity
    gen double u_st = rnormal(0,0.8)          // unit-level shock => cluster(st) is the right level
    gen ttre = .
    replace ttre = 5 if st<=10
    replace ttre = 8 if st>10 & st<=20
    gen wave = cond(ttre<., ttre, 99)
    expand 12
    bys st: gen t = _n
    gen double y = a_st + 0.2*t + `TAU'*(t>=ttre & ttre<.) + u_st + rnormal(0,1)
    gen double uu = runiform()
    drop if uu < 0.08                          // ordinary ragged admin panel
    drop uu
  }
  local a = 1
  foreach spec in "wave t ttre, att cluster(st)" "wave t ttre, att cluster(st) fe(st)" "st t ttre, att cluster(st)" {
    capture quietly wooldid y `spec'
    if _rc == 0 {
      capture quietly matrix M = e(wdidmainresults)
      if _rc == 0 {
        post pf (`a') (`d') (M[1,1]) (M[1,2]) (M[1,5]) (M[1,6])
      }
    }
    local a = `a' + 1
  }
  if mod(`d',25)==0 noisily dis "draw `d' done"
}
postclose pf
use res, clear
gen byte cov = (lb <= `TAU') & (ub >= `TAU')
label define A 1 "A: i=wave cluster(st)  [remedy, no fe]" 2 "B: i=wave cluster(st) fe(st)  [full recipe]" 3 "C: i=st cluster(st)  [at-i control]"
label values arm A
dis "===== RESULTS: true att = `TAU' ====="
table arm, stat(n b) stat(mean b) stat(sd b) stat(mean se) stat(mean cov)
collapse (count) n=b (mean) mb=b (sd) sdb=b (mean) mse=se (mean) covg=cov, by(arm)
gen double sesd = mse/sdb
gen double bias = mb - `TAU'
list arm n mb bias sdb mse sesd covg, noobs abbrev(12)
dis "FEMC_DONE"
