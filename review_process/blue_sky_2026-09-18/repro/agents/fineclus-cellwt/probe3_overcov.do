* probe3_overcov.do -- MECHANISM of the recorded over-coverage of att_i under
* bootstrap(standard) at cluster() == i: 149/150 = .9933 vs nominal .95 (c8c4e2c,
* 2026-09-16).  The record reports the number and attributes it to "at or above the i
* level the covariance stops seeing the treated cells' own variation" -- but that
* effect makes analytic SEs too SMALL (the help's own measured 0.66-0.75 SE/SD), so it
* cannot be why an interval is too WIDE.  No mechanism is on the record.
*
* HYPOTHESIS.  At cluster() == i the re-id leg resamples COHORTS with replacement and
* re-identifies each drawn copy as fresh units, so the bootstrap distribution of att_i
* carries a cohort-COMPOSITION term Var_s(tau_s)/S_treated.  The Monte Carlo that
* measured the coverage holds tau_s FIXED (in its DGP tau is a deterministic function
* of the state index), so the sampling distribution it measures carries no such term.
* The bootstrap is estimating the variance of a superpopulation-of-cohorts estimand;
* the MC is measuring the variance of a fixed-cohorts one.  Hence over-coverage.
*
* FALSIFIABLE PREDICTION.  Hold the truth of att_i fixed at 0.8 and switch the
* ACROSS-COHORT variance of tau_s on and off:
*   HET  tau = 0.4 (10 states) / 1.2 (10 states), Var_s = 0.16 -> boot SE >> MC SD
*   HOM  tau = 0.8 (all 20 treated states),       Var_s = 0    -> boot SE ~= MC SD
* att_i = 0.8 exactly in both; same noise, same cluster count, same everything else.
* PART 0 makes it stark by turning the noise down to sd 0.02, so the MC sampling SD is
* ~0 in BOTH arms while the composition term is untouched: any bootstrap SE the HET arm
* reports there is composition variance and nothing else.
*
* The design re-implements the published dgpE (40 states, T=6, unit RE + AR(1)) from
* its description so this file stands alone.  Expected composition SE in the HET arm:
* sqrt(Var_s(tau)/S_treated) = sqrt(0.16/20) = 0.0894.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

prog def dgp
  args seed mode noisesd
  clear
  qui set obs 40
  gen long state = _n
  gen int nu = 4
  qui replace nu = 2 if mod(state,4)==1
  qui replace nu = 6 if mod(state,4)==3
  qui expand nu
  bysort state: gen byte uw = _n
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 3 if mod(state,4)==1
  qui replace ttre = 4 if mod(state,4)==3
  if "`mode'" == "HET" gen double tau = 0.2*nu if mod(state,2)==1
  else                 gen double tau = 0.8    if mod(state,2)==1
  qui replace tau = 0 if missing(tau)
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  qui by unit: gen double mu = rnormal(0,0.5*`noisesd') if _n==1
  qui by unit: replace mu = mu[1]
  gen double eps = rnormal(0,0.7*`noisesd')
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.05*t + mu + eps
end

tempname fh
file open `fh' using "p3.csv", write replace
file write `fh' "kind,mode,m,rc,b_atti,se_atti,b_att,se_att" _n
file close `fh'

prog def onerun
  args fh kind mode m cv reps rs
  local rc = _rc
  local bi = .
  local si = .
  local ba = .
  local sa = .
  if "`cv'" == "none" {
    cap qui wooldid y state t ttre, att att_i
    local rc = _rc
  }
  else {
    cap qui wooldid y state t ttre, att att_i cluster(`cv') bootstrap(standard) reps(`reps') rseed(`rs') dots(0)
    local rc = _rc
  }
  if `rc' == 0 {
    mat b = e(b)
    mat V = e(V)
    local bi = b[1,colnumb(b,"main_att_i")]
    local ba = b[1,colnumb(b,"main_att")]
    if "`cv'" != "none" {
      local si = sqrt(V[colnumb(b,"main_att_i"),colnumb(b,"main_att_i")])
      local sa = sqrt(V[colnumb(b,"main_att"),colnumb(b,"main_att")])
    }
  }
  file open `fh' using "p3.csv", write append
  file write `fh' "`kind',`mode',`m',`rc',`bi',`si',`ba',`sa'" _n
  file close `fh'
end

*==================== PART 0: noise scaled to 2% -- composition variance ALONE.
* With the noise off, the sampling SD of att_i around its truth of 0.8 is ~0.003.
* Anything the bootstrap reports here is composition variance.
foreach mode in HET HOM {
  qui dgp 300001 `mode' 0.02
  onerun `fh' "LOWN_NEST" "`mode'" 1 state 99 8001
  qui dgp 300001 `mode' 0.02
  onerun `fh' "LOWN_FINE" "`mode'" 1 unit 99 8001
  qui dgp 300001 `mode' 0.02
  onerun `fh' "LOWN_PT"   "`mode'" 1 none . .
  dis "PART0 `mode' DONE"
}

*==================== PART 1: full noise, MC sampling SD of att_i (fast fits)
foreach mode in HET HOM {
  forval m = 1/70 {
    qui dgp `=300000 + `m'' `mode' 1
    onerun `fh' "MC" "`mode'" `m' none . .
  }
  dis "MC `mode' DONE"
}

*==================== PART 2: full noise, bootstrap SE at cluster(state) == i (NEST)
foreach mode in HET HOM {
  forval m = 1/5 {
    qui dgp `=300000 + `m'' `mode' 1
    onerun `fh' "NEST" "`mode'" `m' state 99 `=8000+`m''
  }
  dis "NEST `mode' DONE"
}

*==================== PART 3: full noise, bootstrap SE at cluster(unit) (keep-i)
forval m = 1/5 {
  qui dgp `=300000 + `m'' HET 1
  onerun `fh' "FINE" "HET" `m' unit 99 `=8000+`m''
}
dis "FINE DONE"

import delimited "p3.csv", clear varnames(1) case(preserve)
bysort kind mode: egen double sd_atti = sd(b_atti)
bysort kind mode: egen double mn_seatti = mean(se_atti)
bysort kind mode: egen double mn_batti = mean(b_atti)
bysort kind mode: gen long nn = _N
bysort kind mode: keep if _n == 1
keep kind mode nn mn_batti sd_atti mn_seatti
list, sep(0) noobs
