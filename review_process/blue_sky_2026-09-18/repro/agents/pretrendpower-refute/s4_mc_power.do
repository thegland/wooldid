* s4_mc_power.do -- SELF-CONTAINED (generates its own null AND power arms).
* Does pretrendpower()'s reported mdt_wald(80) deliver 80% power to the joint
* pre-trend test AT THE NOMINAL 5% CUTOFF wooldid itself uses?  And how does the
* reported MDT compare with the design's true minimum detectable trend?
* Own 16-state panel, own seeds, independent of any other probe in this directory.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop skpanel
prog def skpanel
  args delta rho
  clear
  set obs 16
  gen long st = _n
  gen byte coh = 1 + floor((st-1)/4)
  gen int ttre = .
  replace ttre = 6  if coh == 1
  replace ttre = 8  if coh == 2
  replace ttre = 10 if coh == 3
  gen double alpha_s = rnormal()*0.8
  expand 12
  bysort st: gen int t = _n
  sort st t
  gen double yeareff = 0.12*t - 0.02*t*t + 0.25*cos(t)
  gen byte D = (ttre < . & t >= ttre)
  gen double relyr = t - ttre
  gen double trendterm = cond(ttre < ., relyr, 0)
  gen double e = rnormal()
  forvalues k = 2/12 {
    qui replace e = `rho'*e[_n-1] + sqrt(1-`rho'*`rho')*rnormal() if t == `k'
  }
  gen double y = 2 + alpha_s + yeareff + 0.8*D + `delta'*trendterm + e
end

cap prog drop pstat
prog def pstat, rclass
  args lam
  tempname B V S bp Si m w
  matrix `B'  = e(b_es_att)
  matrix `V'  = e(V_es_att)
  matrix `S'  = `V'[1..3, 1..3]
  matrix `S'  = (`S' + `S'')/2
  matrix `bp' = (`B'[1,1] \ `B'[1,2] \ `B'[1,3])
  matrix `Si' = invsym(`S')
  matrix `m'  = (-3 \ -2 \ -1)
  matrix `w'  = `bp''*`Si'*`bp'
  return scalar W  = `w'[1,1]
  matrix `w'  = `m''*`Si'*`m'
  return scalar c2 = `w'[1,1]
  matrix `w'  = `m''*`Si'*`bp'
  return scalar slope = `w'[1,1]/return(c2)
  return scalar mdt80 = sqrt(`lam'/return(c2))
  return scalar b1 = `bp'[1,1]
  return scalar b2 = `bp'[2,1]
  return scalar b3 = `bp'[3,1]
end

local lam80 = npnF(3, 15, invF(3, 15, .95), .20)
local cutW  = 3*invF(3, 15, .95)
dis "lam80 = " `lam80' "    wooldid W cutoff = " `cutW'

* ---------------- arm 0: null ----------------
local R0 250
tempfile n0
postfile P0 int rep double(W mdt80 b1 b2 b3) using `n0', replace
forvalues r = 1/`R0' {
  set seed `= 310000 + 419*`r''
  set sortseed `= 310000 + 419*`r''
  skpanel 0 0.5
  cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
  if _rc != 0 continue
  pstat `lam80'
  post P0 (`r') (r(W)) (r(mdt80)) (r(b1)) (r(b2)) (r(b3))
}
postclose P0
use `n0', clear
export delimited using "s4_null.csv", replace
gen byte rej = (W > `cutW')
qui sum rej
dis "ARM0 null: size of wooldid's joint pre-trend test = " %8.4f r(mean) " (n=" r(N) ")"
_pctile W, p(95)
local cutEmp = r(r1)
dis "ARM0 empirical 95th percentile of W = " `cutEmp' "   (wooldid uses " `cutW' ")"
qui sum mdt80, detail
local dMed = r(p50)
dis "ARM0 median reported mdt_wald(80) = " %12.8f `dMed'
corr b1 b2 b3, cov
matrix SIG = r(C)
matrix MM = (-3 \ -2 \ -1)
matrix QQ = MM'*invsym(SIG)*MM
local c2t = QQ[1,1]
local dTrue = sqrt(`lam80'/`c2t')
dis "ARM0 c2_true = " `c2t' "   mdt_true (wooldid's own lambda, true Sigma) = " %12.8f `dTrue'
dis "ARM0 median mdt_hat / mdt_true = " %8.4f `dMed'/`dTrue'

* ---------------- arms 1 and 2: power ----------------
foreach arm in 1 2 {
  if `arm' == 1 {
    local delta `dMed'
    local R 400
    local base 620000
    local lbl "at the MEDIAN REPORTED mdt_wald(80)"
  }
  if `arm' == 2 {
    local delta `dTrue'
    local R 250
    local base 740000
    local lbl "at mdt_TRUE (the same lambda on the true covariance)"
  }
  tempfile n1
  postfile P1 int rep double(W) using `n1', replace
  forvalues r = 1/`R' {
    set seed `= `base' + 307*`r''
    set sortseed `= `base' + 307*`r''
    skpanel `delta' 0.5
    cap qui wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2)
    if _rc != 0 continue
    pstat `lam80'
    post P1 (`r') (r(W))
  }
  postclose P1
  preserve
  use `n1', clear
  export delimited using "s4_power_arm`arm'.csv", replace
  gen byte rnom = (W > `cutW')
  gen byte remp = (W > `cutEmp')
  qui sum rnom
  local pn = r(mean)
  local nn = r(N)
  qui sum remp
  local pe = r(mean)
  dis "ARM`arm' delta = " %12.8f `delta' "  `lbl'"
  dis "ARM`arm'  power at wooldid's OWN nominal cutoff = " %8.4f `pn' "  (n=" `nn' ", se " %6.4f sqrt(`pn'*(1-`pn')/`nn') ")"
  dis "ARM`arm'  power at the size-corrected empirical cutoff = " %8.4f `pe'
  restore
}
dis "DONE"
