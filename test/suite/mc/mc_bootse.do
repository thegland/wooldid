* mc_bootse.do - Monte Carlo calibration of the bootstrap standard errors.
* Feature-scoped (NOT part of run_full.sh's components; see mc/README).
* Run from test/suite/mc/ : "<stata>" -e do mc_bootse.do
* Measures, per design and estimator: mean reported SE / true sampling SD of the
* att estimate ("ratio", 1.0 = calibrated), and empirical 95% CI coverage.
* Design A (primary) clusters strictly ABOVE the unit; design B is the cluster-at-i
* contrast that exposes the bootalpha(size) shortfall. Authored 2026-08-30.
clear all
set more off
set linesize 200
set graphics off

local ADO `1'
if `"`ADO'"' == "" local ADO ../../../wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
qui do `"`ADO'"'

local DESIGN `2'
if "`DESIGN'" == "" local DESIGN A
local FIRST `3'
if "`FIRST'" == "" local FIRST 1
local LAST `4'
if "`LAST'" == "" local LAST 200
local OUTDIR `5'
if `"`OUTDIR'"' == "" local OUTDIR .
local UA `6'
if "`UA'" == "" local UA 2
local REPS `7'
if "`REPS'" == "" local REPS 199
local DEBUG `8'
if "`DEBUG'" == "" local DEBUG 0
local PLAINONLY `9'
if "`PLAINONLY'" == "" local PLAINONLY 0
local G `10'
if "`G'" == "" local G 24

local TRUEATT = 2
local T = 6
if "`DESIGN'" == "B" local UA = 1
local CVAR = cond("`DESIGN'" == "A", "cl", "unit")
local CLSHOCK = cond("`DESIGN'" == "A", 1, 0)
local SUBEVERY = 5
local TAG = cond(`PLAINONLY', "_plain", "")
if `G' != 24 local TAG `TAG'_G`G'
local OUTFILE `"`OUTDIR'/mc_bootse_`DESIGN'`TAG'.csv"'

di "MCB_CONFIG design=`DESIGN' G=`G' U=`UA' T=`T' cvar=`CVAR' reps=`REPS' draws=`FIRST'..`LAST'"
di "MCB_OUTFILE `OUTFILE'"

*-------------------------------------------------------------------- DGP
cap prog drop mcdgp
prog def mcdgp
  args seed G U T clshock
  clear
  qui set obs `=`G'*`U''
  gen long unit = _n
  gen long cl   = ceil(unit/`U')
  gen byte uw   = unit - (cl-1)*`U'
  gen int ttre = .
  qui replace ttre = 3 if mod(cl,3)==1
  qui replace ttre = 5 if mod(cl,3)==2
  qui expand `T'
  bysort unit: gen int t = _n
  set seed `seed'
  * unit effect (absorbed by the unit FEs; keeps the panel realistic)
  qui by unit: gen double mu = rnormal(0,0.5) if _n==1
  qui by unit: replace mu = mu[1]
  * within-unit AR(1), stationary start
  gen double eps = rnormal(0,0.7)
  qui by unit: replace eps = eps/sqrt(1-0.5^2) if _n==1
  qui by unit: replace eps = 0.5*eps[_n-1] + eps if _n>1
  * cluster-by-period shock: the intra-cluster dependence that survives unit FEs
  gen double s = 0
  if `clshock' {
    tempvar clt
    qui egen long `clt' = group(cl t)
    sort `clt'
    qui by `clt': replace s = rnormal(0,0.5) if _n==1
    qui by `clt': replace s = s[1]
    sort unit t
  }
  gen double y = 2*(t>=ttre & ttre<.) + 0.05*t + mu + eps + s
end

* wrapper for Stata's bootstrap prefix: rebuilds a UNIQUE unit id inside each
* replication (idcluster() renames the resampled CLUSTER; the unit ids inside a
* twice-drawn cluster would otherwise collide and be pooled into one unit).
cap prog drop mcwrap
prog def mcwrap, eclass
  syntax , cv(string)
  tempvar uid
  cap confirm variable bsnewid
  if _rc {
    tempvar bsn
    qui gen long `bsn' = `cv'
    local CV `bsn'
  }
  else local CV bsnewid
  qui egen long `uid' = group(`CV' uw)
  qui wooldid y `uid' t ttre, att cluster(`CV') safety(off)
end

*-------------------------------------------------------------- output file
tempname fh M M2 M3 M4 PB PV
cap confirm file `"`OUTFILE'"'
if _rc {
  file open `fh' using `"`OUTFILE'"', write replace
  file write `fh' "design,draw,dgpseed,nobs,b,std_rc,std_se,std_lb,std_ub,std_nv,std_analse," ///
    "bayeq_rc,bayeq_se,bayeq_lb,bayeq_ub,bayeq_nv,baysz_rc,baysz_se,baysz_lb,baysz_ub,baysz_nv," ///
    "plain_rc,plain_b,plain_se,plain_lb,plain_ub,pfx_rc,pfx_b,pfx_se,pfx_lb,pfx_ub,secs" _n
  file close `fh'
}

*-------------------------------------------------------------------- loop
forval m = `FIRST'/`LAST' {
  local dgpseed = 700000 + cond("`DESIGN'"=="A", 0, 300000) + `m'
  mcdgp `dgpseed' `G' `UA' `T' `CLSHOCK'
  qui count
  local nobs = r(N)
  local sub = (mod(`m', `SUBEVERY') == 1)

  timer clear 9
  timer on 9

  foreach a in std bayeq baysz plain pfx {
    local `a'_rc = .
    local `a'_se = .
    local `a'_lb = .
    local `a'_ub = .
    local `a'_nv = .
  }
  local b = .
  local std_analse = .
  local plain_b = .
  local pfx_b = .

  * ---- arm 1: bootstrap(standard); also yields the hidden analytic covariance
  if `PLAINONLY' local std_rc = .
  if !`PLAINONLY' {
  cap matrix drop VA
  cap qui wooldid y unit t ttre, att cluster(`CVAR') bootstrap(standard) ///
      reps(`REPS') rseed(`=800000 + `m'') dots(0) safety(off)
  local std_rc = _rc
  if `std_rc' == 0 {
    * the Wald family (estimate +/- z * bootstrap SE): the arms compare a normal-theory
    * interval with the prefix's; the displayed table carries the percentile family by
    * default since 2026-09-04
    matrix `M' = e(boot_wald_mainresults)
    local b        = el(`M',1,1)
    local std_se   = el(`M',1,2)
    local std_lb   = el(`M',1,5)
    local std_ub   = el(`M',1,6)
    local std_nv   = e(boot_N_V)
    cap matrix VA = e(V_analytic)
    if _rc == 0 local std_analse = sqrt(el(VA,1,1))
    if `DEBUG' == 1 & `m' == `FIRST' {
      di "MCB_DEBUG rownames: " "`: rownames `M''"
      matrix list `M'
    }
  }

  * ---- arm 2: bootstrap(bayes), default bootalpha(equal)
  cap qui wooldid y unit t ttre, att cluster(`CVAR') bootstrap(bayes) ///
      reps(`REPS') rseed(`=810000 + `m'') dots(0) safety(off)
  local bayeq_rc = _rc
  if `bayeq_rc' == 0 {
    matrix `M2' = e(boot_wald_mainresults)
    local bayeq_se = el(`M2',1,2)
    local bayeq_lb = el(`M2',1,5)
    local bayeq_ub = el(`M2',1,6)
    local bayeq_nv = e(boot_N_V)
  }

  * ---- arm 3: bootstrap(bayes) bootalpha(size); every draw on B, subset on A
  if "`DESIGN'" == "B" | `sub' {
    cap qui wooldid y unit t ttre, att cluster(`CVAR') bootstrap(bayes) ///
        bootalpha(size) reps(`REPS') rseed(`=820000 + `m'') dots(0) safety(off)
    local baysz_rc = _rc
    if `baysz_rc' == 0 {
      matrix `M3' = e(boot_wald_mainresults)
      local baysz_se = el(`M3',1,2)
      local baysz_lb = el(`M3',1,5)
      local baysz_ub = el(`M3',1,6)
      local baysz_nv = e(boot_N_V)
    }
  }
  }

  * ---- arm 4: analytic clustered VCE (no bootstrap)
  cap qui wooldid y unit t ttre, att cluster(`CVAR') safety(off)
  local plain_rc = _rc
  if `plain_rc' == 0 {
    matrix `M4' = e(wdidmainresults)
    local plain_b  = el(`M4',1,1)
    local plain_se = el(`M4',1,2)
    local plain_lb = el(`M4',1,5)
    local plain_ub = el(`M4',1,6)
    if missing(`b') local b = `plain_b'
  }

  * ---- arm 5: Stata's official bootstrap prefix (subset draws)
  if `sub' & !`PLAINONLY' {
    cap qui bootstrap _b[main_att], cluster(`CVAR') idcluster(bsnewid) ///
        reps(`REPS') seed(`=830000 + `m'') nodots nowarn: mcwrap, cv(`CVAR')
    local pfx_rc = _rc
    if `pfx_rc' == 0 {
      matrix `PB' = e(b)
      matrix `PV' = e(V)
      local pfx_b  = el(`PB',1,1)
      local pfx_se = sqrt(el(`PV',1,1))
      local pfx_lb = `pfx_b' - 1.959963984540054*`pfx_se'
      local pfx_ub = `pfx_b' + 1.959963984540054*`pfx_se'
      if `DEBUG' == 1 & `m' == `FIRST' {
        di "MCB_DEBUG pfx b colnames: " "`: colnames `PB''"
      }
    }
  }

  timer off 9
  qui timer list 9
  local secs = r(t9)

  local row "`DESIGN',`m',`dgpseed',`nobs'"
  foreach v in b std_rc std_se std_lb std_ub std_nv std_analse ///
               bayeq_rc bayeq_se bayeq_lb bayeq_ub bayeq_nv ///
               baysz_rc baysz_se baysz_lb baysz_ub baysz_nv ///
               plain_rc plain_b plain_se plain_lb plain_ub ///
               pfx_rc pfx_b pfx_se pfx_lb pfx_ub secs {
    local x = strtrim(string(``v'', "%18.12g"))
    if "`x'" == "" local x .
    local row "`row',`x'"
  }
  file open `fh' using `"`OUTFILE'"', write append
  file write `fh' `"`row'"' _n
  file close `fh'

  di "MCB_DRAW `DESIGN' m=`m' b=" %7.4f `b' " std=" %7.4f `std_se' ///
     " bayeq=" %7.4f `bayeq_se' " anal=" %7.4f `plain_se' " secs=" %6.1f `secs'
}

*----------------------------------------------------------------- summary
* Reads this run's CSV and, when it exists, the cheap plain-only CSV on the SAME dgp
* seeds. The plain pass is what makes the RATIO precise: mean(SE)/SD(att-hat) is
* limited mostly by the SD in the denominator, and the plain arm is ~200x cheaper per
* draw, so it can afford 10-20x the draws the bootstrap arms can.
tempfile plaindta
local HASPLAIN 0
cap confirm file `"`OUTDIR'/mc_bootse_`DESIGN'_plain.csv"'
if _rc == 0 {
  qui import delimited using `"`OUTDIR'/mc_bootse_`DESIGN'_plain.csv"', clear varnames(1) case(preserve)
  qui keep draw plain_b plain_se plain_lb plain_ub
  rename (plain_b plain_se plain_lb plain_ub) (pb pse plb pub)
  qui duplicates drop draw, force
  qui save `plaindta'
  local HASPLAIN 1
}
qui import delimited using `"`OUTFILE'"', clear varnames(1) case(preserve)
qui keep if design == "`DESIGN'"
qui duplicates drop draw, force
if `HASPLAIN' {
  qui merge 1:1 draw using `plaindta', nogen
}
else {
  foreach v in pb pse plb pub {
    qui gen double `v' = .
  }
}
* att-hat is the SAME estimate in every arm, so one column of point estimates serves
* as the sampling-SD denominator for all of them
qui gen double bU = b
qui replace bU = pb if missing(bU)
foreach v in se lb ub {
  qui gen double aplain_`v' = plain_`v'
  qui replace aplain_`v' = p`v' if missing(aplain_`v')
}
qui sum bU
local sdall  = r(sd)
local nall   = r(N)
local meanall = r(mean)
qui sum b
local sdb   = r(sd)
local ndraw = r(N)
di ""
di "MCB_SUMMARY design=`DESIGN' G=`G' U=`UA' T=`T' reps=`REPS' true_att=`TRUEATT'"
di "MCB_NDRAWS bootstrap_arms=`ndraw'  sd_denominator=`nall'"
di "MCB_SD_OF_BHAT bootdraws=" %9.6f `sdb' "  alldraws=" %9.6f `sdall' "  <- ratio denominator"
if `nall' > 1 di "MCB_MEAN_BHAT=" %9.6f `meanall' "  bias=" %9.6f `=`meanall'-`TRUEATT'' ///
   "  bias_t=" %6.3f `=(`meanall'-`TRUEATT')/(`sdall'/sqrt(`nall'))'
di "MCB_ARM arm / n / meanSE / ratio / coverage / coverage_mcse_pp"
foreach a in std bayeq baysz plain pfx {
  local sev `a'_se
  local lbv `a'_lb
  local ubv `a'_ub
  if "`a'" == "plain" {
    local sev aplain_se
    local lbv aplain_lb
    local ubv aplain_ub
  }
  qui count if !missing(`sev')
  local na = r(N)
  if `na' > 0 {
    qui sum `sev'
    local mse = r(mean)
    cap drop cov_`a'
    qui gen byte cov_`a' = (`lbv' <= `TRUEATT' & `ubv' >= `TRUEATT') if !missing(`sev') & !missing(`lbv')
    qui sum cov_`a'
    local cv = r(mean)
    local ncv = r(N)
    di "MCB_ARM `a' " %5.0f `na' "  " %9.6f `mse' "  " %7.4f `=`mse'/`sdall'' ///
       "  " %7.4f `cv' "  " %5.2f `=100*sqrt(`cv'*(1-`cv')/`ncv')'
  }
  else di "MCB_ARM `a' 0 - not run"
}
* PAIRED per-draw SE ratios. These share the draw (same data, same att-hat), so they
* are far more precise than comparing two independently noisy mean/SD ratios, and they
* stay informative at the small n the 199-rep arms can afford.
local PAIRS baysz:bayeq pfx:std std:plain bayeq:std
foreach pr of local PAIRS {
  local cut = strpos("`pr'", ":")
  local nu = substr("`pr'", 1, `cut'-1)
  local de = substr("`pr'", `cut'+1, .)
  cap drop prr
  qui gen double prr = `nu'_se/`de'_se if !missing(`nu'_se) & !missing(`de'_se)
  qui count if !missing(prr)
  if r(N) > 0 {
    qui sum prr
    di "MCB_PAIRED `nu'/`de' n=" %5.0f r(N) " mean=" %7.4f r(mean) ///
       " se=" %7.4f `=r(sd)/sqrt(r(N))' " min=" %7.4f r(min) " max=" %7.4f r(max)
  }
  cap drop prr
}
di "MCB_SIZEALPHA_PREDICTION mean_cluster_obs=`=`UA'*`T''  predicted_size/equal=" ///
   %7.4f `=1/sqrt(`UA'*`T')'
* feature claims that must hold EXACTLY, checked on every draw
cap drop dchk
qui gen double dchk = abs(b - plain_b)
qui sum dchk
di "MCB_CHECK point estimates unchanged by bootstrap: maxabsdiff=" %12.4e r(max) " n=" r(N)
cap drop dchk
qui gen double dchk = abs(pfx_b - b)
qui sum dchk
di "MCB_CHECK Stata-prefix point estimate == wooldid: maxabsdiff=" %12.4e r(max) " n=" r(N)
cap drop dchk
qui gen double dchk = abs(std_analse - plain_se)
qui sum dchk
di "MCB_CHECK e(V_analytic) == analytic run SE: maxabsdiff=" %12.4e r(max) " n=" r(N)
foreach v in std_nv bayeq_nv baysz_nv {
  qui sum `v'
  if r(N) > 0 di "MCB_BOOT_NV `v' n=" r(N) " min=" r(min) " mean=" %7.2f r(mean) " (of `REPS' requested)"
}
qui sum secs
if r(N) > 0 di "MCB_SECS_PER_DRAW mean=" %7.1f r(mean) " total_min=" %8.1f `=r(mean)*r(N)/60'
di "MCB_DONE=1"

* ---------------------------------------------------------------------------
* ZERO-COST EXTENSION (not exercised by the run this file documents): every
* bootstrap arm already stores e(boot_ci_percentile) and e(boot_ci_bc) alongside the
* reported normal-based CI, computed from the SAME replication vector. Capturing
* el(e(boot_ci_percentile),1,1..2) and el(e(boot_ci_bc),1,1..2) next to each arm's
* lb/ub would give coverage for all three CI flavours at no extra replication cost -
* the expensive part (199 refits) is already paid. The coverage numbers recorded here
* are for the REPORTED (normal-based) CI only, which is what a user sees by default.
* ---------------------------------------------------------------------------