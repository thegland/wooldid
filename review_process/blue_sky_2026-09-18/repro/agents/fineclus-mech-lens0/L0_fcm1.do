*======================================================================
* LENS 0 -- INDEPENDENT REPRODUCTION of FCM-1, from the claim alone.
* My own DGP: 20 school DISTRICTS x 5 SCHOOLS = 100 schools, T = 2010..2019.
* i = school (unit-level), cluster(district) is COARSER than i.
* The district label is attached AFTER y is generated, so it enters the
* command ONLY through cluster(): no relabelling can move a point estimate.
*
* FOUR ARMS, one bootstrap seed, one reps():
*   A  baseline                                   -> nests
*   B  1 school changes district from 2015 on      -> straddles
*   D  10 schools change district from 2015 on     -> straddles
*   C  the SAME 10 schools, SAME target districts,
*      but for ALL t                               -> still nests   <-- CONTROL
* C vs D is the discriminating pair: same schools, same targets, the only
* difference is whole-panel vs mid-panel.  If C ~ A and D ~ B, the driver is
* the global flag, not the data change.
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

prog def mkdat
  syntax , [seed(int 5150) nd(int 20) ns(int 5) tt(int 10) k(int 0) whole(int 0)]
  clear
  set seed `seed'
  set sortseed `seed'
  set obs `=`nd'*`ns''
  gen long school = _n
  gen long distr  = ceil(_n/`ns')
  gen int onset = .
  replace onset = 2013 if mod(school,4)==1
  replace onset = 2015 if mod(school,4)==2
  replace onset = 2017 if mod(school,4)==3
  by distr, sort: gen double dshk = rnormal()*0.7  if _n==1
  by distr:       gen double dtau = 1 + rnormal()*0.9 if _n==1
  by distr:       replace dshk = dshk[1]
  by distr:       replace dtau = dtau[1]
  gen double stau = dtau + rnormal()*0.35
  gen double as   = rnormal()*0.6
  expand `tt'
  bysort school (distr): gen int year = 2009 + _n
  gen double e = rnormal()*0.35
  gen byte D = (year >= onset) & !missing(onset)
  gen double y = as + dshk + 0.04*(year-2010) + stau*D + e
  * relabel k schools (the FIRST school of districts 1..k) into district+1.
  * whole(0): only from 2015 on (mid-panel administrative redraw) -> straddle
  * whole(1): for every year (the same move, recorded retrospectively) -> nests
  if `k' > 0 {
    forvalues q = 1/`k' {
      local src = `q'
      local tgt = cond(`q'==`nd', 1, `q'+1)
      local sid = `=(`q'-1)*`ns'+1'
      if `whole'==0 qui replace distr = `tgt' if school==`sid' & year>=2015
      if `whole'==1 qui replace distr = `tgt' if school==`sid'
    }
  }
  sort school year
end

prog def regpred
  syntax varlist(min=2 max=2)
  tokenize `varlist'
  tempvar smp mn mx bad tg
  qui gen byte `smp' = !missing(`2')
  qui gegen double `mn' = min(`2') if `smp', by(`1')
  qui gegen double `mx' = max(`2') if `smp', by(`1')
  cap assert `mn' == `mx' if `smp'
  local p = (_rc==0)
  qui gen byte `bad' = (`mn' != `mx')
  qui gegen byte `tg' = tag(`1')
  qui count if `tg'==1 & `bad'==1
  dis "  PREDNEST=`p'  straddlers=" r(N)
  c_local PRED `p'
end

prog def grab
  syntax , tag(string)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  local K = colsof(`b')
  forvalues q=1/`K' {
     local this : word `q' of `nm'
     foreach w in main_att main_att_i main_att_it {
       if "`this'"=="`w'" {
          dis "`tag'|`w'|b=" %14.12f `b'[1,`q'] "|bhex=" %21x `b'[1,`q'] "|se=" %14.12f sqrt(`V'[`q',`q'])
       }
     }
  }
end

*---------------------------------------------------------------- step 0
* what the estimand columns are actually called on this build
mkdat , k(0)
qui wooldid y school year onset , cluster(distr) att att_i att_it
dis "COLNAMES: " _col(12) "`: colnames e(b)'"
grab , tag(ANALYTIC_A)

*---------------------------------------------------------------- step 1
* analytic arm, all four datasets: confirms the estimation data is identical
foreach arm in A B D C {
  if "`arm'"=="A" mkdat , k(0)
  if "`arm'"=="B" mkdat , k(1)  whole(0)
  if "`arm'"=="D" mkdat , k(10) whole(0)
  if "`arm'"=="C" mkdat , k(10) whole(1)
  dis "=== ARM `arm' (analytic) ==="
  regpred school distr
  qui count
  local NROW = r(N)
  qui wooldid y school year onset , cluster(distr) att att_i att_it
  grab , tag(ANA_`arm')
  dis "ANA_`arm'|Nrows=`NROW'|N=" e(N) "|NCLUST=" e(N_clust)
}

*---------------------------------------------------------------- step 2
* bootstrap arm, identical seed / reps across all four
foreach arm in A B D C {
  if "`arm'"=="A" mkdat , k(0)
  if "`arm'"=="B" mkdat , k(1)  whole(0)
  if "`arm'"=="D" mkdat , k(10) whole(0)
  if "`arm'"=="C" mkdat , k(10) whole(1)
  regpred school distr
  local P = `PRED'
  cap noisily wooldid y school year onset , cluster(distr) att att_i att_it ///
      bootstrap(standard) reps(250) rseed(31415) dots(0)
  dis "BOOTRC_`arm'=" _rc
  grab , tag(BOOT_`arm')
  dis "BOOT_`arm'|PREDNEST=`P'|N=" e(N) "|NCLUST=" e(N_clust) "|completed=" e(boot_completed) "|failed=" e(boot_failed)
  * silence check: any e() name that records the regime?
  local mm : all macros
  local nfound 0
  foreach e in `: e(macros)' {
    if strpos("`e'","nest")|strpos("`e'","reid")|strpos("`e'","regime") local ++nfound
  }
  foreach e in `: e(scalars)' {
    if strpos("`e'","nest")|strpos("`e'","reid")|strpos("`e'","regime") local ++nfound
  }
  dis "EREC_`arm'|names_matching_nest_reid_regime=`nfound'"
}

*---------------------------------------------------------------- step 3
* what the data itself says the pairs-cluster SE of att_i ought to be:
* SD over districts of the district mean of the TREATED schools' stau, / sqrt(G)
mkdat , k(0)
qui keep if !missing(onset)
qui bysort distr school: keep if _n==1
qui bysort distr: egen double dm = mean(stau)
qui bysort distr: keep if _n==1
qui count
local G = r(N)
qui summ dm
dis "DGPCHECK|G=`G'|SD_of_district_mean_tau=" %9.6f r(sd) "|between_component=" %9.6f r(sd)/sqrt(`G')

dis "L0_FCM1_DONE"
