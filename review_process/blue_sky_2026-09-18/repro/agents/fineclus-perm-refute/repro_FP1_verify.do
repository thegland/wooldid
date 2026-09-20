*=============================================================================
* repro_FP1_verify.do   SELF-CONTAINED.  Run with:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_FP1_verify.do
*
* WHAT IT SHOWS (independent of the original prober's script and data)
* Under permonsetcohorts the per-replication rebuild of i makes a covariate that is
* CONSTANT WITHIN THE OBSERVED i go live inside every replication, although it is
* provably inert for the reported point estimate.  The reported permutation p-value and
* null distribution therefore come from a different model than the estimate they judge.
*
* Blocks:
*  (1) observed fit only  -> the covariate arms are BITWISE identical to no covariate
*  (2) permutation, cluster(state) STRICTLY FINER than i -> rebuild fires; p and null SD
*      move, and coarsecohortcontrols() moves them enormously
*  (3) permutation, cluster(coh) == i -> NO rebuild; the same covariate is bit-identical
*  (4) contrast: a STATE-level control in coarsecohortcontrols() (varies within the
*      cohort, so it is live in the observed fit too) does NOT blow up -> it is the
*      inert-in-observed / live-in-replication asymmetry, not coarsecohortcontrols()
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*---- a ragged administrative-shaped panel; i = coh IS the onset grouping
clear
set obs 7
gen int state = _n
gen byte coh = 1
replace coh = 2 if inlist(state,4,5)
replace coh = 3 if inlist(state,6,7)
gen int onset = .
replace onset = 2005 if coh==1
replace onset = 2009 if coh==2
gen int t0 = 1998
replace t0 = 2000 if inlist(state,2,6)
replace t0 = 1999 if inlist(state,3,7)
replace t0 = 2001 if state==4
gen int t1 = 2015
replace t1 = 2014 if state==2
replace t1 = 2013 if state==4
replace t1 = 2012 if state==6
gen double a_s = .
replace a_s = 1.20 in 1
replace a_s = 0.40 in 2
replace a_s = 2.10 in 3
replace a_s = -0.60 in 4
replace a_s = 1.70 in 5
replace a_s = 0.10 in 6
replace a_s = -1.30 in 7
gen double g_s = .
replace g_s =  0.030 in 1
replace g_s = -0.015 in 2
replace g_s =  0.055 in 3
replace g_s = -0.040 in 4
replace g_s =  0.010 in 5
replace g_s =  0.025 in 6
replace g_s = -0.005 in 7
expand t1-t0+1
bys state: gen int t = t0 + _n - 1
drop t0 t1
* a genuine COHORT-level control ("the generosity of the reform that wave adopted")
gen double cohbase = cond(coh==1,0.80,cond(coh==2,1.45,2.10))
* an "adoption wave" grouping
gen byte wavegrp = (coh==1)
gen double post = (onset<. & t>=onset)
set seed 771131
gen double y = a_s + g_s*(t-2006) + 0.35*cohbase + 0.012*cohbase*(t-2006) ///
               + 0.45*rnormal() + 0.22*post
compress
tempfile PAN
save `PAN'
dis "N=" _N "   states=7, M = 7!/(3!2!2!) = 210 enumerated assignments"

dis _n "########## (1) OBSERVED FIT ONLY - is the covariate inert? ##########"
foreach A in NONE CTRL FEWAVE FECOH CCC CCCSTATE {
  use `PAN', clear
  local o
  if "`A'"=="CTRL"     local o controls(cohbase)
  if "`A'"=="FEWAVE"   local o fe(wavegrp)
  if "`A'"=="FECOH"    local o fe(coh)
  if "`A'"=="CCC"      local o coarsecohortcontrols(c.cohbase)
  if "`A'"=="CCCSTATE" local o coarsecohortcontrols(c.a_s)
  cap noi qui wooldid y coh t onset, att cluster(state) `o'
  if _rc==0 {
    mat B = e(b)
    dis "OBS_`A'_B %21x = " %21x B[1,1] "     (rc=" _rc ")"
  }
  else dis "OBS_`A'_RC=" _rc
}

cap prog drop arm
prog def arm
  syntax , TAG(string) PAN(string) CLUS(string) [OPTS(string)]
  use `pan', clear
  cap noi qui wooldid y coh t onset, att cluster(`clus') permute(shuffle) permexact permonsetcohorts `opts'
  local rc = _rc
  if `rc'!=0 {
    dis "ARM_`tag'_RC=`rc'"
    exit
  }
  mat B = e(b)
  dis "ARM_`tag': b %21x = " %21x B[1,1] "  p = " %12.8f e(perm_p_att) ///
      "  nullSD = " %12.8f e(perm_sd_att) "  nfail = " e(perm_failed) ///
      "  irebuilt = " e(perm_irebuilt) "  R = " e(perm_reps)
end

dis _n "########## (2) cluster(state) STRICTLY FINER than i - the rebuild fires ##########"
arm , tag(R_BASE)     pan(`PAN') clus(state)
arm , tag(R_CTRL)     pan(`PAN') clus(state) opts(controls(cohbase))
arm , tag(R_FEWAVE)   pan(`PAN') clus(state) opts(fe(wavegrp))
arm , tag(R_FECOH)    pan(`PAN') clus(state) opts(fe(coh))
arm , tag(R_CCC)      pan(`PAN') clus(state) opts(coarsecohortcontrols(c.cohbase))
arm , tag(R_CCCSTATE) pan(`PAN') clus(state) opts(coarsecohortcontrols(c.a_s))

dis _n "########## (3) cluster(coh) == i - NO rebuild (the control arm) ##########"
arm , tag(N_BASE)     pan(`PAN') clus(coh)
arm , tag(N_CTRL)     pan(`PAN') clus(coh) opts(controls(cohbase))
arm , tag(N_CCC)      pan(`PAN') clus(coh) opts(coarsecohortcontrols(c.cohbase))

dis _n "########## (4) paired, assignment by assignment ##########"
use `PAN', clear
qui wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
mat DB = e(permdist_b)
use `PAN', clear
qui wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts controls(cohbase)
mat DC = e(permdist_b)
clear
svmat double DB, names(bb)
svmat double DC, names(cc)
gen double d = cc1 - bb1
gen int _rep = _n
count if abs(d) > 1e-10
dis "REFERENCE STATISTICS THAT CHANGE: " r(N) " of " _N
sort d
list _rep bb1 cc1 d in 1/3, noobs
gen double ad = abs(d)
sort ad
dis "SMALLEST |difference| (the identity assignment) = " %21x ad[1]
dis _n "REPRO DONE"
