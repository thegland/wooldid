* probe02_arms.do -- MY OWN test of the FP1 claim, built from the claim's own terms.
* Claim: under permonsetcohorts a covariate constant within the observed i (= the onset
* cohort) is INERT for the point estimate but changes the permutation reference
* distribution, because repopts passes the ORIGINAL covariate while the varlist gets the
* REBUILT i.  Control: with cluster() == i there is no rebuild and the option must be
* inert everywhere.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
which wooldid

*----- data: same shape as mkdata_R.do but with a MODERATE effect so p is interior
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
gen double cohbase = cond(coh==1,0.80,cond(coh==2,1.45,2.10))
gen byte wavegrp = (coh==1)
gen double post = (onset<. & t>=onset)
set seed 771131
gen double y = a_s + g_s*(t-2006) + 0.35*cohbase + 0.012*cohbase*(t-2006) ///
               + 0.45*rnormal() + 0.22*post
compress
save panelR2.dta, replace
tab coh state

*----- the observed fit alone: is the covariate really inert? (no permute)
foreach A in NONE CTRL FEWAVE FECOH CCC FESTATE {
  use panelR2.dta, clear
  local o
  if "`A'"=="CTRL"    local o controls(cohbase)
  if "`A'"=="FEWAVE"  local o fe(wavegrp)
  if "`A'"=="FECOH"   local o fe(coh)
  if "`A'"=="CCC"     local o coarsecohortcontrols(c.cohbase)
  if "`A'"=="FESTATE" local o fe(state)
  cap noi qui wooldid y coh t onset, att cluster(state) `o'
  dis "OBSONLY_`A'_RC=" _rc
  if _rc==0 {
    mat B = e(b)
    dis "OBSONLY_`A'_B %21x = " %21x B[1,1]
  }
}

*----- the permutation arms
cap prog drop runarm
prog def runarm
  syntax , TAG(string) CLUS(string) [OPTS(string)]
  use panelR2.dta, clear
  cap noi qui wooldid y coh t onset, att cluster(`clus') permute(shuffle) permexact permonsetcohorts `opts' permsave(S_`tag', replace)
  local rc = _rc
  dis "ARM_`tag'_RC=`rc'"
  if `rc'==0 {
    mat B = e(b)
    mat C = e(perm_comparator_b)
    mat D = e(permdist_b)
    dis "ARM_`tag'_B %21x = " %21x B[1,1]
    dis "ARM_`tag'_CMP %21x = " %21x C[1,1]
    dis "ARM_`tag'_P = " %14.10f e(perm_p_att)
    dis "ARM_`tag'_SD = " %14.10f e(perm_sd_att)
    dis "ARM_`tag'_NFAIL = " e(perm_failed)
    dis "ARM_`tag'_REBUILT = " e(perm_irebuilt)
    dis "ARM_`tag'_R = " e(perm_reps)
    preserve
    clear
    svmat double D, names(dd)
    gen int _rep = _n
    rename dd1 b_`tag'
    sort _rep
    save D_`tag'.dta, replace
    restore
  }
end

timer clear
timer on 1
runarm , tag(BASE)   clus(state)
runarm , tag(CTRL)   clus(state) opts(controls(cohbase))
runarm , tag(FEWAVE) clus(state) opts(fe(wavegrp))
runarm , tag(FECOH)  clus(state) opts(fe(coh))
runarm , tag(CCC)    clus(state) opts(coarsecohortcontrols(c.cohbase))
runarm , tag(FESTATE) clus(state) opts(fe(state))
* CONTROL: cluster at i -> no rebuild
runarm , tag(CBASE)  clus(coh)
runarm , tag(CCTRL)  clus(coh) opts(controls(cohbase))
timer off 1
timer list 1

*----- paired per-assignment comparison of the reference distributions
cap noi {
  use D_BASE.dta, clear
  merge 1:1 _rep using D_CTRL.dta, nogen
  merge 1:1 _rep using D_FEWAVE.dta, nogen
  merge 1:1 _rep using D_FECOH.dta, nogen
  merge 1:1 _rep using D_CCC.dta, nogen
  gen double d_ctrl = b_CTRL - b_BASE
  gen double d_wave = b_FEWAVE - b_BASE
  gen double d_coh  = b_FECOH  - b_BASE
  gen double d_ccc  = b_CCC    - b_BASE
  sum d_ctrl d_wave d_coh d_ccc, det
  count if abs(d_ctrl) > 1e-12
  dis "NDIFF_CTRL=" r(N) " of " _N
  count if abs(d_wave) > 1e-12
  dis "NDIFF_WAVE=" r(N) " of " _N
  count if abs(d_coh) > 1e-12
  dis "NDIFF_FECOH=" r(N) " of " _N
  save PAIRED.dta, replace
}
dis "PAIRED_RC=" _rc

*----- do the two arms enumerate the SAME assignments? (proves the pairing)
cap noi {
  use S_BASE.dta, clear
  keep _rep onset_*
  ds onset_*
  local ov `r(varlist)'
  merge 1:1 _rep using S_CTRL.dta, keepusing(_rep `ov') nogen update replace
  dis "ASGMATCH: merged without conflict"
}
dis "ASG_RC=" _rc
cap noi {
  use S_BASE.dta, clear
  keep _rep onset_*
  foreach v of varlist onset_* {
    rename `v' base_`v'
  }
  merge 1:1 _rep using S_CTRL.dta, nogen
  gen byte mismatch = 0
  foreach v of varlist base_onset_* {
    local raw = subinstr("`v'","base_","",1)
    replace mismatch = 1 if `v' != `raw' & !(missing(`v') & missing(`raw'))
  }
  count if mismatch
  dis "ASG_MISMATCH_ROWS=" r(N) " of " _N
}
dis "ASG2_RC=" _rc
dis "PROBE02 DONE"
