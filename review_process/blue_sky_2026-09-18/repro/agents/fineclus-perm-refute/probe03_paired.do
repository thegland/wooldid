* probe03_paired.do -- BASE arm (correct permsave syntax) + the PAIRED per-assignment
* comparison of the reference distributions.  Two panels (PA = probe02's, PB = a second
* draw) so nothing rests on one dataset.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , SD(string) TE(string) OUT(string)
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
  set seed `sd'
  gen double y = a_s + g_s*(t-2006) + 0.35*cohbase + 0.012*cohbase*(t-2006) ///
                 + 0.45*rnormal() + `te'*post
  compress
  save `out', replace
end

cap prog drop runarm
prog def runarm
  syntax , TAG(string) PAN(string) [OPTS(string)]
  use `pan', clear
  cap noi qui wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts `opts' permsave(S_`tag') replace
  local rc = _rc
  dis "ARM_`tag'_RC=`rc'"
  if `rc'==0 {
    mat B = e(b)
    mat C = e(perm_comparator_b)
    mat D = e(permdist_b)
    dis "ARM_`tag'_B %21x = " %21x B[1,1]
    dis "ARM_`tag'_CMP %21x = " %21x C[1,1]
    dis "ARM_`tag'_BDEC = " %18.14f B[1,1]
    dis "ARM_`tag'_P = " %14.10f e(perm_p_att)
    dis "ARM_`tag'_SD = " %14.10f e(perm_sd_att)
    dis "ARM_`tag'_NFAIL = " e(perm_failed)
    dis "ARM_`tag'_REBUILT = " e(perm_irebuilt)
    scalar bo_`tag' = B[1,1]
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

mkpanel , sd(771131) te(0.22) out(panA.dta)
mkpanel , sd(408822) te(0.22) out(panB.dta)

foreach P in A B {
  dis "#################### PANEL `P' ####################"
  runarm , tag(`P'BASE)   pan(pan`P'.dta)
  runarm , tag(`P'CTRL)   pan(pan`P'.dta) opts(controls(cohbase))
  runarm , tag(`P'FEWAVE) pan(pan`P'.dta) opts(fe(wavegrp))
  runarm , tag(`P'FECOH)  pan(pan`P'.dta) opts(fe(coh))
  runarm , tag(`P'CCC)    pan(pan`P'.dta) opts(coarsecohortcontrols(c.cohbase))
}

foreach P in A B {
  cap noi {
    use D_`P'BASE.dta, clear
    foreach A in CTRL FEWAVE FECOH CCC {
      merge 1:1 _rep using D_`P'`A'.dta, nogen
    }
    local bo = bo_`P'BASE
    dis "PANEL`P'_BOBS = " %18.14f `bo'
    foreach A in CTRL FEWAVE FECOH CCC {
      gen double d_`A' = b_`P'`A'' - b_`P'BASE'
      qui count if abs(d_`A') > 1e-10
      dis "PANEL`P'_NDIFF_`A' = " r(N) " of " _N
      qui sum d_`A'
      dis "PANEL`P'_MEANDIFF_`A' = " %14.10f r(mean) "  MAXABS = " %14.10f max(abs(r(min)),abs(r(max)))
      qui sum b_`P'`A''
      dis "PANEL`P'_SD_`A' = " %14.10f r(sd)
      qui count if b_`P'`A'' >= bo_`P'`A''
      dis "PANEL`P'_NGE_`A' = " r(N)
    }
    qui sum b_`P'BASE'
    dis "PANEL`P'_SD_BASE = " %14.10f r(sd)
    qui count if b_`P'BASE' >= bo_`P'BASE'
    dis "PANEL`P'_NGE_BASE = " r(N)
    * identity row: which _rep is the observed assignment?
    gen byte isobs = abs(b_`P'BASE' - bo_`P'BASE') < 1e-12
    qui count if isobs
    dis "PANEL`P'_NIDENTROWS = " r(N)
    list _rep d_CTRL d_FEWAVE d_FECOH if isobs, noobs
    save PAIRED_`P'.dta, replace
  }
  dis "PAIRED_`P'_RC=" _rc
}

*----- do BASE and CTRL enumerate the SAME assignment sequence?
cap noi {
  use S_ABASE.dta, clear
  keep _rep onset_*
  foreach v of varlist onset_* {
    rename `v' zz`v'
  }
  merge 1:1 _rep using S_ACTRL.dta, nogen keepusing(_rep onset_*)
  gen byte mm = 0
  foreach v of varlist zzonset_* {
    local raw = substr("`v'",3,.)
    replace mm = 1 if (`v' != `raw') & !(missing(`v') & missing(`raw'))
  }
  qui count if mm
  dis "ASG_MISMATCH_ROWS = " r(N) " of " _N
}
dis "ASG_RC=" _rc
dis "PROBE03 DONE"
