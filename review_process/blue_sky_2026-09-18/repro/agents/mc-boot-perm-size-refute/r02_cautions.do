*======================================================================
* r02_cautions -- SELF-CONTAINED.  What does wooldid actually SAY on the
* exact design the seeded findings used (G=10 states, 6 treated in 3 cohorts
* of 2, 4 never-treated, i = cluster = state, sharp null)?
*   (a) MCBP-GUARD claims the few-cluster caution keys only on TREATED clusters
*       and that a total-cluster caution does not exist.  G=9 / G=10 / G=12.
*   (b) MCBP-JOINTSIZE claims a user "reads a statistic that fires three times
*       in four" with nothing said.  Run jointtests NON-quietly.
* Run:  stata-mp -e do r02_cautions.do        (~1 min)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpan
prog def mkpan
  syntax , G(integer) T(integer) SEEDX(integer) [KPER(integer 2) RHO(real 0.4) NTRTC(integer 3)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `g'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  gen double ttre = .
  forvalues c = 1/`ntrtc' {
    local lo = (`c'-1)*`kper' + 1
    local hi = `c'*`kper'
    local on = 2003 + 2*`c'
    replace ttre = `on' if rk >= `lo' & rk <= `hi'
  }
  gen double ai = rnormal()*1.2
  drop u rk
  expand `t'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.6*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal() if _n==1
  by st: replace e = `rho'*e[_n-1] + rnormal()*sqrt(1-`rho'^2) if _n>1
  gen double y = ai + ft + e
  keep st t ttre y
end

dis "##### (b) JOINTTESTS AT cluster(st)==i, NOT QUIET #####"
mkpan , g(10) t(10) seedx(31001) kper(2) ntrtc(3)
qui tab ttre, missing
dis "n treated states = " r(N) - r(r)
cap noi wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial
dis "JT_RC=" _rc
dis "JT_pjoint=" e(joint_horiz) "  JT_padv=" e(adversarial_pval_att)

dis "##### (b2) the SAME with permute(), which the caution recommends #####"
mkpan , g(10) t(10) seedx(31001) kper(2) ntrtc(3)
cap noi wooldid y st t ttre, cluster(st) att esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial permute(shuffle) permreps(79) rseed(5)
dis "JTP_RC=" _rc
dis "JTP_pjoint=" e(joint_horiz)

dis "##### (a) FEW-CLUSTER CAUTION CENSUS, bootstrap(standard) #####"
foreach g in 8 9 10 11 12 {
  dis "----- G=`g' -----"
  mkpan , g(`g') t(10) seedx(`=32000+`g'') kper(2) ntrtc(3)
  qui count if ttre < .
  qui gegen byte tg = tag(st) if ttre < .
  qui count if tg==1
  dis "TREATEDSTATES_G`g'=" r(N)
  drop tg
  cap noi wooldid y st t ttre, cluster(st) att bootstrap(standard) reps(99) rseed(7)
  dis "BS_RC_G`g'=" _rc
}
dis "##### (a2) G=10 with only 4 treated states (2 cohorts of 2) #####"
mkpan , g(10) t(10) seedx(33001) kper(2) ntrtc(2)
cap noi wooldid y st t ttre, cluster(st) att bootstrap(standard) reps(99) rseed(7)
dis "BS_RC_G10T4=" _rc
dis "R02_DONE"
