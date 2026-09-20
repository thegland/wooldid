*======================================================================
* d03 -- diagnostic: how does the bootstrap(bayes) replication distribution
* differ from bootstrap(standard) on the SAME data and the same rseed?
* Dirichlet(1,..,1) cluster reweighting should share the pairs cluster
* bootstrap's first two moments (ado:6498 says so in as many words), so
*   SD(bayes)/SD(standard) ~ 1   and   both centred at the observed estimate.
* Self-contained.  stata-mp -e do d03_bayes_vs_standard.do
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cap prog drop mkdgp
prog def mkdgp
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs `ngrp'
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  local nc = round(`ngrp'*0.2)
  gen double ttre = .
  replace ttre = 2004 if rk <= `nc'
  replace ttre = 2005 if rk >    `nc' & rk <= 2*`nc'
  replace ttre = 2006 if rk >  2*`nc' & rk <= 3*`nc'
  gen double ai = rnormal()
  drop u rk
  expand `nper'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75) if _n>1
  gen double y = ai + ft + e
  keep st t ttre y
end

dis "draw bobs  sd_std  sd_bay  ratio  mean_std-bobs  mean_bay-bobs  p_std  p_bay"
forvalues d = 1/5 {
  local sx = 2000000 + 10*1000 + `d'
  local rs = 700000 + `d'*7 + 10
  mkdgp , ngrp(10) nper(8) seedx(`sx')
  tempfile dat
  qui save "`dat'", replace
  foreach m in standard bayes {
    qui use "`dat'", clear
    ereturn clear
    cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(`m') reps(399) rseed(`rs')
    local rc_`m' = _rc
    if `rc_`m'' == 0 {
      mat BD = e(bootdist)
      local c = colnumb(BD, "main_att")
      if `c' >= . local c = 1
      qui svmat double BD, names(bd)
      qui sum bd`c'
      local sd_`m' = r(sd)
      local mn_`m' = r(mean)
      mat PP = e(boot_p_percentile)
      local p_`m' = PP[1, colnumb(PP,"main_att")]
      mat CC = e(boot_comparator)
      local bo_`m' = CC[1, 1]
      mat PC = e(boot_pctl_mainresults)
      local se_`m' = PC[1, 2]
      local lo_`m' = PC[1, 5]
      local hi_`m' = PC[1, 6]
      qui drop bd*
    }
  }
  dis %3.0f `d' "  " %8.4f `bo_standard' "  " %7.4f `sd_standard' "  " %7.4f `sd_bayes' "  " ///
      %6.3f `=`sd_bayes'/`sd_standard'' "  " %8.4f `=`mn_standard'-`bo_standard'' "  " ///
      %8.4f `=`mn_bayes'-`bo_standard'' "  " %6.4f `p_standard' "  " %6.4f `p_bayes'
  dis "      reported se: std=" %7.4f `se_standard' " bay=" %7.4f `se_bayes' ///
      "   CI std=[" %7.4f `lo_standard' "," %7.4f `hi_standard' "]  bay=[" %7.4f `lo_bayes' "," %7.4f `hi_bayes' "]"
}
