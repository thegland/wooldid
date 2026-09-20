*=======================================================================
* repro_perm_exactness.do -- SELF-CONTAINED.  Two THEOREM-BACKED size checks
* of the permutation engine under a TRUE SHARP NULL.
*
*  EX : 8 states, cohort sizes (2 onset 2004, 1 onset 2006, 5 never).
*       M = 8!/(2! 1! 5!) = 168 distinct onset assignments, so `permexact'
*       ENUMERATES the whole reference set.  The observed assignment is drawn
*       uniformly from those 168 and y does not depend on it, so the rank of
*       the observed statistic is uniform on 1..168 and
*           P(default equal-tailed p <= .05) = 2*4/168 = .0476  EXACTLY.
*  PO : 12 states in 4 equal onset groups; i = COHORT and cluster(st) is FINER
*       than i, which only `permonsetcohorts' (shipped 2026-09-16) admits.
*       Sampled permreps(79); the add-one convention is exact too, so
*           P(p <= .05) <= .05.
*  POC: the same data run the ordinary way (i = st = cluster) as a control.
* A rejection rate materially above nominal in ANY of these is a DEFECT.
*   stata-mp -e do repro_perm_exactness.do      (~45 min at NDRAW 40)
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local NDRAW 40

cap prog drop mkEX
prog def mkEX
  syntax [, SEEDX(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs 8
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  gen double ttre = .
  replace ttre = 2004 if rk <= 2
  replace ttre = 2006 if rk == 3
  gen double ai = rnormal()
  drop u rk
  expand 8
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
cap prog drop mkPO
prog def mkPO
  syntax [, SEEDX(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs 12
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  gen double ttre = .
  gen long coh = 4
  replace ttre = 2004 if rk <= 3
  replace coh  = 1    if rk <= 3
  replace ttre = 2005 if rk > 3 & rk <= 6
  replace coh  = 2    if rk > 3 & rk <= 6
  replace ttre = 2006 if rk > 6 & rk <= 9
  replace coh  = 3    if rk > 6 & rk <= 9
  gen double ai = rnormal()
  drop u rk
  expand 8
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75) if _n>1
  gen double y = ai + ft + e
  keep st coh t ttre y
end

tempname P
postfile `P' int draw str8 arm double p double lo double hi int rc double ncomp double npos ///
  using repro_perm_exact_out.dta, replace
forvalues d = 1/`NDRAW' {
  local sx = 3000000 + `d'
  local rs = 800000 + `d'*13
  foreach arm in EX PO POC {
    if "`arm'" == "EX" mkEX , seedx(`sx')
    else               mkPO , seedx(`sx')
    ereturn clear
    if "`arm'" == "EX"  cap noi qui wooldid y st  t ttre, cluster(st) att permute(shuffle) permexact rseed(`rs')
    if "`arm'" == "PO"  cap noi qui wooldid y coh t ttre, cluster(st) att permute(shuffle) permreps(79) rseed(`rs') permonsetcohorts
    if "`arm'" == "POC" cap noi qui wooldid y st  t ttre, cluster(st) att permute(shuffle) permreps(79) rseed(`rs')
    local rc = _rc
    local p = .
    local lo = .
    local hi = .
    if `rc' == 0 {
      local p = e(perm_p_att)
      cap mat drop TT
      cap mat TT = e(wdidmainresults)
      cap confirm matrix TT
      if !_rc {
        local rn : rownames TT
        local re : roweq TT
        local nr : word count `rn'
        local j = 0
        forvalues z = 1/`nr' {
          local nz : word `z' of `rn'
          local ez : word `z' of `re'
          if `j'==0 & "`nz'"=="att" & "`ez'"=="main" local j = `z'
        }
        if `j' > 0 {
          local lo = TT[`j',5]
          local hi = TT[`j',6]
        }
      }
    }
    post `P' (`d') ("`arm'") (`p') (`lo') (`hi') (`rc') (`=e(perm_completed)') (`=e(perm_npossible)')
  }
  if mod(`d',10)==0 dis "... draw `d' done"
}
postclose `P'
use repro_perm_exact_out.dta, clear
gen byte rej  = (p  <= .05 + 1e-9) if p < .
gen byte excl = (lo > 0 | hi < 0)  if lo < . & hi < .
dis "=== rejection of a TRUE sharp null; EX nominal .0476, PO/POC nominal .0500 ==="
table arm, stat(mean rej) stat(mean excl) stat(count p) stat(min ncomp) stat(mean npos)
