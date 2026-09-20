*======================================================================
* mc_alt.do -- theorem-backed exactness probes of the permutation engine
*              under a TRUE SHARP NULL.
*
*  usage: stata-mp -e do mc_alt.do <design> <d0> <d1> <tag>
*
*  designs
*   EX  : G=8 states (2 onset 2004, 1 onset 2006, 5 never), i = st = cluster.
*         M = 8!/(2! 1! 5!) = 168 distinct assignments -> permexact ENUMERATES
*         them.  Under a sharp null with a uniformly-drawn assignment the rank
*         of the observed statistic among the M is uniform, so
*             P(default p <= .05) = 2*4/168 = .0476  EXACTLY,
*         with NO asymptotics, whatever the error structure.  Any materially
*         larger rejection rate is a defect, not a design limit.
*   PO  : 12 states, i = COHORT (3 treated onsets + never), cluster(st) FINER
*         than i -> requires permonsetcohorts (T1, shipped 2026-09-16).
*         Sampled permreps(79); the add-one convention is exact too
*         (Phipson-Smyth 2010), so P(p <= .05) <= .05 is again a theorem.
*   POC : the SAME data run the ordinary way, i = st = cluster(st).  Control.
*   HV  : G=10 state-year, treated states' error sd x3 (MacKinnon-Webb's
*         atypical-variance case the ado cautions about at ado:6056);
*         permute default + permtstat + bootstrap(standard) + bootstrap(bayes)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
args DES D0 D1 TAG
dis "ARGS des=`DES' d0=`D0' d1=`D1' tag=`TAG'"
which wooldid

*---------------- generators ----------------
* EX design: 8 states, cohort sizes (2,1), 5 never; one row per state-year
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

* PO design: 12 states in 4 equal groups (3 onsets + never); i = cohort code
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

* HV design: 10 states, treated states' error sd x3
cap prog drop mkHV
prog def mkHV
  syntax [, SEEDX(integer 1) FIXEDVAR(integer 0)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs 10
  gen long st = _n
  gen double u = runiform()
  sort u
  gen long rk = _n
  gen double ttre = .
  replace ttre = 2004 if rk <= 2
  replace ttre = 2005 if rk > 2 & rk <= 4
  replace ttre = 2006 if rk > 4 & rk <= 6
  gen double ai = rnormal()
  * fixedvar(0): the error sd is inflated ON THE TREATED states -> the SHARP null
  *              is false (treatment changes the variance), so a rejection is power
  *              against a variance effect AND the superpopulation size distortion an
  *              applied user would call a false positive.  This is the F049 setup.
  * fixedvar(1): the error sd is a FIXED state attribute (odd/even original id),
  *              independent of the onset draw -> the sharp null HOLDS and
  *              randomization inference must still be EXACT.  Any over-rejection
  *              here is a defect, not a modelling artefact.
  if `fixedvar' == 0 gen double esd = cond(ttre < ., 3, 1)
  if `fixedvar' == 1 gen double esd = cond(mod(st,2)==0, 3, 1)
  drop u rk
  expand 8
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()*esd
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(.75)*esd if _n>1
  gen double y = ai + ft + e
  keep st t ttre y
end

local csv "mcalt_`TAG'.csv"
cap confirm file "`csv'"
local isnew = _rc
file open fh using "`csv'", write append text
if `isnew' != 0 file write fh "draw,des,arm,est,p,lo,hi,rc,ncomp,nfail,nblank,npos,exact" _n

forvalues d = `D0'/`D1' {
  local sx = 3000000 + `d'
  local rs = 800000 + `d'*13

  if "`DES'" == "EX" {
    mkEX , seedx(`sx')
    tempfile dat
    qui save "`dat'", replace
    foreach a in EXdef EXabs {
      qui use "`dat'", clear
      ereturn clear
      local extra ""
      if "`a'" == "EXabs" local extra permabsp
      * permexact and permreps() are mutually exclusive (ado:5937)
      cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permexact rseed(`rs') `extra'
      local rc = _rc
      foreach ee in att att_i {
        local p = .
        local lo = .
        local hi = .
        if `rc' == 0 {
          local p  = e(perm_p_`ee')
          cap mat drop TP
          cap mat TP = e(wdidmainresults)
          cap confirm matrix TP
          if !_rc {
            local rn : rownames TP
            local re : roweq TP
            local nr : word count `rn'
            local j = 0
            forvalues z = 1/`nr' {
              local nz : word `z' of `rn'
              local ez : word `z' of `re'
              if `j'==0 & "`nz'"=="`ee'" & "`ez'"=="main" local j = `z'
            }
            if `j' > 0 {
              local lo = TP[`j',5]
              local hi = TP[`j',6]
            }
          }
        }
        file write fh "`d',EX,`a',`ee',`p',`lo',`hi',`rc',`=e(perm_completed)',`=e(perm_failed)',`=e(perm_blank)',`=e(perm_npossible)',`=e(perm_exact)'" _n
      }
    }
  }

  if "`DES'" == "PO" {
    mkPO , seedx(`sx')
    tempfile dat
    qui save "`dat'", replace
    * feature arm: i = cohort, cluster = state (finer than i)
    qui use "`dat'", clear
    ereturn clear
    cap noi qui wooldid y coh t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs') permonsetcohorts
    local rc = _rc
    foreach ee in att att_i {
      local p = .
      local lo = .
      local hi = .
      if `rc' == 0 {
        local p  = e(perm_p_`ee')
        cap mat drop TP
        cap mat TP = e(wdidmainresults)
        cap confirm matrix TP
        if !_rc {
          local rn : rownames TP
          local re : roweq TP
          local nr : word count `rn'
          local j = 0
          forvalues z = 1/`nr' {
            local nz : word `z' of `rn'
            local ez : word `z' of `re'
            if `j'==0 & "`nz'"=="`ee'" & "`ez'"=="main" local j = `z'
          }
          if `j' > 0 {
            local lo = TP[`j',5]
            local hi = TP[`j',6]
          }
        }
      }
      file write fh "`d',PO,POonc,`ee',`p',`lo',`hi',`rc',`=e(perm_completed)',`=e(perm_failed)',`=e(perm_blank)',`=e(perm_npossible)',`=e(perm_irebuilt)'" _n
    }
    * control arm: i = state, cluster = state (ordinary route, same data)
    qui use "`dat'", clear
    ereturn clear
    cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs')
    local rc = _rc
    foreach ee in att att_i {
      local p = .
      local lo = .
      local hi = .
      if `rc' == 0 {
        local p  = e(perm_p_`ee')
        cap mat drop TP
        cap mat TP = e(wdidmainresults)
        cap confirm matrix TP
        if !_rc {
          local rn : rownames TP
          local re : roweq TP
          local nr : word count `rn'
          local j = 0
          forvalues z = 1/`nr' {
            local nz : word `z' of `rn'
            local ez : word `z' of `re'
            if `j'==0 & "`nz'"=="`ee'" & "`ez'"=="main" local j = `z'
          }
          if `j' > 0 {
            local lo = TP[`j',5]
            local hi = TP[`j',6]
          }
        }
      }
      file write fh "`d',PO,POctl,`ee',`p',`lo',`hi',`rc',`=e(perm_completed)',`=e(perm_failed)',`=e(perm_blank)',`=e(perm_npossible)',." _n
    }
  }

  if "`DES'" == "HV" | "`DES'" == "HVF" {
    mkHV , seedx(`sx') fixedvar(`=("`DES'" == "HVF")')
    tempfile dat
    qui save "`dat'", replace
    foreach a in HVperm HVpermt HVboot HVbayes HVan {
      qui use "`dat'", clear
      ereturn clear
      if "`a'" == "HVperm"  cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs')
      if "`a'" == "HVpermt" cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs') permtstat
      if "`a'" == "HVboot"  cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs')
      if "`a'" == "HVbayes" cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(bayes) reps(79) rseed(`rs')
      if "`a'" == "HVan"    cap noi qui wooldid y st t ttre, cluster(st) att att_i
      local rc = _rc
      local nco = .
      if inlist("`a'","HVperm","HVpermt") local nco = e(perm_completed)
      if inlist("`a'","HVboot","HVbayes") local nco = e(boot_completed)
      foreach ee in att att_i {
        local p = .
        local lo = .
        local hi = .
        if `rc' == 0 {
          if inlist("`a'","HVperm","HVpermt") {
            local p  = e(perm_p_`ee')
            cap mat drop TP
            cap mat TP = e(wdidmainresults)
            cap confirm matrix TP
            if !_rc {
              local rn : rownames TP
              local re : roweq TP
              local nr : word count `rn'
              local j = 0
              forvalues z = 1/`nr' {
                local nz : word `z' of `rn'
                local ez : word `z' of `re'
                if `j'==0 & "`nz'"=="`ee'" & "`ez'"=="main" local j = `z'
              }
              if `j' > 0 {
                local lo = TP[`j',5]
                local hi = TP[`j',6]
              }
            }
          }
          else {
            cap mat drop TT
            if "`a'" == "HVan"  cap mat TT = e(wdidmainresults)
            else                cap mat TT = e(boot_pctl_mainresults)
            cap confirm matrix TT
            if !_rc {
              local rn : rownames TT
              local re : roweq TT
              local nr : word count `rn'
              local j = 0
              forvalues z = 1/`nr' {
                local nz : word `z' of `rn'
                local ez : word `z' of `re'
                if `j'==0 & "`nz'"=="`ee'" & "`ez'"=="main" local j = `z'
              }
              if `j' > 0 {
                local p  = TT[`j',4]
                local lo = TT[`j',5]
                local hi = TT[`j',6]
              }
            }
          }
        }
        file write fh "`d',`DES',`a',`ee',`p',`lo',`hi',`rc',`nco',.,.,.,." _n
      }
    }
  }
  if mod(`d',10) == 0 dis "PROGRESS `DES' draw `d'"
}
file close fh
dis "ALT_BLOCK_DONE des=`DES' draws `D0'-`D1'"
