*======================================================================
* mc_main.do -- Monte Carlo size study of wooldid's replication engines
*               under a TRUE SHARP NULL (no treatment effect anywhere).
*
*  usage:  stata-mp -e do mc_main.do <G> <d0> <d1> "<arms>" <hetv*100> <tag>
*  e.g.    stata-mp -e do mc_main.do 10 1 50 "A0 A1 A5 A7" 100 g10core
*
*  arms:  A0 plain analytic (+ES, jointtests, adversarial)
*         A1 permute(shuffle) permreps(79)                 [default equal-tailed add-one p]
*         A2 permute(shuffle) permreps(79) permabsp
*         A3 permute(shuffle) permreps(79) permrawp
*         A4 permute(shuffle) permreps(79) permtstat
*         A5 bootstrap(standard) reps(79)  att att_i       [pctl + bc + wald families]
*         A6 bootstrap(standard) reps(79) boottstat        [studentized]
*         A7 bootstrap(bayes)    reps(79)  att only        [pctl + bc + wald]
*         A8 bootstrap(standard) reps(79)  att only        [non-reid comparator for A7]
*
*  output: appended rows to mcres_<tag>.csv, ONE ROW PER (draw, arm, family, estimand)
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
args GG D0 D1 ARMS HV TAG
local hetv = `HV'/100
dis "ARGS G=`GG' d0=`D0' d1=`D1' arms=`ARMS' hetv=`hetv' tag=`TAG'"
which wooldid

*---------------- DGP: SHARP NULL state-year panel (no treatment effect) -----
* ngrp   states (= i = cluster);  nper periods 2000..2000+nper-1
* onsets 2004/2005/2006 to 20% of states each; 40% never treated
* y = a_i + f_t + e_it ,  e_it AR(1) rho=.5 within state (Bertrand-Duflo-
* Mullainathan serial correlation).  hetv scales the error sd of TREATED
* states only (hetv=1 -> homogeneous).
cap prog drop mkdgp
prog def mkdgp
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1) HETV(real 1) NPERCELL(integer 1)]
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
  gen double esd = cond(ttre < ., `hetv', 1)
  drop u rk
  expand `nper'
  bysort st: gen int t = 1999 + _n
  sort t st
  by t: gen double ft = 0.5*rnormal() if _n==1
  by t: replace ft = ft[1]
  sort st t
  by st: gen double e = rnormal()*esd if _n==1
  by st: replace e = 0.5*e[_n-1] + rnormal()*sqrt(1-0.25)*esd if _n>1
  gen double y = ai + ft + e
  if `npercell' > 1 {
    expand `npercell'
    bysort st t: gen int sub = _n
    gen long unit = st*1000 + sub
    replace y = y + rnormal()*0.5
  }
  keep st t ttre y
  order st t ttre y
end

*---------------- open CSV (append; header only if new) ----------------
local csv "mcres_`TAG'.csv"
cap confirm file "`csv'"
local isnew = _rc
file open fh using "`csv'", write append text
if `isnew' != 0 file write fh "draw,G,hetv,arm,fam,est,b,se,p,lo,hi,rc,ncomp,nfail,nblank,extra" _n

forvalues d = `D0'/`D1' {
  local sx = 2000000 + `HV'*100000 + `GG'*1000 + `d'
  local rs = 700000 + `d'*7 + `GG'
  mkdgp , ngrp(`GG') nper(8) seedx(`sx') hetv(`hetv')
  tempfile dat
  qui save "`dat'", replace

  foreach a of local ARMS {
    qui use "`dat'", clear
    ereturn clear
    cap mat drop T1
    cap mat drop T2
    cap mat drop T3
    local rc = 0
    local nc = .
    local nf = .
    local nb = .
    local ex = .
    local fams ""
    local kind ""

    if "`a'" == "A0" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i
      local rc = _rc
      if `rc' == 0 {
        mat T1 = e(wdidmainresults)
        local fams "an"
        local kind "mat"
      }
    }
    if "`a'" == "A9" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial
      local rc = _rc
      if `rc' == 0 {
        mat T1 = e(wdidmainresults)
        local fams "anes"
        local kind "mat"
        local ex = e(joint_horiz)
      }
    }
    if "`a'" == "A10" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial permute(shuffle) permreps(79) rseed(`rs')
      local rc = _rc
      local kind "perm"
      local fams "permes"
    }
    if "`a'" == "A11" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial bootstrap(standard) reps(79) rseed(`rs')
      local rc = _rc
      local kind "boot"
    }
    if "`a'" == "A12" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i esprelength(3) espostlength(3) esfixedbaseperiod jointtests adversarial bootstrap(standard) reps(79) rseed(`rs') boottstat
      local rc = _rc
      local kind "boott"
    }
    if "`a'" == "A1" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs')
      local rc = _rc
      local kind "perm"
      local fams "perm"
    }
    if "`a'" == "A2" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs') permabsp
      local rc = _rc
      local kind "perm"
      local fams "permabs"
    }
    if "`a'" == "A3" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs') permrawp
      local rc = _rc
      local kind "perm"
      local fams "permraw"
    }
    if "`a'" == "A4" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i permute(shuffle) permreps(79) rseed(`rs') permtstat
      local rc = _rc
      local kind "perm"
      local fams "permt"
    }
    if "`a'" == "A5" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs')
      local rc = _rc
      local kind "boot"
    }
    if "`a'" == "A6" {
      cap noi qui wooldid y st t ttre, cluster(st) att att_i bootstrap(standard) reps(79) rseed(`rs') boottstat
      local rc = _rc
      local kind "boott"
    }
    if "`a'" == "A7" {
      cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(bayes) reps(79) rseed(`rs')
      local rc = _rc
      local kind "boot"
    }
    if "`a'" == "A8" {
      cap noi qui wooldid y st t ttre, cluster(st) att bootstrap(standard) reps(79) rseed(`rs')
      local rc = _rc
      local kind "boot"
    }

    *---- harvest ----
    if "`kind'" == "perm" {
      local nc = e(perm_completed)
      local nf = e(perm_failed)
      local nb = e(perm_blank)
      local ex = e(perm_npossible)
      * CI MUST come from the TABLE (e(wdidmainresults) cols 5/6).  The scalars
      * e(perm_lb_<est>)/e(perm_ub_<est>) are the deprecated NULL BAND, NOT an
      * interval for the estimand (sthlp:1353) - verified on d04_perm_dual.do.
      cap mat drop TP
      if `rc' == 0 cap mat TP = e(wdidmainresults)
      foreach ee in att att_i {
        local b  = .
        local se = .
        local p  = .
        local lo = .
        local hi = .
        if `rc' == 0 {
          local p  = e(perm_p_`ee')
          local se = e(perm_sd_`ee')
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
              local b  = TP[`j',1]
              local lo = TP[`j',5]
              local hi = TP[`j',6]
            }
          }
        }
        file write fh "`d',`GG',`hetv',`a',`fams',`ee',`b',`se',`p',`lo',`hi',`rc',`nc',`nf',`nb',`ex'" _n
      }
    }
    else if "`kind'" == "boot" | "`kind'" == "boott" {
      local nc = e(boot_completed)
      local nf = e(boot_failed)
      local nb = e(boot_N_V)
      if `rc' == 0 {
        if "`kind'" == "boot" {
          cap mat T1 = e(boot_pctl_mainresults)
          cap mat T2 = e(boot_bc_mainresults)
          cap mat T3 = e(boot_wald_mainresults)
          local fams "pctl bc wald"
        }
        else {
          cap mat T1 = e(boot_tpctl_mainresults)
          local fams "tpctl"
        }
        local kind "mat"
      }
      else local fams ""
    }

    if "`kind'" == "mat" {
      local q = 0
      foreach f of local fams {
        local ++q
        cap confirm matrix T`q'
        if _rc {
          foreach ee in att att_i {
            file write fh "`d',`GG',`hetv',`a',`f',`ee',.,.,.,.,.,`rc',`nc',`nf',`nb',`ex'" _n
          }
          continue
        }
        local rn : rownames T`q'
        local re : roweq T`q'
        local nr : word count `rn'
        foreach ee in att att_i {
          local b  = .
          local se = .
          local p  = .
          local lo = .
          local hi = .
          local j = 0
          forvalues z = 1/`nr' {
            local nz : word `z' of `rn'
            local ez : word `z' of `re'
            if `j' == 0 & "`nz'" == "`ee'" & "`ez'" == "main" local j = `z'
          }
          if `j' > 0 {
            local b  = T`q'[`j',1]
            local se = T`q'[`j',2]
            local p  = T`q'[`j',4]
            local lo = T`q'[`j',5]
            local hi = T`q'[`j',6]
          }
          file write fh "`d',`GG',`hetv',`a',`f',`ee',`b',`se',`p',`lo',`hi',`rc',`nc',`nf',`nb',`ex'" _n
        }
      }
    }
    else if `rc' != 0 {
      file write fh "`d',`GG',`hetv',`a',NA,att,.,.,.,.,.,`rc',`nc',`nf',`nb',`ex'" _n
    }
    if inlist("`a'","A9","A10","A11","A12") {
      local jp = .
      local ap = .
      if `rc' == 0 {
        local jp = e(joint_horiz)
        local ap = e(adversarial_pval_att)
      }
      file write fh "`d',`GG',`hetv',`a',joint,horiz,.,.,`jp',.,.,`rc',`nc',`nf',`nb',`ex'" _n
      file write fh "`d',`GG',`hetv',`a',adv,att,.,.,`ap',.,.,`rc',`nc',`nf',`nb',`ex'" _n
    }
  }
  if mod(`d',10) == 0 dis "PROGRESS draw `d' done"
}
file close fh
dis "BLOCK_DONE tag=`TAG' G=`GG' draws `D0'-`D1' arms=`ARMS'"
