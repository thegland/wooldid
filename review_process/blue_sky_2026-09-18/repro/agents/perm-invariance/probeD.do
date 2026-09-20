*==========================================================================
* probeD.do -- (1) FIDELITY: is the permsave() assignment block actually
*                  the assignment the replication used?  Reconstructed
*                  independently and re-estimated.
*              (2) permreps() prefix property of the draws.
*              (3) invariance of the att p to the estimand list requested.
*              (4) permstrata(constant) == no permstrata.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

mkpanel
save panelD.dta, replace

*---------------------------------------------------------------- D0
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(40) rseed(20260917) esprelength(2) espostlength(3) dots(0) replace permsave(ps_D0)
dis "RC_D0=" _rc
dmp D0
asgdmp D0

*---------------------------------------------------------------- FIDELITY
foreach r in 1 7 23 40 {
  use "ps_D0.dta", clear
  qui keep if _rep == `r'
  local sv_att  = main_att[1]
  local sv_atti = main_att_i[1]
  qui ds onset_*
  local ov `r(varlist)'
  local map
  foreach v of local ov {
     local id = real(subinstr("`v'", "onset_", "", 1))
     local on = `v'[1]
     local map `map' `id'|`on'
  }
  use panelD.dta, clear
  qui gen double nt = .
  foreach pair of local map {
     local id  = substr("`pair'", 1, strpos("`pair'", "|") - 1)
     local on  = substr("`pair'", strpos("`pair'", "|") + 1, .)
     if "`on'" != "." qui replace nt = `on' if state == `id'
  }
  qui count if !missing(nt)
  local ntreat = r(N)
  cap noi qui wooldid y unit t nt , att att_i fastnose cluster(state) esprelength(2) espostlength(3) esrelativeto(-1)
  local frc = _rc
  dis "FID rep=`r' rc=`frc' treatedrows=`ntreat'"
  if `frc' == 0 {
    dis "FID rep=`r' att   saved=" %21x `sv_att'  " recon=" %21x _b[main_att]   " equal=" (`sv_att'  == _b[main_att])
    dis "FID rep=`r' att_i saved=" %21x `sv_atti' " recon=" %21x _b[main_att_i] " equal=" (`sv_atti' == _b[main_att_i])
  }
}

*---------------------------------------------------------------- D1 permreps prefix
use panelD.dta, clear
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(80) rseed(20260917) esprelength(2) espostlength(3) dots(0) replace permsave(ps_D1r80)
dis "RC_D1r80=" _rc
cap dmp D1r80
cap asgdmp D1r80

*---------------------------------------------------------------- D2 att only
use panelD.dta, clear
cap noi wooldid y unit t ttre , att cluster(state) permute(shuffle) permreps(40) rseed(20260917) esprelength(2) espostlength(3) dots(0) replace permsave(ps_D2attonly)
dis "RC_D2attonly=" _rc
cap dmp D2attonly
cap asgdmp D2attonly

*---------------------------------------------------------------- D3 permstrata(constant)
use panelD.dta, clear
qui gen byte one = 1
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permstrata(one) permreps(40) rseed(20260917) esprelength(2) espostlength(3) dots(0) replace permsave(ps_D3str1)
dis "RC_D3str1=" _rc
cap dmp D3str1
cap asgdmp D3str1

*---------------------------------------------------------------- D4 no ES window
use panelD.dta, clear
cap noi wooldid y unit t ttre , att att_i cluster(state) permute(shuffle) permreps(40) rseed(20260917) dots(0) replace permsave(ps_D4noes)
dis "RC_D4noes=" _rc
cap dmp D4noes
cap asgdmp D4noes

dis "PROBED_DONE"
