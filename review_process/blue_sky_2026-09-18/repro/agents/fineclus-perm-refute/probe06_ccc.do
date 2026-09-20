* probe06_ccc.do -- diagnose the coarsecohortcontrols / interactivecontrols arms, and test
* whether the ONE diagnostic wooldid already owns for "this control has no within-cell
* variation" (the interactivecontrols() caution, sthlp:1068) survives permute()'s
* fast-mode default.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute/../../build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
which wooldid
use panA.dta, clear
* a_s is a STATE-level variable (varies within cohort) -- the contrast arm
desc

dis "=============== (1) OBSERVED FIT ONLY, no permute ==============="
foreach A in NONE CCCCOH CCCSTATE IC ICSTATE {
  use panA.dta, clear
  local o
  if "`A'"=="CCCCOH"   local o coarsecohortcontrols(c.cohbase)
  if "`A'"=="CCCSTATE" local o coarsecohortcontrols(c.a_s)
  if "`A'"=="IC"       local o interactivecontrols(cohbase)
  if "`A'"=="ICSTATE"  local o interactivecontrols(a_s)
  dis "---- OBS arm `A' ----"
  cap noi wooldid y coh t onset, att cluster(state) `o'
  dis "OBS_`A'_RC=" _rc
  if _rc==0 {
    mat B = e(b)
    dis "OBS_`A'_B %21x = " %21x B[1,1]
  }
}

dis "=============== (2) SAME ARMS UNDER permute() -- is the caution still printed? ==============="
foreach A in IC ICSTATE {
  use panA.dta, clear
  local o
  if "`A'"=="IC"      local o interactivecontrols(cohbase)
  if "`A'"=="ICSTATE" local o interactivecontrols(a_s)
  dis "---- PERM-CAUTION arm `A' ----"
  cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts `o'
  dis "PC_`A'_RC=" _rc
}
dis "---- PERM-CAUTION arm IC with permanalyticse (non-fast path) ----"
use panA.dta, clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(20) permonsetcohorts permanalyticse interactivecontrols(cohbase)
dis "PCASE_RC=" _rc

dis "=============== (3) THE CCC BLOW-UP: is it the REBUILD? ==============="
cap prog drop arm
prog def arm
  syntax , TAG(string) CLUS(string) [OPTS(string)]
  use panA.dta, clear
  cap noi qui wooldid y coh t onset, att cluster(`clus') permute(shuffle) permexact permonsetcohorts `opts'
  local rc = _rc
  dis "ARM_`tag'_RC=`rc'"
  if `rc'==0 {
    mat B = e(b)
    mat D = e(permdist_b)
    dis "ARM_`tag'_B %21x = " %21x B[1,1]
    dis "ARM_`tag'_P = " %14.10f e(perm_p_att)
    dis "ARM_`tag'_SD = " %14.10f e(perm_sd_att)
    dis "ARM_`tag'_NFAIL = " e(perm_failed)
    dis "ARM_`tag'_REBUILT = " e(perm_irebuilt)
    preserve
    clear
    svmat double D, names(dd)
    qui sum dd1, det
    dis "ARM_`tag'_DIST p1=" %10.4f r(p1) " p25=" %10.4f r(p25) " p50=" %10.4f r(p50) " p75=" %10.4f r(p75) " p99=" %10.4f r(p99) " min=" %10.4f r(min) " max=" %10.4f r(max)
    restore
  }
end
* rebuild fires
arm , tag(R_BASE)     clus(state)
arm , tag(R_CCCCOH)   clus(state) opts(coarsecohortcontrols(c.cohbase))
arm , tag(R_CCCSTATE) clus(state) opts(coarsecohortcontrols(c.a_s))
arm , tag(R_IC)       clus(state) opts(interactivecontrols(cohbase))
* no rebuild (cluster == i)
arm , tag(N_BASE)     clus(coh)
arm , tag(N_CCCCOH)   clus(coh) opts(coarsecohortcontrols(c.cohbase))
arm , tag(N_IC)       clus(coh) opts(interactivecontrols(cohbase))
dis "PROBE06 DONE"
