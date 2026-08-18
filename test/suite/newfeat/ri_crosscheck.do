* ri_crosscheck.do — external cross-check of F6 randomization inference (exact mode).
* On tiny data (6 clusters, 3 treated at onset 4, 3 never; C(6,3)=20 assignments) run
* wooldid permute(shuffle) permexact, then EXTERNALLY enumerate the same assignment set
* (clusters across ttre levels, level counts held fixed — wooldid's scheme) with plain
* fastnose wooldid calls (the engine estimates draws on the fastnose path), and check:
*   (a) multiset of external coefficients == e(permdist)   [assert 1e-10; bitwise reported]
*   (b) exact p recomputed from the external distribution == e(perm_p_att)
*       (exact mode uses the plain share nex/nuse — the observed assignment is a member
*        of the enumeration; the add-one form is sampled-mode only)
*   (c) count == e(perm_npossible) == e(perm_completed) == 20
* Args: 1 = path to wooldid.ado (file or its directory). No data dir needed (synthetic).
* Run:  "<stata>" -e do ri_crosscheck.do <adopath>     -> expect RI_CROSSCHECK_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 24601
set sortseed 24601

local ADO `1'
if `"`ADO'"' == "" local ADO /tmp/wdsuite_build/newfeat/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'

global pf 0

* ---- tiny data: 6 clusters (i = cluster), t = 1..8, onset 4 for i<=3, never for i>3 ----
clear
qui set obs 6
qui gen i = _n
qui gen ttre = cond(_n <= 3, 4, .)
qui expand 8
qui bysort i: gen t = _n
qui gen y = 3*(t >= ttre & ttre < .) + 0.1*i + rnormal()
sort i t

* ---- the engine run: exact enumeration ----
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
if _rc != 0 {
  dis as error "R0 FAIL: engine run rc=" _rc
  global pf = $pf + 1
  dis as error "RI_CROSSCHECK_FAILS=$pf"
  exit, clear
}

* capture engine results BEFORE any further estimation clobbers e()
tempname PD pengine npos ncomp
matrix `PD' = e(permdist)
scalar `pengine' = e(perm_p_att)
scalar `npos'   = e(perm_npossible)
scalar `ncomp'  = e(perm_completed)
local pmethod `e(perm_method)'
local pstat `e(perm_stat)'

* ---- R1: counts and surface ----
if `npos' != 20 | `ncomp' != 20 | rowsof(`PD') != 20 | colsof(`PD') != 1 {
  dis as error "R1 FAIL: npossible=" `npos' " completed=" `ncomp' " permdist is " rowsof(`PD') "x" colsof(`PD') " (want 20x1)"
  global pf = $pf + 1
}
if "`pmethod'" != "shuffle" | "`pstat'" != "b" {
  dis as error "R1b FAIL: perm_method=[`pmethod'] perm_stat=[`pstat']"
  global pf = $pf + 1
}

* ---- rep-0 observed statistic: fastnose on the TRUE assignment (the engine's own
*      comparison statistic; scalar storage, never a local) ----
sort i t
cap noi wooldid y i t ttre, att cluster(i) fastnose
if _rc != 0 {
  dis as error "R2 FAIL: observed fastnose run rc=" _rc
  global pf = $pf + 1
}
tempname OBS
scalar `OBS' = _b[main_att]

* ---- external enumeration: every way to place the 3 onset-4 labels on 6 clusters ----
tempname EXT
matrix `EXT' = J(20, 1, .)
local k 0
local extfail 0
forval a = 1/4 {
  forval b = `=`a'+1'/5 {
    forval c = `=`b'+1'/6 {
      local ++k
      cap drop ttre2
      qui gen double ttre2 = cond(inlist(i, `a', `b', `c'), 4, .)
      sort i t
      cap qui wooldid y i t ttre2, att cluster(i) fastnose
      if _rc != 0 {
        dis as error "R3 FAIL: external assignment `k' (treated `a' `b' `c') rc=" _rc
        local ++extfail
        global pf = $pf + 1
      }
      else matrix `EXT'[`k', 1] = _b[main_att]
    }
  }
}
if `k' != 20 {
  dis as error "R3b FAIL: enumerated `k' assignments (want 20)"
  global pf = $pf + 1
}

* ---- (a) multiset comparison: sorted external coefs vs sorted e(permdist) ----
if `extfail' == 0 {
  mata: pd = st_matrix(st_local("PD")); ex = st_matrix(st_local("EXT"))
  mata: pds = sort(pd, 1); exs = sort(ex, 1)
  mata: st_numscalar("__maxrd", max(reldif(pds, exs)))
  mata: st_numscalar("__nbit", sum(pds :== exs))
  if __maxrd > 1e-10 {
    dis as error "RA FAIL: sorted external vs e(permdist) max reldif = " %12.4e __maxrd " (> 1e-10)"
    global pf = $pf + 1
  }
  dis "RA info: max reldif (sorted) = " %12.4e __maxrd "; bitwise-equal entries = " __nbit " of 20"

  * membership: the observed rep-0 statistic is a member of the external distribution
  mata: st_numscalar("__member", anyof(st_matrix(st_local("EXT")), st_numscalar(st_local("OBS"))))
  mata: st_numscalar("__minod", min(abs(st_matrix(st_local("EXT")) :- st_numscalar(st_local("OBS")))))
  if __member != 1 {
    dis as error "RM FAIL: observed fastnose stat is not a bitwise member of the external set; min |diff| = " %12.4e __minod
    dis as error "         obs = " %21x `OBS'
    global pf = $pf + 1
  }

  * ---- (b) exact p recomputed from the external distribution (plain share) ----
  mata: st_numscalar("__nex", sum(abs(st_matrix(st_local("EXT"))) :>= abs(st_numscalar(st_local("OBS")))))
  tempname PEXT
  scalar `PEXT' = __nex/20
  if reldif(`PEXT', `pengine') > 1e-12 {
    dis as error "RB FAIL: external exact p = " `PEXT' " (nex=" __nex "/20) vs e(perm_p_att) = " `pengine'
    dis as error "         obs = " %21x `OBS'
    * ulp-boundary triage: print the external draw closest to |obs|
    mata: dd = abs(abs(st_matrix(st_local("EXT"))) :- abs(st_numscalar(st_local("OBS")))); w = selectindex(dd :== min(dd)); st_numscalar("__near", st_matrix(st_local("EXT"))[w[1], 1])
    dis as error "         nearest external draw = " %21x __near
    global pf = $pf + 1
  }
  dis "RB info: exact p = " `PEXT' " (engine " `pengine' ")"
}
else dis as error "RA/RB skipped: `extfail' external run(s) failed"

* ---- (c) count identity already asserted in R1; re-assert against enumeration ----
if `npos' != `k' {
  dis as error "RC FAIL: e(perm_npossible)=" `npos' " != external count `k'"
  global pf = $pf + 1
}

if $pf == 0 dis "RI_CROSSCHECK_ALL_PASS"
else dis as error "RI_CROSSCHECK_FAILS=$pf"
exit, clear
