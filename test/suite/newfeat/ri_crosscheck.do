* ri_crosscheck.do — external cross-check of F6 randomization inference.
* Block R (exact mode, original): on tiny data (6 clusters, 3 treated at onset 4, 3
* never; C(6,3)=20 assignments) run wooldid permute(shuffle) permexact, then EXTERNALLY
* enumerate the same assignment set (clusters across ttre levels, level counts held
* fixed — wooldid's scheme) with plain fastnose wooldid calls (the engine estimates
* draws on the fastnose path), and check:
*   (a) multiset of external coefficients == e(permdist)   [assert 1e-10; bitwise reported]
*   (b) exact p recomputed from the external distribution == e(perm_p_att)
*       (exact mode uses the plain share nex/nuse — the observed assignment is a member
*        of the enumeration; the add-one form is sampled-mode only)
*   (c) count == e(perm_npossible) == e(perm_completed) == 20
* Block S (sampled mode, review F3): permreps(200) on the same panel; membership,
*   binomial frequency bands, counts, the add-one p, and a rep-by-rep comparison with
*   an externally replicated draw stream (the only check with power against a single
*   reused/skipped draw).
* Block U (review F5): a second exact panel with UNEQUAL level counts {3,3,5,.,.,.}
*   (60 assignments) — separates the level-counts-fixed scheme from iid shuffles,
*   pooled-treated-levels, and status-only schemes that the 20-assignment design
*   cannot distinguish.
* Block T (review F4): permtstat VALUE check on a heteroskedastic-cluster panel where
*   the |t| ranking differs from the |b| ranking: e(permdist) must equal the external
*   b/se multiset, not the b multiset, and the exact p must use the t ranking.
* Block P (review F6): permute(probabilities) draw-frequency check against the
*   specified asymmetric design (per-cluster marginal bands, membership, accounting).
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
* (TRI additionally records each assignment's treated triple, for block S's draw map)
tempname EXT TRI
matrix `EXT' = J(20, 1, .)
matrix `TRI' = J(20, 3, .)
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
      matrix `TRI'[`k', 1] = `a'
      matrix `TRI'[`k', 2] = `b'
      matrix `TRI'[`k', 3] = `c'
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
  * equal-tailed convention (2026-09-05): p = min(1, 2 min(#{s <= obs}, #{s >= obs})/M); __nex
  * reports the smaller tail count for the diagnostics below
  mata: st_numscalar("__nex", min((sum(st_matrix(st_local("EXT")) :<= st_numscalar(st_local("OBS"))), sum(st_matrix(st_local("EXT")) :>= st_numscalar(st_local("OBS"))))))
  tempname PEXT
  scalar `PEXT' = min(1, 2*__nex/20)
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

**************************** S: SAMPLED-mode external check (review F3) ****************************
* permute(shuffle) permreps(200), fixed seed, same panel and assignment set. Checks:
*  (Sa) every e(permdist) value is a bitwise member of the 20 external coefficients
*       (catches wrong-scheme draws);
*  (Sb) per-assignment frequencies within a generous binomial band: >= 1 and within
*       5 sd of 200/20 (catches systematic reuse / degenerate mixing);
*  (Sc) completed==200, failed==0, npossible still 20, permdist 200x1;
*  (Sd) rep-by-rep equality with an EXTERNALLY REPLICATED draw stream: with the same
*       seed, the sampled branch of wooldid_permx (per-rep jumble of the cluster-level
*       onset vector, frame rows in ascending cluster order, one stratum here) is
*       reproduced in mata and each rep's predicted coefficient is compared to
*       e(permdist)[rep]. This is the only check with power against a SINGLE reused or
*       skipped draw (the mutant class the review demonstrated): such a draw is itself
*       a valid member and moves each frequency count by only 1, so (Sa)/(Sb)/(Sc) are
*       structurally blind to it. (Sd) deliberately pins the current drawing scheme
*       and RNG consumption under this seed; if the drawing scheme legitimately
*       changes, update the replication loop to match.
*  (Se) the sampled add-one p: e(perm_p_att) == (1+nex)/(nuse+1) recomputed from
*       e(permdist) and the rep-0 observed statistic.
if `extfail' == 0 {
  cap drop ttre2
  sort i t
  * external draw-stream replication (must IMMEDIATELY follow the seed set). The
  * sampled branch of wooldid_permx draws, per rep, jumble(ons[idx]) with idx = the
  * single stratum's full row set, i.e. one jumble() call per rep on the same
  * 6-element onset vector — bitwise-identical stream consumption to the line below.
  set seed 5150
  mata: ons = (4\4\4\.\.\.)
  mata: A = J(6, 200, .)
  mata: for (r=1; r<=200; r++) A[., r] = jumble(ons)
  * map each replicated draw to its enumerated assignment via a treated-set bitmask
  mata: codeA = (2:^(0::5))' * (A :== 4)
  mata: T_ = st_matrix(st_local("TRI")); codeT = 2:^(T_[.,1]:-1) + 2:^(T_[.,2]:-1) + 2:^(T_[.,3]:-1)
  mata: E_ = st_matrix(st_local("EXT")); LU = J(63, 1, .)
  mata: for (kk=1; kk<=20; kk++) LU[codeT[kk,1], 1] = E_[kk,1]
  mata: PRED = LU[codeA', 1]
  mata: st_numscalar("__sunmap", sum(PRED :>= .))
  mata: st_numscalar("__sguard", PRED[11,1] != PRED[12,1])
  if __sunmap != 0 {
    dis as error "S1 FAIL: " __sunmap " externally replicated draws map to no enumerated assignment (replication broken)"
    global pf = $pf + 1
  }
  * design guard: under THIS seed the review's escaped-mutant site (rep 12 reusing rep
  * 11's assignment) produces distinct predicted coefficients, so (Sd) retains power
  * against that demonstrated escape class
  if __sguard != 1 {
    dis as error "S1b FAIL: seed guard - predicted draws 11 and 12 coincide; pick a different seed"
    global pf = $pf + 1
  }

  set seed 5150
  cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(200)
  if _rc != 0 {
    dis as error "S2 FAIL: sampled engine run rc=" _rc
    global pf = $pf + 1
  }
  else {
    tempname PD2 peng2 npos2 ncomp2 nfail2
    matrix `PD2' = e(permdist)
    scalar `peng2' = e(perm_p_att)
    scalar `npos2' = e(perm_npossible)
    scalar `ncomp2' = e(perm_completed)
    scalar `nfail2' = e(perm_failed)
    * (Sc) accounting
    if `ncomp2' != 200 | `nfail2' != 0 | rowsof(`PD2') != 200 | colsof(`PD2') != 1 | `npos2' != 20 {
      dis as error "S3 FAIL: completed=" `ncomp2' " failed=" `nfail2' " permdist " rowsof(`PD2') "x" colsof(`PD2') " npossible=" `npos2' " (want 200/0/200x1/20)"
      global pf = $pf + 1
    }
    mata: pd = st_matrix(st_local("PD2")); E_ = st_matrix(st_local("EXT"))
    mata: comp = ((pd*J(1, rows(E_), 1)) :== (J(rows(pd), 1, 1)*E_'))
    mata: st_numscalar("__snonm", sum(rowsum(comp) :== 0))
    mata: cnt = colsum(comp)'
    mata: st_numscalar("__scmin", min(cnt)); st_numscalar("__scmax", max(cnt))
    mata: Es = sort(E_, 1); st_numscalar("__semind", min(Es[|2,1 \ 20,1|] - Es[|1,1 \ 19,1|]))
    mata: st_numscalar("__smaxrd", (rows(pd) == rows(PRED)) ? max(reldif(pd, PRED)) : .)
    mata: st_numscalar("__snbit", (rows(pd) == rows(PRED)) ? sum(pd :== PRED) : 0)
    mata: st_numscalar("__snex", min((sum(pd :<= st_numscalar(st_local("OBS"))), sum(pd :>= st_numscalar(st_local("OBS"))))))
    mata: st_numscalar("__snuse", rows(pd))
    * distinctness guard: the by-value frequency keying needs 20 distinct coefficients
    if __semind <= 0 {
      dis as error "S4a FAIL: external coefficient set is not distinct (min sorted gap = " %12.4e __semind ")"
      global pf = $pf + 1
    }
    * (Sa) bitwise membership of every sampled draw
    if __snonm != 0 {
      dis as error "S4 FAIL: " __snonm " of 200 sampled e(permdist) values are not bitwise members of the external 20-coefficient set"
    global pf = $pf + 1
    }
    * (Sb) frequency band: 200/20 = 10 per assignment, sd = sqrt(200*.05*.95) = 3.082
    if __scmin < 1 | __scmax > 10 + 5*sqrt(200*0.05*0.95) {
      dis as error "S5 FAIL: per-assignment frequency out of band: min=" __scmin " max=" __scmax " (want within [1, " %6.2f 10 + 5*sqrt(200*0.05*0.95) "])"
      global pf = $pf + 1
    }
    dis "S5 info: per-assignment frequencies span [" __scmin ", " __scmax "] over 200 draws (uniform mean 10)"
    * (Sd) rep-by-rep vs the externally replicated stream
    if __smaxrd > 1e-10 {
      dis as error "S6 FAIL: rep-by-rep sampled draws vs external replication: max reldif = " %12.4e __smaxrd " (> 1e-10); bitwise-equal reps = " __snbit " of 200"
      global pf = $pf + 1
    }
    dis "S6 info: rep-by-rep max reldif = " %12.4e __smaxrd "; bitwise-equal reps = " __snbit " of 200"
    * (Se) sampled-mode add-one p
    tempname PS
    scalar `PS' = min(1, 2*(1 + __snex)/(__snuse + 1))
    if reldif(`PS', `peng2') > 1e-12 {
      dis as error "S7 FAIL: add-one p recomputed = " `PS' " (nex=" __snex "/" __snuse ") vs e(perm_p_att) = " `peng2'
      global pf = $pf + 1
    }
    dis "S7 info: sampled add-one p = " `PS' " (engine " `peng2' ")"
  }
}
else dis as error "S block skipped: `extfail' external run(s) failed (already counted)"

**************************** U: unequal-level-counts exact design (review F5) ****************************
* onsets {3, 3, 5, ., ., .} across 6 clusters: 6!/(2! 1! 3!) = 60 distinct assignments.
* The equal-counts 20-assignment design cannot distinguish wooldid's level-counts-fixed
* scheme from an iid label shuffle, a pooled-treated-levels scheme, or a treated-status
* permutation that scrambles onset labels; with two unequal treated levels those
* schemes produce different supports/counts (pooled-treated gives 20; status-only and
* iid schemes leave the {3,5} split unfixed), so multiset equality with the externally
* enumerated 60-assignment set + npossible==60 kills them.
clear
qui set obs 6
qui gen i = _n
qui gen ttre = cond(_n <= 2, 3, cond(_n == 3, 5, .))
qui expand 8
qui bysort i: gen t = _n
set seed 7311
qui gen y = 2.5*(t >= ttre & ttre < .) + 0.1*i + rnormal()
sort i t

cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
if _rc != 0 {
  dis as error "U0 FAIL: unequal-counts engine run rc=" _rc
  global pf = $pf + 1
}
else {
  tempname UPD upeng unpos uncomp
  matrix `UPD' = e(permdist)
  scalar `upeng' = e(perm_p_att)
  scalar `unpos' = e(perm_npossible)
  scalar `uncomp' = e(perm_completed)
  if `unpos' != 60 | `uncomp' != 60 | rowsof(`UPD') != 60 | colsof(`UPD') != 1 {
    dis as error "U1 FAIL: npossible=" `unpos' " completed=" `uncomp' " permdist is " rowsof(`UPD') "x" colsof(`UPD') " (want 60x1)"
    global pf = $pf + 1
  }

  * observed statistic: fastnose on the true assignment
  sort i t
  cap noi wooldid y i t ttre, att cluster(i) fastnose
  if _rc != 0 {
    dis as error "U2 FAIL: observed fastnose run rc=" _rc
    global pf = $pf + 1
  }
  tempname UOBS
  scalar `UOBS' = _b[main_att]

  * external enumeration of the level-counts-fixed set: pick {a<b} for onset 3 and
  * c (distinct from a,b) for onset 5; the remaining three clusters stay never-treated
  tempname UEXT
  matrix `UEXT' = J(60, 1, .)
  local uk 0
  local uextfail 0
  forval a = 1/5 {
    forval b = `=`a'+1'/6 {
      forval c = 1/6 {
        if `c' != `a' & `c' != `b' {
          local ++uk
          cap drop ttre2
          qui gen double ttre2 = cond(inlist(i, `a', `b'), 3, cond(i == `c', 5, .))
          sort i t
          cap qui wooldid y i t ttre2, att cluster(i) fastnose
          if _rc != 0 {
            dis as error "U3 FAIL: external assignment `uk' (onset3: `a' `b'; onset5: `c') rc=" _rc
            local ++uextfail
            global pf = $pf + 1
          }
          else matrix `UEXT'[`uk', 1] = _b[main_att]
        }
      }
    }
  }
  cap drop ttre2
  if `uk' != 60 {
    dis as error "U3b FAIL: enumerated `uk' assignments (want 60)"
    global pf = $pf + 1
  }
  if `unpos' != `uk' {
    dis as error "U3c FAIL: e(perm_npossible)=" `unpos' " != external count `uk'"
    global pf = $pf + 1
  }

  if `uextfail' == 0 {
    * multiset equality (the separating assertion)
    mata: upd = sort(st_matrix(st_local("UPD")), 1); uex = sort(st_matrix(st_local("UEXT")), 1)
    mata: st_numscalar("__umaxrd", rows(upd) == rows(uex) ? max(reldif(upd, uex)) : .)
    mata: st_numscalar("__unbit", rows(upd) == rows(uex) ? sum(upd :== uex) : 0)
    if __umaxrd > 1e-10 {
      dis as error "UA FAIL: sorted external vs e(permdist) max reldif = " %12.4e __umaxrd " (> 1e-10)"
      global pf = $pf + 1
    }
    dis "UA info: max reldif (sorted) = " %12.4e __umaxrd "; bitwise-equal entries = " __unbit " of 60"

    * membership + exact p (plain share)
    mata: st_numscalar("__umember", anyof(st_matrix(st_local("UEXT")), st_numscalar(st_local("UOBS"))))
    if __umember != 1 {
      dis as error "UM FAIL: observed fastnose stat is not a bitwise member of the external 60-assignment set"
      global pf = $pf + 1
    }
    mata: st_numscalar("__unex", min((sum(st_matrix(st_local("UEXT")) :<= st_numscalar(st_local("UOBS"))), sum(st_matrix(st_local("UEXT")) :>= st_numscalar(st_local("UOBS"))))))
    tempname UP
    scalar `UP' = min(1, 2*__unex/60)
    if reldif(`UP', `upeng') > 1e-12 {
      dis as error "UB FAIL: external exact p = " `UP' " (nex=" __unex "/60) vs e(perm_p_att) = " `upeng'
      global pf = $pf + 1
    }
    dis "UB info: exact p = " `UP' " (engine " `upeng' ")"
  }
  else dis as error "UA/UB skipped: `uextfail' external run(s) failed (already counted)"
}

**************************** T: permtstat value check (review F4) ****************************
* Heteroskedastic clusters (sd 0.3 for i<=3, sd 3 for i>3) make the |t| ranking differ
* from the |b| ranking across the 20 assignments (assignments that treat the noisy
* clusters get large |b| but larger se). The external enumeration uses FULL
* (non-fastnose) wooldid runs — permtstat rep runs go through the full inference path
* (fastnose is refused alongside permtstat) — capturing _b and _se of att. Asserts:
*   e(permdist) multiset == external b/se multiset (NOT the b multiset), the ranking
*   flip is present (nex_t != nex_b, so the two exact p-values differ), and the exact
*   p uses the t ranking (plain share; floor 1/M under permtstat).
clear
qui set obs 6
qui gen i = _n
qui gen ttre = cond(_n <= 3, 4, .)
qui expand 8
qui bysort i: gen t = _n
set seed 31415
qui gen y = 1.5*(t >= ttre & ttre < .) + 0.2*i + rnormal()*cond(i <= 3, 0.3, 3)
sort i t

* external enumeration first (full runs; observed = assignment {1,2,3} = entry 1)
tempname TB TT
matrix `TB' = J(20, 1, .)
matrix `TT' = J(20, 1, .)
local tk 0
local textfail 0
forval a = 1/4 {
  forval b = `=`a'+1'/5 {
    forval c = `=`b'+1'/6 {
      local ++tk
      cap drop ttre2
      qui gen double ttre2 = cond(inlist(i, `a', `b', `c'), 4, .)
      sort i t
      cap qui wooldid y i t ttre2, att cluster(i)
      if _rc != 0 {
        dis as error "T0 FAIL: external full run `tk' (treated `a' `b' `c') rc=" _rc
        local ++textfail
        global pf = $pf + 1
      }
      else {
        matrix `TB'[`tk', 1] = _b[main_att]
        matrix `TT'[`tk', 1] = _b[main_att]/_se[main_att]
      }
    }
  }
}
cap drop ttre2

sort i t
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permtstat
if _rc != 0 {
  dis as error "T1 FAIL: permtstat engine run rc=" _rc
  global pf = $pf + 1
}
else if `textfail' == 0 {
  tempname TPD tpeng
  matrix `TPD' = e(permdist)
  scalar `tpeng' = e(perm_p_att)
  if "`e(perm_stat)'" != "t" | e(perm_completed) != 20 | e(perm_npossible) != 20 | rowsof(`TPD') != 20 {
    dis as error "T2 FAIL: perm_stat=[`e(perm_stat)'] completed=" e(perm_completed) " npossible=" e(perm_npossible) " permdist rows=" rowsof(`TPD') " (want t/20/20/20)"
    global pf = $pf + 1
  }
  mata: tpd = sort(st_matrix(st_local("TPD")), 1); tt = sort(st_matrix(st_local("TT")), 1); tb = sort(st_matrix(st_local("TB")), 1)
  mata: st_numscalar("__tmaxrd", max(reldif(tpd, tt))); st_numscalar("__tnbit", sum(tpd :== tt))
  mata: st_numscalar("__tbrd", max(reldif(tpd, tb)))
  mata: ttf = st_matrix(st_local("TT")); tbf = st_matrix(st_local("TB"))
  mata: st_numscalar("__tnext", sum(abs(ttf) :>= abs(ttf[1,1]))); st_numscalar("__tnexb", sum(abs(tbf) :>= abs(tbf[1,1])))
  * equal-tailed convention (2026-09-05): the engine's exact p is min(1, 2 min(#{t <= obs}, #{t >= obs})/M) on
  * the SIGNED t-ratios; the |t| counts above stay as the ranking-flip design guard (TC)
  mata: st_numscalar("__tnexeq", min((sum(ttf :<= ttf[1,1]), sum(ttf :>= ttf[1,1]))))
  * (i) the value check: e(permdist) is the t multiset ...
  if __tmaxrd > 1e-10 {
    dis as error "TA FAIL: sorted e(permdist) vs external t multiset (b/se): max reldif = " %12.4e __tmaxrd " (> 1e-10)"
    global pf = $pf + 1
  }
  dis "TA info: t-multiset max reldif = " %12.4e __tmaxrd "; bitwise-equal entries = " __tnbit " of 20"
  * ... and NOT the b multiset (design guarantees a wide gap)
  if __tbrd <= 1e-2 {
    dis as error "TB FAIL: e(permdist) is indistinguishable from the b multiset (max reldif " %12.4e __tbrd "); the design lost its separating power"
    global pf = $pf + 1
  }
  dis "TB info: permdist vs b multiset max reldif = " %12.4e __tbrd " (b-impostor would show ~0)"
  * (ii) ranking flip present: the exact p differs between the t and b rankings
  if __tnext == __tnexb {
    dis as error "TC FAIL: nex_t = nex_b = " __tnext "; the design no longer flips the ranking"
    global pf = $pf + 1
  }
  dis "TC info: nex_t = " __tnext ", nex_b = " __tnexb " of 20"
  * (iii) the engine's exact p uses the t ranking (plain share)
  tempname TP
  scalar `TP' = min(1, 2*__tnexeq/20)
  if reldif(`TP', `tpeng') > 1e-12 {
    dis as error "TD FAIL: external exact equal-tailed t-p = " `TP' " (smaller tail " __tnexeq "/20) vs e(perm_p_att) = " `tpeng' " (|t| ranking would give " __tnext "/20, |b| " __tnexb "/20)"
    global pf = $pf + 1
  }
  dis "TD info: exact t-p = " `TP' " (engine " `tpeng' ")"
}
else dis as error "TA-TD skipped: `textfail' external run(s) failed (already counted)"

**************************** P: permute(probabilities) frequency check (review F6) ****************************
* 6 clusters; clusters 1-4 draw onset 4 vs never with deliberately ASYMMETRIC
* probabilities (.8/.6/.4/.2), clusters 5-6 are anchored never-treated (p=1 on never,
* keeping every non-empty draw estimable). permreps(400), fixed seed. The support has
* 2^4 = 16 assignments; the all-never draw (probability q = .2*.4*.6*.8 = .0384) fails
* to estimate and is excluded, so per-cluster marginals conditional on completion are
* p_c/(1-q). Externally enumerate the 15 estimable assignments (fastnose), map each
* e(permdist) value to its treated subset by bitwise membership, and assert each
* cluster's realized treatment frequency within +-5 sd binomial bands, plus membership,
* accounting (npossible = product of support sizes = 16; completed+failed = 400;
* failed within its own band), and the sampled add-one p.
clear
qui set obs 6
qui gen i = _n
qui gen ttre = cond(_n <= 3, 4, .)
qui expand 8
qui bysort i: gen t = _n
set seed 2718
qui gen y = 2*(t >= ttre & ttre < .) + 0.15*i + rnormal()
sort i t

* external enumeration: the 15 non-empty subsets of {1,2,3,4} (bitmask order)
tempname QEXT QSUB
matrix `QEXT' = J(15, 1, .)
matrix `QSUB' = J(15, 4, 0)
local qextfail 0
forval m = 1/15 {
  cap drop ttre2
  qui gen double ttre2 = .
  forval c = 1/4 {
    if mod(floor(`m'/2^(`c'-1)), 2) == 1 {
      qui replace ttre2 = 4 if i == `c'
      matrix `QSUB'[`m', `c'] = 1
    }
  }
  sort i t
  cap qui wooldid y i t ttre2, att cluster(i) fastnose
  if _rc != 0 {
    dis as error "Q0 FAIL: external subset `m' rc=" _rc
    local ++qextfail
    global pf = $pf + 1
  }
  else matrix `QEXT'[`m', 1] = _b[main_att]
}
cap drop ttre2

* observed statistic (rep-0 fastnose on the true assignment = subset {1,2,3})
sort i t
cap noi wooldid y i t ttre, att cluster(i) fastnose
if _rc != 0 {
  dis as error "Q1 FAIL: observed fastnose run rc=" _rc
  global pf = $pf + 1
}
tempname QOBS
scalar `QOBS' = _b[main_att]

sort i t
set seed 424242
cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(400) permprobs((ttre 4 .) (1 .8 .2) (2 .6 .4) (3 .4 .6) (4 .2 .8) (5 0 1) (6 0 1))
if _rc != 0 {
  dis as error "Q2 FAIL: probabilities engine run rc=" _rc
  global pf = $pf + 1
}
else if `qextfail' == 0 {
  tempname QPD qpeng qnpos qncomp qnfail
  matrix `QPD' = e(permdist)
  scalar `qpeng' = e(perm_p_att)
  scalar `qnpos' = e(perm_npossible)
  scalar `qncomp' = e(perm_completed)
  scalar `qnfail' = e(perm_failed)
  * accounting: npossible = 2*2*2*2*1*1 = 16; completed + failed = 400; failed draws
  * are the all-never assignments: E = 400q = 15.36, sd = sqrt(400q(1-q)) = 3.84
  if `qnpos' != 16 | `qncomp' + `qnfail' != 400 | rowsof(`QPD') != `qncomp' {
    dis as error "Q3 FAIL: npossible=" `qnpos' " completed=" `qncomp' " failed=" `qnfail' " permdist rows=" rowsof(`QPD') " (want 16 / sum 400 / rows==completed)"
    global pf = $pf + 1
  }
  if `qnfail' > 400*.0384 + 5*sqrt(400*.0384*(1 - .0384)) {
    dis as error "Q4 FAIL: failed=" `qnfail' " exceeds the 5-sd band around 400*q=15.36"
    global pf = $pf + 1
  }
  dis "Q4 info: failed draws = " `qnfail' " (expected ~15.4, 5-sd cap ~34.6)"
  mata: qpd = st_matrix(st_local("QPD")); QE = st_matrix(st_local("QEXT")); QS = st_matrix(st_local("QSUB"))
  mata: QEs = sort(QE, 1); st_numscalar("__qmind", min(QEs[|2,1 \ 15,1|] - QEs[|1,1 \ 14,1|]))
  mata: qcomp = ((qpd*J(1, rows(QE), 1)) :== (J(rows(qpd), 1, 1)*QE'))
  mata: st_numscalar("__qnonm", sum(rowsum(qcomp) :== 0))
  mata: ccnt = (colsum(qcomp) * QS)'
  mata: st_matrix("__qccnt", ccnt)
  * equal-tailed convention (2026-09-05): the smaller signed tail count
  mata: st_numscalar("__qnex", min((sum(qpd :<= st_numscalar(st_local("QOBS"))), sum(qpd :>= st_numscalar(st_local("QOBS")))))); st_numscalar("__qnuse", rows(qpd))
  if __qmind <= 0 {
    dis as error "Q5a FAIL: external subset coefficients are not distinct (min sorted gap = " %12.4e __qmind "); the frequency mapping is ill-defined"
    global pf = $pf + 1
  }
  * membership: every completed draw's coefficient belongs to the 15-subset support
  * (a draw treating an anchored cluster or an off-design onset would land outside)
  if __qnonm != 0 {
    dis as error "Q5 FAIL: " __qnonm " of " __qnuse " e(permdist) values are not bitwise members of the enumerated design support"
    global pf = $pf + 1
  }
  * per-cluster marginal frequency bands (conditional on completion): p_c/(1-q)
  local qq = .0384
  local qpc1 = .8/(1 - `qq')
  local qpc2 = .6/(1 - `qq')
  local qpc3 = .4/(1 - `qq')
  local qpc4 = .2/(1 - `qq')
  forval c = 1/4 {
    local pc = `qpc`c''
    local mu = __qnuse*`pc'
    local sd = sqrt(__qnuse*`pc'*(1 - `pc'))
    local ob = __qccnt[`c', 1]
    if abs(`ob' - `mu') > 5*`sd' {
      dis as error "Q6 FAIL: cluster `c' treated in " `ob' " of " __qnuse " completed draws; expected " %8.2f `mu' " +- " %8.2f 5*`sd'
      global pf = $pf + 1
    }
    else dis "Q6 info: cluster `c' treated " `ob' " of " __qnuse " (expected " %8.2f `mu' " +- " %8.2f 5*`sd' ")"
  }
  * sampled add-one p from the design distribution
  tempname QP
  scalar `QP' = min(1, 2*(1 + __qnex)/(__qnuse + 1))
  if reldif(`QP', `qpeng') > 1e-12 {
    dis as error "Q7 FAIL: add-one p recomputed = " `QP' " (nex=" __qnex "/" __qnuse ") vs e(perm_p_att) = " `qpeng'
    global pf = $pf + 1
  }
  dis "Q7 info: add-one p = " `QP' " (engine " `qpeng' ")"
}
else dis as error "Q5-Q7 skipped: `qextfail' external run(s) failed (already counted)"

if $pf == 0 dis "RI_CROSSCHECK_ALL_PASS"
else dis as error "RI_CROSSCHECK_FAILS=$pf"
exit, clear
