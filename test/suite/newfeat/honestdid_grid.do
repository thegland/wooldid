* honestdid_grid.do — preregistered M5 grid (test/suite/PREREG_NEWFEAT_2026-08-24.md):
* wooldid -> HonestDiD positional handoff validated across 11 cases spanning dta_binary
* (baseline, anticipation(1), asymmetric window, aw weights, per-subgroup ES), the
* realistic admin panel dta_s_admin (incl. aw=pop), the rich suite datasets binary_fine
* and binary_mixed_imbalanced, and dta_cont (att side). Every case satisfies the
* wooldid.sthlp e(b_es) handoff geometry: esfixedbaseperiod at the default
* esrelativeto(-1), >= 1 pre + >= 1 post horizon, consecutive horizons (label
* signature: rightmost pre column pre2, post0 present, no other gaps), no degenerate
* VCV rows (all probed clean 2026-08-24). Per case honestdid runs twice — delta(rm)
* and delta(sd), each with mvec(0.5(0.5)2). Preregistered rules per case:
*   - honestdid rc 0 (task carve-out: a delta(sd) refusal with a clear rc on a case
*     where delta(rm) ran is recorded as an info FINDING line, not a script failure);
*   - the Original row reproduces the wooldid post0 CI (est +/- invnormal(.975)*se
*     from e(b_es_*)/e(V_es_*)) at reldif <= 1e-6;
*   - every robust row finite with lb < ub;
*   - delta(rm): CI length monotone nondecreasing in M (jitter tolerance 1e-8).
* The breakdown M (largest mvec M whose robust CI excludes 0) is recorded per case
* and delta as an info line (a record for the report, not pass/fail).
* Layout note (verified empirically in scratch, honestdid v1.1.2): delta(sd) returns
* the SAME 5x3 s(HonestEventStudy).CI layout as delta(rm) — row 1 Original with a
* missing M entry, rows 2..5 carrying the mvec M values — so both deltas share one
* extraction/check path. delta(sd) runs in ~1s per case here (k <= 6), so the full
* mvec is kept everywhere (no runtime reduction needed).
* Args: 1 = path to wooldid.ado (file or dir); 2 = data dir = <repo>/test/data
*       (the suite's own data dir is derived as <arg2>/../suite/data).
* Run:  "<stata>" -e do honestdid_grid.do <adopath> <datadir>
*    -> expect HONESTDID_GRID_ALL_PASS
*
* ADJUDICATED PREREG DEVIATION (2026-08-24, documented in the campaign report): the
* authoring smoke runs surfaced that on data with sharply significant pre-period
* coefficients, honestdid's C-LF robust confidence sets come back EMPTY — (., .)
* rows. That is a documented possibility of the Rambachan-Roth procedure (the data
* reject the imposed restriction class at that M), NOT a handoff defect: on the same
* inputs the Original row reproduces the wooldid post0 CI at reldif 0.0 (22/22
* extractions) and delta(sd) runs clean on the same b/V. The preregistered
* "every robust row finite" rule did not anticipate legitimately-empty sets, so
* empty rows are accepted ONLY on the explicit allowlist below, verified case by
* case under honestdid v1.1.2 (pre2 t = -4.1 on g09_fine, pre t = -5.5 on dta_cont;
* g10_mixed_imb rm M=0.5 empty with the finite rows above it flagged open at the
* upper endpoint by honestdid itself). An empty row NOT on the allowlist FAILs; a
* half-missing row FAILs everywhere; an allowlisted row turning finite prints an
* info line (likely an upstream honestdid change — retire the entry then).
* Known upstream package bug, recorded as a FINDING (not a FAIL): g10_mixed_imb
* delta(sd) crashes honestdid v1.1.2 with mata assert rc 3498 in
* _flciFindWorstCaseBiasGivenH (parallel and sequential paths both assert).

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 8642
set sortseed 8642

local ADO `1'
if `"`ADO'"' == "" local ADO /tmp/wdsuite_build/newfeat/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local DATA `2'
if `"`DATA'"' == "" local DATA /Users/hegland/Documents/Research/wooldid_dev/test/data
local SDATA `DATA'/../suite/data

global pf 0
* verified legitimately-empty C-LF rows under honestdid v1.1.2 (case|delta|row; see header)
global HDG_EMPTYOK "g09_fine|rm|2 g09_fine|rm|3 g09_fine|rm|4 g09_fine|rm|5 g11_cont_att|rm|2 g11_cont_att|rm|3 g11_cont_att|rm|4 g11_cont_att|rm|5 g10_mixed_imb|rm|2"
global ncase 0

* ============================================================================
* gridcase — one grid case: geometry guard, then honestdid delta(rm) and
* delta(sd) each at mvec(0.5(0.5)2) with the preregistered checks.
* Expects GRID_B (1 x k event-study b) and GRID_V (k x k VCV) already set.
* ============================================================================
cap program drop gridcase
program define gridcase
  syntax , case(string) npre(string)
  global ncase = $ncase + 1

  * ---- inputs available? (graceful-continue guards: a lost case FAILs, the grid goes on)
  cap confirm matrix GRID_B
  local rcb = _rc
  cap confirm matrix GRID_V
  if `rcb' != 0 | _rc != 0 {
    dis as error "`case' FAIL: event-study matrices unavailable (GRID_B/GRID_V not set)"
    global pf = $pf + 1
    exit
  }
  local np = real("`npre'")
  if `np' == . | `np' < 1 | `np' != int(`np') {
    dis as error "`case' FAIL: invalid e(numpre) [`npre']"
    global pf = $pf + 1
    exit
  }
  local k = colsof(GRID_B)
  if `k' < `np' + 1 {
    dis as error "`case' FAIL: geometry — numpre=`np' with only `k' horizons (need >= 1 pre and >= 1 post)"
    global pf = $pf + 1
    exit
  }

  * ---- geometry guard: the full sthlp handoff label signature —
  * pre(np+1)..pre2 then post0..post(k-np-1), i.e. consecutive horizons apart from
  * the intended jump over the omitted -1 reference (pre1 absent, post0 present)
  local lbl : colnames GRID_B
  forval j = 1/`k' {
    local lj : word `j' of `lbl'
    if `j' <= `np' local expnum = `np' + 2 - `j'
    else local expnum = `j' - `np' - 1
    local expect = cond(`j' <= `np', "pre`expnum'", "post`expnum'")
    if "`lj'" != "`expect'" {
      dis as error "`case' FAIL: geometry — b_es entry `j' is [`lj'], expected [`expect'] under the esfixedbaseperiod/esrelativeto(-1) handoff signature"
      global pf = $pf + 1
    }
  }

  * ---- wooldid post0 CI target (est +/- invnormal(.975)*se; the crosscheck's H3 convention)
  tempname EST SE LB UB
  scalar `EST' = GRID_B[1, `np'+1]
  scalar `SE'  = sqrt(GRID_V[`np'+1, `np'+1])
  scalar `LB'  = `EST' - invnormal(.975)*`SE'
  scalar `UB'  = `EST' + invnormal(.975)*`SE'

  matrix b = GRID_B
  matrix V = GRID_V
  local rmok 0
  foreach dl in rm sd {
    cap noi honestdid, numpre(`np') b(b) vcov(V) delta(`dl') mvec(0.5(0.5)2)
    local hrc = _rc
    if `hrc' != 0 {
      if "`dl'" == "rm" {
        dis as error "`case' rm FAIL: honestdid rc=`hrc'"
        global pf = $pf + 1
      }
      else if `rmok' {
        * preregistered task carve-out: an sd refusal where rm ran clean is a
        * recordable finding for the author, not a script failure
        dis "`case' sd FINDING: honestdid delta(sd) refused with rc=`hrc' on a case where delta(rm) ran clean"
      }
      else {
        dis as error "`case' sd FAIL: honestdid rc=`hrc' (delta(rm) also failed on this case)"
        global pf = $pf + 1
      }
      continue
    }
    if "`dl'" == "rm" local rmok 1

    * ---- H2b-style guarded CI extraction from s(HonestEventStudy)
    local hobj `s(HonestEventStudy)'
    if `"`hobj'"' == "" {
      dis as error "`case' `dl' FAIL: honestdid returned no s(HonestEventStudy) object"
      global pf = $pf + 1
      continue
    }
    tempname HCI
    cap noi mata: st_matrix(st_local("HCI"), `hobj'.CI)
    local mrc = _rc
    cap confirm matrix `HCI'
    local crc = _rc
    if `mrc' != 0 | `crc' != 0 {
      dis as error "`case' `dl' FAIL: could not extract the CI matrix (mata rc=`mrc', confirm rc=`crc')"
      global pf = $pf + 1
      continue
    }
    if rowsof(`HCI') != 5 | colsof(`HCI') != 3 {
      dis as error "`case' `dl' FAIL: CI matrix is " rowsof(`HCI') "x" colsof(`HCI') " (want 5x3: Original + mvec 0.5(0.5)2)"
      global pf = $pf + 1
      continue
    }

    * ---- Original row (row 1, missing M entry) reproduces the wooldid post0 CI at 1e-6
    if !missing(`HCI'[1, 1]) {
      dis as error "`case' `dl' FAIL: row 1 M = " `HCI'[1,1] " (expected missing, marking the Original row)"
      global pf = $pf + 1
    }
    if reldif(`HCI'[1, 2], `LB') > 1e-6 | reldif(`HCI'[1, 3], `UB') > 1e-6 {
      dis as error "`case' `dl' FAIL: Original CI (" `HCI'[1,2] ", " `HCI'[1,3] ") vs wooldid post0 CI (" `LB' ", " `UB' ")"
      dis as error "         reldif lb = " %12.4e reldif(`HCI'[1,2], `LB') ", ub = " %12.4e reldif(`HCI'[1,3], `UB')
      global pf = $pf + 1
    }
    else dis "`case' `dl' info: Original CI reldif lb = " %12.4e reldif(`HCI'[1,2], `LB') ", ub = " %12.4e reldif(`HCI'[1,3], `UB')

    * ---- robust rows: M values match the mvec; every row finite with lb < ub,
    * EXCEPT rows on the verified legitimately-empty allowlist (see header): there a
    * fully-empty (., .) row is a documented Rambachan-Roth outcome, a half-missing
    * row still FAILs, and a finite row prints an upstream-change info line.
    forval r = 2/5 {
      if reldif(`HCI'[`r', 1], (`r'-1)*0.5) > 1e-12 {
        dis as error "`case' `dl' FAIL: row `r' M = " `HCI'[`r',1] " (expected " (`r'-1)*0.5 ")"
        global pf = $pf + 1
      }
      local allowempty = strpos(" $HDG_EMPTYOK ", " `case'|`dl'|`r' ") > 0
      if missing(`HCI'[`r', 2]) & missing(`HCI'[`r', 3]) & `allowempty' {
        dis "`case' `dl' info: row `r' (M=" (`r'-1)*0.5 ") is a documented legitimately-empty C-LF set (significant pre-trends; see header)"
        continue
      }
      if missing(`HCI'[`r', 2]) | missing(`HCI'[`r', 3]) | `HCI'[`r', 2] >= `HCI'[`r', 3] {
        dis as error "`case' `dl' FAIL: row `r' (M=" (`r'-1)*0.5 ") CI (" `HCI'[`r',2] ", " `HCI'[`r',3] ") is nonfinite or inverted"
        global pf = $pf + 1
      }
      else if `allowempty' {
        dis "`case' `dl' info: allowlisted row `r' came back FINITE (" `HCI'[`r',2] ", " `HCI'[`r',3] ") — likely an upstream honestdid change; consider retiring the allowlist entry"
      }
    }

    * ---- delta(rm) only: CI length monotone nondecreasing in M (jitter tolerance 1e-8).
    * A comparison against a row without finite bounds is undefined and that row has
    * already FAILed the finiteness rule, so only pairs of finite lengths are compared.
    if "`dl'" == "rm" {
      forval r = 3/5 {
        local lenprev = `HCI'[`r'-1, 3] - `HCI'[`r'-1, 2]
        local lencur  = `HCI'[`r', 3] - `HCI'[`r', 2]
        if !missing(`lenprev') & !missing(`lencur') & !(`lencur' >= `lenprev' - 1e-8) {
          dis as error "`case' rm FAIL: CI length not monotone in M — len(M=" (`r'-1)*0.5 ") = " %21.0g `lencur' " < len(M=" (`r'-2)*0.5 ") = " %21.0g `lenprev'
          global pf = $pf + 1
        }
      }
    }

    * ---- breakdown M (info, preregistered record; not pass/fail)
    local bkm .
    forval r = 2/5 {
      if !missing(`HCI'[`r',2]) & !missing(`HCI'[`r',3]) & (`HCI'[`r',2] > 0 | `HCI'[`r',3] < 0) local bkm = (`r'-1)*0.5
    }
    if `bkm' == . dis "`case' `dl' info: breakdown M — no mvec M has a robust CI excluding 0 (already includes 0 at M=0.5)"
    else if `bkm' == 2 dis "`case' `dl' info: breakdown M > 2 (the largest mvec M, 2, still has a robust CI excluding 0; breakdown beyond the grid)"
    else dis "`case' `dl' info: breakdown M = `bkm' (largest mvec M with a robust CI excluding 0)"
  }
end

* feedcase — pull e(b_es_att`sfx')/e(V_es_att`sfx')/e(numpre_att`sfx') from the
* wooldid run just completed (BEFORE any rclass/sclass package call) and dispatch.
cap program drop feedcase
program define feedcase
  syntax , case(string) [sfx(string)]
  cap matrix drop GRID_B
  cap matrix drop GRID_V
  cap noi matrix GRID_B = e(b_es_att`sfx')
  cap noi matrix GRID_V = e(V_es_att`sfx')
  gridcase, case(`case') npre(`=e(numpre_att`sfx')')
end

* ============================================================================
* the 11 grid cases (a failing wooldid run FAILs its case(s) and the grid goes on)
* ============================================================================

* ---- g01_bin_base: dta_binary, the suite's basic ES window under esfixedbaseperiod
* (the existing honestdid_crosscheck.do covers esprelength(4) espostlength(4); the
* grid baseline deliberately uses the 3/3 window so the single case is not duplicated)
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod
if _rc != 0 {
  dis as error "g01_bin_base FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g01_bin_base)

* ---- g02_bin_antic: anticipation(1) (help-file example 9 configuration, sans pretrendpower)
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod anticipation(1)
if _rc != 0 {
  dis as error "g02_bin_antic FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g02_bin_antic)

* ---- g03_bin_asym: asymmetric window esprelength(2) espostlength(4)
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(4) esfixedbaseperiod
if _rc != 0 {
  dis as error "g03_bin_asym FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g03_bin_asym)

* ---- g04_bin_aw: analytic weights
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre [aw=w], att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod
if _rc != 0 {
  dis as error "g04_bin_aw FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g04_bin_aw)

* ---- g05_bin_sg1 / g06_bin_sg2: per-subgroup handoff via e(b_es_att_s#)/e(V_es_att_s#).
* Both subgroups' objects are captured BEFORE the first honestdid call (honestdid is
* sclass and each case must be fed from a pristine wooldid capture).
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod subgroups(sg)
if _rc != 0 {
  dis as error "g05_bin_sg1 FAIL: wooldid subgroup run rc=" _rc
  dis as error "g06_bin_sg2 FAIL: wooldid subgroup run rc=" _rc
  global pf = $pf + 2
}
else {
  local np2 = e(numpre_att_s2)
  cap matrix drop GRID_B2
  cap matrix drop GRID_V2
  cap noi matrix GRID_B2 = e(b_es_att_s2)
  cap noi matrix GRID_V2 = e(V_es_att_s2)
  feedcase, case(g05_bin_sg1) sfx(_s1)
  cap matrix drop GRID_B
  cap matrix drop GRID_V
  cap noi matrix GRID_B = GRID_B2
  cap noi matrix GRID_V = GRID_V2
  gridcase, case(g06_bin_sg2) npre(`np2')
}

* ---- g07_adm_asym: realistic admin panel, cohort-level i, asymmetric window
use `"`SDATA'/dta_s_admin.dta"', clear
cap noi wooldid y cohid year ttre, att cluster(state) esprelength(2) espostlength(3) esfixedbaseperiod
if _rc != 0 {
  dis as error "g07_adm_asym FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g07_adm_asym)

* ---- g08_adm_aw: population-weighted admin panel (configs.csv adm_aw window)
use `"`SDATA'/dta_s_admin.dta"', clear
cap noi wooldid y cohid year ttre [aw=pop], att cluster(state) esprelength(2) espostlength(2) esfixedbaseperiod
if _rc != 0 {
  dis as error "g08_adm_aw FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g08_adm_aw)

* ---- g09_fine: rich suite dataset, fine-grained i (configs_ext ext_bin_ext_fine
* varlist/cluster; safety(off) per that config)
use `"`SDATA'/binary_fine.dta"', clear
cap noi wooldid y cohort t ttre, att cluster(cohort) esprelength(3) espostlength(3) esfixedbaseperiod safety(off)
if _rc != 0 {
  dis as error "g09_fine FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g09_fine)

* ---- g10_mixed_imb: rich suite dataset, mixed-grain i on an imbalanced panel
* (configs_ext ext_bin_ext_mixed_imbalanced varlist/cluster, but at the default
* esrelativeto(-1) — that config's esrelativeto(-2) would break the handoff geometry)
use `"`SDATA'/binary_mixed_imbalanced.dta"', clear
cap noi wooldid y mixed_i t ttre, att cluster(cohort) esprelength(3) espostlength(3) esfixedbaseperiod safety(off)
if _rc != 0 {
  dis as error "g10_mixed_imb FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g10_mixed_imb)

* ---- g11_cont_att: dta_cont att side (10 clusters; dataset diversity)
use `"`DATA'/dta_cont.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod
if _rc != 0 {
  dis as error "g11_cont_att FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else feedcase, case(g11_cont_att)

* ============================================================================
dis "honestdid_grid: $ncase case dispatches completed (11 expected), each under delta(rm) and delta(sd)"
if $ncase != 11 {
  dis as error "HONESTDID_GRID FAIL: only $ncase of 11 cases were dispatched"
  global pf = $pf + 1
}
if $pf == 0 dis "HONESTDID_GRID_ALL_PASS"
else dis as error "HONESTDID_GRID_FAILS=$pf"
exit, clear
