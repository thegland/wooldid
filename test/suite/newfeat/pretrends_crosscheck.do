* pretrends_crosscheck.do — preregistered M3a (test/suite/PREREG_NEWFEAT_2026-08-24.md):
* feed the SAME e(b_es_att)/e(V_es_att)/e(numpre_att) that wooldid's pretrendpower()
* diagnostics were computed from to Roth's pretrends package (v0.5.1, Caceres Bravo
* translation) on 4 designs: dta_binary at esprelength(4)/espostlength(4), an
* asymmetric window (2/4), a 3-pre-period window (3/3), and the realistic admin panel
* dta_s_admin (3/3). All designs keep numpre <= 4 (wooldid's exact mdt_nis quadrature
* is computed only for <= 4 pre-periods) and use esfixedbaseperiod + cluster.
* Preregistered rules per design:
*   (ii)  the package's slope-for-power-50/-80 (pretrends power #.#, numpre() b() v())
*         vs wooldid's mdt_nis rows of e(roth_power_att): |reldif| <= 5e-2 (wooldid =
*         exact quadrature; the package inverts its own power curve for the same
*         any-individual-pre-period-significant convention);
*   (iii) the package's power evaluated AT wooldid's mdt_nis(80) slope
*         (pretrends, numpre() b() v() slope()) must land in [.78, .82].
* Rule (i) — package fitted GLS slope vs e(roth_slope_att) — turned out to be
* structurally N/A: pretrends v0.5.1 reports NO data-fitted linear-trend slope.
* Its r(slope) only echoes the hypothesized slope (slope mode) or the power-solved
* slope (power mode), and the mata PreTrendsResults.slope field does the same
* (verified in the package sthlp r() list and empirically in scratch 2026-08-24).
* That structural difference is reported as a FINDING line, wooldid's GLS slope is
* recorded per design, and no substitute tolerance check is manufactured.
* Additional structural gate per design (pass/fail): the package's r(results)
* time/coefficient columns must reproduce wooldid's event-time geometry exactly —
* pre rows at e(escomp)'s relyr values, the omitted reference inserted at t = -1
* with a zero coefficient, post rows at their relyr values, and betahat equal to
* e(b_es_att) positionally. If that gate fails the design's tolerance checks (ii)
* and (iii) are SKIPPED and the structural difference reported (per the prereg
* directive: a convention mismatch is a finding, not something to absorb).
* Recorded info (not pass/fail): the package's implied pre-test size (its power at
* slope 0) vs e(roth_nissize_att); its Bayes factor and likelihood ratio at the
* mdt_nis(80) slope; its power at the mdt_nis(50) slope.
* State management: pretrends is rclass — everything needed from wooldid is captured
* into locals/matrices BEFORE the first package call; each design refits wooldid.
* Args: 1 = path to wooldid.ado (file or dir); 2 = data dir = <repo>/test/data
*       (the suite's own data dir is derived as <arg2>/../suite/data).
* Run:  "<stata>" -e do pretrends_crosscheck.do <adopath> <datadir>
*    -> expect PRETRENDS_CROSSCHECK_ALL_PASS

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
global ncase 0

dis "PRETRENDS_CROSSCHECK FINDING (prereg M3a check i): pretrends v0.5.1 reports no data-fitted"
dis "  linear-trend slope — r(slope) echoes the hypothesized or power-solved slope and mata"
dis "  PreTrendsResults.slope does the same — so check (i) vs e(roth_slope_att) has no package-side"
dis "  counterpart and is recorded as structurally N/A; wooldid's GLS slope is recorded per design."

* ============================================================================
* ptcase — one design: wooldid (with pretrendpower(50 80)) has just run with rc 0;
* capture its results, then run the pretrends package on the same b/V/numpre.
* ============================================================================
cap program drop ptcase
program define ptcase
  syntax , case(string)
  global ncase = $ncase + 1

  * ---- capture EVERYTHING from wooldid before any rclass package call ----
  local np     = e(numpre_att)
  local rnp    = e(roth_npre_att)
  local rok    = e(roth_ok)
  local wslope = e(roth_slope_att)
  local wslse  = e(roth_slope_se_att)
  local wnis   = e(roth_nissize_att)
  cap matrix drop PT_B
  cap matrix drop PT_V
  cap matrix drop PT_RP
  cap matrix drop PT_EC
  cap noi matrix PT_B  = e(b_es_att)
  cap noi matrix PT_V  = e(V_es_att)
  cap noi matrix PT_RP = e(roth_power_att)
  cap noi matrix PT_EC = e(escomp)

  * ---- W1: wooldid-side prerequisites ----
  cap confirm matrix PT_B
  local r1 = _rc
  cap confirm matrix PT_V
  local r2 = _rc
  cap confirm matrix PT_RP
  local r3 = _rc
  cap confirm matrix PT_EC
  local r4 = _rc
  if `r1' | `r2' | `r3' | `r4' | `np' == . | `np' < 1 | `rok' != 1 {
    dis as error "`case' W1 FAIL: wooldid outputs incomplete (b_es rc=`r1', V_es rc=`r2', roth_power rc=`r3', escomp rc=`r4', numpre=`np', roth_ok=`rok')"
    global pf = $pf + 1
    exit
  }
  if `rnp' != `np' {
    dis as error "`case' W1 FAIL: e(roth_npre_att)=`rnp' != e(numpre_att)=`np'"
    global pf = $pf + 1
  }
  if rowsof(PT_RP) != 2 | colsof(PT_RP) != 5 {
    dis as error "`case' W1 FAIL: e(roth_power_att) is " rowsof(PT_RP) "x" colsof(PT_RP) " (want 2x5 for pretrendpower(50 80))"
    global pf = $pf + 1
    exit
  }
  local cpow = colnumb(PT_RP, "power")
  local cnis = colnumb(PT_RP, "mdt_nis")
  if `cpow' == . | `cnis' == . {
    local rplbl : colnames PT_RP
    dis as error "`case' W1 FAIL: e(roth_power_att) lacks power/mdt_nis columns (got: `rplbl')"
    global pf = $pf + 1
    exit
  }
  if PT_RP[1, `cpow'] != 50 | PT_RP[2, `cpow'] != 80 {
    dis as error "`case' W1 FAIL: e(roth_power_att) power column is (" PT_RP[1,`cpow'] ", " PT_RP[2,`cpow'] ") (want 50, 80)"
    global pf = $pf + 1
    exit
  }
  local nis50 = PT_RP[1, `cnis']
  local nis80 = PT_RP[2, `cnis']
  if `nis50' == . | `nis80' == . | `nis50' <= 0 | `nis80' <= 0 {
    dis as error "`case' W1 FAIL: mdt_nis rows are not positive finite (nis50=`nis50', nis80=`nis80')"
    global pf = $pf + 1
    exit
  }
  local k = colsof(PT_B)
  dis "`case' info: numpre=`np', k=`k'; wooldid GLS pre-trend slope = " %12.0g `wslope' " (se " %12.0g `wslse' ") [no package counterpart — see the check-(i) FINDING above]"

  * ---- S1: structural geometry gate — the package must reproduce wooldid's
  * event-time layout from numpre() alone (slope(0) call; also yields the implied size)
  matrix b = PT_B
  matrix V = PT_V
  cap noi pretrends, numpre(`np') b(b) v(V) slope(0)
  if _rc != 0 {
    dis as error "`case' S1 FAIL: pretrends slope(0) rc=" _rc
    global pf = $pf + 1
    exit
  }
  local psize = r(Power)
  tempname RES
  cap noi matrix `RES' = r(results)
  local mrc = _rc
  cap confirm matrix `RES'
  local crc = _rc
  if `mrc' != 0 | `crc' != 0 {
    dis as error "`case' S1 FAIL: could not extract r(results) (copy rc=`mrc', confirm rc=`crc')"
    global pf = $pf + 1
    exit
  }
  local structok 1
  if rowsof(`RES') != `k' + 1 | colsof(`RES') < 2 {
    dis as error "`case' S1 FAIL: r(results) is " rowsof(`RES') "x" colsof(`RES') " (want " `k'+1 " rows = k horizons + inserted reference)"
    global pf = $pf + 1
    local structok 0
  }
  else {
    local cry = colnumb(PT_EC, "relyr")
    if `cry' == . | rowsof(PT_EC) != `k' {
      dis as error "`case' S1 FAIL: e(escomp) lacks a relyr column or has " rowsof(PT_EC) " rows (want `k')"
      global pf = $pf + 1
      local structok 0
    }
    else {
      * row alignment guard (H1-style): escomp rows are in b_es order by name
      local eslbl : colnames PT_B
      local eclbl : rownames PT_EC
      forval j = 1/`k' {
        local esl : word `j' of `eslbl'
        local ecl : word `j' of `eclbl'
        if "`esl'" != "`ecl'" {
          dis as error "`case' S1 FAIL: row alignment broken at `j': b_es [`esl'] vs escomp [`ecl']"
          global pf = $pf + 1
          local structok 0
          continue, break
        }
        * package row index: pre rows map 1:1; post rows sit after the inserted reference
        local pr = cond(`j' <= `np', `j', `j' + 1)
        if reldif(`RES'[`pr', 1], PT_EC[`j', `cry']) > 1e-12 {
          dis as error "`case' S1 FAIL: structural time mismatch at [`esl']: package t = " `RES'[`pr',1] " vs wooldid relyr = " PT_EC[`j',`cry'] " — the package's event-time convention differs from wooldid's"
          global pf = $pf + 1
          local structok 0
        }
        if reldif(`RES'[`pr', 2], PT_B[1, `j']) > 1e-12 {
          dis as error "`case' S1 FAIL: coefficient mismatch at [`esl']: package betahat = " %21.0g `RES'[`pr',2] " vs b_es = " %21.0g PT_B[1,`j']
          global pf = $pf + 1
          local structok 0
        }
      }
      * the inserted reference row: t = -1 (wooldid's fixed reference) with a zero coefficient
      if `RES'[`np'+1, 1] != -1 | `RES'[`np'+1, 2] != 0 {
        dis as error "`case' S1 FAIL: package reference row is (t=" `RES'[`np'+1,1] ", beta=" `RES'[`np'+1,2] ") (want t=-1, beta=0) — reference-period convention differs"
        global pf = $pf + 1
        local structok 0
      }
    }
  }
  dis "`case' info: package implied pre-test size (power at slope 0) = " %9.6f `psize' " vs wooldid e(roth_nissize_att) = " %9.6f `wnis' " (reldif " %12.4e reldif(`psize', `wnis') "; record only)"
  if `structok' == 0 {
    dis as error "`case': structural gate failed — tolerance checks (ii)/(iii) SKIPPED for this design; the package's convention difference reported above is the finding"
    exit
  }

  * ---- C2 (prereg rule ii): package slope-for-power vs wooldid mdt_nis, |reldif| <= 5e-2
  foreach pw in 50 80 {
    local target = cond(`pw' == 50, "0.5", "0.8")
    cap noi pretrends power `target', numpre(`np') b(b) v(V)
    if _rc != 0 {
      dis as error "`case' C2 FAIL: pretrends power `target' rc=" _rc
      global pf = $pf + 1
      continue
    }
    local pslope = r(slope)
    if reldif(`pslope', `nis`pw'') > 5e-2 {
      dis as error "`case' C2 FAIL: package slope-for-power-`pw' = " %12.0g `pslope' " vs wooldid mdt_nis(`pw') = " %12.0g `nis`pw'' " (reldif " %12.4e reldif(`pslope', `nis`pw'') " > 5e-2)"
      global pf = $pf + 1
    }
    else dis "`case' C2 info: slope-for-power-`pw' = " %12.0g `pslope' " vs mdt_nis(`pw') = " %12.0g `nis`pw'' " (reldif " %12.4e reldif(`pslope', `nis`pw'') ")"
  }

  * ---- C3 (prereg rule iii): package power AT wooldid's mdt_nis(80) in [.78, .82]
  cap noi pretrends, numpre(`np') b(b) v(V) slope(`nis80')
  if _rc != 0 {
    dis as error "`case' C3 FAIL: pretrends slope(mdt_nis80) rc=" _rc
    global pf = $pf + 1
  }
  else {
    local ppow = r(Power)
    local pbay = r(Bayes)
    local plr  = r(LR)
    if !(`ppow' >= .78 & `ppow' <= .82) {
      dis as error "`case' C3 FAIL: package power at wooldid mdt_nis(80) = " %9.6f `ppow' " outside [.78, .82]"
      global pf = $pf + 1
    }
    else dis "`case' C3 info: package power at wooldid mdt_nis(80) = " %9.6f `ppow' " (in [.78, .82])"
    dis "`case' info: at the mdt_nis(80) slope the package reports Bayes factor = " %12.0g `pbay' ", likelihood ratio = " %12.0g `plr' " (record only)"
  }

  * ---- info: package power at wooldid's mdt_nis(50) (record only; no preregistered corridor)
  cap noi pretrends, numpre(`np') b(b) v(V) slope(`nis50')
  if _rc == 0 dis "`case' info: package power at wooldid mdt_nis(50) = " %9.6f r(Power) " (record only)"
  else dis "`case' info: pretrends slope(mdt_nis50) rc=" _rc " (record only)"
end

* ============================================================================
* the 4 designs (a failing wooldid run FAILs its design and the script goes on)
* ============================================================================

* ---- d1_bin_p4: dta_binary baseline, 4 pre / 4 post (the quadrature cap numpre=4)
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod pretrendpower(50 80)
if _rc != 0 {
  dis as error "d1_bin_p4 FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else ptcase, case(d1_bin_p4)

* ---- d2_bin_asym: asymmetric window, 2 pre / 4 post
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(4) esfixedbaseperiod pretrendpower(50 80)
if _rc != 0 {
  dis as error "d2_bin_asym FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else ptcase, case(d2_bin_asym)

* ---- d3_bin_p3: 3 pre-period variant, 3 pre / 3 post
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod pretrendpower(50 80)
if _rc != 0 {
  dis as error "d3_bin_p3 FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else ptcase, case(d3_bin_p3)

* ---- d4_adm_p3: realistic admin panel, cohort-level i, 48 state clusters, 3 pre / 3 post
use `"`SDATA'/dta_s_admin.dta"', clear
cap noi wooldid y cohid year ttre, att cluster(state) esprelength(3) espostlength(3) esfixedbaseperiod pretrendpower(50 80)
if _rc != 0 {
  dis as error "d4_adm_p3 FAIL: wooldid rc=" _rc
  global pf = $pf + 1
}
else ptcase, case(d4_adm_p3)

* ============================================================================
dis "pretrends_crosscheck: $ncase design dispatches completed (4 expected)"
if $ncase != 4 {
  dis as error "PRETRENDS_CROSSCHECK FAIL: only $ncase of 4 designs were dispatched"
  global pf = $pf + 1
}
if $pf == 0 dis "PRETRENDS_CROSSCHECK_ALL_PASS"
else dis as error "PRETRENDS_CROSSCHECK_FAILS=$pf"
exit, clear
