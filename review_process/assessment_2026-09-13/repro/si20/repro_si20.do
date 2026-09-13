* repro_si20.do - SI-20 reproduction (stage-B proof agent, 2026-09-12)
* Usage (one build per Stata process, from its own cwd so the logs do not collide):
*   cd <workdir>/run_old && "$STATA" -e do <repro>/si20/repro_si20.do old <repo> <workdir>
*   cd <workdir>/run_new && "$STATA" -e do <repro>/si20/repro_si20.do new <repo> <workdir>
* ARGUMENTS (positional, all three required):
*   `1' build   - "old" or "new"; sources <workdir>/<build>_i/wooldid.ado (the INSTRUMENTED copy)
*   `2' repo    - the wooldid_dev checkout root (supplies test/suite/data/*.dta)
*   `3' workdir - scratch directory holding old_i/ and new_i/ (see instrument.py in repro/si20/)
* Builds: old = v1.3 (main 301b706), new = HEAD 08ec0b4; each an INSTRUMENTED COPY (old_i/, new_i/)
* carrying only the read-only REPRO_SI20 probe (see instrument.py): the joint-test margins call
* re-run without contrast(overjoint), so the reference-cell (_t == 0) level margin is visible.
* Sections: fix 1 on si20_lat_pw and on si20_lat_aw_bt (the ticket's fixture of record);
*           fix 2 on si20_ic_coh (the flagship). Per-config scrub block mirrors test/suite/mkrun.py.
args build repo W_
clear all
set more off
set maxvar 20000
set graphics off
local W "`W_'"
if "`build'" == "" | "`repo'" == "" | "`W'" == "" {
  di as error "usage: do repro_si20.do <old|new> <repo-root> <workdir>"
  exit 198
}
local D "`repo'/test/suite/data"
do "`W'/`build'_i/wooldid.ado"
dis "BUILD|`build'"

* ---- dumpres: every e() scalar/macro/matrix cell in %21x (harness_template.do pattern) ----
cap program drop dumpres
program define dumpres
  args cfg rc
  tempname fh
  file open `fh' using "`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    if _rc != 0 continue
    local R = rowsof(`M')
    local C = colsof(`M')
    local rn : rowfullnames `M'
    local cn : colfullnames `M'
    file write `fh' `"matnames|`m'|`R'x`C'|`rn'|`cn'"' _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  file close `fh'
end

* ---- printres: the decisive numbers as short tagged lines (batch logs wrap at 80 cols) ----
cap program drop printres
program define printres
  args cfg rc
  dis "RC|`cfg'|`rc'"
  foreach tag in F1 F2 {
    cap confirm matrix REPRO_ANCH_`tag'
    if _rc == 0 {
      local cn : colnames REPRO_ANCH_`tag'
      local eq : coleq REPRO_ANCH_`tag'
      local k 0
      foreach c of local cn {
        local ++k
        local e : word `k' of `eq'
        dis "ANCH|`cfg'|`tag'|`e':`c'|" strtrim(string(el(REPRO_ANCH_`tag',1,`k'), "%21.15g"))
        dis "ANCHV|`cfg'|`tag'|`e':`c'|" strtrim(string(el(REPRO_ANCHV_`tag',`k',`k'), "%21.15g"))
      }
    }
    else dis "ANCH|`cfg'|`tag'|ABSENT"
  }
  if `rc' != 0 exit
  foreach s in N df_r F joint_horiz joint_Choriz ll rss {
    dis "SCAL|`cfg'|`s'|" strtrim(string(e(`s'), "%21.15g"))
  }
  tempname B M
  matrix `B' = e(b)
  local cn : colnames `B'
  local k 0
  foreach c of local cn {
    local ++k
    dis "B|`cfg'|`c'|" strtrim(string(el(`B',1,`k'), "%21.15g"))
  }
  cap matrix `M' = e(wdidmainresults)
  if _rc == 0 {
    local rn : rowfullnames `M'
    local r 0
    foreach nm of local rn {
      local ++r
      dis "W|`cfg'|`nm'|se|" strtrim(string(el(`M',`r',2), "%21.15g"))
      dis "W|`cfg'|`nm'|p|" strtrim(string(el(`M',`r',4), "%21.15g"))
    }
  }
end

* ---- scrub: the suite's per-config state reset (test/suite/mkrun.py block()) ----
cap program drop scrub
program define scrub
  set seed 98765
  set sortseed 98765
  ereturn clear
  cap macro drop WOOLDID_LASTRC WOOLDID_CLEANUP
  cap matrix drop _all
  cap scalar drop _all
  cap estimates drop _all
end

* ======== section 1a: fix 1 on si20_lat_pw (lattice(all 3) + basic, pweights, jointtests) ========
scrub
use "`D'/continuous_control_dense.dta", clear
cap noi wooldid y_cont_common cohort t ttre [pw=probability_weight], att contreat(dose_common) contreatcontrols(basic) lattice(all 3) esprelength(2) espostlength(2) jointtests safety(off)
local rc = _rc
dumpres si20_lat_pw_`build' `rc'
printres si20_lat_pw `rc'
dis "REPRO_SI20F1_si20_lat_pw_DONE"

* ======== section 1b: fix 1 on si20_lat_aw_bt (the ticket's fixture of record: anchor 1.0198262516390375) ========
scrub
use "`D'/continuous_control_dense.dta", clear
cap noi wooldid y_cont_common cohort t ttre [aw=analytic_weight], att contreat(dose_common) contreatcontrols(basic) lattice(baseTreated 2) esprelength(2) espostlength(2) robust jointtests safety(off)
local rc = _rc
dumpres si20_lat_aw_bt_`build' `rc'
printres si20_lat_aw_bt `rc'
dis "REPRO_SI20F1_si20_lat_aw_bt_DONE"

* ======== section 1c: fix 1 on the ticket's UNWEIGHTED fixture of record (Findings block, Mechanism 2: anchor 1.0198262516390375, V 1.784e-3) ========
scrub
use "`D'/continuous_control_dense.dta", clear
cap noi wooldid y_cont_common cohort t ttre, att contreat(dose_common) contreatcontrols(basic) lattice(baseTreated 2) esprelength(2) espostlength(2) robust jointtests safety(off)
local rc = _rc
dumpres ticket_bt2_`build' `rc'
printres ticket_bt2 `rc'
dis "REPRO_SI20F1_ticket_bt2_DONE"

* ======== section 2: fix 2 on si20_ic_coh (cohort-level i via prep covars_coh + fe(i) + time-varying x1) ========
scrub
use "`D'/dta_s_covars.dta", clear
qui gen coh = cond(ttre==5,1,cond(ttre==9,2,3))
cap noi wooldid y coh t ttre, att cluster(i) fe(i) esprelength(2) espostlength(3) interactivecontrols(x1) jointtests
local rc = _rc
dumpres si20_ic_coh_`build' `rc'
printres si20_ic_coh `rc'
dis "REPRO_SI20F2_si20_ic_coh_DONE"

dis "REPRO_SI20_DONE"
exit, clear
