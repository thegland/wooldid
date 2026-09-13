* REPRO A12-ttb - tcen f7019c9: v1.3 timetrends numerically INERT at t ~ 5e5
* Suite config `trig_tt_big' (dta_s_ttbig: t = 500000..500012, onsets 500006/500009).
* One build per Stata process.
* Legs: A = the suite config (timetrends, raw t); B = same without timetrends;
*       C = timetrends on the recoded panel t-499999 (t = 1..13, onsets 7/10), the
*           small-code benchmark every build fits correctly.
* Old-build defect signature: A == B (rc 0, trends inert), A != C.
* New-build fix signature:    A == C (centered trend is shift-invariant), A != B.
*
* ARGUMENTS (positional, all three required):
*   `1' build  - "old" or "new"; selects the ado in <workdir>/old or <workdir>/new
*   `2' repo   - the wooldid_dev checkout root (supplies test/suite/data/*.dta)
*   `3' workdir- scratch directory holding old/, new/; logs and dumps_<build>/ are written here
* EXACT INVOCATION (from any cwd):
*   stata-mp -e do repro_A12TTB.do old /path/to/wooldid_dev /path/to/scratch/legacy
*   stata-mp -e do repro_A12TTB.do new /path/to/wooldid_dev /path/to/scratch/legacy
args build repo WD_
clear all
set more off
set maxvar 20000
set graphics off
local WD "`WD_'"
if "`build'" == "" | "`repo'" == "" | "`WD'" == "" {
  di as error "usage: do repro_A12TTB.do <old|new> <repo-root> <workdir>"
  exit 198
}
cap log close
log using "`WD'/repro_A12TTB_`build'.log", replace text
adopath ++ "`WD'/`build'"
cap mkdir "`WD'/dumps_`build'"

cap program drop dumpres
program define dumpres
  args cfg rc dir
  tempname fh
  file open `fh' using "`dir'/`cfg'.csv", write replace
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

cap program drop showfit
program define showfit
  args build leg
  tempname b
  cap matrix `b' = e(b)
  if _rc != 0 {
    dis "REPRO_A12TTB_`build' B_`leg'=NOMATRIX"
    exit
  }
  local v : di %21x `b'[1,1]
  dis "REPRO_A12TTB_`build' B_`leg'=`v'"
  local v : di %21x e(rss)
  dis "REPRO_A12TTB_`build' RSS_`leg'=`v'"
  local v : di %21x e(r2)
  dis "REPRO_A12TTB_`build' R2_`leg'=`v'"
  dis "REPRO_A12TTB_`build' PROPS_`leg'=`e(properties)'"
  dis "REPRO_A12TTB_`build' dec `leg' att=" %14.0g `b'[1,1] " rss=" %14.0g e(rss) " r2=" %10.0g e(r2)
end

* ---- leg A: trig_tt_big exactly as mkrun.py emits it ----
set seed 98765
set sortseed 98765
ereturn clear
use "`repo'/test/suite/data/dta_s_ttbig.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
local rc = _rc
dis "REPRO_A12TTB_`build'_A_RC=`rc'"
dumpres trig_tt_big `rc' "`WD'/dumps_`build'"
showfit `build' TT_RAW

* ---- leg B: same data, no timetrends ----
set seed 98765
set sortseed 98765
ereturn clear
use "`repo'/test/suite/data/dta_s_ttbig.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dis "REPRO_A12TTB_`build'_B_RC=`rc'"
showfit `build' NOTT_RAW

* ---- leg C: recoded small codes (t 1..13, onsets 7/10), timetrends ----
set seed 98765
set sortseed 98765
ereturn clear
use "`repo'/test/suite/data/dta_s_ttbig.dta", clear
qui replace t = t - 499999
qui replace ttre = ttre - 499999
qui sum t
dis "REPRO_A12TTB_`build' recoded t range " r(min) " .. " r(max)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends
local rc = _rc
dis "REPRO_A12TTB_`build'_C_RC=`rc'"
showfit `build' TT_RECODED
dis "REPRO_A12TTB_DONE"
log close
exit, clear
