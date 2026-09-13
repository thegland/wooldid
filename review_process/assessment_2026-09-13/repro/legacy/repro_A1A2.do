* REPRO A1/A2 - B1 (control rows = pretreat rows) + B2 (sub-sample IQR = full IQR)
* Suite config `basic' (test/suite/configs.csv), one build per Stata process:
* Builds: old = v1.3 301b706 (<workdir>/old/wooldid.ado); new = 08ec0b4 (<workdir>/new/wooldid.ado).
* Prints e(summarystats) cells and the hand truth computed on the SAME e(sample):
* the user ttre (2005/2008) never exceeds max t (2011), so user ttre == internal ttre.
*
* ARGUMENTS (positional, all three required):
*   `1' build  - "old" or "new"; selects the ado in <workdir>/old or <workdir>/new
*   `2' repo   - the wooldid_dev checkout root (supplies test/data/*.dta)
*   `3' workdir- scratch directory holding old/, new/; logs and dumps_<build>/ are written here
* EXACT INVOCATION (from any cwd; see repro/README.md for the two build extractions):
*   stata-mp -e do repro_A1A2.do old /path/to/wooldid_dev /path/to/scratch/legacy
*   stata-mp -e do repro_A1A2.do new /path/to/wooldid_dev /path/to/scratch/legacy
args build repo WD_
clear all
set more off
set maxvar 20000
set graphics off
local WD "`WD_'"
if "`build'" == "" | "`repo'" == "" | "`WD'" == "" {
  di as error "usage: do repro_A1A2.do <old|new> <repo-root> <workdir>"
  exit 198
}
cap log close
log using "`WD'/repro_A1A2_`build'.log", replace text
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

* ---- config `basic', exactly as mkrun.py emits it ----
set seed 98765
set sortseed 98765
ereturn clear
use "`repo'/test/data/dta_binary.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats jointtests safety(off)
local rc = _rc
dis "REPRO_A1A2_`build'_RC=`rc'"
dumpres basic `rc' "`WD'/dumps_`build'"

tempname S
matrix `S' = e(summarystats)
local rn : rownames `S'
dis "REPRO_A1A2_`build' rows=`rn'"
local v : di %21x `S'[2,1]
dis "REPRO_A1A2_`build' SS_pretreat_mean=`v'"
local v : di %21x `S'[4,1]
dis "REPRO_A1A2_`build' SS_control_mean=`v'"
local v : di %21x `S'[1,4]
dis "REPRO_A1A2_`build' SS_insample_iqr=`v'"
local v : di %21x `S'[2,4]
dis "REPRO_A1A2_`build' SS_pretreat_iqr=`v'"
local v : di %21x `S'[3,4]
dis "REPRO_A1A2_`build' SS_posttreat_iqr=`v'"
local v : di %21x `S'[4,4]
dis "REPRO_A1A2_`build' SS_control_iqr=`v'"
local v : di %21x e(ypre_iqr)
dis "REPRO_A1A2_`build' E_ypre_iqr=`v'"
dis "REPRO_A1A2_`build' dec pretreat_mean=" %12.0g `S'[2,1] " control_mean=" %12.0g `S'[4,1]
dis "REPRO_A1A2_`build' dec insample_iqr=" %12.0g `S'[1,4] " pre_iqr=" %12.0g `S'[2,4]
dis "REPRO_A1A2_`build' dec post_iqr=" %12.0g `S'[3,4] " control_iqr=" %12.0g `S'[4,4]

* ---- hand truth on the realized e(sample) ----
qui gstats sum y if e(sample) & ttre == ., detail
local v : di %21x r(mean)
dis "REPRO_A1A2_`build' TRUTH_control_mean=`v'"
local v : di %21x (r(p75) - r(p25))
dis "REPRO_A1A2_`build' TRUTH_control_iqr=`v'"
dis "REPRO_A1A2_`build' dec TRUTH control mean=" %12.0g r(mean) " N=" r(N)
qui gstats sum y if e(sample) & ttre < . & t < ttre, detail
local v : di %21x (r(p75) - r(p25))
dis "REPRO_A1A2_`build' TRUTH_pre_iqr=`v'"
dis "REPRO_A1A2_`build' dec TRUTH pre iqr=" %12.0g (r(p75) - r(p25)) " mean=" %12.0g r(mean)
qui gstats sum y if e(sample) & ttre < . & t >= ttre, detail
local v : di %21x (r(p75) - r(p25))
dis "REPRO_A1A2_`build' TRUTH_post_iqr=`v'"
dis "REPRO_A1A2_`build' dec TRUTH post iqr=" %12.0g (r(p75) - r(p25))
qui gstats sum y if e(sample), detail
local v : di %21x (r(p75) - r(p25))
dis "REPRO_A1A2_`build' TRUTH_insample_iqr=`v'"
dis "REPRO_A1A2_`build' dec TRUTH insample iqr=" %12.0g (r(p75) - r(p25))
dis "REPRO_A1A2_DONE"
log close
exit, clear
