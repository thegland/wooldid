* ============================================================================
* wooldid testing harness template (known-good patterns from the B1-B23 campaign)
* Copy this into YOUR working directory and adapt. Do not run it as-is.
*
* Golden rules baked in below:
*  - one command per line (Stata does NOT treat ";" as a separator)
*  - cap noi + immediate sentinel for every command whose rc you care about
*  - deterministic data: set BOTH seed and sortseed
*  - pinned ado copies, never the repo working tree
*  - full e() dumps in %21x for exact comparison (compare with test/compare.py)
* ============================================================================
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345

* ---- load pinned program versions (extract BEFORE this script, in shell:) ----
*   git -C /Users/hegland/Documents/Research/wooldid_dev show <SHA_NEW>:wooldidinstall/wooldid.ado > wooldid_new.ado
*   git -C /Users/hegland/Documents/Research/wooldid_dev show <SHA_OLD>:wooldidinstall/wooldid.ado | sed 's/^cap prog drop wooldid$/cap prog drop wooldid_old/; s/^prog def wooldid, eclass$/prog def wooldid_old, eclass/' > wooldid_old.ado
do "wooldid_new.ado"
do "wooldid_old.ado"
* now `wooldid` = new build, `wooldid_old` = old build, in one session

* ---- dumpres: write EVERY e() scalar/macro/matrix cell in exact %21x hex ----
cap program drop dumpres
program define dumpres
  args cfg rc
  tempname fh
  file open `fh' using "`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
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
* compare two dump directories:  python3 <repo>/test/compare.py DIR_A DIR_B [--rtol 0]

* ---- freshstart: reset state between configs (stale e()/data leaks otherwise) ----
cap program drop freshstart
program define freshstart
  args dtafile
  ereturn clear
  use "`dtafile'", clear
end

* ============================ example test blocks ============================

* Pattern 1: rc sentinel (works in -e batch mode because cap prevents the abort)
freshstart "/Users/hegland/Documents/Research/wooldid_dev/test/data/dta_binary"
cap noi wooldid y i t ttre, att cluster(i)
dis "T1_RC=" _rc
dumpres t1_new `=_rc'

* Pattern 2: old-vs-new bit comparison on identical data
freshstart "/Users/hegland/Documents/Research/wooldid_dev/test/data/dta_binary"
cap noi wooldid_old y i t ttre, att cluster(i)
dumpres t1_old `=_rc'
* -> shell: python3 <repo>/test/compare.py . .   (put old/new dumps in separate dirs)

* Pattern 3: in-session matrix comparison with a guard against missing matrices
freshstart "/Users/hegland/Documents/Research/wooldid_dev/test/data/dta_binary"
cap noi wooldid y i t ttre, att cluster(i)
cap confirm matrix e(wdidmainresults)
if _rc == 0 {
  matrix NEWM = e(wdidmainresults)
}
freshstart "/Users/hegland/Documents/Research/wooldid_dev/test/data/dta_binary"
cap noi wooldid_old y i t ttre, att cluster(i)
matrix OLDM = e(wdidmainresults)
dis "T2_mreldif=" mreldif(NEWM, OLDM)
if mreldif(NEWM, OLDM) < 1e-14 dis "T2_PASS"
else dis "T2_FAIL"

* Pattern 4: expected-error check (error rc AND message both matter)
freshstart "/Users/hegland/Documents/Research/wooldid_dev/test/data/dta_binary"
cap noi wooldid y i t ttre, att unconditionalse
dis "T3_RC=" _rc
* grep the .log afterward for the expected message text

* Always end with an unambiguous completion marker (lets the caller detect a hung run):
dis "HARNESS_DONE"
exit, clear
