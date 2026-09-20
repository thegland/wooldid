*! SELF-CONTAINED REPRO -- state-leak agent, finding 2
*!
*! `verbose' adds a variable to the USER'S DATASET and leaves it there.
*! The help documents that verbose retains the stored model ___wooldidfullmodel
*! "for inspection" (sthlp:712 for the poisson twin); it does not say that Stata's
*! estimates-store bookkeeping variable `_est____wooldidfullmodel' is left in the
*! data with it.  A user who saves, exports, or loops `foreach v of varlist _all'
*! after a verbose run gets a column they did not create.
*!
*! Run:  cd <this directory>
*!       "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_verbose_leaves_variable.do
*!       then read repro_verbose_leaves_variable.log
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop mkpanel
prog def mkpanel
  qui {
    clear
    set obs 60
    gen int unit  = _n
    gen int state = mod(_n-1, 12) + 1
    gen int ttre = .
    replace ttre = 2005 if mod(_n,4) == 1
    replace ttre = 2007 if mod(_n,4) == 2
    replace ttre = 2009 if mod(_n,4) == 3
    gen double ufe = rnormal()
    expand 12
    bys unit: gen int year = 2000 + _n
    gen byte post = (year >= ttre) & !missing(ttre)
    gen double y = 1 + ufe + 0.03*(year-2000) + 0.5*post + rnormal()
  }
end

* ---------------- control: a plain run adds nothing
mkpanel
qui ds
local before `r(varlist)'
dis as text "k before: " c(k)
capture noisily wooldid y unit year ttre, cluster(state) att
dis "PLAIN_RC=" _rc
dis as text "k after plain run: " c(k)
qui ds
dis as text "varlist after plain run: `r(varlist)'"
capture confirm variable _est____wooldidfullmodel
dis "PLAIN_LEFT_A_VARIABLE=" (_rc == 0)

* ---------------- the case: verbose
mkpanel
dis as text "k before: " c(k)
capture noisily wooldid y unit year ttre, cluster(state) att verbose
dis "VERBOSE_RC=" _rc
dis as text "k after verbose run: " c(k)
qui ds
dis as text "varlist after verbose run: `r(varlist)'"
capture confirm variable _est____wooldidfullmodel
dis "VERBOSE_LEFT_A_VARIABLE=" (_rc == 0)
if _rc == 0 {
  local lb : variable label _est____wooldidfullmodel
  local ty : type _est____wooldidfullmodel
  dis as text "left-behind variable: type `ty', label: `lb'"
  qui sum _est____wooldidfullmodel, meanonly
  dis as text "  N nonmissing = " r(N) ", mean = " r(mean)
}
* what a user would then ship
qui describe, varlist
dis as text "describe varlist: `r(varlist)'"

* ==================================================================
* MEASURED (2026-09-17, pinned build c1d8979), via the full session snapshot in
* p16_snapdiff.do arm s3 -- `diff snap_s3_pre.txt snap_s3_post.txt':
*   [C] k          12 -> 13
*   [D] nvars      12 -> 13
*   + [VAR] _est____wooldidfullmodel | byte | %8.0g | vl= |
*           lbl=esample() from estimates store
*   [D] datasignature 720:12(95947):... -> 720:13(90315):...
* and on the plain (non-verbose) arm s1 the same diff shows NO variable added,
* because ado:4539 (`if "`verbose'" != "verbose" ... qui est drop `wdfm'') removes
* the stored estimate, and dropping it removes the marker variable with it.
dis "REPRO_DONE"
