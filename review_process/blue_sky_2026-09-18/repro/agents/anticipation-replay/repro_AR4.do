*==============================================================================
* repro_AR4.do -- SELF-CONTAINED.
* Claim: wooldid refuses 1152 of the 9000 legal two-decimal confidence levels,
* including the upper bound 99.99 that sthlp:540 names and the ordinary 99.9 --
* AND it refuses them even when the user never types level(), because a session
* `set level 99.9' is resolved into the same code path (wooldid_validate
* ado:5507 `if `level' == 0 local level = c(level)', then wooldid_worker
* ado:574-576 `set level `level'').
* The error message, "level() can have at most two digits after the decimal
* point", is false about the value typed (99.9 has ONE decimal) and names an
* option the user may never have used.
* Mechanism: `syntax , Level(real 0)' renders the double into the macro with 16
* significant digits, so 99.9 arrives at `set level' as 99.90000000000001.
* The refused set is exactly {L : %.16g(L) != %.15g(L)} and starts at 64.01.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

quietly {
    set obs 60
    gen long i = _n
    gen byte cohort = mod(_n,3)
    gen double ttre = .
    replace ttre = 2004 if cohort==1
    replace ttre = 2006 if cohort==2
    gen double a_i = rnormal()
    expand 10
    bysort i: gen int t = 1999 + _n
    gen double y = a_i + 0.1*(t-2000) + 0.8*(ttre<. & t>=ttre) + 0.3*rnormal()
}

dis as result "==== A. Stata itself accepts every one of these ===="
foreach L in 99.9 99.99 98.9 79.9 99.98 92.5 {
  cap noi set level `L'
  dis "  set level `L'  rc=" _rc "  c(level)=" c(level)
  set level 95
}

dis as result "==== B. wooldid with level() refuses some of them ===="
foreach L in 99.9 99.99 98.9 79.9 99.98 92.5 {
  cap noi wooldid y i t ttre, att robust level(`L')
  dis "  wooldid level(`L')  rc=" _rc
}

dis as result "==== C. and refuses them with NO level() option at all ===="
foreach L in 99.9 99.99 99.98 {
  set level `L'
  cap noi wooldid y i t ttre, att robust
  dis "  session set level `L', plain wooldid  rc=" _rc "  c(level) after=" c(level)
  dis "     leftover globals: LEVEL0=[$WOOLDID_LEVEL0] LASTRC=[$WOOLDID_LASTRC] CLEANUP=[$WOOLDID_CLEANUP]"
  cap est dir
  set level 95
}
dis "REPRO_AR4_DONE"
