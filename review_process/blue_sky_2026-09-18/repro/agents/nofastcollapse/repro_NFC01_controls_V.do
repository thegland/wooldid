*==========================================================================
* REPRO NFC-01
*   RUN:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_NFC01_controls_V.do
*   Reads only the pinned build and the blue-sky recon data generator
*   (scratchpad/bluesky/recon/gen/gen_core.do, which the brief offers for reuse).
*   No repo test-suite file and no repo .dta is touched.
*
* CLAIM  (all four numbers below come out of this file)
*   With an ordinary continuous covariate in controls(), the default cell collapse
*   and nofastcollapse disagree by ~1e-6 RELATIVE on the reported variance, and the
*   disagreement reaches the digits wooldid PRINTS.  The project's indifference
*   threshold is rtol 1e-10.  wooldid.sthlp says the opposite four times:
*     sthlp:332  bounded by about ONE ULP of the linear prediction, "every difference
*                measured was a small integer multiple of that one quantum"
*     sthlp:338  "Standard errors and variances do not amplify at all: they were
*                measured at about 1e-16 relative"
*     sthlp:350  "no digit wooldid displays is affected"
*     sthlp:357  "merely permuting the order of the rows ... moves the same near-zero
*                coefficient about 19 times MORE than the collapse does"
*   NOTE sthlp:357 is CORRECT as written: on the near-zero COEFFICIENT the measured
*   ratio is 15.5.  Coefficients are fine everywhere (1e-16..1e-15).  It is the
*   VARIANCE that moves, by 1e-7 to 5e-6 depending on the cell, and arms F and G
*   below show the second-order consequence: with the collapse OFF a full row
*   permutation moves e(V) by ~7e-16, with it ON by ~9e-8 - the collapse makes the
*   reported standard error depend on the order your observations are stored in.
*
* THE ARMS
*   N  noise floor, the same arm twice, control present     -> must be exactly 0
*   A  no controls                                          -> clean  (~1e-14)
*   C  one continuous control varying within the cell       -> FAILS (1e-7..5e-6,
*                                                                cell-dependent)
*   D  the same control made CONSTANT within the cell       -> clean
*   E  the control absorbed by fe() so it has no coefficient-> clean
*   F  a plain ROW PERMUTATION with the collapse OFF        -> clean  (~1e-15 on V)
*   G  the same row permutation with the collapse ON        -> FAILS  (~9e-8 on V)
*   Arms D and E isolate the trigger: a CONTINUOUS control that carries a coefficient
*   and varies WITHIN the collapse's own (cohort, relative-period, targetflag) cell.
*   Arm G is the consequence: under the collapse the reported variance depends on the
*   order the user's observations happen to be stored in - by 8 orders of magnitude
*   more than it does with the collapse off - because the collapse keeps the cell's
*   FIRST row, so re-sorting changes which control value survives.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*---- the panel: 150 units, 3 treated cohorts + never-treated, 10 periods, i = cohort
local gN 150
local gT 10
local gonsets 2004 2006 2008
local gnever 0.35
local gysd 0.30
local gdose unitconst
local gwtype unit
local gwspread 2
local gpois 1
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen/gen_core.do"
gen double xpct = 50 + 8*x1
bysort cohort t (unit): gen double xcellconst = xpct[1]
label var xpct "an ordinary percentage-scale covariate"
save "_nfc01.dta", replace
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "_nfc01_shuffled.dta", replace

local ES att att_i att_it esprelength(3) espostlength(3) jointtests
tempname VA VB
dis _n as text "{hline 86}"
dis as text "arm                                              fcol A/B | max rel e(V) | max rel e(b)"
dis as text "{hline 86}"
foreach arm in N A C D E F G {
  local ctrl controls(xpct)
  local dtaA "_nfc01.dta"
  local dtaB "_nfc01.dta"
  local nofcA
  local nofcB nofastcollapse
  local lbl
  if "`arm'" == "N" {
    local lbl "N  noise floor: collapse arm run TWICE"
    local nofcB
  }
  if "`arm'" == "A" {
    local lbl "A  no controls"
    local ctrl
  }
  if "`arm'" == "C" local lbl "C  controls(xpct)  <- varies within the cell"
  if "`arm'" == "D" {
    local lbl "D  controls(xcellconst)  <- constant in the cell"
    local ctrl controls(xcellconst)
  }
  if "`arm'" == "E" {
    local lbl "E  fe(region)  <- absorbed, no coefficient"
    local ctrl fe(region)
  }
  if "`arm'" == "F" {
    local lbl "F  ROW SHUFFLE, collapse OFF in both"
    local nofcA nofastcollapse
    local dtaB "_nfc01_shuffled.dta"
  }
  if "`arm'" == "G" {
    local lbl "G  ROW SHUFFLE, collapse ON in both"
    local nofcB
    local dtaB "_nfc01_shuffled.dta"
  }
  use "`dtaA'", clear
  cap noi qui wooldid y cohort t ttre, `ES' `ctrl' `nofcA'
  local rcA = _rc
  if `rcA' == 0 {
    matrix VVA = e(V)
    matrix BBA = e(b)
    local fA `e(wdid_fcol)'
  }
  use "`dtaB'", clear
  cap noi qui wooldid y cohort t ttre, `ES' `ctrl' `nofcB'
  local rcB = _rc
  if `rcB' == 0 {
    matrix VVB = e(V)
    matrix BBB = e(b)
    local fB `e(wdid_fcol)'
  }
  if `rcA' != 0 | `rcB' != 0 {
    dis as error "`lbl' : rcA=`rcA' rcB=`rcB'"
    continue
  }
  mata: VA=st_matrix("VVA"); VB=st_matrix("VVB"); bA=st_matrix("BBA"); bB=st_matrix("BBB"); sV=abs(VA); sV=sV:*(sV:>=abs(VB))+abs(VB):*(sV:<abs(VB)); st_numscalar("mxV", max(abs(VA-VB):/(sV:+1e-300))); sB=abs(bA); sB=sB:*(sB:>=abs(bB))+abs(bB):*(sB:<abs(bB)); st_numscalar("mxb", max(abs(bA-bB):/(sB:+1e-300)))
  dis as result %-48s "`lbl'" "   `fA'/`fB'    | " %9.2e mxV "    | " %9.2e mxb
}
dis as text "{hline 86}"
dis as text "PASS = at or below 1e-10.  Arm N must be exactly 0 or the harness, not wooldid, is at fault."

dis _n as text "THE PRINTED TABLE - diff the two blocks.  Two continuous controls."
dis "<<<<< PRINTED_A  (collapse on, the default) >>>>>"
use "_nfc01.dta", clear
cap noi wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xpct x2)
dis "<<<<< PRINTED_B  (nofastcollapse) >>>>>"
use "_nfc01.dta", clear
cap noi wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xpct x2) nofastcollapse
dis "<<<<< PRINTED_END >>>>>"
erase "_nfc01.dta"
erase "_nfc01_shuffled.dta"
dis "REPRO_NFC01_DONE"
