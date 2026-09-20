*==========================================================================
* REPRO NFC-01   (self-contained; generates its own data; uses no repo file
*                 except the pinned build)
*
* CLAIM
*   With an ordinary continuous covariate in controls(), the default cell
*   collapse and nofastcollapse disagree by up to ~5e-6 in relative terms on
*   the reported VARIANCE / standard error, and the disagreement reaches the
*   digits wooldid PRINTS.  wooldid.sthlp states the opposite twice:
*     sthlp:338    "Standard errors and variances do not amplify at all:
*                   they were measured at about 1e-16 relative regardless of
*                   the scale of the predictions."
*     sthlp:350    "For almost all work, not at all: no digit wooldid
*                   displays is affected"
*   The project's own indifference threshold is rtol 1e-10; this is 4-5
*   orders past it.  Point estimates are fine (~1e-15) - it is the delta-
*   method variance that moves, and the CI/p columns with it.
*
* MECHANISM.  What the four arms below ESTABLISH is the trigger.  The numerical-
*   analysis story underneath (last paragraph) is consistent with them but is not
*   itself measured by this file.
*   The disagreement appears ONLY when a control is (a) continuous and
*   (b) varies WITHIN the collapse's (cohort, period, targetflag) cell.
*   Make the same control constant within the cell -> clean.  Absorb it with
*   fe() so it has no coefficient -> clean.  Use a factor control -> clean.
*   The collapsed row carries ONE row's value of the control, and that value
*   reaches margins' delta-method Jacobian, whose entry for a control
*   coefficient is analytically zero but is formed as a difference of two
*   nearly-equal predictions and so survives as cancellation residue.
*
* RUN
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_NFC01_controls_V.do
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 8675309
set sortseed 8675309

*---- an ordinary staggered panel: 120 firms in 3 treated cohorts + never-treated,
*---- 10 years, i = the cohort id (wooldid's own guidance for a coarser i)
set obs 120
gen long unit   = _n
gen byte cohort = mod(_n-1,4)+1
gen double ttre = .
replace ttre = 2004 if cohort==1
replace ttre = 2006 if cohort==2
replace ttre = 2008 if cohort==3
gen double a_i = rnormal()
expand 10
bysort unit: gen int t = 2000 + _n - 1
* xvary varies WITHIN a (cohort,t) cell.   xcellconst does not.
gen double xvary = rnormal()
bysort cohort t (unit): gen double xcellconst = xvary[1]
* the same covariate on an administrative scale (e.g. a dollar amount)
gen double xdollar = 52000 + 9000*xvary
gen byte   region  = mod(unit,5) + 1
gen double y = a_i + 0.4*(t-2000) + 0.5*xvary + 1.2*(ttre<. & t>=ttre) + rnormal()*0.3
compress
save "_nfc01.dta", replace

cap prog drop nfcrun
prog def nfcrun
  syntax , TAG(string) [CTRL(string) SECOND]
  * arm A = collapse (default);  arm B = nofastcollapse, or the SAME arm again
  * when -second- is asked for (the noise floor)
  foreach arm in A B {
    use "_nfc01.dta", clear
    local nofc nofastcollapse
    if "`arm'" == "A" local nofc
    if "`second'" != "" local nofc
    qui wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests `ctrl' `nofc'
    matrix V`arm' = e(V)
    matrix b`arm' = e(b)
    scalar p`arm' = e(joint_horiz)
    scalar fc`arm' = ("`e(wdid_fcol)'" == "1")
  }
  mata: VA=st_matrix("VA"); VB=st_matrix("VB"); bA=st_matrix("bA"); bB=st_matrix("bB"); dv=abs(VA-VB):/(((abs(VA):>abs(VB)):*abs(VA)+(abs(VA):<=abs(VB)):*abs(VB)):+1e-300); db=abs(bA-bB):/(((abs(bA):>abs(bB)):*abs(bA)+(abs(bA):<=abs(bB)):*abs(bB)):+1e-300); st_numscalar("mxV", max(dv)); st_numscalar("mxb", max(db))
  local pp = abs(pA-pB)/max(abs(pA),abs(pB))
  dis as result "NFC01| " %-22s "`tag'" " |fcolA=" fcA " fcolB=" fcB ///
     " |max rel e(V)=" %9.2e mxV " |max rel e(b)=" %9.2e mxb " |rel joint p=" %9.2e `pp'
end

dis " "
dis as text "=== 0. NOISE FLOOR: the SAME arm twice, with the control.  Must be exactly 0. ==="
nfcrun , tag(noisefloor_controls) ctrl(controls(xvary)) second

dis " "
dis as text "=== 1. no controls -- the cell the project has already measured ==="
nfcrun , tag(no controls)

dis " "
dis as text "=== 2. ONE ordinary continuous control that varies within the cell ==="
nfcrun , tag(controls(xvary))          ctrl(controls(xvary))

dis " "
dis as text "=== 3. the same control on a dollar scale ==="
nfcrun , tag(controls(xdollar))        ctrl(controls(xdollar))

dis " "
dis as text "=== 4. MECHANISM: the same control made constant within the cell ==="
nfcrun , tag(controls(xcellconst))     ctrl(controls(xcellconst))

dis " "
dis as text "=== 5. MECHANISM: absorbed instead (no coefficient) ==="
nfcrun , tag(fe(region))               ctrl(fe(region))

dis " "
dis as text "=== 6. MECHANISM: a FACTOR control ==="
nfcrun , tag(controls(i.region))       ctrl(controls(i.region))

dis " "
dis as text "=============================================================="
dis as text " THE PRINTED TABLE.  Two continuous controls.  Diff the two blocks."
dis as text "=============================================================="
dis "<<<<< PRINTED_A  (collapse on, the default) >>>>>"
use "_nfc01.dta", clear
noisily wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xdollar xvary)
dis "<<<<< PRINTED_B  (nofastcollapse) >>>>>"
use "_nfc01.dta", clear
noisily wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xdollar xvary) nofastcollapse
dis "<<<<< PRINTED_END >>>>>"
erase "_nfc01.dta"
dis "REPRO_NFC01_DONE"
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"
local OO "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/out"
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"
local BASE wooldid y cohort t ttre
local ES att esprelength(3) espostlength(3) jointtests

*==================== SECTION 2: p12 framing + dropped control sets ====================
local BASE wooldid y cohort t ttre
local ES att esprelength(3) espostlength(3) jointtests

* ===== THE FRAMING CELL =====
* Is the V movement specific to the COLLAPSE, or is margins' delta-method V simply
* row-order sensitive at 1e-6 once a continuous control is present?
* nofastcollapse on the ORIGINAL row order vs a shuffled copy of the same data.
use "`DD'/dCc.dta", clear
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "`DD'/dCc_shuf.dta", replace

qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/y01_orig_nofc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/y01_shuf_nofc.txt")
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/y01_orig_fc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/y01_shuf_fc.txt")

* ===== the two control sets that were dropped on my syntax error =====
wdcell , dta("`DD'/dCc.dta") tag(y02_ccc)     base(`BASE') out("`OO'") opts(`ES' coarsecohortcontrols(c.xdollar))
wdcell , dta("`DD'/dCc.dta") tag(y03_cccabs)  base(`BASE') out("`OO'") opts(`ES' ccc_absorb(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(y04_ccci)    base(`BASE') out("`OO'") opts(`ES' coarsecohortcontrols(i.region))
wdcell , dta("`DD'/dCc.dta") tag(y05_inter)   base(`BASE') out("`OO'") opts(`ES' interactivecontrols(xdollar))
dis "P12_DONE"

*==================== SECTION 3: p04 state probes ====================
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"

*==================== A. retired option names ====================
use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att fastcollapse
dis "A1_fastcollapse_RC=" _rc
cap noisily wooldid y cohort t ttre, att fullprecisionslow
dis "A2_fullprecisionslow_RC=" _rc
cap noisily wooldid y cohort t ttre, att nofastcollapse
dis "A3_nofastcollapse_RC=" _rc

*==================== B. hidden results contract ====================
use "`DD'/dC.dta", clear
qui wooldid y cohort t ttre, att esprelength(2) espostlength(2)
dis "B_fcol=" e(wdid_fcol) " fcncell=" e(wdid_fcncell) " fcskip=<" e(wdid_fcskip) "> fcguard=<" e(wdid_fcguard) ">"
qui wooldid y cohort t ttre, att esprelength(2) espostlength(2) nofastcollapse
dis "B_NOFC_fcol=" e(wdid_fcol) " fcncell=" e(wdid_fcncell) " fcskip=<" e(wdid_fcskip) "> fcguard=<" e(wdid_fcguard) ">"
* ineligible skip should CAUTION and set fcskip
qui wooldid y cohort t ttre, att esprelength(2) espostlength(2) unconditionalse
dis "B_UNC_fcol=" e(wdid_fcol) " fcskip=<" e(wdid_fcskip) ">"

*==================== C. data state: signature / sortorder / N / k ====================
foreach arm in "" "nofastcollapse" {
  use "`DD'/dC.dta", clear
  qui datasignature set, reset
  local s0 = r(datasignature)
  local n0 = _N
  local k0 = c(k)
  local sb0 : sortedby
  qui wooldid y cohort t ttre, att att_i att_it att_itime esprelength(3) espostlength(3) jointtests `arm'
  qui datasignature
  local s1 = r(datasignature)
  local sb1 : sortedby
  dis "C|arm=<`arm'>|sigsame=" ("`s0'"=="`s1'") "|N=`n0'->" _N "|k=`k0'->" c(k) "|sortedby=<`sb0'>-><`sb1'>"
}

*==================== D. user's own preserve around wooldid ====================
use "`DD'/dC.dta", clear
local n0 = _N
preserve
qui drop if t > 2005
qui count
dis "D_inner_N=" r(N)
cap noisily wooldid y cohort t ttre, att esprelength(2) espostlength(2)
dis "D1_RC=" _rc
cap noisily restore
dis "D1_restore_RC=" _rc
dis "D1_after_N=" _N " expected=`n0'"

*==================== E. error raised INSIDE the collapse window ====================
* graph save to a nonexistent directory, from wooldid_hist, which runs inside the window
use "`DD'/dC.dta", clear
qui datasignature set, reset
local s0 = r(datasignature)
local n0 = _N
cap noisily wooldid y cohort t ttre, att att_i esprelength(2) espostlength(2) histogramcohorteffects saveplots(/no/such/dir/zz)
dis "E1_RC=" _rc
qui datasignature
dis "E1_sigsame=" ("`s0'"=="`=r(datasignature)'") " N=" _N " expected=`n0'" " k=" c(k)
qui describe, short
cap confirm variable y
dis "E1_y_present_RC=" _rc

*==================== F. same, but the user had preserved first ====================
use "`DD'/dC.dta", clear
local n0 = _N
preserve
qui keep if cohort != 2
cap noisily wooldid y cohort t ttre, att esprelength(2) espostlength(2) histogramcohorteffects saveplots(/no/such/dir/zz)
dis "F1_RC=" _rc
dis "F1_N_now=" _N
cap noisily restore
dis "F1_restore_RC=" _rc
dis "F1_after_N=" _N " expected=`n0'"

*==================== G. wooldid twice in a row after an error ====================
use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att histogramcohorteffects saveplots(/no/such/dir/zz)
cap noisily wooldid y cohort t ttre, att esprelength(2) espostlength(2)
dis "G1_RC=" _rc
cap dis "G1_b_att=" %21x el(e(b),1,1)
dis "G1_disp_rc=" _rc
dis "P04_DONE"

*==================== SECTION 4: p10 rows-per-cell ladder ====================
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"

* ---- how does the controls() disagreement scale with ROWS PER CELL? ----
* same design, same DGP, only the number of units per cohort changes, so the
* number of rows the collapse folds into one grows 4 -> 8 -> 32 -> 128 -> 512.
foreach u in 4 8 32 128 512 {
  clear
  set seed 20260917
  set sortseed 20260917
  qui set obs `=4*`u''
  qui gen long unit   = _n
  qui gen byte cohort = mod(_n-1,4)+1
  qui gen double ttre = .
  qui replace ttre = 2004 if cohort==1
  qui replace ttre = 2006 if cohort==2
  qui replace ttre = 2008 if cohort==3
  qui gen double a_i = rnormal()
  qui expand 10
  qui bysort unit: gen int t = 2000 + _n - 1
  qui gen double x1 = rnormal()
  qui gen double y = a_i + 0.4*(t-2000) + 0.5*x1 + 1.2*(ttre<. & t>=ttre) + rnormal()*0.3
  qui save "`DD'/cs`u'.dta", replace
  wdcell , dta("`DD'/cs`u'.dta") tag(cs`u') base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests controls(x1))
  wdcell , dta("`DD'/cs`u'.dta") tag(cs`u'_noc) base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests)
}
dis "P10_DONE"
dis "PALL_DONE"
