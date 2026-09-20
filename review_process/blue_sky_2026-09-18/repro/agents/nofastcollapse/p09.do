*==========================================================================
* REPRO NFC-01      (self-contained: generates its own data; touches no repo
*                    file except the pinned build; writes one temp .dta in
*                    the current directory and erases it)
*
* RUN:
*   cd <this directory>
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_NFC01_controls_V.do
*   then read repro_NFC01_controls_V.log
*
*--------------------------------------------------------------------------
* CLAIM
*   With an ordinary continuous covariate in controls(), the default cell
*   collapse and nofastcollapse disagree by up to ~5e-6 in RELATIVE terms on
*   the reported variance / standard error, and the disagreement reaches the
*   digits wooldid PRINTS.  The point estimates are fine (~1e-15).
*   The project's indifference threshold is rtol 1e-10.
*
* wooldid.sthlp says the opposite, four times (line numbers on the pinned build):
*   sthlp:332  "The difference is bounded, in absolute terms, by about one unit
*               in the last place of the linear prediction ... and every
*               difference measured was a small integer multiple of that one
*               quantum.  On a panel whose predictions are of order 10, that is
*               about 2e-15."
*   sthlp:338  "Standard errors and variances do not amplify at all: they were
*               measured at about 1e-16 relative regardless of the scale of the
*               predictions."
*   sthlp:344  (continuous treatment) "no cell moved by more than 1e-9"
*   sthlp:350  "For almost all work, not at all: no digit wooldid displays is
*               affected"
*
* WHAT THE ARMS BELOW SHOW
*   0. noise floor: the SAME arm run twice, with the control present -> 0
*   1. no controls                                       -> ~1e-14   (clean)
*   2. one continuous control varying within the cell    -> 1e-7 .. 5e-6 (FAIL)
*   3. the same control on a dollar scale                -> FAIL
*   4. the same control made CONSTANT within the cell    -> ~1e-14   (clean)
*   5. the control absorbed by fe() (no coefficient)     -> ~1e-14   (clean)
*   6. a FACTOR control                                  -> ~1e-14   (clean)
*   then the two printed tables, to diff by eye.
*
*   So the trigger is precisely: a CONTINUOUS control that carries a
*   coefficient and VARIES WITHIN the collapse's own (cohort, relative-period,
*   targetflag) cell.  The collapsed row carries one arbitrary row's value of
*   it.  The gate's argument (wooldid.ado:5588-5610) is that such a control
*   cancels between xb(targetflag=1) and xb(targetflag=0) - true for the point
*   estimate, as arm 2's e(b) column confirms - but the argument is never
*   extended to margins' delta-method Jacobian, where the same cancellation is
*   between two nearly equal quantities and leaves residue instead of zero.
*
*   Corroboration from the wider sweep this file is distilled from: every
*   impvalues run is clean because impvalues sets `nose' and so has no
*   delta-method VCE at all; and in every failing cell the REGRESSION-stage
*   scalars (F N df_r ll rmse rss r2) are BIT-IDENTICAL between the arms, so
*   the divergence is confined to the margins stage.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 8675309
set sortseed 8675309

*---- an ordinary staggered panel: 120 firms in 3 treated cohorts + never-treated,
*---- 10 years, i = the cohort id (which is wooldid's own guidance for a coarser i)
set obs 120
gen long unit   = _n
gen byte cohort = mod(_n-1,4)+1
gen double ttre = .
replace ttre = 2004 if cohort==1
replace ttre = 2006 if cohort==2
replace ttre = 2008 if cohort==3
gen double a_i  = rnormal()
gen byte region = mod(unit,5) + 1
expand 10
bysort unit: gen int t = 2000 + _n - 1
* xvary varies WITHIN a (cohort,t) cell;  xcellconst does not
gen double xvary = rnormal()
bysort cohort t (unit): gen double xcellconst = xvary[1]
* the same covariate on an administrative scale (a dollar amount)
gen double xdollar = 52000 + 9000*xvary
gen double y = a_i + 0.4*(t-2000) + 0.5*xvary + 1.2*(ttre<. & t>=ttre) + rnormal()*0.3
compress
save "_nfc01.dta", replace

dis _n as text "{hline 78}"
dis as text "arm                        collapse? | max rel e(V) | max rel e(b) | rel joint p"
dis as text "{hline 78}"

foreach case in NOISEFLOOR NOCTRL CTRLVARY CTRLDOLLAR CTRLCELLCONST FEABSORB CTRLFACTOR {
    local ctrl
    if "`case'" == "NOISEFLOOR"    local ctrl controls(xvary)
    if "`case'" == "CTRLVARY"      local ctrl controls(xvary)
    if "`case'" == "CTRLDOLLAR"    local ctrl controls(xdollar)
    if "`case'" == "CTRLCELLCONST" local ctrl controls(xcellconst)
    if "`case'" == "FEABSORB"      local ctrl fe(region)
    if "`case'" == "CTRLFACTOR"    local ctrl controls(i.region)

    * ---- arm A : the collapse, which is the default ----
    use "_nfc01.dta", clear
    cap noi qui wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests `ctrl'
    local rcA = _rc
    if `rcA' == 0 {
        matrix VA = e(V)
        matrix bA = e(b)
        scalar pA = e(joint_horiz)
        local fcA `e(wdid_fcol)'
    }

    * ---- arm B : nofastcollapse.  For the NOISEFLOOR case, arm A again instead. ----
    local nofc nofastcollapse
    if "`case'" == "NOISEFLOOR" local nofc
    use "_nfc01.dta", clear
    cap noi qui wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests `ctrl' `nofc'
    local rcB = _rc
    if `rcB' == 0 {
        matrix VB = e(V)
        matrix bB = e(b)
        scalar pB = e(joint_horiz)
        local fcB `e(wdid_fcol)'
    }

    if `rcA' != 0 | `rcB' != 0 {
        dis as error "`case': rcA=`rcA' rcB=`rcB' - arm failed, skipping"
        continue
    }
    mata: VA=st_matrix("VA"); VB=st_matrix("VB"); bA=st_matrix("bA"); bB=st_matrix("bB"); dV=abs(VA-VB); sV=abs(VA); sV=sV:*(sV:>=abs(VB))+abs(VB):*(sV:<abs(VB)); st_numscalar("mxV", max(dV:/(sV:+1e-300))); dB=abs(bA-bB); sB=abs(bA); sB=sB:*(sB:>=abs(bB))+abs(bB):*(sB:<abs(bB)); st_numscalar("mxb", max(dB:/(sB:+1e-300)))
    local pp = abs(pA-pB)/max(abs(pA),abs(pB))
    dis as result %-26s "`case'" "  A=`fcA' B=`fcB'  | " %9.2e mxV " | " %9.2e mxb " | " %9.2e `pp'
}
dis as text "{hline 78}"
dis as text "PASS means <= 1e-10 (the project's indifference threshold)."
dis as text "NOISEFLOOR must be exactly 0, or this harness - not wooldid - is at fault."

dis _n as text "{hline 78}"
dis as text "THE PRINTED TABLE.  Two continuous controls.  Diff the two blocks below."
dis as text "{hline 78}"
dis "<<<<< PRINTED_A  (collapse on, the default) >>>>>"
use "_nfc01.dta", clear
cap noi wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xdollar xvary)
dis "<<<<< PRINTED_B  (nofastcollapse) >>>>>"
use "_nfc01.dta", clear
cap noi wooldid y cohort t ttre, att esprelength(3) espostlength(3) jointtests controls(xdollar xvary) nofastcollapse
dis "<<<<< PRINTED_END >>>>>"
erase "_nfc01.dta"
dis "REPRO_NFC01_DONE"
