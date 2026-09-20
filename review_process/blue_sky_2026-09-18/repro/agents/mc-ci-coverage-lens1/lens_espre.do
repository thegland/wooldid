*=========================================================================
* lens_espre.do -- DOCUMENTED-BEHAVIOUR lens on mcci-03-esprelength-moves-att
*
* The behaviour (esprelength() > 0 re-bases the headline att) is already
* documented at sthlp:1493-1497 and was ruled as item 13 / T-BJS1 of the
* 2026-09-13 slate.  What this probe tests is the ONE sentence in that
* shipped passage that may not be true:
*
*   "requesting pre-periods removes them from that pool and shifts every
*    post horizon BY A CONSTANT (measured 4.8e-02 on a common-onset panel)"
*
* The project's own measurement file (held_for_ruling_2026-09-13.md:735)
* records the shift as constant on common-onset fixtures but NOT constant
* on staggered ones (costag 1.17e-2 "not constant", dta_binary post0 shift
* 3.71e-2 "not constant").  So: measure the per-horizon shift on BOTH a
* common-onset and a staggered panel and see which way it falls.
*
* Also re-checks, independently, the skeptic's correction that
* espostlength() alone moves nothing.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do lens_espre.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*-------------------------------------------------------------------------
* helper: build a panel.  mode = "stag" (3 onsets) or "hom" (common onset)
*-------------------------------------------------------------------------
capture program drop mkpanel
program define mkpanel
    args mode
    clear
    set obs 80
    gen long unit = _n
    gen byte cg = mod(_n-1,4)+1
    gen double ttre = .
    if "`mode'" == "stag" {
        replace ttre = 2004 if cg==1
        replace ttre = 2006 if cg==2
        replace ttre = 2008 if cg==3
    }
    else {
        replace ttre = 2006 if cg<=3
    }
    gen double ui = rnormal(0,0.70)
    expand 12
    bys unit: gen int t = 2000 + _n - 1
    gen double e = .
    bys unit (t): replace e = rnormal(0,0.55) if _n==1
    bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
    gen int rel = t - ttre
    gen double truetau = 0
    replace truetau = 1.00 + 0.10*rel if !missing(ttre) & rel>=0
    gen double y = ui + 0.25*(t-2000) + truetau + e
end

*-------------------------------------------------------------------------
foreach mode in hom stag {

    set seed 771103
    set sortseed 771103
    mkpanel `mode'
    tempfile d
    qui save "`d'"

    di as txt _n(2) "==================== PANEL = `mode' ===================="

    * --- headline att under four specs -----------------------------------
    di as txt _n "HEADLINE_`mode'  spec            b                se              N   df_r  refN"
    foreach spec in none post3 pre3 both3 {
        use "`d'", clear
        local opt ""
        if "`spec'"=="post3" local opt "espostlength(3)"
        if "`spec'"=="pre3"  local opt "esprelength(3) espostlength(-999)"
        if "`spec'"=="both3" local opt "esprelength(3) espostlength(3)"
        capture noisily qui wooldid y cg t ttre, att cluster(unit) `opt'
        di as txt "RC_`mode'_`spec'=" _rc
        if _rc == 0 {
            mat T = e(wdidmainresults)
            local rN = "."
            capture local rN = e(escomp_refN)
            di as txt "HEADLINE_`mode' " %-8s "`spec'" "  " %16.12f T[1,1] "  " %16.12f T[1,2] "  " ///
                %5.0f e(N) " " %5.0f e(df_r) "  `rN'"
            di as txt "   hex b = " %21x T[1,1] "   hex se = " %21x T[1,2]
        }
    }

    * --- per-horizon post coefficients: esprelength(0) vs esprelength(3) --
    use "`d'", clear
    capture noisily qui wooldid y cg t ttre, att cluster(unit) esprelength(0) espostlength(3)
    di as txt "RC_`mode'_A=" _rc
    mat A = e(wdidmainresults)

    use "`d'", clear
    capture noisily qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3)
    di as txt "RC_`mode'_B=" _rc
    mat B = e(wdidmainresults)

    di as txt _n "SHIFT_`mode'   row                          b(pre0)          b(pre3)          shift"
    local rnB : rownames B
    local nb = rowsof(B)
    forvalues r = 1/`nb' {
        local nm : word `r' of `rnB'
        local ra = rownumb(A, "`nm'")
        if `ra' < . {
            local bA = A[`ra',1]
            local bB = B[`r',1]
            local dd = `bB' - `bA'
            di as txt "SHIFT_`mode' " %-24s "`nm'" " " %16.12f `bA' " " %16.12f `bB' " " %16.12f `dd'
        }
    }
}

di as txt _n "ALLDONE_LENS1"
