*=====================================================================
* repro_lattice_maxvar.do   -- SELF-CONTAINED.  Run as-is, ~2 minutes.
*
*   cd <this directory>
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_lattice_maxvar.do
*
* SHOWS: wooldid's maxvar preflight (ado:1342-1374) never reads lattice(),
* so on an ordinary continuous-treatment DiD at Stata's DEFAULT maxvar the
* guard predicts 1320 variables, stays silent, and the run dies on a raw
* r(907) carrying no wooldid message.  The identical spec without the
* lattice completes in well under a minute.
*
* DATA: 50 states, 45 with staggered onsets 2005-2012, 15 counties each,
* 20 years.  Continuous county-level EXPOSURE, defined on every row and
* time-invariant (Bartik / shift-share shape).  Nothing pathological.
*
* NOTE on the counting: `fvexpand' cannot be used to count the lattice
* block, because it materialises the full factor grid BEFORE applying
* `set emptycells drop' and so hits r(907) itself at the default maxvar.
* The identity below counts the NONEMPTY cells directly.  It is calibrated:
* on the smaller 25-state panel it reproduces fvexpand's answer (2579)
* exactly.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"

dis as text "Stata " c(flavor) " " c(stata_version) ",  DEFAULT c(maxvar) = " c(maxvar) ///
            ",  c(max_k_theory) = " c(max_k_theory)
which wooldid

clear
set seed 424242
set sortseed 424242
set obs 50
gen long state = _n
gen int ttre = cond(_n <= 45, 2005 + mod(_n-1,8), .)
expand 15
sort state
by state: gen int cnum = _n
gen long county = state*1000 + cnum
gen double exposure = exp(rnormal()*0.7)
gen double cfe      = rnormal()
expand 20
sort county
by county: gen int t = 1999 + _n
gen double yfe  = 0.30*sin(t/3)
gen byte  post  = (t >= ttre) & ttre < .
gen double y    = cfe + yfe + 0.25*exposure + 0.80*exposure*post + rnormal()*0.40
tempfile D
save `D'

*--- what si8demand (ado:1371) computes: cells*(1+contreat) + ilev + tlev + 200 ---
gen int h = t - ttre
tempvar tg
qui gegen byte `tg' = tag(state h) if post
qui count if `tg' == 1
local cells = r(N)
qui levelsof state, local(L)
local ilev : word count `L'
qui levelsof t, local(L2)
local tlev : word count `L2'
local pred = `cells'*2 + `ilev' + `tlev' + 200
dis as result "PREFLIGHT PREDICTS si8demand = " `pred' ///
              "   (cells=" `cells' ", ilev=" `ilev' ", tlev=" `tlev' ")"
dis as result "c(maxvar) = " c(maxvar) "  ->  " `pred' " < " c(maxvar) ", so the preflight will NOT fire"

*--- what the regression actually demands (nonempty-cell identity) -------------
gen double _ii = cond(post, state, 0)
qui egen long _hh = group(h) if post
qui replace _hh = 0 if missing(_hh)
gen byte _tf = post
qui gquantiles _lat = exposure , xtile nquantiles(5)
tempvar c1 c2
qui gegen byte `c1' = tag(_ii _hh _tf)
qui count if `c1' == 1
local ncell = r(N)
qui gegen byte `c2' = tag(_ii _hh _tf _lat)
qui count if `c2' == 1
local ncelllat = r(N)
dis as result "ACTUAL columns, no lattice      = " 2*`ncell' + 1
dis as result "ACTUAL columns, lattice(all 5)  = " 2*(`ncell' + `ncelllat') + 2*5 + 1 ///
              "   <-- exceeds c(maxvar)"

*--- RUN 1: no lattice -> completes --------------------------------------------
use `D', clear
dis as text "=== RUN 1: contreat(), NO lattice, default maxvar ==="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
dis as result "RUN1_RC = " _rc

*--- RUN 2: add lattice(all 5) -> r(907), preflight silent ----------------------
use `D', clear
dis as text "=== RUN 2: same spec + lattice(all 5), default maxvar ==="
dis as text "--- everything below is what the user sees ---"
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
local rc = _rc
dis as text "--- end of what the user sees ---"
dis as result "RUN2_RC = " `rc'
if `rc' == 907 dis as result "CONFIRMED: preflight predicted `pred' and said nothing; the run died r(907)."
dis "REPRO_DONE"
