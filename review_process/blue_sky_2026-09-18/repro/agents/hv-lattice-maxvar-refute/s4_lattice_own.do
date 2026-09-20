*=====================================================================
* s4_lattice_own.do  -- SELF-CONTAINED.  ~6 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s4_lattice_own.do
*
* MY OWN repro of seeded claim hv-lat-01, on MY OWN data shape, at Stata's
* UNTOUCHED DEFAULT c(maxvar).  Deliberately a different design from the
* prober's (theirs: 50 states x 15 counties x 20 years, lattice(all 5); mine:
* 40 states x 10 counties x 18 years, lattice(all 4) and (all 5)) so that
* nothing rests on their script.
*
* Q1: does si8demand stay BELOW maxvar while the lattice run dies r(907)?
* Q2: does the identical spec WITHOUT lattice() complete?
* Q3: what exactly does the user see, verbatim?           (hv-lat-03)
* Q4: does the dataset survive, and is WOOLDID_LASTRC set? (hv-lat-03)
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
which wooldid
dis as text "UNTOUCHED default c(maxvar) = " c(maxvar)

*---------------- data: ordinary shift-share style exposure panel -----
clear
set obs 40
gen long state = _n
gen int ttre = cond(_n <= 36, 2004 + mod(_n-1,8), .)
gen double sfe = rnormal()
expand 10
sort state
by state: gen int cnum = _n
gen long county = state*100 + cnum
gen double exposure = exp(rnormal()*0.6)
gen double cfe = rnormal()*0.5
expand 18
sort county
by county: gen int t = 2000 + _n - 1
gen double yfe = 0.25*cos(t/4)
gen byte post = (t >= ttre) & ttre < .
gen double y = sfe + cfe + yfe + 0.20*exposure + 0.70*exposure*post + rnormal()*0.45
dis as result "PANEL rows=" _N "  states=40 (36 treated)  counties=400  years=18"
tempfile D
qui save `D'

*---------------- what si8demand (ado:1371) will compute ---------------
* si8cells = distinct TREATED (cohort, relative-period) cells; si8ilev = distinct i;
* si8tlev = distinct t; si8demand = cells*(1+contreat) + ilev + tlev + 200
gen int rel = t - ttre
tempvar tg
qui egen byte `tg' = tag(state rel) if post
qui count if `tg'==1
local cells = r(N)
qui levelsof state, local(LS)
local ilev : word count `LS'
qui levelsof t, local(LT)
local tlev : word count `LT'
local pred = `cells'*2 + `ilev' + `tlev' + 200
dis as result "PREDICTED si8demand = " `pred' "  (cells=" `cells' ", ilev=" `ilev' ", tlev=" `tlev' ")"
dis as result "c(maxvar)=" c(maxvar) " -> guard fires? " (`pred' > c(maxvar))

*---------------- independent count of the ACTUAL lattice block --------
* fvexpand itself dies r(907) at this size, so count nonempty cells directly.
* The reghdfe design (ado:1476) is
*    i._i#i._t#i.targetflag ## i.lat ## c.dose
* whose nonempty-column count is (cells_no_lat + cells_with_lat)*2 + 2*L + 1.
qui gen long _i = cond(post, state, 0)
qui egen long _t = group(rel) if post
qui replace _t = 0 if missing(_t)
qui gen byte tf = post
foreach L in 4 5 {
  cap drop lat`L'
  qui gquantiles lat`L' = exposure , xtile nquantiles(`L')
  tempvar a b
  qui egen byte `a' = tag(_i _t tf)
  qui count if `a'==1
  local nc = r(N)
  qui egen byte `b' = tag(_i _t tf lat`L')
  qui count if `b'==1
  local ncl = r(N)
  dis as result "ACTUAL columns, no lattice        ~ " 2*`nc' + 1
  dis as result "ACTUAL columns, lattice(all `L')   ~ " 2*(`nc'+`ncl') + 2*`L' + 1
  drop `a' `b'
}

*---------------- RUN 1: contreat, NO lattice --------------------------
use "`D'", clear
dis as text "=== RUN 1: att contreat(exposure) contreatcontrols(basic), NO lattice, DEFAULT maxvar ==="
timer clear 1
timer on 1
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
local rc1 = _rc
timer off 1
dis as result "RUN1_RC=" `rc1'
timer list 1

*---------------- RUN 2: same + lattice(all 4) -------------------------
use "`D'", clear
dis as text "=== RUN 2: same + lattice(all 4), DEFAULT maxvar ==="
dis as text "vvvvvv EVERYTHING BELOW IS WHAT THE USER SEES vvvvvv"
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 4)
local rc2 = _rc
dis as text "^^^^^^ END OF WHAT THE USER SEES ^^^^^^"
dis as result "RUN2_RC=" `rc2'
dis as result "dataset still in memory?  _N=" _N "  vars=" c(k)
dis as result "WOOLDID_LASTRC = " "$WOOLDID_LASTRC"

*---------------- RUN 3: same + lattice(all 5) -------------------------
use "`D'", clear
dis as text "=== RUN 3: same + lattice(all 5), DEFAULT maxvar ==="
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
dis as result "RUN3_RC=" _rc

*---------------- RUN 4: lattice(all 4) at maxvar high is not tested here
* (a separate session is needed to change maxvar; see s5).

dis "S4_DONE"
