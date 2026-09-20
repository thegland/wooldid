*=====================================================================
* s6_lattice_big.do  -- SELF-CONTAINED.  ~8 minutes.
*   "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do s6_lattice_big.do
*
* MY OWN repro of seeded claim hv-lat-01, scaled so the actual column count
* crosses Stata's UNTOUCHED DEFAULT c(maxvar) = 5000 while si8demand stays far
* below it.  Own data shape (50 states x 8 counties x 22 years, lattice(all 5)),
* not the prober's.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424299
set sortseed 424299
which wooldid
dis as text "UNTOUCHED default c(maxvar) = " c(maxvar)

clear
set obs 50
gen long state = _n
gen int ttre = cond(_n <= 46, 2004 + mod(_n-1,9), .)
gen double sfe = rnormal()
expand 8
sort state
by state: gen int cnum = _n
gen long county = state*100 + cnum
gen double exposure = exp(rnormal()*0.6)
gen double cfe = rnormal()*0.5
expand 22
sort county
by county: gen int t = 2000 + _n - 1
gen double yfe = 0.25*cos(t/4)
gen byte post = (t >= ttre) & ttre < .
gen double y = sfe + cfe + yfe + 0.20*exposure + 0.70*exposure*post + rnormal()*0.45
dis as result "PANEL rows=" _N "  states=50 (46 treated)  counties=400  years=22"
tempfile D
qui save `D'

*--- what si8demand (ado:1371) will compute ---
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
dis as result "c(maxvar)=" c(maxvar) " -> preflight fires? " (`pred' > c(maxvar))

*--- independent nonempty-cell count of the actual lattice block ---
qui gen long _i = cond(post, state, 0)
qui egen long _t = group(rel) if post
qui replace _t = 0 if missing(_t)
qui gen byte tf = post
qui gquantiles lat5 = exposure , xtile nquantiles(5)
tempvar a b
qui egen byte `a' = tag(_i _t tf)
qui count if `a'==1
local nc = r(N)
qui egen byte `b' = tag(_i _t tf lat5)
qui count if `b'==1
local ncl = r(N)
dis as result "ACTUAL columns, no lattice       ~ " 2*`nc' + 1
dis as result "ACTUAL columns, lattice(all 5)   ~ " 2*(`nc'+`ncl') + 2*5 + 1

*--- RUN 1: no lattice ---
use "`D'", clear
dis as text "=== RUN 1: att contreat(exposure) contreatcontrols(basic), NO lattice ==="
timer clear 1
timer on 1
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic)
local rc1 = _rc
timer off 1
dis as result "RUN1_RC=" `rc1'
timer list 1

*--- RUN 2: + lattice(all 5) ---
use "`D'", clear
dis as text "=== RUN 2: same + lattice(all 5), DEFAULT maxvar ==="
dis as text "vvvvvv EVERYTHING BELOW IS WHAT THE USER SEES vvvvvv"
timer clear 2
timer on 2
cap noi wooldid y state t ttre , att contreat(exposure) contreatcontrols(basic) lattice(all 5)
local rc2 = _rc
timer off 2
dis as text "^^^^^^ END OF WHAT THE USER SEES ^^^^^^"
dis as result "RUN2_RC=" `rc2'
timer list 2
dis as result "dataset still in memory?  _N=" _N "  vars=" c(k)
dis as result "WOOLDID_LASTRC = " "$WOOLDID_LASTRC"

*--- RUN 3: bare replay after the failure (hv-lat-03 side claim) ---
dis as text "=== RUN 3: bare -wooldid- replay after the r(907) ==="
cap noi wooldid
dis as result "RUN3_RC=" _rc

dis "S6_DONE"
