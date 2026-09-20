clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- self-contained staggered-DiD panel (60 units, T=10, 3 onsets) ----
clear
set seed 4242
set sortseed 4242
qui set obs 60
gen long unit = _n
gen double a_i = rnormal(0,1)
gen double u = runiform()
gen double ttre = .
qui replace ttre = 2003 if u > 0.35 & u <= 0.5667
qui replace ttre = 2005 if u > 0.5667 & u <= 0.7833
qui replace ttre = 2007 if u > 0.7833
qui expand 10
bysort unit: gen int t = 1999 + _n
gen byte D = (ttre < . & t >= ttre)
gen double x1 = rnormal()
gen double y = a_i + 0.4*(t-2000) + 2*D + 0.2*x1 + 0.3*rnormal()
gen double dose = D*(0.5 + runiform())
gen long i = unit
egen long icoh = group(ttre)
qui replace icoh = 99 if ttre == .
drop u a_i

* ============ FUZZ-07: `update' on a full command line ============
* WARNING: THIS PERFORMS A NETWORK PACKAGE INSTALL INTO THE USER'S PLUS
* DIRECTORY (ado:559 -- net install wooldid, from(github .../wooldidinstall/) replace)
* and it is DELIBERATELY LEFT COMMENTED OUT here.  The mechanism is visible in
* the source without running it: wooldid_worker parses its syntax line, then
*
*     if "`update'" == "update" {
*         net install wooldid, from(...) replace
*         dis as text "Update complete; please restart stata ..."
*         exit
*     }
*
* runs BEFORE wooldid_validate and exits unconditionally, so the estimation
* request on the same line is discarded and nothing is posted to e().
* Uncomment the two lines below only if you accept the network install.
* cap noi wooldid y i t ttre , att robust update
* dis "UPDATE rc=" _rc "  e(cmd)=[" e(cmd) "]"
dis "see the source at ado:555-563; not executed (network install)"
