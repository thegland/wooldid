*=======================================================================
* SELF-CONTAINED REPRO -- fuzz agent
* CLAIM: `estimates replay <name>' after a successful wooldid run fails
*        r(198) and prints a message that blames the USER'S COMMAND LINE
*        ("command line input should feature all of y i t ttre; too few
*        variables specified"), although the user typed no wooldid
*        command line at all.  Bare `wooldid' replay works, and
*        `estimates restore' + bare `wooldid' works, so the stored
*        results are fine; only the standard `estimates replay' entry
*        point is broken, and its diagnostic is actively misleading.
* RUN:  cd <somewhere writable>
*       /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do REPRO_fuzz_estreplay.do
*=======================================================================
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

* ---------------- a plain staggered panel, built from scratch ----------------
set seed 20260917
set sortseed 20260917
clear
set obs 60
gen long unit = _n
gen byte region = 1 + mod(_n,5)
gen double a_i  = rnormal()
gen ttre = .
replace ttre = 2003 if unit <= 12
replace ttre = 2005 if unit > 12 & unit <= 24
expand 8
bysort unit: gen int t = 2000 + _n - 1
gen double g_t = 0.3*(t-2000)
gen byte post = !missing(ttre) & t >= ttre
gen double y = a_i + g_t + 1.2*post + 0.3*rnormal()

* ---------------- estimate ----------------
wooldid y unit t ttre, att att_i cluster(region) esprelength(2) espostlength(2)
di "ESTIMATION_RC = " _rc

di _n "=== (1) bare `wooldid' replay: WORKS ==="
capture noisily wooldid
di "BARE_RC = " _rc

di _n "=== (2) estimates store: WORKS ==="
capture noisily estimates store m1
di "STORE_RC = " _rc

di _n "=== (3) estimates replay m1: FAILS r(198) with a misleading message ==="
capture noisily estimates replay m1
di "REPLAY_RC = " _rc

di _n "=== (4) estimates restore m1 + bare wooldid: WORKS (so e() is intact) ==="
qui wooldid y unit t ttre, att robust
capture noisily estimates restore m1
di "RESTORE_RC = " _rc
capture noisily wooldid
di "BARE2_RC = " _rc

di _n "EXPECTED: REPLAY_RC = 0 and the stored table redisplayed (or, at minimum,"
di "          a refusal that does not blame a command line the user never typed)."
di "OBSERVED: REPLAY_RC = 198 with 'command line input should feature all of y i t ttre'."
