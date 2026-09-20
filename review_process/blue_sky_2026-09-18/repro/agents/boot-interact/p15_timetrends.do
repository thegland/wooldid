*=============================================================================
* p15 — SELF-CONTAINED REPRO.
* CLAIM: with timetrends + cluster() at the level of i, bootstrap() RESCUES the printed
* table (every se/t/p/CI is a bootstrap quantity) but `ereturn repost V =` fails rc 503
* because the estimate was posted without a covariance, so e(V) is never the bootstrap
* covariance, e(wdid_vcvstate) reads "absent", e(V_analytic) does not exist, and every
* downstream reader (test, lincom, esttab, coefplot, honestdid) sees nothing on a run whose
* own table is fully populated.  sthlp:1059 describes this design without carving out
* bootstrap(): "every standard error, test statistic, p-value and interval in the run is
* missing" - which is false under bootstrap().
* Run: stata-mp -e do p15_timetrends.do
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
clear
set obs 40
gen int state = _n
gen int onset = .
replace onset = 5 if state <= 10
replace onset = 7 if state > 10 & state <= 20
replace onset = 9 if state > 20 & state <= 30
gen double sfe = rnormal()
expand 12
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + sfe + 0.05*t + te + rnormal()*0.8
local O cluster(state) esfixedbaseperiod esprelength(3) espostlength(3) timetrends
dis "##### analytic #####"
cap noi wooldid y state t onset, `O'
dis "AN_RC=" _rc " VCV=[`e(wdid_vcvstate)']"
cap confirm matrix e(V)
dis "AN_HASV_RC=" _rc
dis "##### bootstrap #####"
cap noi wooldid y state t onset, `O' bootstrap(standard) reps(60) rseed(4242) dots(0)
dis "BO_RC=" _rc " VCV=[`e(wdid_vcvstate)'] VPOSTED=" e(boot_v_posted) " COMPLETED=" e(boot_completed)
cap confirm matrix e(V)
dis "BO_HASV_RC=" _rc
cap confirm matrix e(V_boot)
dis "BO_HASVBOOT_RC=" _rc
cap confirm matrix e(V_analytic)
dis "BO_HASVANALYTIC_RC=" _rc
matrix M = e(wdidmainresults)
matrix list M, format(%12.6g)
cap noi test main_att = 0
dis "BO_TEST_RC=" _rc " p=" r(p)
cap noi lincom main_att
dis "BO_LINCOM_RC=" _rc " se=" r(se)
dis "P15_DONE"
