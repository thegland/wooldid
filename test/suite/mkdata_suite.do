* mkdata_suite.do — deterministic builder for the NEW wooldid-suite datasets.
* Usage:  do mkdata_suite.do <outdir>     (default outdir: suite_data)
* Every dataset block sets its own seed AND sortseed; two runs are value-identical
* (verify with cf _all / datasignature — .dta headers embed a timestamp, so file
* hashes will differ). A hard self-test assert block at the end reloads every
* dataset and checks N, unit counts, and the trigger-relevant properties.
*
* Schema convention (battery style, where sensible): i t ttre y [ycount w sg]
* plus dataset-specific probe columns documented per block.
* Ledger references (divergence_ledger_draft.md): A6-B9, A14, A12, A10+A11,
* A16, A15/B-F13, A8(b), B-clus, A7-B10.
* Datasets-audit fixtures (2026-09-13): dta_s_stateyear (P1), dta_s_biennial (P2).

version 19
clear all
set more off
set linesize 255

args outdir
local outdir = subinstr(`"`outdir'"', "\", "/", .)
if `"`outdir'"' == "" local outdir "suite_data"
capture mkdir `"`outdir'"'
local OUT `"`outdir'"'

* =====================================================================
* dta_s_extmiss — trigger A6-B9 (extended missings; v1.3's `== .` misses .a)
* Faithful clone of test/data/dta_binary.dta (mkdata.do recipe, same seeds,
* DEFAULT FLOAT types so the pre-injection clone is cf-identical to the repo
* file), then .a injected into y, t, w on deterministic mod()-selected rows.
* 80 units (ttre: 2005 for i<=25, 2008 for 26..50, never for 51..80) x
* t=2000..2011 -> N=960. Injected rows (computed analytically, disjoint
* selections except one w&t overlap): y=.a on 24 rows [mod(i+t,37)==0],
* t=.a on 20 rows [mod(3i+t,53)==0], w=.a on 23 rows [mod(i+2t,41)==0].
* Hits land in both treated cohorts and the never pool.
* Cleaned-reference twin for the comparator: drop if y>=. | t>=. | w>=.
* =====================================================================
clear
set type float
set seed 20260719
set sortseed 20260719
set obs 80
gen i = _n
gen ttre = 2005 if i <= 25
replace ttre = 2008 if i > 25 & i <= 50
gen sg = 1 + mod(i,2)
gen ife = rnormal()
expand 12
bysort i: gen t = 1999 + _n
gen w = 0.25 + rchi2(3)
gen eff = cond(ttre==2005, 1.0, 0.6)
gen y = 2 + ife + 0.05*(t-2000) + eff*(t >= ttre & ttre != .) + rnormal(0, 0.5)
gen ycount = round(exp(0.2*y))
* --- injection: flag first (so selections use the ORIGINAL t), then replace
gen byte _hity = mod(i + t, 37) == 0
gen byte _hitt = mod(3*i + t, 53) == 0
gen byte _hitw = mod(i + 2*t, 41) == 0
replace y = .a if _hity
replace t = .a if _hitt
replace w = .a if _hitw
drop _hity _hitt _hitw
save `"`OUT'/dta_s_extmiss.dta"', replace

* =====================================================================
* dta_s_deeppre — trigger A14 (SG adversarial esf skiplocation: v1.3 r(2000),
* HEAD rc 0 with finite per-SG p). Recorded campaign recipe: 60 units,
* t=1..14, cohorts ttre in {8,10} (20 units each) + 20 never-treated,
* sg = mod(i,2), y = 2*(t>=ttre & ttre<.) + 0.05*i + rnormal().
* NOTE: deliberately NO time-trend term in y — a growing outcome trips an
* unrelated r(459) estimability guard (campaign note).  N=840.
* =====================================================================
clear
set type double
set seed 60814
set sortseed 60814
set obs 60
gen i = _n
gen ttre = cond(i <= 20, 8, cond(i <= 40, 10, .))
gen sg = mod(i,2)
expand 14
bysort i: gen t = _n
gen y = 2*(t >= ttre & ttre < .) + 0.05*i + rnormal()
gen w = 0.5 + runiform()
save `"`OUT'/dta_s_deeppre.dta"', replace

* =====================================================================
* dta_s_ttyear — trigger A12 small (timetrends at 4-digit calendar years:
* ~1e-8-class carve-out old vs new). 40 units x t=2000..2012 (N=520),
* cohorts 2006 (i<=12) / 2009 (13..24) / never (25..40). Mild common trend
* + small deterministic unit-trend heterogeneity (0.01*mod(i,5) per year)
* so timetrends has real signal. ycount included for the poisson+timetrends
* variant (tt_pois rss/rmse pin).
* =====================================================================
clear
set type double
set seed 20120
set sortseed 20120
set obs 40
gen i = _n
gen ttre = cond(i <= 12, 2006, cond(i <= 24, 2009, .))
gen sg = 1 + mod(i,2)
expand 13
bysort i: gen t = 1999 + _n
gen w = 0.5 + runiform()
gen y = 1 + 0.05*i + 0.04*(t-2000) + 0.01*mod(i,5)*(t-2000) + 1.5*(t >= ttre & ttre < .) + rnormal(0, 0.4)
gen ycount = max(1, round(exp(0.15*y)))
save `"`OUT'/dta_s_ttyear.dta"', replace

* =====================================================================
* dta_s_ttbig — trigger A12 large (t ~ 5e5: v1.3 timetrends numerically
* INERT — rc 0, results == no-timetrends run; HEAD applies real trends).
* Same shape as dta_s_ttyear but t=500000..500012, cohorts 500006/500009.
* =====================================================================
clear
set type double
set seed 50812
set sortseed 50812
set obs 40
gen i = _n
gen ttre = cond(i <= 12, 500006, cond(i <= 24, 500009, .))
gen sg = 1 + mod(i,2)
expand 13
bysort i: gen t = 499999 + _n
gen w = 0.5 + runiform()
gen y = 1 + 0.05*i + 0.04*(t-500000) + 0.01*mod(i,5)*(t-500000) + 1.5*(t >= ttre & ttre < .) + rnormal(0, 0.4)
gen ycount = max(1, round(exp(0.15*y)))
save `"`OUT'/dta_s_ttbig.dta"', replace

* =====================================================================
* dta_s_bigcodes — triggers A10 (unit IDs in [2^24,2^31): v1.3 float _i
* pools distinct units) + A11 (YYYYMMDD onset/time codes: v1.3 float gen
* sites misplace horizons / r(461)). 40 units, IDs 30000001..30000040
* (CONSECUTIVE — float spacing at 3e7 is 2, so odd IDs round onto even
* neighbors: pooling is provable in-data). Time grid: 12 consecutive
* month-starts 20200101..20201201 (uniform spacing 100, all odd values).
* Cohorts ttre = 20200601 (i_small<=13) / 20200901 (14..26) / never.
* Small-recode twin columns for the HEAD == recode-twin pin:
*   i_small = rank of i (1..40), t_small = month index (1..12),
*   ttre_small = 6/9/. — y is built from the SMALL codes, so runs under
* either coding target identical estimands. N=480.
* =====================================================================
clear
set type double
set seed 31710
set sortseed 31710
set obs 40
gen double i_small = _n
gen double i = 30000000 + _n
gen double ttre_small = cond(_n <= 13, 6, cond(_n <= 26, 9, .))
gen double ttre = cond(_n <= 13, 20200601, cond(_n <= 26, 20200901, .))
gen sg = 1 + mod(_n,2)
expand 12
bysort i: gen double t_small = _n
gen double t = 20200001 + t_small*100
gen w = 0.5 + runiform()
gen y = 2*(t_small >= ttre_small & ttre_small < .) + 0.05*i_small + 0.03*t_small + rnormal(0, 0.5)
save `"`OUT'/dta_s_bigcodes.dta"', replace

* =====================================================================
* dta_s_imdegen — trigger A16 (IM degeneracy under impvalues: per-cluster
* estimates numerically identical -> HEAD sets IM se/t/p/CI missing +
* Caution; v1.3 posts zero-diagonal e(V)). NOISELESS outcomes: exact
* functions of (cohort group, t, horizon) — absorbed i/t FE and the cell
* dummies fit them with zero residuals, so every per-cluster estimate
* equals the truth exactly. 300 units x t=1..12 (N=3600),
* cohgrp = mod(i-1,3): 0=never, 1=ttre 5, 2=ttre 8 (100 units each).
* clus = 1+mod(floor((i-1)/3),12): 12 clusters, each built from unit
* TRIPLES spanning all three cohort groups (cluster orthogonal to cohort).
* Two outcome columns (verified on HEAD with
* `att cluster(clus) esprelength(2) espostlength(2) impvalues`):
*   y     — FULL degeneracy: an exact 0.5 bump on relative years -1/-2
*           makes every reported cell (incl. the pre window) NONZERO with
*           per-cluster-identical estimates, so the strict-relative guard
*           flags ALL cells -> e(V) ABSENT on HEAD, e(properties) "b".
*   ypart — PARTIAL degeneracy: no pre bump; pre cells are structural
*           zeros (b ~1e-16 FP noise, never flagged by the strict-relative
*           test) -> post cells missing, pre cells finite: A17 sentinel
*           path, e(V) PRESENT with placeholder diagonal, "3 of 6" Caution.
* No RNG used (seeds set for policy only).
* =====================================================================
clear
set type double
set seed 16263
set sortseed 16263
set obs 300
gen i = _n
gen cohgrp = mod(i-1, 3)
gen ttre = cond(cohgrp == 1, 5, cond(cohgrp == 2, 8, .))
gen clus = 1 + mod(floor((i-1)/3), 12)
gen sg = mod(i,2)
expand 12
bysort i: gen t = _n
gen ypart = 1 + 0.5*cohgrp + 0.1*t + 2*(t >= ttre & ttre < .)
gen y = ypart + 0.5*((t == ttre - 1 | t == ttre - 2) & ttre < .)
save `"`OUT'/dta_s_imdegen.dta"', replace

* =====================================================================
* dta_s_mintsent — triggers A15 + B-F13 (onset varying within unit).
* 30 units x t=1..10 (N=300); true onset: ttre_good = 5 for i<=10, never
* for 11..30. Three onset codings frozen side by side:
*   ttre      — onset recorded ONLY from the onset period onward (missing
*               pre-onset, 5 from t=5 on) for ALL treated units: the B-F13
*               probe (NEW default safety r(459); safety(off) runs).
*   ttre_sent — sentinel coding: correct everywhere except unit 1's
*               pre-onset rows carry 0 (mkmint-"one" style): the A15 probe
*               under safety(off) (v1.3 biased; HEAD == unit-1-deleted bench).
*   ttre_good — fully correct column, for benchmark/twin runs.
* =====================================================================
clear
set type double
set seed 51015
set sortseed 51015
set obs 30
gen i = _n
gen ttre_good = cond(i <= 10, 5, .)
gen sg = mod(i,2)
expand 10
bysort i: gen t = _n
gen y = 2*(t >= ttre_good & ttre_good < .) + 0.05*i + rnormal()
gen w = 0.5 + runiform()
gen ttre = ttre_good
replace ttre = . if ttre_good == 5 & t < 5
gen ttre_sent = ttre_good
replace ttre_sent = 0 if i == 1 & t < 5
save `"`OUT'/dta_s_mintsent.dta"', replace

* =====================================================================
* dta_s_covars — trigger A8(b) + zero-coverage closure for controls()/
* fe()/interactivecontrols()/coarsecohortcontrols()/customweightmultiplier.
* 150 units x t=1..12 (N=1800); cohorts by mod(i,3): 1 -> ttre 5,
* 2 -> ttre 9, 0 -> never (50 units each).
*   x1 x2    — continuous controls; y built from COMPLETE x's, then missing
*              injected on known rows: x1=. on 60 rows [mod(i+t,29)==0],
*              x2=. on 56 rows [mod(2i+t,31)==0].
*   fecat    — categorical fe() variable, levels mod(i,10) PLUS three
*              singleton levels on single rows: 101 @(i=7,t=3, treated pre),
*              102 @(i=22,t=8, treated POST), 103 @(i=149,t=6, treated pre)
*              -> reghdfe singleton drops (A7-adjacent) when fe(fecat) used.
*   icvar    — well-behaved unit-level numeric for interactivecontrols().
*   icbad    — A8(b) probe: unit-constant-plus-t DOUBLE control,
*              icbad = 1e6 + 137*i + t (level ~1e6, within-unit time
*              variation ~1): v1.3's float copy destroys it -> rc 0 garbage,
*              no e(V); HEAD keeps double.
*   cccvar   — coarse categorical (2 unit-level classes cutting across
*              cohorts) for coarsecohortcontrols().
*   cwmvar   — positive multiplier for customweightmultiplier().
* =====================================================================
clear
set type double
set seed 15012
set sortseed 15012
set obs 150
gen i = _n
gen ttre = cond(mod(i,3) == 1, 5, cond(mod(i,3) == 2, 9, .))
gen sg = 1 + mod(i,2)
gen icvar = 1 + mod(i,7) + 0.5*mod(i,3)
gen cccvar = 1 + (mod(i,5) < 2)
expand 12
bysort i: gen t = _n
gen x1 = 0.3*mod(i,5) + 0.1*t + 0.2*sin(i + t)
gen x2 = cos(0.5*i + 0.3*t)
gen fecat = mod(i,10)
replace fecat = 101 if i == 7   & t == 3
replace fecat = 102 if i == 22  & t == 8
replace fecat = 103 if i == 149 & t == 6
gen cwmvar = 0.5 + 0.1*mod(i + t, 5)
gen double icbad = 1000000 + 137*i + t
gen w = 0.5 + runiform()
gen y = 2 + 0.05*i + 0.1*t + 0.5*x1 - 0.3*x2 + 0.2*icvar + cond(ttre == 5, 1.5, 1.0)*(t >= ttre & ttre < .) + rnormal(0, 0.5)
gen ycount = max(1, round(exp(0.15*y)))
replace x1 = . if mod(i + t, 29) == 0
replace x2 = . if mod(2*i + t, 31) == 0
save `"`OUT'/dta_s_covars.dta"', replace

* =====================================================================
* dta_s_strclus — trigger B-clus (string cluster() variable: v1.3 accepts
* the varlist and passes it through; HEAD refuses at parse, r(109) from
* CLUSter(varlist numeric)). 24 units x t=1..10 (N=240); ttre = 5 for
* i<=8, never otherwise. strc = "c"+string(g) with g = 1+mod(i,5) the
* numeric twin (5 cluster groups).
* =====================================================================
clear
set type double
set seed 10924
set sortseed 10924
set obs 24
gen i = _n
gen ttre = cond(i <= 8, 5, .)
gen g = 1 + mod(i,5)
gen sg = mod(i,2)
expand 10
bysort i: gen t = _n
gen y = 1.5*(t >= ttre & ttre < .) + 0.05*i + rnormal()
gen w = 0.5 + runiform()
gen strc = "c" + string(g)
save `"`OUT'/dta_s_strclus.dta"', replace

* =====================================================================
* dta_s_singleton — trigger A7-B10 (e(sample) < screened sample via
* reghdfe singleton drops). 40 units x t=1..10 base; cohorts ttre = 4
* (i<=12) / 7 (13..24) / never (25..40) — two treated cohorts with long
* post windows so att_it is well-defined. Imbalance:
*   - scattered mod() row drops [mod(7i+t,23)==0, t>2 only, so every
*     treated unit keeps >=2 pre rows]: 14 rows;
*   - never-treated units 37-40 observed ONCE (t=5 only): singletons in
*     the ALWAYS-absorbed i.i FE -> reghdfe drops them on every config
*     (makes e(sample) < screened visible);
*   - fevar: base 1+mod(i,4) plus three singleton levels each on ONE row:
*     91 @(i=2,t=6, cohort-4 POST row), 92 @(i=15,t=9, cohort-7 POST row),
*     93 @(i=26,t=4, never row) -> with fe(fevar), reghdfe strands TREATED
*     POST rows, which is what actually moves the att_it weights (A7).
* N=352.
* =====================================================================
clear
set type double
set seed 40710
set sortseed 40710
set obs 40
gen i = _n
gen ttre = cond(i <= 12, 4, cond(i <= 24, 7, .))
gen sg = 1 + mod(i,2)
expand 10
bysort i: gen t = _n
gen w = 0.5 + runiform()
gen y = cond(ttre == 4, 1.2, 0.8)*(t >= ttre & ttre < .) + 0.05*i + 0.02*t + rnormal(0, 0.5)
gen ycount = max(1, round(exp(0.2*y)))
drop if mod(7*i + t, 23) == 0 & !inlist(t, 1, 2)
drop if i >= 37 & t != 5
gen fevar = 1 + mod(i,4)
replace fevar = 91 if i == 2  & t == 6
replace fevar = 92 if i == 15 & t == 9
replace fevar = 93 if i == 26 & t == 4
save `"`OUT'/dta_s_singleton.dta"', replace

* =====================================================================
* dta_s_stateyear — P1 of the datasets audit (2026-09-12/13): realistic
* county-in-state policy panel. 12 states of 5-12 counties (106 units) x
* year 1995-2014; onset is a STATE attribute: 1999 (st 1-2), 2002 (3-4),
* 2005 (5), 2008 (6-7), 2011 (8), 2013 (9), never (10-12) -> six cohorts
* of 17/21/6/19/5/10 counties plus 28 never-treated counties. Errors:
* AR(1) .6 within county, county-level heteroskedastic sig, plus a
* state-by-year shock nu (sd .15) that sits INSIDE the st x year cell.
* Outcomes:
*   y      - shock included: for i = unit configs with cluster(st);
*   y_nosh - y minus the shock: for i = st configs with cluster(unit)
*            (with i = st the shock lives inside the i x t cell and makes
*            cluster(unit) invalid - audit probe run2 joint pre p = .003);
*   y_dose - continuous-dose outcome (.4*dose_it*post), no shock, for
*            i = st contreat(dose_it) cluster(unit) runs (NO lattice: on
*            5-12-unit cells lattice(all 3) zeroes the VCV, audit 3.12);
*   ycount - count outcome (~22% zeros) for poisson at i = st: HEAD's
*            separation screen refuses any singleton treated cell whose
*            count is zero, so a rarer count (the audit's intercept -.8,
*            58% zeros) is refused at unit-level i and leaves all-zero
*            st x year cells at i = st (patch Run 3, 2026-09-13).
* Also: sg = 20/80 county subgroup (mod(cty,4)==0), lognormal pop (aw/pw),
* region = 1+mod(st-1,3): three regions of four states, EACH holding one
* never-treated state, so fe(i.year#i.region) is identified in every
* region-year and permstrata(region) keeps a control state per stratum
* (the audit's 1+mod(st,4) left region 2 without one: real_fe r322);
* dose_it varies within every st x year cell.
* Imbalance: late entry / early exit drawn per COUNTY at
* the one-row-per-unit stage after `sort unit` (a unique key), so the
* draws cannot depend on row order (critique A-2); ~2% internal holes
* [mod(unit*7+year,41)==0]; ~1.5% missing y [mod(unit*3+year,67)==0].
* y_nosh / y_dose are built before the y-missing injection and carry no
* missings (intended: missing-y robustness is y's job, i = unit).
* Scale: nu sd .3 -> .15 and sig .6 -> .4 (audit 5(d)). Every constant
* the self-test pins was measured on this exact block (patch Run 1).
* =====================================================================
clear
set type double
set seed 24118
set sortseed 24118
set obs 12
gen st = _n
gen ncty = 5 + mod(5*st, 8)
gen region = 1 + mod(st - 1, 3)
gen ttre = cond(st <= 2, 1999, cond(st <= 4, 2002, cond(st == 5, 2005, cond(st <= 7, 2008, cond(st == 8, 2011, cond(st == 9, 2013, .))))))
expand ncty
bysort st: gen cty = _n
gen unit = st*100 + cty
gen sg = 1 + (mod(cty, 4) == 0)
gen alpha = rnormal()
gen sig = 0.4*exp(0.4*rnormal())
gen dose = exp(0.5*rnormal())
gen pop = exp(9 + 1.2*rnormal())
replace pop = 500 if pop < 500
* entry/exit draws at the unit level on a unique sort key (critique A-2):
* independent of whatever row order an upstream edit leaves behind
sort unit
gen ulate = runiform()
gen uexit = runiform()
expand 20
bysort unit: gen year = 1994 + _n
bysort st year: gen nu = rnormal(0, 0.15) if _n == 1
bysort st year: replace nu = nu[1]
bysort unit (year): gen e = sig*rnormal()/sqrt(1 - .36) if _n == 1
bysort unit (year): replace e = .6*e[_n-1] + sig*rnormal() if _n > 1
gen post = (ttre < . & year >= ttre)
gen rel = year - ttre
gen beta = cond(ttre == 1999, .8, cond(ttre == 2002, .5, cond(ttre == 2005, 1.0, cond(ttre == 2008, .3, cond(ttre == 2011, .6, cond(ttre == 2013, .4, 0))))))
gen tau = post*(beta + .05*min(rel, 6) + .3*(sg == 2))
gen dose_it = dose*exp(.1*rnormal())
gen y = 5 + alpha + .03*(year - 1995) + .2*sin(year/2) + tau + nu + e
gen y_nosh = y - nu
gen y_dose = 5 + alpha + .03*(year - 1995) + .4*dose_it*post + .1*dose_it + e
gen ycount = rpoisson(exp(.4 + .25*post + .1*(sg == 2) + .2*alpha))
* pre-imbalance grid: every st x year cell holds its state's 5-12 counties
* (every random draw is complete here, so this by-sort moves no value)
bysort st year: gen long __ncell = _N
assert __ncell >= 5 & __ncell <= 12
drop __ncell
drop if ulate < .05 & year < 1995 + ceil(ulate*100)
drop if uexit > .95 & year > 2014 - ceil((1 - uexit)*100)
drop if mod(unit*7 + year, 41) == 0
replace y = . if mod(unit*3 + year, 67) == 0
drop ulate uexit alpha sig e nu post rel beta tau
sort unit year
save `"`OUT'/dta_s_stateyear.dta"', replace

* =====================================================================
* dta_s_biennial — P2 of the datasets audit (2026-09-12/13): biennial
* survey grid with a GLOBAL missing wave. 120 units x t = 1997, 1999, ...
* 2019 (spacing 2) with the 2009 wave dropped for every unit -> 11 waves,
* N = 1320. Onsets 2003 (i 1-30) and 2011 (i 31-60) + 60 never; sg 40/80
* [mod(i,3)==0]; w = .5 + gamma(2,.5) analytic weight; y = alpha_i +
* .04(t-1997) + {1.2, .7}*post + N(0,.6); ycount = 1 + poisson(exp(-.5 +
* .3 post)) - zero-free ON PURPOSE: with unit-level i every treated cell
* is a single observation and HEAD's separation screen refuses a zero
* count there (patch Run 3, 2026-09-13); the count's job here is the
* poisson path on a spacing-2 calendar, not count realism (see P1).
* Facts the bien_* configs pin (audit probe_p1p2 / probe2, HEAD rc 0):
* relative year -1 exists for NO cohort (the default reference is
* off-grid; only even horizons exist, so esprelength(2) espostlength(2)
* shows pre2 + post0 only); cohort 2011's relative year -2 is the 2009
* hole, so esfixedbaseperiod esrelativeto(-2) legitimately drops that
* cohort (N 990); anticipation(1) lands off-grid (caution and run);
* att-only is unbiased across seeds. Block verbatim from the probe.
* =====================================================================
clear
set type double
set seed 19972019
set sortseed 19972019
set obs 120
gen i = _n
gen ttre = cond(i <= 30, 2003, cond(i <= 60, 2011, .))
gen sg = 1 + (mod(i, 3) == 0)
gen alpha = rnormal()
expand 12
bysort i: gen t = 1995 + 2*_n
drop if t == 2009
gen w = .5 + rgamma(2, .5)
gen post = (t >= ttre & ttre < .)
gen y = alpha + .04*(t - 1997) + cond(ttre == 2003, 1.2, .7)*post + rnormal(0, .6)
gen ycount = 1 + rpoisson(exp(-.5 + .3*post))
drop post alpha
save `"`OUT'/dta_s_biennial.dta"', replace

* =====================================================================
* SELF-TEST — hard deterministic asserts (no tolerances). Reloads every
* dataset and checks N, unit counts, and the trigger-relevant properties.
* =====================================================================
di as result "=== mkdata_suite self-test ==="

* ---- dta_s_extmiss
use `"`OUT'/dta_s_extmiss.dta"', clear
assert _N == 960
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 80
drop __f
qui count if y == .a
assert r(N) == 24
qui count if t == .a
assert r(N) == 20
qui count if w == .a
assert r(N) == 23
* injections hit both cohorts and the never pool
foreach v in y t w {
  qui count if `v' == .a & ttre == 2005
  assert r(N) > 0
  qui count if `v' == .a & ttre == 2008
  assert r(N) > 0
  qui count if `v' == .a & ttre >= .
  assert r(N) > 0
}
bys i: gen byte __f2 = _n == 1
qui count if __f2 & ttre == 2005
assert r(N) == 25
qui count if __f2 & ttre == 2008
assert r(N) == 25
drop __f2
di as result "extmiss OK"

* ---- dta_s_deeppre
use `"`OUT'/dta_s_deeppre.dta"', clear
assert _N == 840
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 60
qui count if __f & ttre == 8
assert r(N) == 20
qui count if __f & ttre == 10
assert r(N) == 20
qui count if __f & ttre >= .
assert r(N) == 20
drop __f
qui sum t
assert r(min) == 1 & r(max) == 14
assert inlist(sg, 0, 1)
assert !missing(y)
di as result "deeppre OK"

* ---- dta_s_ttyear
use `"`OUT'/dta_s_ttyear.dta"', clear
assert _N == 520
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 40
qui count if __f & ttre == 2006
assert r(N) == 12
qui count if __f & ttre == 2009
assert r(N) == 12
qui count if __f & ttre >= .
assert r(N) == 16
drop __f
qui sum t
assert r(min) == 2000 & r(max) == 2012
assert ycount >= 1 & ycount == round(ycount)
di as result "ttyear OK"

* ---- dta_s_ttbig
use `"`OUT'/dta_s_ttbig.dta"', clear
assert _N == 520
qui sum t
assert r(min) == 500000 & r(max) == 500012
assert inlist(ttre, 500006, 500009) | ttre >= .
bys i: gen byte __f = _n == 1
qui count if __f & ttre == 500006
assert r(N) == 12
drop __f
di as result "ttbig OK"

* ---- dta_s_bigcodes
use `"`OUT'/dta_s_bigcodes.dta"', clear
assert _N == 480
assert i == round(i) & i > 2^24 & i < 2^31
* float-pooling premise: casting IDs to float collapses distinct units
gen double __fi = float(i)
bys __fi: gen byte __ff = _n == 1
qui count if __ff
assert r(N) < 40
drop __fi __ff
* rank-recode bijection: ascending i <-> i_small
sort i t
egen double __gi = group(i)
assert __gi == i_small
drop __gi
* dense uniform YYYYMMDD grid, on-grid onsets
assert mod(t - 20200101, 100) == 0 & t >= 20200101 & t <= 20201201
assert inlist(ttre, 20200601, 20200901) | ttre >= .
assert (ttre >= .) == (ttre_small >= .)
* identical treatment classification under both codings
assert (t >= ttre & ttre < .) == (t_small >= ttre_small & ttre_small < .)
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 40
drop __f
di as result "bigcodes OK"

* ---- dta_s_imdegen
use `"`OUT'/dta_s_imdegen.dta"', clear
assert _N == 3600
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 300
drop __f
assert !missing(y) & !missing(ypart)
* outcomes are EXACTLY the documented functions (bit-exact double recompute)
assert ypart == 1 + 0.5*cohgrp + 0.1*t + 2*(t >= ttre & ttre < .)
assert y == ypart + 0.5*((t == ttre - 1 | t == ttre - 2) & ttre < .)
qui count if y != ypart
assert r(N) == 400
* 12 clusters, each containing all three cohort groups
bys clus cohgrp: gen byte __f1 = _n == 1
bys clus: egen __ng = total(__f1)
assert __ng == 3
bys clus: gen byte __f2 = _n == 1
qui count if __f2
assert r(N) == 12
drop __f1 __ng __f2
di as result "imdegen OK"

* ---- dta_s_mintsent
use `"`OUT'/dta_s_mintsent.dta"', clear
assert _N == 300
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 30
qui count if __f & ttre_good == 5
assert r(N) == 10
drop __f
* onset-only-post coding on ALL treated units
assert ttre >= . if ttre_good == 5 & t < 5
assert ttre == 5 if ttre_good == 5 & t >= 5
assert ttre >= . if ttre_good >= .
* sentinel coding on unit 1 only
qui count if ttre_sent == 0
assert r(N) == 4
assert ttre_sent == ttre_good if !(i == 1 & t < 5)
di as result "mintsent OK"

* ---- dta_s_covars
use `"`OUT'/dta_s_covars.dta"', clear
assert _N == 1800
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 150
drop __f
qui count if x1 >= .
assert r(N) == 60
qui count if x2 >= .
assert r(N) == 56
assert !missing(y)
foreach l in 101 102 103 {
  qui count if fecat == `l'
  assert r(N) == 1
}
* singleton fe level 102 sits on a treated POST row
assert ttre < . & t >= ttre if fecat == 102
assert icbad == 1000000 + 137*i + t
assert icbad > 1e6
assert cwmvar > 0
assert inlist(cccvar, 1, 2)
di as result "covars OK"

* ---- dta_s_strclus
use `"`OUT'/dta_s_strclus.dta"', clear
assert _N == 240
confirm string variable strc
assert strc == "c" + string(g)
bys i: gen byte __f = _n == 1
qui count if __f & ttre == 5
assert r(N) == 8
drop __f
di as result "strclus OK"

* ---- dta_s_singleton
use `"`OUT'/dta_s_singleton.dta"', clear
assert _N == 352
* exactly four single-observation units, all never-treated
bys i: gen __cnt = _N
bys i: gen byte __f = _n == 1
qui count if __f & __cnt == 1
assert r(N) == 4
assert ttre >= . & inlist(i, 37, 38, 39, 40) & t == 5 if __cnt == 1
* every treated unit keeps at least 2 pre rows
bys i: egen __npre = total(t < ttre)
assert __npre >= 2 if ttre < .
drop __cnt __f __npre
* the three singleton fe levels: one row each; 91/92 on treated POST rows
foreach l in 91 92 93 {
  qui count if fevar == `l'
  assert r(N) == 1
}
assert ttre < . & t >= ttre if inlist(fevar, 91, 92)
assert ttre >= . if fevar == 93
di as result "singleton OK"

* ---- dta_s_stateyear (P1, datasets audit; every constant measured on the
* generator block as shipped, 2026-09-13)
use `"`OUT'/dta_s_stateyear.dta"', clear
assert _N == 2050
bys unit: gen byte __f = _n == 1
qui count if __f
assert r(N) == 106
* six state-level onset cohorts with fixed county counts, 28 never-treated
qui count if __f & ttre == 1999
assert r(N) == 17
qui count if __f & ttre == 2002
assert r(N) == 21
qui count if __f & ttre == 2005
assert r(N) == 6
qui count if __f & ttre == 2008
assert r(N) == 19
qui count if __f & ttre == 2011
assert r(N) == 5
qui count if __f & ttre == 2013
assert r(N) == 10
qui count if __f & ttre >= .
assert r(N) == 28
* onset is a state attribute; unit = st*100 + cty; region = 1 + mod(st - 1, 3)
* and every region holds exactly one never-treated state
bys st (ttre): gen byte __same = ttre[1] == ttre[_N]
assert __same
assert unit == st*100 + cty
assert region == 1 + mod(st - 1, 3)
forval r = 1/3 {
  qui count if __f & ttre >= . & region == `r'
  assert r(N) > 0
}
* one row per unit-year on the 1995-2014 calendar
bys unit year: gen long __c = _N
assert __c == 1
qui sum year
assert r(min) == 1995 & r(max) == 2014
* 20/80 subgroup: 22 of the 106 counties are sg == 2
assert inlist(sg, 1, 2)
qui count if __f & sg == 2
assert r(N) == 22
* imbalance: 8 late-entry and 4 early-exit counties; 33 missing y (mod()
* injection); y_nosh / y_dose / ycount / pop / dose_it are complete
bys unit (year): gen byte __late = _n == 1 & year > 1995
bys unit (year): gen byte __early = _n == _N & year < 2014
qui count if __late
assert r(N) == 8
qui count if __early
assert r(N) == 4
qui count if y == .
assert r(N) == 33
assert !missing(y_nosh) & !missing(y_dose) & !missing(ycount) & !missing(pop) & !missing(dose_it)
* count outcome: integer-valued with 403 zeros (19.7%) and no all-zero
* st x year cell (poisson at i = st runs under default safety)
assert ycount >= 0 & ycount == round(ycount)
qui count if ycount == 0
assert r(N) == 403
bys st year: egen double __mx = max(ycount)
assert __mx > 0
* weights bounded below; post-drop st x year cells hold 4-12 counties and
* dose_it varies within every one of them (the contreat identification bed)
assert pop >= 500
bys st year: gen long __nc = _N
assert __nc >= 4 & __nc <= 12
bys st year: egen double __sd = sd(dose_it)
assert __sd > 0 & __sd < .
drop __f __same __c __late __early __nc __sd __mx
di as result "stateyear OK"

* ---- dta_s_biennial (P2, datasets audit 2026-09-13)
use `"`OUT'/dta_s_biennial.dta"', clear
assert _N == 1320
bys i: gen byte __f = _n == 1
qui count if __f
assert r(N) == 120
qui count if __f & ttre == 2003
assert r(N) == 30
qui count if __f & ttre == 2011
assert r(N) == 30
qui count if __f & ttre >= .
assert r(N) == 60
qui count if __f & sg == 2
assert r(N) == 40
drop __f
* biennial calendar 1997..2019 (odd years only): 11 waves, 2009 missing for all
qui sum t
assert r(min) == 1997 & r(max) == 2019
assert mod(t, 2) == 1
qui tab t
assert r(r) == 11
qui count if t == 2009
assert r(N) == 0
bys i t: gen long __c = _N
assert __c == 1
drop __c
* on-grid onsets: every relative year is even, relative year -1 exists for
* no cohort, and cohort 2011's relative year -2 (2009) is the hole
assert mod(t - ttre, 2) == 0 if ttre < .
qui count if t - ttre == -1
assert r(N) == 0
qui count if ttre == 2011 & t == ttre - 2
assert r(N) == 0
qui count if ttre == 2003 & t == ttre - 2
assert r(N) == 30
assert w > .5 & !missing(y)
* the count is zero-free by construction (unit-level i + a zero singleton
* treated cell = separation refusal on HEAD; see the generator comment)
assert ycount >= 1 & ycount == round(ycount)
di as result "biennial OK"

di as result "MKDATA_SUITE_SELFTEST_PASS"
exit, clear
