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

di as result "MKDATA_SUITE_SELFTEST_PASS"
exit, clear
