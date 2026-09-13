* xtevent_identity.do - external oracle for wooldid's FIXED-BASE-PERIOD event study,
*   checked against xtevent (Freyaldenhoven-Hansen-Perez Perez-Shapiro) run in
*   impute() mode.  xtevent shares no code path with wooldid beyond the data and
*   Stata's own OLS: wooldid fits a saturated cohort-by-period cell regression and
*   aggregates cells to relative horizons, while xtevent fits ONE pooled coefficient
*   per relative period on a lead/lag dummy set with unit and time fixed effects.  The
*   agreement below is therefore an AGGREGATION result, not a same-regression result,
*   and it pins cell construction, base-period logic, window/labelling and aggregation
*   at once.
*
* SCOPE - COMMON ONSET ONLY, on a balanced panel with no covariates, no weights,
*   cluster(i), one observation per i-t cell.  With a SINGLE treated cohort calendar
*   time and relative time are in bijection on the treated units, so each wooldid
*   horizon is one cell and comes out equal to xtevent's pooled per-relative-period
*   coefficient.  With TWO cohorts xtevent's OLS pools the cells under its own weights
*   rather than wooldid's and the two part company: measured on this very fixture left
*   UNDROPPED (onsets 2005 and 2008), mreldif 6.2e-02 at window(-4 5) and 6.3e-02 at
*   window(-3 4), and an earlier probe on a two-cohort panel whose cohorts share a
*   CONSTANT homogeneous effect still measured a 6.0e-02 ABSOLUTE gap - so the
*   divergence is estimator structure, not treatment-effect heterogeneity.  Staggered
*   adoption is consequently NOT tested here.  Neither are covariates, unbalanced
*   panels or weights: none of the three was ever measured (an 8%-thinned common-onset
*   panel came apart at a 2.4e-01 absolute gap).
*   xtevent's Sun-Abraham mode (cohort()/control_cohort()) requires the avar package,
*   which is NOT installed on this machine, so that mode cannot be tested here at all.
*
* STANDARD ERRORS ARE NOT GATED: the two packages' variance formulas genuinely differ
*   (measured SE ratio 1.045 at post0), so only the coefficients are an identity.  That
*   is a scope statement, not a loosened rule.
*
* THE MAPPING:
*     wooldid y i t ttre, att cluster(i) esprelength(P) espostlength(Q+1)
*                         esfixedbaseperiod              (default esrelativeto(-1))
*     gen byte z = (ttre < . & t >= ttre)
*     xtevent y, policyvar(z) panelvar(i) timevar(t) window(-P Q) impute(stag)
*                vce(cluster i)
*   NAME MAP: xtevent _k_eq_m<h> <-> wooldid pre<h>, _k_eq_p<h> <-> post<h>.  Both
*   sides normalise at k = -1 (wooldid's omitted "pre1"), so the labels line up with no
*   offset.  impute(stag) is what keeps the sample whole: in plain TWFE mode xtevent
*   collapses a common-onset panel to a two-period sample and drops every lead and lag
*   but _k_eq_p0 for collinearity.  (impute(nuchange) and impute(instag) returned
*   bit-identical vectors when probed, so the choice among the three is immaterial.)
*
* THE BINNED ENDPOINTS - the one trap in this mapping.  window(-P Q) returns
*   _k_eq_m(P+1).._k_eq_m2 and _k_eq_p0.._k_eq_p(Q+1), and the OUTERMOST lead and lag
*   are BINS, not horizons: _k_eq_m(P+1) means "rel <= -(P+1)" and _k_eq_p(Q+1) means
*   "rel >= Q+1".  wooldid's pre<h>/post<h> are always single horizons, so the two
*   endpoints are NOT comparable and are excluded from every gate below.  The exclusion
*   is ASYMMETRIC on purpose: espostlength(Q+1) stops wooldid at post<Q>, so xtevent's
*   _k_eq_p(Q+1) simply has no wooldid counterpart, whereas wooldid's pre(P+1) does
*   exist and is deliberately dropped.  Do not "fix" that asymmetry.
*   MEASURED SIGNATURE of getting this wrong, on this fixture at window(-3 4): the
*   binned _k_eq_m4 is 0.0133446 while wooldid pre4 is 0.0835880, a 7.0e-02 gap at
*   exactly the endpoint, while all seven interior horizons agree to 7.3e-15.  A
*   failure at exactly the two endpoints and nowhere else is this bug, not a wooldid
*   defect.  Both excluded endpoints are printed INFO and UNCOUNTED by every window
*   below so that signature stays visible in the log.
*
* THE TWO WINDOWS.  Relative time on this fixture runs -5..+6, so window(-4 5) is the
*   widest xtevent accepts, and it is widest on BOTH sides: window(-5 5), window(-4 6)
*   and window(-5 6) each ask for an endpoint bin no observation can fall in, and each
*   was MEASURED here to error r(198) "Window outside event-time range".  Anything
*   wider than these two windows is therefore not available, not merely unused.
*   At window(-4 5) BOTH bins hold a single
*   relative period, so the excluded lead coincidentally agrees to 2.1e-17; at
*   window(-3 4) each bin pools two periods and the excluded lead diverges by 7.0e-02.
*   One window where the exclusion is load-bearing and one where it is not.
*
* CHECKS: X-rc and X-lab are wooldid-only structure (EXACT); X-align and X-b are the
*   oracle (set equality, then TOLERANCE 1e-10).  Each rule is stated where it is made.
*
* OPTIONAL DEPENDENCY.  If xtevent is absent, X-align and X-b are SKIPPED and
*   UNCOUNTED, X-rc and X-lab still run and are counted, and the file STILL reaches
*   ALL_PASS (suite_log_pass grades a SKIPPED-only exit as FAIL).  Only wooldid's own
*   calls are trapped: if xtevent is present but errors on a frozen fixture the file
*   aborts with its r() and no ALL_PASS line, which the harness grades FAIL.
*
* Args: 1 = path to wooldid.ado (file or its directory); 2 = <repo>/test/data.
* Run:  "<stata>" -e do xtevent_identity.do <adopath> <datadir>
*    -> expect XTEVENT_IDENTITY_ALL_PASS
* Sentinels: RES|<leg>|<item>|<value>|<PASS|FAIL|INFO|SKIPPED>

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 20260913
set sortseed 20260913
timer clear
timer on 9
global pf 0
global nchk 0
* The equivalence gate, read as mreldif = max |a-b|/(1+|b|) - see mrchk.  The
* preregistered rule is the mixed form |b_w - b_xt| <= atol + rtol*|b_xt| with atol
* 1e-12 and rtol 1e-10; mreldif collapses the two constants into one, and 1e-10 (the
* rtol) is the one taken, which makes this gate identical to the rule on the relative
* term and looser by at most 1e-10 on the absolute term.  1e-12 was NOT taken: it is
* TIGHTER than the preregistered rule and the only thing it could buy is a false
* failure.  Measured worst over the gated legs here is mreldif 4.7e-15, i.e. ~4.3
* orders of headroom, while the two structural divergences measured on this fixture
* (6.2e-02 staggered, 7.0e-02 pooled lead bin) sit 8-9 orders ABOVE the gate.
global TOL 1e-10

if `"`1'"' == "" | `"`2'"' == "" {
  dis "RES|X0|args|arg1 wooldid.ado and arg2 <repo>/test/data are required|FAIL"
  dis "XTEVENT_IDENTITY_FAILS=1"
  exit, clear
}
local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
qui do `"`ADO'"'
cap which xtevent
local HASXT = (_rc == 0)
if !`HASXT' {
  dis "XTEVENT_IDENTITY_SKIPPED: xtevent is not installed; the oracle legs X-align and"
  dis "  X-b are SKIPPED and UNCOUNTED; X-rc and X-lab still run and are counted."
}

* The one fixture is the common-onset panel left by dropping the 2008 cohort from the
* repo's dta_binary (which is 80 units, t=2000..2011, onsets 2005 and 2008, 30
* never-treated, balanced, one observation per i-t cell).  What remains: 55 units, 660
* observations, a single onset at 2005, relative time -5..+6.  z is xtevent's policy
* variable.  The path travels as a GLOBAL so a directory with a space in it survives
* the trip through args.
use `"`2'/dta_binary.dta"', clear
qui drop if ttre == 2008
qui gen byte z = (ttre < . & t >= ttre)
tempfile fco
qui save `fco'
global FX `fco'

* ---- helpers ----
cap program drop bfail
program define bfail
  args id item msg
  dis "RES|`id'|`item'|`msg'|FAIL"
  global pf = $pf + 1
  global nchk = $nchk + 1
end

* pairup: copy named columns of AM and BM into conformable 1xK vectors PA and PB,
* paired BY NAME, never by position.  A name absent on either side comes back in
* r(miss) so the caller fails loudly - a bare _b[_k_eq_p4] on a missing column is
* r(111).
cap program drop pairup
program define pairup, rclass
  args AM ANM BM BNM
  local K : word count `ANM'
  matrix PA = J(1, `K', .)
  matrix PB = J(1, `K', .)
  local miss
  forval j = 1/`K' {
    local an : word `j' of `ANM'
    local bn : word `j' of `BNM'
    local ca = colnumb(`AM', "`an'")
    local cb = colnumb(`BM', "`bn'")
    if `ca' >= . | `cb' >= . local miss `miss' `an'/`bn'
    else {
      matrix PA[1,`j'] = `AM'[1,`ca']
      matrix PB[1,`j'] = `BM'[1,`cb']
    }
  }
  return local miss `miss'
end

* mrchk: one counted check on two conformable matrices.  mreldif(A,B) is
* max |a-b|/(1+|b|), i.e. exactly the mixed gate |a-b| <= tol + tol*|b| elementwise -
* the right form because pre-period coefficients pass through zero, where a bare
* relative difference is meaningless (the worst RELATIVE gap measured here is 1.3e-12,
* on a pre coefficient of 0.0032, whose absolute gap is only 4.2e-15).  A failure lists
* both matrices.
cap program drop mrchk
program define mrchk
  args id item A B tol
  global nchk = $nchk + 1
  local d = mreldif(`A', `B')
  local ok = (`d' <= `tol')
  dis "RES|`id'|`item'|mreldif " %9.2e `d' " vs gate " %9.2e `tol' "|" cond(`ok', "PASS", "FAIL")
  if !`ok' global pf = $pf + 1
  matrix MRD = `A' \ `B'
  if !`ok' matlist MRD, format(%21.0g) title("`id' `item': `A' over `B'")
end

* ---- one window: wooldid structure, then the xtevent oracle ----
* args: id P Q hasxt    (P/Q are xtevent's window(-P Q); wooldid gets esprelength(P)
* and espostlength(Q+1), which is what makes wooldid's top posted horizon post<Q>
* coincide with xtevent's top NON-binned lag _k_eq_p<Q>.)
cap program drop do_window
program define do_window
  args id P Q hasxt
  local EQ = `Q' + 1
  use "$FX", clear
  cap noi wooldid y i t ttre, att cluster(i) esprelength(`P') espostlength(`EQ') esfixedbaseperiod
  local rc = _rc
  if `rc' bfail "`id'" X-rc "wooldid rc=`rc'"
  if `rc' exit
  cap confirm matrix e(b_es_att)
  local rc = _rc
  if `rc' bfail "`id'" X-rc "e(b_es_att) not posted"
  if `rc' exit
  * X-rc: EXACT.  A nonzero return code, or a missing e(b_es_att), on a frozen fixture
  * is never acceptable - and e(b_es_att) is what the oracle comparison below reads.
  global nchk = $nchk + 1
  dis "RES|`id'|X-rc|rc 0, e(b_es_att) posted, N=" e(N) "|PASS"
  matrix WB = e(b_es_att)
  local np = e(numpre_att)
  local got : colnames WB
  * X-lab: EXACT string/integer equality.  esprelength(P) must report pre(P+1)..pre2
  * (the reference pre1 is omitted), espostlength(Q+1) must report post0..post(Q), and
  * e(numpre_att) must be P.  The whole _k_eq_m/_k_eq_p name map AND the endpoint
  * exclusion are built on this contract, so it is asserted before it is used.
  local want
  forval h = `=`P'+1'(-1)2 {
    local want `want' pre`h'
  }
  forval h = 0/`Q' {
    local want `want' post`h'
  }
  global nchk = $nchk + 1
  if "`got'" != "`want'" | "`np'" != "`P'" global pf = $pf + 1
  dis "RES|`id'|X-lab|[`got'] numpre_att=`np'|" cond("`got'" == "`want'" & "`np'" == "`P'", "PASS", "FAIL")
  if "`got'" != "`want'" exit
  if !`hasxt' dis "RES|`id'|X-align,X-b|xtevent absent|SKIPPED (uncounted)"
  if !`hasxt' exit
  use "$FX", clear
  qui xtset i t
  xtevent y, policyvar(z) panelvar(i) timevar(t) window(-`P' `Q') impute(stag) vce(cluster i)
  matrix XB = e(b)
  * The compared set: pre(P)..pre2 and post0..post(Q), i.e. every horizon EXCEPT the two
  * binned endpoints _k_eq_m(P+1) and _k_eq_p(Q+1).  K = P + Q.
  local wn
  local xn
  forval h = `P'(-1)2 {
    local wn `wn' pre`h'
    local xn `xn' _k_eq_m`h'
  }
  forval h = 0/`Q' {
    local wn `wn' post`h'
    local xn `xn' _k_eq_p`h'
  }
  local K : word count `wn'
  local maxm = 0
  local maxp = -1
  local cl : colnames XB
  foreach nm of local cl {
    if substr("`nm'", 1, 7) == "_k_eq_m" & real(substr("`nm'", 8, .)) > `maxm' local maxm = real(substr("`nm'", 8, .))
    if substr("`nm'", 1, 7) == "_k_eq_p" & real(substr("`nm'", 8, .)) > `maxp' local maxp = real(substr("`nm'", 8, .))
  }
  pairup WB "`wn'" XB "`xn'"
  local miss `r(miss)'
  * X-align: EXACT.  All K paired names must exist on both sides, AND xtevent's
  * outermost lead must be exactly _k_eq_m(P+1) and its outermost lag exactly
  * _k_eq_p(Q+1).  That second half is what protects the exclusion: if xtevent ever
  * stopped binning its endpoints, or shifted them, this fails with a legible message
  * instead of X-b quietly comparing a bin against a horizon.
  local aok = ("`miss'" == "" & `maxm' == `=`P'+1' & `maxp' == `=`Q'+1')
  global nchk = $nchk + 1
  if !`aok' global pf = $pf + 1
  dis "RES|`id'|X-align|`K' paired, xtevent outermost m`maxm'/p`maxp' (want m`=`P'+1'/p`=`Q'+1'), missing [`miss']|" cond(`aok', "PASS", "FAIL")
  if !`aok' exit
  * X-b: TOLERANCE, mreldif <= 1e-10 over every non-binned horizon.  On a single treated
  * cohort the two estimators are algebraically the same object, so the only gap is
  * floating-point ordering: measured 3.3e-15 at window(-4 5) and 4.7e-15 at
  * window(-3 4) on this fixture (worst ABSOLUTE gap 7.3e-15), and worst absolute gaps
  * of 1.1e-14 and 2.1e-14 on two other balanced common-onset fixtures probed.  See
  * $TOL above for why 1e-10 and not 1e-12.
  mrchk "`id'" X-b PA PB $TOL
  * The two excluded endpoints, printed UNCOUNTED, so the binned-endpoint signature
  * described in the header stays visible whether or not anything failed.
  local cw = colnumb(WB, "pre`=`P'+1'")
  local cm = colnumb(XB, "_k_eq_m`=`P'+1'")
  if `cw' < . & `cm' < . dis "RES|`id'|X-bin-lead|wooldid pre`=`P'+1' vs binned _k_eq_m`=`P'+1', |d|=" %9.2e abs(WB[1,`cw'] - XB[1,`cm']) " (excluded)|INFO"
  local cp = colnumb(XB, "_k_eq_p`=`Q'+1'")
  if `cp' < . dis "RES|`id'|X-bin-lag|binned _k_eq_p`=`Q'+1'=" %9.2e XB[1,`cp'] ", no wooldid counterpart (excluded)|INFO"
end

* window(-4 5) is the widest this fixture supports and puts a single period in each
* bin; window(-3 4) pools two periods into each bin.  9 and 7 gated horizons.
do_window X-w45 4 5 `HASXT'
do_window X-w34 3 4 `HASXT'

timer off 9
qui timer list 9
dis "RES|ALL|counted|$nchk checks, $pf failing, " %6.1f r(t9) " s, xtevent " cond(`HASXT', "present", "ABSENT") "|INFO"
if $pf == 0 dis "XTEVENT_IDENTITY_ALL_PASS"
else dis "XTEVENT_IDENTITY_FAILS=$pf"
exit, clear
