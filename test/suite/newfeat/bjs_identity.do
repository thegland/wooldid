* bjs_identity.do - external oracle for wooldid's POOLED-REFERENCE event study,
*   checked against did_imputation (Borusyak-Jaravel-Spiess).  Scope is deliberately
*   narrow: ONLY the configuration in which equivalence can reasonably be expected -
*   esprelength(0) (a pooled reference period, no pre horizons) on a BALANCED panel,
*   one observation per i-t cell, binary staggered treatment, no covariates, no
*   weights, cluster(i).  One common-onset fixture, one staggered fixture, no more.
*
* THE MAPPING (wooldid.sthlp, "Comparing Output to did_imputation"):
*     wooldid y i t ttre, att cluster(i) esprelength(0) espostlength(H+1)
*     did_imputation y i t ttre, horizons(0/H) minn(0) cluster(i)
*   BJS tau<h> <-> wooldid ES_att_post<h>; BJS's aggregate tau (from a run with NO
*   horizons()) <-> wooldid main_att.  Paired BY NAME, never by position.
*
* CHECKS: B1-b/B1-att the identity (TOLERANCE 1e-6), B4-neg the negative control that
*   keeps B1 from going vacuous (INVERTED, must exceed 1e-4), B2-lab/B2-b/B2-V the
*   anticipation identity (EXACT, gate 0).  Every rule and the measurement behind it
*   is stated where the check is made.
*
* OPTIONAL DEPENDENCY.  If did_imputation is absent, B1 and B4 are SKIPPED and
*   UNCOUNTED, the wooldid-internal B2 legs still run and are counted, and the file
*   STILL reaches ALL_PASS (suite_log_pass grades a SKIPPED-only exit as FAIL).  Only
*   wooldid's own calls are trapped: if did_imputation is present but errors on a
*   frozen fixture the file aborts with its r() and no ALL_PASS line, which the
*   harness grades FAIL.
*
* Args: 1 = path to wooldid.ado (file or its directory); 2 = <repo>/test/data.
* Run:  "<stata>" -e do bjs_identity.do <adopath> <datadir>
*    -> expect BJS_IDENTITY_ALL_PASS
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
* TOL is the equivalence gate and NEGFLOOR the floor B4 must exceed; both are read as
* mreldif, i.e. max |a-b|/(1+|b|) - see mrchk.
global TOL 1e-6
global NEGFLOOR 1e-4

if `"`1'"' == "" | `"`2'"' == "" {
  dis "RES|B0|args|arg1 wooldid.ado and arg2 <repo>/test/data are required|FAIL"
  dis "BJS_IDENTITY_FAILS=1"
  exit, clear
}
local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
qui do `"`ADO'"'
cap which did_imputation
local HASBJS = (_rc == 0)
if !`HASBJS' {
  dis "BJS_IDENTITY_SKIPPED: did_imputation is not installed; B1 (the identity) and"
  dis "  B4 (the negative control) are SKIPPED and UNCOUNTED; B2 still runs."
}

* Both fixtures come from the repo's dta_binary (80 units, t=2000..2011, onsets 2005
* and 2008, 30 never-treated, balanced, one observation per i-t cell): the staggered
* panel itself, and the common-onset panel left by dropping the 2008 cohort.  Top post
* horizon is 6 on both (2011-2005).  Paths travel as GLOBALS so a directory with a
* space in it survives the trip through args.
global FX_stag `"`2'/dta_binary.dta"'
use `"$FX_stag"', clear
qui drop if ttre == 2008
tempfile fco
qui save `fco'
global FX_common `fco'

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
* r(miss) so the caller fails loudly - a bare _b[tau6] on a missing column is r(111).
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
* the right form because event-study coefficients pass through zero, where a bare
* relative difference is meaningless.  mode tol: mreldif <= tol (tol 0 means EXACT).
* mode neg: mreldif must EXCEED tol.  A failure lists both matrices.
cap program drop mrchk
program define mrchk
  args id item mode A B tol
  global nchk = $nchk + 1
  local d = mreldif(`A', `B')
  local ok = cond("`mode'" == "neg", `d' > `tol', `d' <= `tol')
  dis "RES|`id'|`item'|mreldif " %9.2e `d' " vs gate " %9.2e `tol' " (`mode')|" cond(`ok', "PASS", "FAIL")
  if !`ok' global pf = $pf + 1
  matrix MRD = `A' \ `B'
  if !`ok' matlist MRD, format(%21.0g) title("`id' `item': `A' over `B'")
end

* ---- B1 (the identity) and B4 (the negative control), per fixture ----
* args: id fx H hasbjs neg    (neg=1 also runs B4 on this fixture)
cap program drop do_fixture
program define do_fixture
  args id fx H hasbjs neg
  local Q = `H' + 1
  local wn
  local bn
  forval h = 0/`H' {
    local wn `wn' ES_att_post`h'
    local bn `bn' tau`h'
  }
  use "${FX_`fx'}", clear
  qui sum y
  local ymax = max(abs(r(min)), abs(r(max)))
  cap noi wooldid y i t ttre, att cluster(i) esprelength(0) espostlength(`Q')
  local rc = _rc
  if `rc' bfail "`id'" wooldid-rc "esprelength(0) espostlength(`Q') rc=`rc'"
  if `rc' exit
  matrix W0 = e(b)
  dis "RES|`id'|fit|N=" e(N) " max|y|=" %6.2f `ymax' " feps*max|y|=" %9.2e 1.1921e-07*`ymax' "|INFO"
  if !`hasbjs' dis "RES|`id'|B1-b,B1-att,B4-neg|did_imputation absent|SKIPPED (uncounted)"
  if !`hasbjs' exit
  use "${FX_`fx'}", clear
  did_imputation y i t ttre, horizons(0/`H') minn(0) cluster(i)
  matrix BH = e(b)
  use "${FX_`fx'}", clear
  did_imputation y i t ttre, minn(0) cluster(i)
  matrix BA = e(b)
  * B1-b: TOLERANCE, mreldif <= 1e-6 over every post horizon (|w-bjs| <= 1e-6 +
  * 1e-6*|bjs| elementwise).  NOT a solver tolerance - wooldid at regtol(13) against
  * did_imputation at tol(1e-13) maxit(10000) leaves the gap unmoved.  It is
  * did_imputation's storage type: it builds the imputed counterfactual and the
  * per-observation effect with untyped gen's (did_imputation.ado:354, :393), single
  * precision under the default `set type float'.  That rounding lands on the
  * COUNTERFACTUAL, magnitude ~|y| and not ~|b|, so the per-observation ceiling is
  * float-eps*max|y| (printed on the fit line above) and averaging those effects into
  * a coefficient buys the headroom: measured worst 4.5e-08, 22x inside the gate.  An
  * outcome in the hundreds would need the gate rescaled - a fixture-scale change, not
  * a wooldid regression.  Do NOT tighten this under `set type double': the first
  * double-precision call in a session differs bitwise from every later one on
  * identical data, so that gate would not be deterministic.
  pairup W0 "`wn'" BH "`bn'"
  if "`r(miss)'" != "" bfail "`id'" B1-b "columns absent: `r(miss)'"
  else mrchk "`id'" B1-b tol PA PB $TOL
  * B1-att: TOLERANCE 1e-6, the same gate on the aggregate.
  pairup W0 main_att BA tau
  if "`r(miss)'" != "" bfail "`id'" B1-att "columns absent: `r(miss)'"
  else mrchk "`id'" B1-att tol PA PB $TOL
  if !`neg' exit
  use "${FX_`fx'}", clear
  cap noi wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(`Q')
  local rc = _rc
  if `rc' bfail "`id'" B4-neg "wooldid esprelength(4) rc=`rc' (control not evaluable)"
  if `rc' exit
  matrix W4 = e(b)
  * B4-neg: INVERTED, mreldif must EXCEED 1e-4.  The help's condition list is
  * incomplete: a pooled reference period is necessary but not sufficient, because the
  * reference pool is the set of horizons the run ASKS FOR.  At esprelength(4) - the
  * natural way to ask for an event study - the pre periods leave the pool and every
  * post coefficient shifts, while did_imputation pretrends(4) moves its taus by
  * exactly 0; without this a future change could quietly make B1 vacuous.  The floor
  * is 100x above the equivalence gate and ~460x below the divergence measured here
  * (mreldif 4.65e-02).  Do NOT loosen it - re-read the help's "uses a pooled
  * reference period" sentence instead.
  pairup W4 "`wn'" BH "`bn'"
  if "`r(miss)'" != "" bfail "`id'" B4-neg "columns absent: `r(miss)'"
  else mrchk "`id'" B4-neg neg PA PB $NEGFLOOR
end

do_fixture B1-common common 6 `HASBJS' 1
do_fixture B1-stag   stag   6 `HASBJS' 0

* ---- B2: the anticipation identity (wooldid only, so it runs on the skip path) ----
* EXACT, gate 0.  anticipation(#) is documented as shifting the treatment date by #
* periods, so it must be the SAME regression as wooldid on a manually shifted ttre-#:
* identical labels, coefficients and variances, bit for bit (measured equal at 21
* significant figures).  A tolerance would hide a real behavioural break, so there is
* none.  B2-lab runs first: equal labels are what make B2-b and B2-V conformable.
cap program drop ant_exact
program define ant_exact
  args id fx K Q
  use "${FX_`fx'}", clear
  cap noi wooldid y i t ttre, att cluster(i) esprelength(0) espostlength(`Q') anticipation(`K')
  local rc = _rc
  if `rc' bfail "`id'" wooldid-rc "anticipation(`K') rc=`rc'"
  if `rc' exit
  matrix AB = e(b)
  matrix AV = e(V)
  use "${FX_`fx'}", clear
  qui gen double ttre_sh = ttre - `K'
  cap noi wooldid y i t ttre_sh, att cluster(i) esprelength(0) espostlength(`Q')
  local rc = _rc
  if `rc' bfail "`id'" wooldid-rc "manual ttre-`K' shift rc=`rc'"
  if `rc' exit
  matrix MB = e(b)
  matrix MV = e(V)
  local la : colnames AB
  local lm : colnames MB
  global nchk = $nchk + 1
  if "`la'" != "`lm'" global pf = $pf + 1
  dis "RES|`id'|B2-lab|`:word count `la'' vs `:word count `lm'' posted columns|" cond("`la'" == "`lm'", "PASS", "FAIL")
  if "`la'" != "`lm'" exit
  mrchk "`id'" B2-b tol AB MB 0
  mrchk "`id'" B2-V tol AV MV 0
end

ant_exact B2-ant1 common 1 7

timer off 9
qui timer list 9
dis "RES|ALL|counted|$nchk checks, $pf failing, " %6.1f r(t9) " s, did_imputation " cond(`HASBJS', "present", "ABSENT") "|INFO"
if $pf == 0 dis "BJS_IDENTITY_ALL_PASS"
else dis "BJS_IDENTITY_FAILS=$pf"
exit, clear
