* csdid_identity.do - external oracle for wooldid's FIXED-BASE-PERIOD event study,
*   checked against Callaway-Sant'Anna as implemented by Rios-Avila's csdid.  csdid
*   shares no code path with wooldid beyond the data and Stata's own OLS, so agreement
*   pins cell construction, control-group selection, base-period logic, aggregation
*   and ES labelling at once.
*
* SCOPE.  The identity holds ONLY on BALANCED panels with NO covariates, one
*   observation per i-t cell, binary staggered treatment, no weights, cluster(i): one
*   common-onset fixture and one staggered fixture, nothing else.  On a realistic
*   IMBALANCED administrative panel the same comparison degrades to qualitative
*   agreement (measured: up to 0.42 of a standard error apart on a near-zero pre
*   coefficient), so that case is deliberately not tested here.
*
* THE MAPPING:
*     wooldid y i t ttre, att cluster(i) esprelength(P) espostlength(Q)
*                         esfixedbaseperiod esrelativeto(-1)
*     csdid  y, ivar(i) time(t) gvar(g) long2   with g = cond(ttre<., ttre, 0)
*            then "estat event" and "estat simple"
*   csdid's long2 makes every 2x2 a comparison against g-1, which is what
*   esfixedbaseperiod with esrelativeto(-1) does - hence an identity, not a
*   similarity.  NAME MAP: csdid Tm<h> <-> wooldid pre<h>, Tp<h> <-> post<h>, "estat
*   simple" ATT <-> wooldid main_att.  Neither side emits a rel = -1 record (it is the
*   reference), so the label sets line up with no offset.  Horizons are paired BY NAME
*   and C-align asserts the two sets are equal - misaligning the index is the one way
*   to manufacture a false failure here.  No "notyet" option: on these fixtures
*   never-treated and not-yet-treated controls coincide numerically.
*
* STANDARD ERRORS ARE NOT GATED: the two tools' variance formulas genuinely differ on
*   panels this small (measured cross-tool SE ratios 0.80-1.74 at 40 units), so only
*   the coefficients and the ATT are an identity.  A scope statement, not a loosened
*   rule.
*
* CHECKS: S-rc and S-lab are wooldid-only structure (EXACT); C-align, C-b and C-att
*   are the oracle (set equality, then TOLERANCE 1e-10).  Each rule is stated where
*   the check is made.
*
* OPTIONAL DEPENDENCY.  If csdid or csdid_estat is absent, C-align/C-b/C-att are
*   SKIPPED and UNCOUNTED, S-rc and S-lab still run and are counted, and the file
*   STILL reaches ALL_PASS (suite_log_pass grades a SKIPPED-only exit as FAIL).  Only
*   wooldid's own calls are trapped: if csdid is present but errors on a frozen
*   fixture the file aborts with its r() and no ALL_PASS line, which the harness
*   grades FAIL.
*
* Args: 1 = path to wooldid.ado (file or its directory); 2 = <repo>/test/data.
* Run:  "<stata>" -e do csdid_identity.do <adopath> <datadir>
*    -> expect CSDID_IDENTITY_ALL_PASS
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
* the equivalence gate, read as mreldif = max |a-b|/(1+|b|) - see mrchk
global TOL 1e-10

if `"`1'"' == "" | `"`2'"' == "" {
  dis "RES|C0|args|arg1 wooldid.ado and arg2 <repo>/test/data are required|FAIL"
  dis "CSDID_IDENTITY_FAILS=1"
  exit, clear
}
local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
qui do `"`ADO'"'
local HASCS = 1
cap which csdid
if _rc local HASCS = 0
cap which csdid_estat
if _rc local HASCS = 0
if !`HASCS' {
  dis "CSDID_IDENTITY_SKIPPED: csdid or csdid_estat is not installed; the oracle legs"
  dis "  C-align, C-b and C-att are SKIPPED and UNCOUNTED; S-rc and S-lab still run."
}

* Both fixtures come from the repo's dta_binary (80 units, t=2000..2011, onsets 2005
* and 2008, 30 never-treated, balanced, one observation per i-t cell): the staggered
* panel itself, and the common-onset panel left by dropping the 2008 cohort.  g is
* csdid's cohort variable (0 = never treated).  Paths travel as GLOBALS so a directory
* with a space in it survives the trip through args.
use `"`2'/dta_binary.dta"', clear
qui gen double g = cond(ttre < ., ttre, 0)
tempfile fstag fcommon
qui save `fstag'
global FX_stag `fstag'
qui drop if ttre == 2008
qui save `fcommon'
global FX_common `fcommon'

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
* r(miss) so the caller fails loudly - a bare _b[Tp6] on a missing column is r(111).
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
* relative difference is meaningless.  A failure lists both matrices.
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

* ---- one fixture: wooldid structure, then the csdid oracle ----
* args: id fx P Q hascs
cap program drop do_fixture
program define do_fixture
  args id fx P Q hascs
  use "${FX_`fx'}", clear
  cap noi wooldid y i t ttre, att cluster(i) esprelength(`P') espostlength(`Q') esfixedbaseperiod esrelativeto(-1)
  local rc = _rc
  if `rc' bfail "`id'" S-rc "wooldid rc=`rc'"
  if `rc' exit
  cap confirm matrix e(b_es_att)
  local rc = _rc
  if `rc' bfail "`id'" S-rc "e(b_es_att) not posted"
  if `rc' exit
  * S-rc: EXACT.  A nonzero return code, or a missing e(b_es_att), on a frozen fixture
  * is never acceptable - and it is what the oracle comparison below reads.
  global nchk = $nchk + 1
  dis "RES|`id'|S-rc|rc 0, e(b_es_att) posted, N=" e(N) "|PASS"
  matrix WB = e(b_es_att)
  matrix WEB = e(b)
  local np = e(numpre_att)
  local got : colnames WB
  * S-lab: EXACT string/integer equality.  esprelength(P) must report pre(P+1)..pre2
  * (the reference pre1 is omitted), espostlength(Q) must report post0..post(Q-1), and
  * e(numpre_att) must be P.  The whole Tm/Tp name map is built on this contract, so it
  * is asserted before it is used.
  local want
  forval h = `=`P'+1'(-1)2 {
    local want `want' pre`h'
  }
  forval h = 0/`=`Q'-1' {
    local want `want' post`h'
  }
  global nchk = $nchk + 1
  if "`got'" != "`want'" | "`np'" != "`P'" global pf = $pf + 1
  dis "RES|`id'|S-lab|[`got'] numpre_att=`np'|" cond("`got'" == "`want'" & "`np'" == "`P'", "PASS", "FAIL")
  if "`got'" != "`want'" exit
  if !`hascs' dis "RES|`id'|C-align,C-b,C-att|csdid absent|SKIPPED (uncounted)"
  if !`hascs' exit
  use "${FX_`fx'}", clear
  csdid y, ivar(i) time(t) gvar(g) long2
  qui estimates store __CS
  estat event, post
  matrix CB0 = e(b)
  qui estimates restore __CS
  estat simple, post
  matrix CSB = e(b)
  qui estimates drop __CS
  local cn
  foreach nm of local got {
    if substr("`nm'", 1, 3) == "pre" local cn `cn' Tm`=substr("`nm'", 4, .)'
    else local cn `cn' Tp`=substr("`nm'", 5, .)'
  }
  local ncs = 0
  local cl : colnames CB0
  foreach nm of local cl {
    if inlist(substr("`nm'", 1, 2), "Tm", "Tp") local ncs = `ncs' + 1
  }
  pairup WB "`got'" CB0 "`cn'"
  local miss `r(miss)'
  local K : word count `got'
  * C-align: EXACT set equality.  Every wooldid horizon must exist in csdid's estat
  * event vector AND csdid must report no Tm/Tp horizon the window did not ask for, so
  * a silent offset or a fixture drift cannot hide inside the comparison below.
  global nchk = $nchk + 1
  if "`miss'" != "" | `ncs' != `K' global pf = $pf + 1
  dis "RES|`id'|C-align|`K' wooldid horizons, `ncs' csdid Tm/Tp, missing [`miss']|" cond("`miss'" == "" & `ncs' == `K', "PASS", "FAIL")
  if "`miss'" != "" exit
  * C-b: TOLERANCE, mreldif <= 1e-10 over every horizon.  The two estimators are
  * algebraically the same object here, so the only gap is floating-point ordering:
  * measured worst 1.02e-14 over 100 Monte Carlo draws x 2 DGPs x 14 horizons, and
  * 3.8e-15 on these two fits.  The gate keeps ~4 orders of headroom and still sits
  * ~4 orders below any gap that would mean the estimators had parted company.
  mrchk "`id'" C-b PA PB $TOL
  * C-att: TOLERANCE, the same gate on wooldid main_att vs csdid "estat simple"
  * (measured 2.9e-16 here, i.e. at the double-precision floor).
  local cw = colnumb(WEB, "main_att")
  if `cw' >= . bfail "`id'" C-att "main_att column absent from e(b)"
  else {
    matrix PA = WEB[1, `cw']
    matrix PB = CSB[1, 1]
    mrchk "`id'" C-att PA PB $TOL
  }
end

* Windows are the exact support each fixture can report: the common-onset panel runs
* rel -5..+6 (onset 2005, t 2000..2011), the staggered panel rel -8..+6 (the 2008
* cohort reaches back eight periods).
do_fixture C-common common 4 7 `HASCS'
do_fixture C-stag   stag   7 7 `HASCS'

timer off 9
qui timer list 9
dis "RES|ALL|counted|$nchk checks, $pf failing, " %6.1f r(t9) " s, csdid " cond(`HASCS', "present", "ABSENT") "|INFO"
if $pf == 0 dis "CSDID_IDENTITY_ALL_PASS"
else dis "CSDID_IDENTITY_FAILS=$pf"
exit, clear
