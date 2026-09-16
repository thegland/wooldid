* bootprefix_oracle.do — EXTERNAL ORACLE for wooldid's bootstrap(standard) engine,
*   checked against Stata's OWN `bootstrap` prefix driving wooldid as the command.
*
* WHY THIS IS AN EXACTNESS TEST AND NOT A MONTE-CARLO TEST.  Run at a common seed,
*   the prefix and wooldid's engine make the SAME `bsample` calls on the same RNG
*   stream, so they draw the IDENTICAL resamples — not two independent bootstrap
*   samples.  Every comparison below is therefore a floating-point identity, and
*   the pass rules are 1e-12 (or exact 0), not an MC band.  For reference, the MC
*   band the loose alternative would use is 1/sqrt(2R) = 5 PERCENT in relative terms
*   at R=199 (se/sqrt(2R) = 3.9e-3 on this fixture's scale) - thirteen orders of
*   magnitude looser than what these comparisons actually achieve, so it would pass a
*   thoroughly broken engine.  Asserting it would throw the point of the test away.
*   NO check in this file is MC-ERROR; each is labelled [EXACT], [EXACT-CONVENTION] (an identity
*   that holds only after the documented convention difference is encoded),
*   [SCOPE] (a boundary assertion, no oracle) or [RECORDED] (printed, uncounted).
*   Measured agreements from the 2026-09-13 feasibility probe, for calibration:
*   draws 8.1e-15 (dta_binary R=199), 4.3e-14 (admin R=99), 0 (boottstat leg);
*   percentile/BC endpoints 0 against Stata's own `bstat` on identical draws.
*
* THE TWO CONVENTION DIFFERENCES, encoded rather than compared across:
*   (1) REFERENCE DISTRIBUTION.  The prefix always reports z.  wooldid defaults to
*       bootdf(clusters) = t(G-1) when cluster() is given.  So the Wald/normal
*       interval is NOT compared head-on under the default: block A asserts the
*       exact t(G-1)/z critical-value RATIO instead, and block B re-runs wooldid
*       with bootdf(normal) and then asserts the endpoints head-on at 1e-12.
*   (2) PRINTED p.  The prefix prints a normal p on the bootstrap SE; wooldid's
*       default bootdisplay(percentile) prints a rank p with no add-one.  The
*       prefix has no percentile p at all, so NO p-value comparison is made.
*
* Blocks / fixtures (every leg is one same-seed pair):
*   A  12-cluster subset of dta_binary (G=12, N=144), R=120, seed 31415.
*      R divisible by 40 is deliberate: at level(95), 120*.025 = 3 is an INTEGER, the
*      branch of the percentile rule that averages v[3] and v[4].  R=99/199/499/999
*      never reach it.  The spec measured this pair at R=200; the agreement is an
*      identity, not an MC statement, so it does not depend on R, and the draws at a
*      common seed are the same stream (R=120 is the first 120 of them).  120 is the
*      cheapest multiple of 40 that still compares several hundred draw cells.
*      Covers C-BOOT-1 (averaging branch), C-BOOT-3, C-BOOT-4, C-BOOT-5, C-BOOT-7(a).
*   B  dta_binary entire (G=80, N=960), R=99, seed 20260912, bootdf(normal).
*      Covers C-BOOT-1 (ceil branch) and C-BOOT-2 (the Wald convention pin).
*   C  dta_s_admin (G=48, N=29,067 of 29,667; cluster(state) is NOT the i variable,
*      cluster sizes unequal, 600 rows leave the sample), R=99, seed 4242.
*      Covers C-BOOT-6, C-BOOT-7(b), and C-BOOT-8 on its state<=8 subset.
*      SKIPPED, uncounted, if the suite fixture is absent — the way reghdfe_oracle.do
*      skips its ppmlhdfe variant.  The boottest half of C-BOOT-8 is skipped the
*      same way if `boottest` is not installed.
*
* NOT IMPLEMENTED (stated, not silently dropped): C-BOOT-7(c), the divergence where
*   a replication that loses a horizon is DROPPED by the prefix (e(N_misreps)>0) but
*   KEPT with missing cells by wooldid (e(boot_N_V) < reps()).  No fixture in the
*   suite produces a shrunk replication — N_misreps was 0 on every run of the probe
*   and on every run here — and engineering one is a separate job.  The precondition
*   asserts below (A3, B7, C7) pin N_misreps == 0, so if a future fixture change
*   ever does produce one, the draw-by-draw comparisons fail loudly rather than
*   silently comparing misaligned rows.
*
* The prefix calls all carry `fastnose` when only _b is needed: the prefix throws the
*   per-replication analytic variance away, and computing it costs 1.6x (measured).
*   The one prefix leg that does NOT is the re-identification leg, which exists
*   precisely to harvest each replication's own standard error.
*
* Args: 1 = path to wooldid.ado (file or dir); 2 = path to <repo>/test/data.
*       The suite fixtures are derived as `2'/../suite/data.  No absolute paths.
* Run:  "<stata>" -e do bootprefix_oracle.do <adopath> <datadir>
*    -> expect BOOTPREFIX_ORACLE_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 20260913
set sortseed 20260913

local ADO `1'
if `"`ADO'"' == "" {
  dis as error "BOOTPREFIX_ORACLE_FAILS=1 (arg 1: path to wooldid.ado is required)"
  exit 198
}
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local DATA `2'
if `"`DATA'"' == "" {
  dis as error "BOOTPREFIX_ORACLE_FAILS=1 (arg 2: path to test/data is required)"
  exit 198
}
local SDATA `DATA'/../suite/data

global pf 0
global bn

***************************** helpers *****************************
* counted matrix comparison.  Prints the measured max reldif on every line so a
* later reader can see the headroom, the way reghdfe_oracle.do does.
cap program drop _mchk
program define _mchk
  args tag tol A B lab
  cap confirm matrix `A'
  local ra = _rc
  cap confirm matrix `B'
  local rb = _rc
  if `ra' | `rb' {
    dis as error "`tag' FAIL: matrix missing (`A' rc=`ra', `B' rc=`rb') — `lab'"
    global pf = $pf + 1
    exit
  }
  if rowsof(`A') != rowsof(`B') | colsof(`A') != colsof(`B') {
    dis as error "`tag' FAIL: nonconformable `A' " rowsof(`A') "x" colsof(`A') " vs `B' " rowsof(`B') "x" colsof(`B') " — `lab'"
    global pf = $pf + 1
    exit
  }
  local d = mreldif(`A', `B')
  local n = rowsof(`A')*colsof(`A')
  if `d' > `tol' | `d' >= . {
    dis as error "`tag' FAIL: max reldif = " %12.4e `d' " > `tol' over `n' cells — `lab'"
    global pf = $pf + 1
  }
  else dis "`tag' info: max reldif = " %12.4e `d' " (rule <= `tol', `n' cells) — `lab'"
end

* counted scalar comparison
cap program drop _schk
program define _schk
  args tag tol a b lab
  local d = reldif(`a', `b')
  if `d' > `tol' | `d' >= . {
    dis as error "`tag' FAIL: " %21.0g `a' " vs " %21.0g `b' " reldif = " %12.4e `d' " > `tol' — `lab'"
    global pf = $pf + 1
  }
  else dis "`tag' info: reldif = " %12.4e `d' " (rule <= `tol') — `lab'"
end

* counted equality of an integer-valued scalar
cap program drop _echk
program define _echk
  args tag a b lab
  if `a' != `b' {
    dis as error "`tag' FAIL: " `a' " != " `b' " — `lab'"
    global pf = $pf + 1
  }
  else dis "`tag' info: " `a' " == " `b' " — `lab'"
end

* build a draw matrix (reps x K) out of a `bootstrap` prefix saving() file, taking
* the columns BY NAME from $bn (the prefix names them _b_<e(b) name>) rather than by
* position, so a reordering cannot pass silently
cap program drop _pfmat
program define _pfmat
  args fn out
  preserve
  qui use `"`fn'"', clear
  local want
  foreach c of global bn {
    local want `want' _b_`c'
  }
  cap confirm variable `want'
  if _rc {
    dis as error "  _pfmat: saving file `fn' lacks one of [`want']"
    global pf = $pf + 1
    restore
    exit
  }
  qui mkmat `want', matrix(`out')
  restore
end

* build the same matrix out of a wooldid bootsave() file (columns are the plain e(b)
* names; row _rep == 0 is the OBSERVED estimate and is dropped here)
cap program drop _bsmat
program define _bsmat
  args fn out
  preserve
  qui use `"`fn'"', clear
  cap confirm variable $bn
  if _rc {
    dis as error "  _bsmat: bootsave file `fn' lacks one of [$bn]"
    global pf = $pf + 1
    restore
    exit
  }
  qui drop if _rep == 0
  sort _rep
  qui mkmat $bn, matrix(`out')
  restore
end

* posted bootstrap SE = sqrt of the diagonal of wooldid's reposted e(V)
cap program drop _sediag
program define _sediag
  args V out
  local k = colsof(`V')
  matrix `out' = J(1, `k', .)
  forval j = 1/`k' {
    matrix `out'[1, `j'] = sqrt(`V'[`j', `j'])
  }
end

timer clear

****************************************************************************
* BLOCK A — 12-cluster fixture, R = 120, seed 31415
*   the integer-n*alpha percentile branch, the draw-by-draw identity, the
*   studentized/re-identification oracle, and the count-vs-physical-resample
*   design claim.  Cheapest fixture per certified claim (measured 0.23-0.26 s/rep).
****************************************************************************
local S1 31415
local R1 120

cap confirm file `"`DATA'/dta_binary.dta"'
if _rc {
  dis as error "A0 FAIL: `DATA'/dta_binary.dta not found (arg 2 must be <repo>/test/data)"
  global pf = $pf + 1
  dis as error "BOOTPREFIX_ORACLE_FAILS=$pf"
  exit, clear
}
use `"`DATA'/dta_binary.dta"', clear
keep if inlist(i,1,2,3,26,27,28) | (i >= 51 & i <= 56)
tempfile f12
qui save `"`f12'"', replace

* A1 [EXACT] preconditions: the fixture is the one the pass rules were calibrated on
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1)
if _rc != 0 {
  dis as error "A1 FAIL: analytic run rc=" _rc
  global pf = $pf + 1
  dis as error "BOOTPREFIX_ORACLE_FAILS=$pf"
  exit, clear
}
_echk A1a e(N) 144 "12-cluster fixture N"
_echk A1b e(N_clust) 12 "12-cluster fixture G"
scalar a_att = _b[main_att]
scalar a_se  = _se[main_att]
matrix AB0 = e(b)
global bn : colnames AB0
* A1c [EXACT] the e(b) names the prefix has to address by hand (help defect D1: the
* help never spells them out; `bootstrap b=_b[att]` is r(111)).  Pinning them here
* means a rename cannot silently turn every oracle leg below into a skipped test.
if "$bn" != "main_att pretreatment_att ES_att_pre1 ES_att_post0" {
  dis as error "A1c FAIL: e(b) colnames [$bn] != [main_att pretreatment_att ES_att_pre1 ES_att_post0]"
  global pf = $pf + 1
}
else dis "A1c info: e(b) colnames [$bn]"
dis "A1 info: att = " %21.0g a_att " analytic se = " %21.0g a_se

* ---- wooldid's own engine, count leg, default bootdf(clusters) -> t(G-1)
use `"`f12'"', clear
cap erase wdbpA.dta
timer on 1
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) ///
    bootstrap(standard) reps(`R1') rseed(`S1') bootsave(wdbpA) replace
local rcA = _rc
timer off 1
if `rcA' != 0 {
  dis as error "A2 FAIL: wooldid bootstrap(standard) rc=`rcA'"
  global pf = $pf + 1
  dis as error "BOOTPREFIX_ORACLE_FAILS=$pf"
  exit, clear
}
matrix WB   = e(b)
matrix WD   = e(bootdist)
matrix WCP  = e(boot_ci_percentile)
matrix WCB  = e(boot_ci_bc)
matrix WZ0  = e(boot_z0)
matrix WCN  = e(boot_ci_normal)
matrix WV   = e(V)
scalar w_NV = e(boot_N_V)
scalar w_df = e(df_r)
scalar w_lv = cond(e(level) < ., e(level), c(level))
_sediag WV WSEP
matrix colnames WSEP = $bn

* A2 [EXACT] alignment preconditions.  If a replication ever shrank, the prefix
* would DROP it and wooldid would KEEP it with missing cells, and the row-by-row
* comparisons below would be comparing different draws.  See NOT IMPLEMENTED above.
_echk A2a w_NV `R1' "wooldid completed replications e(boot_N_V)"
_echk A2b rowsof(WD) `R1' "rows of e(bootdist)"
_echk A2c w_df 11 "e(df_r) under the default bootdf(clusters) = G-1"

* A3 [EXACT] bootsave() writes the same draws that e(bootdist) carries, in the same
* order, under the plain e(b) names.  Rule: max reldif exactly 0 (one accumulator,
* two writers).  This certifies the bootsave layout the C-BOOT-4 leg depends on.
_bsmat wdbpA BSA
_mchk A3 0 WD BSA "bootsave(wdbpA) draws vs e(bootdist)"

* ---- Stata's bootstrap prefix, same seed, count leg
use `"`f12'"', clear
tempfile pfA
timer on 2
cap noi bootstrap _b, reps(`R1') seed(`S1') cluster(i) idcluster(ni) nowarn nodots ///
    saving(`"`pfA'"', double every(1) replace): ///
    wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) fastnose
local rcP = _rc
timer off 2
if `rcP' != 0 {
  dis as error "A4 FAIL: bootstrap prefix rc=`rcP'"
  global pf = $pf + 1
  dis as error "BOOTPREFIX_ORACLE_FAILS=$pf"
  exit, clear
}
matrix PB  = e(b)
matrix PSE = e(se)
scalar p_mis = e(N_misreps)
scalar p_rep = e(N_reps)
cap noi estat bootstrap, all
matrix PCP = e(ci_percentile)
matrix PCB = e(ci_bc)
matrix PZ0 = e(z0)
matrix PCN = e(ci_normal)
_echk A4a p_rep `R1' "prefix e(N_reps)"
_echk A4b p_mis 0 "prefix e(N_misreps) — zero is the precondition for row alignment"
_pfmat `"`pfA'"' PFA

* ======================= C-BOOT-1 on the 12-cluster fixture =======================
* every rule 1e-12, i.e. EXACT.  Probe calibration on this shape: draws 2.1e-15.
_mchk A5 1e-12 WB PB   "[EXACT] point estimates e(b)"
_mchk A6 1e-12 WD PFA  "[EXACT] EVERY replication draw: e(bootdist) vs the prefix saving() file"
_mchk A7 1e-12 WSEP PSE "[EXACT] bootstrap SE: sqrt(diag(wooldid e(V))) vs prefix e(se)"
* A8: the same SE recomputed as the plain SD of the draws (n-1 divisor), which is
* what the prefix's e(se) is by construction
preserve
qui use wdbpA, clear
qui drop if _rep == 0
matrix WSD = J(1, colsof(WB), .)
local k = 0
foreach c of global bn {
  local ++k
  qui summarize `c'
  matrix WSD[1, `k'] = r(sd)
}
matrix colnames WSD = $bn
restore
_mchk A8 1e-12 WSD PSE "[EXACT] bootstrap SE recomputed as SD(draws) vs prefix e(se)"
* A9/A10/A11: the interval families.  R=120 => n*alpha = 3 is an INTEGER, so the
* percentile rule takes its AVERAGING branch (mean of v[3] and v[4]); R=99/199/499/999
* never exercise it.  Probe measurement on this branch (at R=200): 1.4e-15.
_mchk A9  1e-12 WCP PCP "[EXACT] percentile CI — INTEGER n*alpha averaging branch (R divisible by 40)"
_mchk A10 1e-12 WCB PCB "[EXACT] bias-corrected CI"
_mchk A11 1e-12 WZ0 PZ0 "[EXACT] bias-correction z0"

* ======================= C-BOOT-2, guard half =======================
* A12 [EXACT-CONVENTION].  The prefix's normal interval uses z; wooldid's default
* bootdf(clusters) uses t(G-1).  They are NOT compared head-on (measured gap 2.5e-3
* — a convention, not a mismatch).  What IS asserted is that the ONLY difference is
* the critical value: (ul_w - b)/(ul_p - b) == invttail(G-1, alpha/2)/invnormal(1-alpha/2)
* for every coefficient, at 1e-12.  Head-on endpoint equality is block B's job.
local alpha = (100 - w_lv)/100
local cvr = invttail(w_df, `alpha'/2) / invnormal(1 - `alpha'/2)
dis "A12 info: level = " w_lv ", t(" w_df ")/z critical-value ratio = " %18.12f `cvr'
local k = 0
local a12bad = 0
foreach c of global bn {
  local ++k
  local rw = (WCN[2, `k'] - WB[1, `k']) / (PCN[2, `k'] - PB[1, `k'])
  local rl = (WCN[1, `k'] - WB[1, `k']) / (PCN[1, `k'] - PB[1, `k'])
  if reldif(`rw', `cvr') > 1e-12 | reldif(`rl', `cvr') > 1e-12 {
    dis as error "A12 FAIL: `c' half-width ratios " %18.12f `rl' " / " %18.12f `rw' " vs t/z " %18.12f `cvr'
    local a12bad = 1
  }
}
if `a12bad' {
  global pf = $pf + 1
}
else dis "A12 info: all " `k' " coefficients' normal-interval endpoints differ from the prefix's by exactly the t(" w_df ")/z ratio (rule <= 1e-12)"

* ======================= C-BOOT-3, averaging branch =======================
* Feed wooldid's OWN draws to Stata's bstat and ask estat bootstrap for the same
* interval families.  No second estimator is involved: identical input, so the rule
* is max reldif EXACTLY 0 — the quantile rule, the BC z0 count and the order-statistic
* lookup are either the same function or they are not.  (A 1e-15 answer here would
* mean a different order statistic was picked, not rounding.)
preserve
qui use wdbpA, clear
qui drop if _rep == 0
cap noi bstat $bn, stat(WB) n(144)
local rcbs = _rc
if `rcbs' != 0 {
  dis as error "A13 FAIL: bstat on wooldid's own draws rc=`rcbs'"
  global pf = $pf + 1
}
else {
  cap noi estat bootstrap, all
  matrix DCP = e(ci_percentile)
  matrix DCB = e(ci_bc)
  matrix DZ0 = e(z0)
}
restore
if `rcbs' == 0 {
  _mchk A13 0 WCP DCP "[EXACT] percentile CI vs bstat on IDENTICAL draws (averaging branch)"
  _mchk A14 0 WCB DCB "[EXACT] bias-corrected CI vs bstat on IDENTICAL draws"
  _mchk A15 0 WZ0 DZ0 "[EXACT] z0 vs bstat on IDENTICAL draws"
}

* ======================= C-BOOT-4: the studentized / re-identification oracle ====
* wooldid's boottstat leg physically re-identifies each drawn cluster, so each
* replication has its own standard error — the studentization's denominator.  The
* prefix reproduces that leg exactly by being driven as `wooldid y ni t ttre,
* cluster(ni)`: idcluster() decides only which id the estimator is GIVEN.
use `"`f12'"', clear
cap erase wdbpB.dta
timer on 3
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) ///
    bootstrap(standard) reps(`R1') rseed(`S1') boottstat bootsave(wdbpB) replace
local rcT = _rc
timer off 3
if `rcT' != 0 {
  dis as error "A16 FAIL: wooldid boottstat rc=`rcT'"
  global pf = $pf + 1
}
cap matrix TB = e(b)
cap matrix WTS = e(boot_tstud)
cap scalar t_NV = e(boot_N_V)

use `"`f12'"', clear
tempfile pfB
timer on 4
cap noi bootstrap b=_b[main_att] se=_se[main_att], reps(`R1') seed(`S1') ///
    cluster(i) idcluster(ni) nowarn nodots saving(`"`pfB'"', double every(1) replace): ///
    wooldid y ni t ttre, att cluster(ni) esprelength(1) espostlength(1)
local rcR = _rc
timer off 4
if `rcR' != 0 {
  dis as error "A17 FAIL: prefix re-identification leg rc=`rcR'"
  global pf = $pf + 1
}
cap scalar r_mis = e(N_misreps)

if `rcT' == 0 & `rcR' == 0 {
  _echk A16a t_NV `R1' "boottstat completed replications"
  _echk A17a r_mis 0 "prefix re-id leg e(N_misreps)"
  preserve
  qui use `"`pfB'"', clear
  qui gen long _rep = _n
  keep _rep b se
  sort _rep
  tempfile reid
  qui save `"`reid'"', replace
  qui use wdbpB, clear
  qui keep if _rep > 0
  cap confirm variable se_main_att
  local rcse = _rc
  if `rcse' {
    dis as error "A18 FAIL: bootsave file under boottstat has no se_main_att column"
    global pf = $pf + 1
  }
  else {
    keep _rep main_att se_main_att
    sort _rep
    qui merge 1:1 _rep using `"`reid'"'
    qui count if _merge == 3
    local nm = r(N)
    qui gen double db = reldif(main_att, b)
    qui gen double ds = reldif(se_main_att, se)
    qui summarize db, meanonly
    local mxb = r(max)
    qui summarize ds, meanonly
    local mxs = r(max)
  }
  restore
  if !`rcse' {
    _echk A18a `nm' `R1' "replications matched 1:1 between bootsave() and the prefix saving() file"
    * A18/A19 [EXACT], rule 1e-12.  Probe measurement: exactly 0 on both, 200/200.
    if `mxb' > 1e-12 | `mxb' >= . {
      dis as error "A18 FAIL: per-replication main_att max reldif = " %12.4e `mxb' " > 1e-12"
      global pf = $pf + 1
    }
    else dis "A18 info: per-replication main_att max reldif = " %12.4e `mxb' " (rule <= 1e-12, `R1' reps) — [EXACT] boottstat draws vs prefix re-id leg"
    if `mxs' > 1e-12 | `mxs' >= . {
      dis as error "A19 FAIL: per-replication se_main_att max reldif = " %12.4e `mxs' " > 1e-12"
      global pf = $pf + 1
    }
    else dis "A19 info: per-replication se_main_att max reldif = " %12.4e `mxs' " (rule <= 1e-12, `R1' reps) — [EXACT] the studentization denominators"
  }
  * A20 [EXACT] bootsave()'s row _rep == 0 is the OBSERVED estimate, not a draw
  preserve
  qui use wdbpB, clear
  qui keep if _rep == 0
  local obs0 = main_att[1]
  local obss0 = se_main_att[1]
  restore
  _schk A20a 1e-12 `obs0' a_att "[EXACT] bootsave row _rep==0 carries the observed att"
  _schk A20b 1e-12 `obss0' a_se "[EXACT] bootsave row _rep==0 carries the observed analytic se"
}

* ======================= the studentized STATISTIC, not just its inputs ========
* A26 [EXACT].  A18/A19 certify the two ingredients of the percentile-t kernel (each
* replication's estimate and its own standard error).  A26 certifies the statistic
* wooldid actually reports from them: build the signed t-draws t_r = (b_r - b_obs)/se_r
* from the PREFIX's columns, hand them to Stata's bstat, and compare Stata's percentile
* endpoints against wooldid's e(boot_tstud) rows qlo_pctl / qhi_pctl.  Both sides use
* the same order-statistic rule (A13 proves wooldid's IS Stata's, bit for bit), so the
* only thing being compared is the studentized kernel itself.  Rule 1e-12.
* NOTE the p-value in that matrix is NOT compared: it is an ADD-ONE two-sided rank p
* (2*(1+min(m0,m1))/(n+1)), a wooldid convention with no counterpart in the prefix.
cap confirm matrix WTS
local hasts = _rc == 0
if `rcT' == 0 & !`hasts' {
  dis as error "A26 FAIL: the boottstat run posted no e(boot_tstud)"
  global pf = $pf + 1
}
if `rcT' == 0 & `rcR' == 0 & `hasts' {
  preserve
  qui use `"`pfB'"', clear
  qui gen double tdraw = (b - a_att)/se
  keep tdraw
  matrix TOBS = J(1, 1, a_att/a_se)
  matrix colnames TOBS = tdraw
  cap noi bstat tdraw, stat(TOBS) n(144)
  local rcbt = _rc
  if `rcbt' == 0 {
    cap noi estat bootstrap, percentile
    matrix TCI = e(ci_percentile)
  }
  restore
  if `rcbt' != 0 {
    dis as error "A26 FAIL: bstat on the prefix-built studentized draws rc=`rcbt'"
    global pf = $pf + 1
  }
  else {
    _schk A26a 1e-12 TCI[1,1] WTS[4,1] "[EXACT] studentized lower quantile qlo_pctl vs bstat on prefix-built t-draws"
    _schk A26b 1e-12 TCI[2,1] WTS[5,1] "[EXACT] studentized upper quantile qhi_pctl vs bstat on prefix-built t-draws"
    _schk A26c 1e-12 `=a_att/a_se' WTS[2,1] "[EXACT] t_obs is the ANALYTIC t (the studentization's reference statistic)"
    _echk A26d WTS[1,1] `R1' "studentized draws used (row n of e(boot_tstud))"
  }
}

* ======================= C-BOOT-5: count leg == physical resample ===============
* A21 [EXACT].  Drive the prefix twice at ONE seed — once with the copies pooled back
* into their original cohort (the count leg's statistic) and once with every copy its
* own cohort (physical re-identification) — and the att draws must agree.  This is
* the design claim that licenses wooldid running the cheap count leg for att/ame and
* the expensive re-id leg only for the cell-weighted estimands; no other test in the
* suite would notice if a change to the margins weights broke it.  Probe: 2.1e-15.
if `rcP' == 0 & `rcR' == 0 {
  preserve
  qui use `"`pfA'"', clear
  keep _b_main_att
  qui gen long _rep = _n
  sort _rep
  tempfile poolm
  qui save `"`poolm'"', replace
  qui use `"`pfB'"', clear
  keep b
  qui gen long _rep = _n
  sort _rep
  qui merge 1:1 _rep using `"`poolm'"', nogen
  qui gen double d5 = reldif(_b_main_att, b)
  qui summarize d5, meanonly
  local mx5 = r(max)
  restore
  if `mx5' > 1e-12 | `mx5' >= . {
    dis as error "A21 FAIL: pooled-id vs fresh-id prefix draws max reldif = " %12.4e `mx5' " > 1e-12"
    global pf = $pf + 1
  }
  else dis "A21 info: pooled-id vs fresh-id prefix draws max reldif = " %12.4e `mx5' " (rule <= 1e-12, `R1' reps) — [EXACT] count leg == physical resample for att"
}

* A22 [EXACT] wooldid's own two legs, same seed: the count leg's draws (block A's
* first run) against the boottstat run's re-identified draws.  Same claim as A21 from
* inside wooldid rather than from the prefix.
if `rcT' == 0 {
  preserve
  qui use wdbpB, clear
  qui keep if _rep > 0
  keep _rep main_att
  rename main_att ts_att
  sort _rep
  tempfile tsm
  qui save `"`tsm'"', replace
  qui use wdbpA, clear
  qui keep if _rep > 0
  keep _rep main_att
  sort _rep
  qui merge 1:1 _rep using `"`tsm'"', nogen
  qui gen double d22 = reldif(main_att, ts_att)
  qui summarize d22, meanonly
  local mx22 = r(max)
  restore
  if `mx22' > 1e-12 | `mx22' >= . {
    dis as error "A22 FAIL: wooldid count-leg vs re-id-leg draws max reldif = " %12.4e `mx22' " > 1e-12"
    global pf = $pf + 1
  }
  else dis "A22 info: wooldid count-leg vs re-id-leg draws max reldif = " %12.4e `mx22' " (rule <= 1e-12) — [EXACT] wooldid's own two legs at one seed"
}

* ======================= C-BOOT-7(a): the weights scope boundary ================
* [SCOPE], no oracle.  The prefix refuses a weighted command outright, so there is
* no external oracle on the weight axis at all; wooldid's own engine supports aw/pw.
* A test that states the boundary stops a future agent "fixing" the disagreement.
* wbp, not w: dta_binary already carries a variable named w, and `gen w' is r(110)
use `"`f12'"', clear
qui gen double wbp = 0.5 + mod(i, 4)/4
tempfile f12w
qui save `"`f12w'"', replace
cap noi bootstrap _b, reps(2) seed(1) cluster(i) idcluster(ni) nowarn nodots: ///
    wooldid y i t ttre [aw=wbp], att cluster(i) esprelength(1) espostlength(1) fastnose
local rcW = _rc
if `rcW' == 0 {
  dis as error "A23 FAIL: the bootstrap prefix ACCEPTED a weighted wooldid call (rc 0); the documented scope boundary (prefix refuses weights) no longer holds"
  global pf = $pf + 1
}
else dis "A23 info: prefix refuses the weighted command, rc=" `rcW' " — [SCOPE] no external oracle exists on the weight axis"
use `"`f12w'"', clear
cap noi wooldid y i t ttre [aw=wbp], att cluster(i) esprelength(1) espostlength(1) ///
    bootstrap(standard) reps(19) rseed(24)
local rcWW = _rc
if `rcWW' != 0 {
  dis as error "A24 FAIL: wooldid's own engine refused [aw=w] bootstrap, rc=`rcWW'"
  global pf = $pf + 1
}
else {
  _echk A24 e(boot_N_V) 19 "[SCOPE] wooldid's own engine completes a weighted bootstrap the prefix cannot drive"
}

* A25 [RECORDED, uncounted] probe defect: bootsave() takes a plain file name and its
* name check has no comma handling, so `bootsave(nm, replace)' is swallowed whole as
* the filename rather than parsed as an option.  Recorded, not asserted: it is a
* cosmetic parsing defect, and asserting a buggy behaviour would pin it in place.
use `"`f12'"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) ///
    bootstrap(standard) reps(5) rseed(1) bootsave(wdbpZ, replace)
local rcZ = _rc
cap confirm file "wdbpZ, replace.dta"
local rcZf = _rc
dis "A25 info (uncounted): bootsave(wdbpZ, replace) rc=" `rcZ' "; file {wdbpZ, replace}.dta exists rc=" `rcZf' " (0 = the comma form was swallowed as a filename)"
cap erase "wdbpZ, replace.dta"
cap erase wdbpZ.dta

****************************************************************************
* BLOCK B — dta_binary entire (G=80, N=960), R = 99, seed 20260912.
*   A second, larger fixture for the C-BOOT-1 asserts (99*.025 = 2.475, so the
*   percentile rule takes its CEIL branch here, the complement of block A), and the
*   asserted half of C-BOOT-2: with bootdf(normal) wooldid must reproduce the
*   prefix's normal interval endpoint for endpoint.
****************************************************************************
local S2 20260912
local R2 99

use `"`DATA'/dta_binary.dta"', clear
timer on 5
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) ///
    bootstrap(standard) reps(`R2') rseed(`S2') bootdf(normal)
local rcB = _rc
timer off 5
if `rcB' != 0 {
  dis as error "B1 FAIL: wooldid bootstrap(standard) bootdf(normal) rc=`rcB'"
  global pf = $pf + 1
}
cap matrix XB  = e(b)
cap matrix XD  = e(bootdist)
cap matrix XCP = e(boot_ci_percentile)
cap matrix XCB = e(boot_ci_bc)
cap matrix XZ0 = e(boot_z0)
cap matrix XCN = e(boot_ci_normal)
cap matrix XV  = e(V)
cap scalar x_NV = e(boot_N_V)
cap scalar x_df = e(df_r)
cap _sediag XV XSEP
cap matrix colnames XSEP = $bn

use `"`DATA'/dta_binary.dta"', clear
tempfile pfX
timer on 6
cap noi bootstrap _b, reps(`R2') seed(`S2') cluster(i) idcluster(ni) nowarn nodots ///
    saving(`"`pfX'"', double every(1) replace): ///
    wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) fastnose
local rcBP = _rc
timer off 6
if `rcBP' != 0 {
  dis as error "B2 FAIL: bootstrap prefix rc=`rcBP'"
  global pf = $pf + 1
}
cap matrix YB  = e(b)
cap matrix YSE = e(se)
cap scalar y_mis = e(N_misreps)
cap noi estat bootstrap, all
cap matrix YCP = e(ci_percentile)
cap matrix YCB = e(ci_bc)
cap matrix YZ0 = e(z0)
cap matrix YCN = e(ci_normal)

if `rcB' == 0 & `rcBP' == 0 {
  _pfmat `"`pfX'"' PFX
  _echk B3 x_NV `R2' "wooldid completed replications"
  _echk B4 y_mis 0 "prefix e(N_misreps) — row-alignment precondition"
  _mchk B5 1e-12 XB YB   "[EXACT] point estimates e(b), G=80"
  _mchk B6 1e-12 XD PFX  "[EXACT] EVERY replication draw, G=80 (`R2' x 4 cells)"
  _mchk B7 1e-12 XSEP YSE "[EXACT] bootstrap SE: sqrt(diag(e(V))) vs prefix e(se)"
  _mchk B8 1e-12 XCP YCP "[EXACT] percentile CI — CEIL branch (n*alpha = 2.475)"
  _mchk B9 1e-12 XCB YCB "[EXACT] bias-corrected CI"
  _mchk B10 1e-12 XZ0 YZ0 "[EXACT] bias-correction z0"
  * B11 [EXACT-CONVENTION].  THE convention pin: under bootdf(normal) wooldid adopts
  * the prefix's reference distribution and the Wald family must then agree endpoint
  * for endpoint.  Under the default it must NOT (block A's A12 asserts the ratio).
  * A failure here means the SD definition or the interval construction moved, not
  * that the two packages disagree about degrees of freedom.
  _mchk B11 1e-12 XCN YCN "[EXACT-CONVENTION] normal/Wald CI under bootdf(normal) vs prefix e(ci_normal)"
  * B12 [EXACT-CONVENTION] bootdf(normal) blanks e(df_r) so a hand-run test/lincom on
  * the bootstrap covariance reproduces the printed z statistics rather than a t
  if x_df < . {
    dis as error "B12 FAIL: e(df_r) = " x_df " under bootdf(normal); expected missing (z reference)"
    global pf = $pf + 1
  }
  else dis "B12 info: e(df_r) is missing under bootdf(normal) — [EXACT-CONVENTION] z reference, matching the prefix"
}

****************************************************************************
* BLOCK C — dta_s_admin (SKIPPED, uncounted, if the suite fixture is absent).
*   The shape blocks A and B do not cover: the cluster variable is NOT the i
*   variable, cluster sizes are unequal, and 600 rows leave the estimation sample.
*   Probe measurement here is one decade looser (4.3e-14) than on dta_binary purely
*   because the panel is 30x larger and the aggregation sums more terms; 1e-12 still
*   leaves two decades of headroom, and it is still exactness, not MC agreement.
****************************************************************************
local S3 4242
local R3 99

cap confirm file `"`SDATA'/dta_s_admin.dta"'
if _rc {
  dis "BOOTPREFIX_ORACLE_SKIPPED_ADMIN: dta_s_admin.dta not found under `SDATA'"
  dis "  (block C — C-BOOT-6, C-BOOT-7b, C-BOOT-8 — skipped, uncounted)"
}
else {
  use `"`SDATA'/dta_s_admin.dta"', clear
  cap erase wdbpC.dta
  timer on 7
  cap noi wooldid y cohid year ttre, att cluster(state) esprelength(1) espostlength(1) ///
      bootstrap(standard) reps(`R3') rseed(`S3') bootsave(wdbpC) replace
  local rcC = _rc
  timer off 7
  if `rcC' != 0 {
    dis as error "C1 FAIL: wooldid bootstrap(standard) on the admin panel rc=`rcC'"
    global pf = $pf + 1
  }
  cap matrix AB  = e(b)
  cap matrix AD  = e(bootdist)
  cap matrix ACP = e(boot_ci_percentile)
  cap matrix ACB = e(boot_ci_bc)
  cap matrix AZ0 = e(boot_z0)
  cap matrix ACN = e(boot_ci_normal)
  cap matrix AV  = e(V)
  cap scalar c_NV = e(boot_N_V)
  cap scalar c_df = e(df_r)
  cap scalar c_lv = cond(e(level) < ., e(level), c(level))
  cap scalar c_G  = e(N_clust)
  cap _sediag AV ASEP
  cap matrix colnames ASEP = $bn

  use `"`SDATA'/dta_s_admin.dta"', clear
  tempfile pfC
  timer on 8
  cap noi bootstrap _b, reps(`R3') seed(`S3') cluster(state) idcluster(ns) nowarn nodots ///
      saving(`"`pfC'"', double every(1) replace): ///
      wooldid y cohid year ttre, att cluster(state) esprelength(1) espostlength(1) fastnose
  local rcCP = _rc
  timer off 8
  if `rcCP' != 0 {
    dis as error "C2 FAIL: bootstrap prefix on the admin panel rc=`rcCP'"
    global pf = $pf + 1
  }
  cap matrix CB  = e(b)
  cap matrix CSE = e(se)
  cap scalar c_mis = e(N_misreps)
  cap noi estat bootstrap, all
  cap matrix CCP = e(ci_percentile)
  cap matrix CCB = e(ci_bc)
  cap matrix CZ0 = e(z0)
  cap matrix CCN = e(ci_normal)

  if `rcC' == 0 & `rcCP' == 0 {
    _pfmat `"`pfC'"' PFC
    _echk C3 c_NV `R3' "wooldid completed replications, admin panel"
    _echk C4 c_mis 0 "prefix e(N_misreps), admin panel"
    _echk C5 c_G 48 "admin panel G (cluster(state), not the i variable)"
    _bsmat wdbpC BSC
    _mchk C5b 0 AD BSC "[EXACT] bootsave(wdbpC) draws vs e(bootdist), admin panel"
    _mchk C6  1e-12 AB CB   "[EXACT] point estimates e(b), admin panel"
    _mchk C7  1e-12 AD PFC  "[EXACT] EVERY replication draw, admin panel (`R3' x 4 cells)"
    _mchk C8  1e-12 ASEP CSE "[EXACT] bootstrap SE: sqrt(diag(e(V))) vs prefix e(se)"
    _mchk C9  1e-12 ACP CCP "[EXACT] percentile CI, admin panel (ceil branch)"
    _mchk C10 1e-12 ACB CCB "[EXACT] bias-corrected CI, admin panel"
    _mchk C11 1e-12 AZ0 CZ0 "[EXACT] bias-correction z0, admin panel"
    * C12 [EXACT-CONVENTION]: the t(G-1)/z ratio again, on a panel where G-1 = 47 and
    * the cluster variable is not i.  Same encoding as A12 — never compared head-on.
    local alpha3 = (100 - c_lv)/100
    local cvr3 = invttail(c_df, `alpha3'/2) / invnormal(1 - `alpha3'/2)
    local k = 0
    local c12bad = 0
    foreach c of global bn {
      local ++k
      local rw = (ACN[2, `k'] - AB[1, `k']) / (CCN[2, `k'] - CB[1, `k'])
      local rl = (ACN[1, `k'] - AB[1, `k']) / (CCN[1, `k'] - CB[1, `k'])
      if reldif(`rw', `cvr3') > 1e-12 | reldif(`rl', `cvr3') > 1e-12 {
        dis as error "C12 FAIL: `c' half-width ratios " %18.12f `rl' " / " %18.12f `rw' " vs t/z " %18.12f `cvr3'
        local c12bad = 1
      }
    }
    if `c12bad' {
      global pf = $pf + 1
    }
    else dis "C12 info: all " `k' " normal-interval endpoints differ from the prefix's by exactly the t(" c_df ")/z ratio " %18.12f `cvr3' " (rule <= 1e-12) — [EXACT-CONVENTION]"
    * C13 [EXACT] C-BOOT-3 on the ceil branch: wooldid's own admin draws through bstat
    preserve
    qui use wdbpC, clear
    qui drop if _rep == 0
    cap noi bstat $bn, stat(AB) n(29067)
    local rcbs3 = _rc
    if `rcbs3' == 0 {
      cap noi estat bootstrap, all
      matrix ECP = e(ci_percentile)
      matrix ECB = e(ci_bc)
      matrix EZ0 = e(z0)
    }
    restore
    if `rcbs3' != 0 {
      dis as error "C13 FAIL: bstat on wooldid's own admin draws rc=`rcbs3'"
      global pf = $pf + 1
    }
    else {
      _mchk C13 0 ACP ECP "[EXACT] percentile CI vs bstat on IDENTICAL draws (ceil branch)"
      _mchk C14 0 ACB ECB "[EXACT] bias-corrected CI vs bstat on IDENTICAL draws"
      _mchk C15 0 AZ0 EZ0 "[EXACT] z0 vs bstat on IDENTICAL draws"
    }
  }

  * ============ C-BOOT-7(b): the re-identification leg under a cross-cutting cluster ==
  * [SCOPE].  On this panel cluster(state) does not nest i = cohid (a cohort spans
  * states).  Until T3 (2026-09-16) the physical re-identification leg refused that
  * with r(459); it now has a rule for it - keep i, and let a twice-drawn cluster
  * arrive as duplicate rows inside its unit's cell - so what this leg pins is that
  * the request RUNS and completes its replications.  Re-identifying here would split
  * a cohort that straddles two states into two cohorts in EVERY replication, drawn or
  * not, which is the failure the keep-i rule exists to avoid.
  use `"`SDATA'/dta_s_admin.dta"', clear
  cap noi wooldid y cohid year ttre, att cluster(state) esprelength(1) espostlength(1) ///
      bootstrap(standard) reps(5) rseed(1) boottstat
  local rcN = _rc
  local ndN = e(boot_completed)
  _echk C16 `rcN' 0 "[SCOPE] boottstat RUNS where cluster() does not nest i (was r(459) before T3)"
  _echk C16b `ndN' 5 "[SCOPE] every replication of that run completed"

  * ============ C-BOOT-8: hand-rebuilt stage 1 (and boottest, recorded only) ======
  * [EXACT re-derivation, rule 1e-10].  An independent reconstruction of wooldid's
  * att and its clustered analytic standard error: reghdfe with ONE absorbed set
  * (i.cell + i.year, absorb(cohid)) plus a count-weighted lincom over the post
  * cells.  This is the only leg in this file that touches the ANALYTIC standard
  * error; every other leg is about the bootstrap engine.  Probe: 1.4e-13 / 5.9e-13
  * at G=8 (and 6.3e-14 / 8.4e-12 on the full 48-state panel).
  use `"`SDATA'/dta_s_admin.dta"', clear
  keep if state <= 8
  cap noi wooldid y cohid year ttre, att cluster(state) esprelength(1) espostlength(1)
  local rcH = _rc
  if `rcH' != 0 {
    dis as error "C17 FAIL: analytic wooldid on the 8-state subset rc=`rcH'"
    global pf = $pf + 1
  }
  else {
    scalar h_att = _b[main_att]
    scalar h_se  = _se[main_att]
    qui gen relidx = year - ttre + 2 if ttre < . & year >= ttre - 1
    qui gen long cell = 0
    qui replace cell = cohid*100 + relidx if relidx < .
    qui levelsof cell if cell > 0, local(hc)
    cap noi qui reghdfe y i.cell i.year, absorb(cohid) cluster(state) tol(1e-9)
    local rcRH = _rc
    if `rcRH' != 0 {
      dis as error "C17 FAIL: hand-rebuilt reghdfe rc=`rcRH'"
      global pf = $pf + 1
    }
    else {
      local hexpr 0
      local hwtot 0
      foreach c of local hc {
        qui count if cell == `c' & e(sample)
        local wc = r(N)
        if mod(`c', 100) >= 2 {
          local hexpr `hexpr' + `wc'*_b[`c'.cell]
          local hwtot = `hwtot' + `wc'
        }
      }
      cap noi qui lincom (`hexpr')/`hwtot'
      local rcL = _rc
      if `rcL' != 0 {
        dis as error "C17 FAIL: count-weighted lincom rc=`rcL'"
        global pf = $pf + 1
      }
      else {
        local hb = r(estimate)
        local hs = r(se)
        _schk C17 1e-10 `hb' h_att "[EXACT] hand-rebuilt single-absorb att vs wooldid's"
        _schk C18 1e-10 `hs' h_se "[EXACT] hand-rebuilt clustered analytic se vs wooldid's"
        * the boottest half is QUALITATIVE and mostly RECORDED: boottest is a wild
        * cluster bootstrap-t with the null imposed, wooldid's engine a pairs bootstrap
        * without it, so their p-values are not expected to agree to any tolerance —
        * and on every fixture in the suite both sit at their resolution floor, so no
        * numeric rule for a p comparison can be calibrated at all.  What IS asserted
        * is that boottest runs (rc 0) and that the t-statistic it reports is wooldid's
        * analytic t, which pins the rebuild onto wooldid's estimand.
        cap which boottest
        if _rc {
          dis "BOOTPREFIX_ORACLE_SKIPPED_BOOTTEST: boottest not installed"
          dis "  (C-BOOT-8's wild-cluster half skipped, uncounted; C17/C18 above still ran)"
        }
        else {
          cap noi boottest ((`hexpr')/`hwtot' = 0), cluster(state) reps(499) seed(11) nograph ptype(equaltail)
          local rcBT = _rc
          if `rcBT' != 0 {
            dis as error "C19 FAIL: boottest on the hand-rebuilt model rc=`rcBT'"
            global pf = $pf + 1
          }
          else {
            local bt_t = r(t)
            local bt_p = r(p)
            _schk C19 1e-6 `bt_t' `=h_att/h_se' "[EXACT] boottest's t vs wooldid's analytic t (rule 1e-6: boottest prints 8 significant digits)"
            dis "C19 info (uncounted, RECORDED not asserted): boottest equal-tail p = " %10.6f `bt_p' " at 499 wild replications; wooldid's own percentile p on this fixture is equally saturated, so no p comparison is asserted anywhere in this file"
            cap noi matrix list r(CI), format(%12.0g)
          }
        }
      }
    }
  }
}

cap erase wdbpA.dta
cap erase wdbpB.dta
cap erase wdbpC.dta
timer list
if $pf == 0 dis "BOOTPREFIX_ORACLE_ALL_PASS"
else dis as error "BOOTPREFIX_ORACLE_FAILS=$pf"
exit, clear
