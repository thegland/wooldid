* ritest_crosscheck.do - EXTERNAL cross-check of wooldid's BUILT-IN randomization
*   inference against two outside engines that randomize the DATA and then run an
*   ordinary wooldid call on each randomization:
*     - the user-written ritest package (ritest 1.1.7; an OPTIONAL dependency), and
*     - Stata's own built-in permute prefix (always available).
*   The question answered is the author's: "compare randomization inference built in
*   to what happens if we externally randomize the data and execute regular wooldid
*   on those randomizations."  Every block below does exactly that - the external
*   engine owns the assignment draw, an ordinary "wooldid ... fastnose" call (fastnose
*   is wooldid's own fast-mode per-draw engine under permute()) is the per-draw
*   command, and the resulting reference distribution and p-value are compared with
*   what wooldid's internal permute() produced.
*
* WHAT THE EXTERNAL ENGINES ARE.  ritest's cluster() kernel and wooldid's
*   permute(shuffle) both draw a uniform permutation of the SAME cluster-level onset
*   multiset within each stratum (wooldid.ado Mata wooldid_permx vs ritest.ado
*   permute_simple), so they induce the same uniform law over the same assignment set.
*   permute has no cluster() option, so it is run on a CLUSTER-LEVEL FRAME through an
*   eclass sortpreserve wrapper (see wd_cl below) - there it permutes the identical
*   one-row-per-cluster onset vector, within strata(), and its enumerate mode walks
*   the same distinct-assignment set that wooldid's permexact enumerates.
*
* CONVENTIONS (source-verified; the numeric mapping the checks rely on).
*   ritest : p = c/n with c = #{|T*| >= |T_obs| - eps}, eps 1e-7, NO add-one.
*   permute: c_lower = #{T <= T_obs + eps}, c_upper = #{T >= T_obs - eps},
*            p_twosided = min(n, 2*min(c_lo,c_hi))/n, NO add-one.
*   wooldid: exact comparison (no eps); sampled p's carry the add-one (1+c)/(n+1);
*            permexact drops it.  So the ONLY wooldid convention that can equal an
*            external p is permrawp (= c_abs/n, which wooldid.sthlp:590 says in as
*            many words), or wooldid's plain tails recomputed offline.  Under
*            permexact wooldid already uses plain shares, so permute's p_twosided
*            equals e(perm_p_att) there exactly.  Block X6 pins all of this.
*   Two numerical traps, both measured and both handled below:
*     (a) ritest's and permute's saving() files are FLOAT unless double is given -
*         without it a bitwise 199/199 becomes 0/199 at a 5.6e-08 spread.  Every
*         saving() here passes double.
*     (b) BOTH prefixes report the observed statistic r(b) through a ~15-significant-
*         digit decimal round trip (ritest.ado:607-608, permute.ado:802/:942), so r(b)
*         matches e(perm_comparator_b) to ~1e-15 RELATIVE, not bitwise.  Observed
*         statistics are therefore compared at reldif <= 1e-13, never bitwise.
*
* BLOCKS (each pass rule is numeric and stated beside the check):
*   X1  EXACT: permute with enumerate on the cluster frame vs wooldid permexact,
*       tiny 6-cluster panel (C(6,3) = 20 assignments).  Deterministic on both
*       sides: draw sets at reldif <= 1e-10 (bitwise count reported), counts, both
*       one-sided p's and the two-sided p at 1e-12, observed statistic at 1e-13.
*   X2  permute SAMPLED on the same panel: every draw a member of the 20-value exact
*       set at reldif <= 1e-10; exactly 20 distinct cluster patterns; 0 multiset
*       violations; chi2(19) of pattern frequencies < 43.82 (p .001).  Tails compared
*       against the exact-set counts RECOMPUTED UNDER permute's own eps rule, with the
*       z-test skipped (equality asserted instead) for any tail whose exact share is
*       0 or 1 - a blind 4*sqrt(p(1-p)/n) rule is a zero-tolerance test there.
*   X3  EXACT, and the highest-value check: ritest REPLAYS wooldid's own saved draws
*       through samplingsourcefile(), so there is no Monte-Carlo tolerance at all.
*       Rep-by-rep equality at reldif <= 1e-10 (bitwise count reported), SD at 1e-10,
*       observed statistic at 1e-13, and ritest's p identical to wooldid's own
*       permrawp engine number at 1e-12.  On dta_binary with a NULL outcome.
*   X4  MONTE-CARLO: independent draws, ritest cluster(i) and permute (wrapper), on
*       the same null outcome.  ONE-SIDED tails only (never the two-sided p: the
*       min-and-double rule amplifies a tail difference near 1/2) and a NULL outcome
*       only (on the effect-carrying outcome every engine sits at its floor and the
*       test is vacuous - measured: wooldid .04 = its 2/(n+1) floor, ritest 0).
*   X5  permstrata() vs ritest strata(), and permute(probabilities) vs a ritest
*       samplingprogram() that walks the same per-cluster probability table.
*   X6  Convention identities from wooldid's own draws (exact and sampled), which is
*       what licenses X3/X4 to compare like with like.
*
* OPTIONAL DEPENDENCY.  ritest is user-written.  If it is absent, the ritest legs
*   (X3's replay, X4's ritest leg, and all of X5) are skipped UNCOUNTED and a SKIPPED
*   line is printed; X1/X2/X6 and X4's permute leg still run (permute is built in) and
*   the file still ends in its ALL_PASS sentinel, the way reghdfe_oracle.do treats a
*   missing ppmlhdfe.  The suite grader requires an ALL_PASS line and no FAILS= line,
*   so a missing optional dependency must not exit without the sentinel.
*
* Args: 1 = path to wooldid.ado (file or its directory)
*       2 = path to <repo>/test/data  (dta_binary.dta lives there; a suite fixture
*           would be `2'/../suite/data - none is needed here)
*       3 = optional, "fast" for a reduced-replication smoke pass (development only;
*           the suite harness never passes a third argument)
* Run:  "<stata>" -e do ritest_crosscheck.do <adopath> <datadir>
*    -> expect RITEST_CROSSCHECK_ALL_PASS
* Cost: ~14 min at the shipped counts (850 s measured on 2026-09-13 with 4-5 other
*   Stata jobs running).  The newfeat harness timeout is 1800 s, which is why the
*   dta_binary blocks run at reps 99 rather than the 199 the feasibility study
*   measured; every tolerance below is computed from the realised n, so the counts can
*   be raised or lowered without touching a pass rule.

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 20260913
set sortseed 20260913

global pf 0

local ADO `1'
if `"`ADO'"' == "" {
  dis as error "X0 FAIL: arg 1 (path to wooldid.ado or its directory) is required"
  dis as error "RITEST_CROSSCHECK_FAILS=1"
  exit, clear
}
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'

local DATA `2'
if `"`DATA'"' == "" {
  dis as error "X0 FAIL: arg 2 (path to <repo>/test/data) is required"
  dis as error "RITEST_CROSSCHECK_FAILS=1"
  exit, clear
}

* replication counts.  Shipped values reproduce the feasibility study's designs at
* counts the 1800 s newfeat timeout can carry; "fast" is for development only.
* RB   = dta_binary reps (X3, X4); the study measured at 199
* RMC  = permute sampled reps on the 20-assignment panel (X2)
* RSMP = sampled reps for the X6 convention identities (permreps() minimum is 20)
* RSTR = ritest reps on the 24-assignment stratified design (X5)
* RPRB = reps on the probabilities design (X5)
local RB   99
local RMC  500
local RSMP 60
local RSTR 720
local RPRB 400
if "`3'" == "fast" {
  local RB   20
  local RMC  60
  local RSMP 20
  local RSTR 48
  local RPRB 60
  dis "X0 info: FAST mode - development counts, NOT the shipped ones"
}
dis "X0 info: reps  dta_binary=`RB'  permute-mc=`RMC'  sampled=`RSMP'  strata=`RSTR'  probs=`RPRB'"

* ---- optional dependency ----
cap which ritest
local HASRI = (_rc == 0)
if !`HASRI' {
  dis "RITEST_CROSSCHECK_SKIPPED: ritest is not installed; the ritest legs of X3/X4 and all of X5 are skipped and UNCOUNTED."
  dis "X0 info: the built-in permute blocks X1/X2/X4 and the convention block X6 still run."
}

***************************** shared helpers *****************************
* rxc_rec records the CLUSTER-LEVEL onset vector the external engine handed to
* wooldid on this call, one column per call, so the ASSIGNMENT SPACE itself (not just
* the coefficient it produced) can be checked.  zeromiss=1 for the permute cluster
* frame, where never-treated is carried as the sentinel level 0.
mata:
void rxc_reset()
{
    external real matrix RXCREC
    RXCREC = J(0, 0, .)
}
void rxc_rec(string scalar iv, string scalar ov, real scalar zeromiss)
{
    external real matrix RXCREC
    real matrix D
    real colvector v
    real scalar r
    D = sort(st_data(., (iv, ov)), 1)
    if (zeromiss == 1) {
        D[., 2] = editvalue(D[., 2], 0, .)
    }
    v = D[1, 2]
    for (r = 2; r <= rows(D); r++) {
        if (D[r, 1] != D[r - 1, 1]) {
            v = v \ D[r, 2]
        }
    }
    if (rows(RXCREC) == 0) {
        RXCREC = v
    }
    else {
        RXCREC = RXCREC, v
    }
}
// rxc_pat takes the LAST nrep recorded columns - the lead-in columns are the prefix's
// own observed-data calls (1 for permute, 2 for ritest) - and reports the lead-in
// count, so a change in either prefix surfaces as a failure rather than as a silent
// misalignment of the pattern analysis.  M > 0 asks for a uniformity chi2 over M cells.
void rxc_pat(string scalar tag, real scalar nrep, real colvector obs, real colvector str, real scalar M)
{
    external real matrix RXCREC
    real matrix P, P2
    real scalar R, c, nd, si, ok, e
    real colvector o2, idx, freq
    string colvector keys, u
    st_numscalar("__" + tag + "_lead", cols(RXCREC) - nrep)
    if (cols(RXCREC) < nrep | nrep < 1) {
        st_numscalar("__" + tag + "_R", cols(RXCREC))
        st_numscalar("__" + tag + "_nd", -1)
        st_numscalar("__" + tag + "_bad", -1)
        st_numscalar("__" + tag + "_minf", -1)
        st_numscalar("__" + tag + "_maxf", -1)
        st_numscalar("__" + tag + "_chi2", .)
        return
    }
    P = RXCREC[., (cols(RXCREC) - nrep + 1)::cols(RXCREC)]
    R = cols(P)
    P2 = editmissing(P, -1)
    o2 = editmissing(obs, -1)
    ok = 0
    for (c = 1; c <= R; c++) {
        for (si = 1; si <= rows(uniqrows(str)); si++) {
            idx = selectindex(str :== uniqrows(str)[si])
            if (sort(P2[idx, c], 1) != sort(o2[idx], 1)) {
                ok++
            }
        }
    }
    keys = J(R, 1, "")
    for (c = 1; c <= R; c++) {
        keys[c] = invtokens(strofreal(P2[., c]'))
    }
    u = uniqrows(keys)
    nd = rows(u)
    freq = J(nd, 1, 0)
    for (c = 1; c <= nd; c++) {
        freq[c] = sum(keys :== u[c])
    }
    st_numscalar("__" + tag + "_R", R)
    st_numscalar("__" + tag + "_nd", nd)
    st_numscalar("__" + tag + "_bad", ok)
    st_numscalar("__" + tag + "_minf", min(freq))
    st_numscalar("__" + tag + "_maxf", max(freq))
    if (M > 0) {
        e = R / M
        st_numscalar("__" + tag + "_chi2", sum((freq :- e):^2 :/ e) + (M - nd) * e)
    }
    else {
        st_numscalar("__" + tag + "_chi2", .)
    }
}
// membership of every draw in an exact reference set, at a RELATIVE tolerance.  The
// bitwise count is reported alongside but never asserted: float saving() files and
// decimal round trips make a bitwise assertion a flake (header traps a and b).
void rxc_member(string scalar dname, string scalar rname, real scalar tol)
{
    real colvector d, r, m
    real matrix A
    d = st_matrix(dname)
    r = st_matrix(rname)
    if (cols(d) > 1) {
        d = d[., 1]
    }
    if (cols(r) > 1) {
        r = r[., 1]
    }
    d = select(d, d :< .)
    A = abs(d * J(1, rows(r), 1) - J(rows(d), 1, 1) * r')
    m = rowmin(A :/ (J(rows(d), 1, 1) * abs(r)' :+ 1))
    st_numscalar("__mem_n", rows(d))
    st_numscalar("__mem_bad", sum(m :> tol))
    st_numscalar("__mem_maxrel", max(m))
    st_numscalar("__mem_bit", sum(rowmin(A) :== 0))
}
// rep-by-rep comparison of the external engine's saved draws (variable _pm_1 in the
// dataset in memory) with wooldid's own e(permdist_b) column, in replication order
void rxc_repcmp(string scalar mname)
{
    real colvector a, b, rd
    a = st_data(., "_pm_1")
    b = st_matrix(mname)
    if (cols(b) > 1) {
        b = b[., 1]
    }
    st_numscalar("__rc_na", rows(a))
    st_numscalar("__rc_nb", rows(b))
    if (rows(a) != rows(b)) {
        st_numscalar("__rc_n", -1)
        return
    }
    rd = abs(a - b) :/ (abs(b) :+ 1)
    st_numscalar("__rc_n", rows(a))
    st_numscalar("__rc_bit", sum(a :== b))
    st_numscalar("__rc_maxrel", max(rd))
    st_numscalar("__rc_bad", sum(rd :> 1e-10))
    st_numscalar("__rc_sd", sqrt(variance(a)))
}
// the tail counts and the ddof=1 SD of a wooldid draw set, recomputed OFFLINE from
// e(permdist_b) and the exact-bit comparator, under BOTH wooldid's exact rule and the
// external engines' eps rule.  __t_mintie is the closest any draw comes to |s_obs|,
// which is the precondition that makes the two tie conventions coincide.
void rxc_tails(string scalar mname, real scalar obs, real scalar eps)
{
    real colvector b
    b = st_matrix(mname)
    if (cols(b) > 1) {
        b = b[., 1]
    }
    b = select(b, b :< .)
    st_numscalar("__t_n", rows(b))
    st_numscalar("__t_clo", sum(b :<= obs))
    st_numscalar("__t_chi", sum(b :>= obs))
    st_numscalar("__t_cabs", sum(abs(b) :>= abs(obs)))
    st_numscalar("__t_cloe", sum(b :<= obs + eps))
    st_numscalar("__t_chie", sum(b :>= obs - eps))
    st_numscalar("__t_cabse", sum(abs(b) :>= abs(obs) - eps))
    st_numscalar("__t_mintie", min(abs(abs(b) :- abs(obs))))
    st_numscalar("__t_sd", sqrt(variance(b)))
}
// bitwise identity of two draw matrices (X6 and X3: same rseed, different convention)
void rxc_bitcmp(string scalar m1, string scalar m2)
{
    real colvector a, b
    a = st_matrix(m1)
    b = st_matrix(m2)
    if (cols(a) > 1) {
        a = a[., 1]
    }
    if (cols(b) > 1) {
        b = b[., 1]
    }
    st_numscalar("__bc_n", rows(a))
    if (rows(a) != rows(b)) {
        st_numscalar("__bc_same", -1)
        return
    }
    st_numscalar("__bc_same", sum(a :== b))
}
// sorted elementwise comparison of two equal-length draw sets
void rxc_setcmp(string scalar m1, string scalar m2)
{
    real colvector a, b, rd
    a = st_matrix(m1)
    b = st_matrix(m2)
    if (cols(a) > 1) {
        a = a[., 1]
    }
    if (cols(b) > 1) {
        b = b[., 1]
    }
    a = sort(a, 1)
    b = sort(b, 1)
    rd = abs(a - b) :/ (abs(b) :+ 1)
    st_numscalar("__sc_n", rows(a))
    st_numscalar("__sc_bad", sum(rd :> 1e-10))
    st_numscalar("__sc_bit", sum(a :== b))
    st_numscalar("__sc_maxrel", max(rd))
}
end
mata: RXCREC = J(0, 0, .)

* The permute wrapper.  Requirements, each learned from a failure: sortpreserve,
* because frlink re-sorts the linked cluster frame and permute then aborts r(5)
* "changed sort order"; nodrop on the prefix, because e(sample) lives in the panel
* frame, not in the frame permute is permuting; never-treated recoded to the
* non-missing sentinel level 0 (permute's permvar cannot be missing) and mapped back
* here.  e() results are global across frames, so the wrapper's e(b) is wooldid's.
cap prog drop wd_cl
prog def wd_cl, eclass sortpreserve
  syntax varname
  mata: rxc_rec("i", "`varlist'", 1)
  frame rxpnl {
    cap drop __lk
    cap drop __ons
    cap drop ttre_p
    qui frlink m:1 i, frame(rxcl) generate(__lk)
    qui frget __ons = `varlist', from(__lk)
    qui gen double ttre_p = cond(__ons == 0, ., __ons)
    sort i t
    wooldid $RXCY i t ttre_p, att cluster(i) fastnose
  }
end

* per-draw command for ritest: an ordinary wooldid call on whatever assignment ritest
* wrote into the treatment variable, plus the recording hook
* ritest's permute_simple leaves the panel sorted by (strata, cluster-marker, rorder),
* a different row order on every draw; restoring the canonical i-t order before the
* estimate keeps the accumulation order the same as wooldid's own engine uses, so the
* comparison is numerical noise from the assignment and nothing else.  It changes no
* pass rule (membership is asserted at reldif 1e-10 either way), only how many draws
* land bitwise.
cap prog drop wd_rec
prog def wd_rec, eclass
  syntax varlist(min=4 max=4) [, *]
  tokenize `varlist'
  mata: rxc_rec("`2'", "`4'", 0)
  sort `2' `3'
  wooldid `1' `2' `3' `4', `options'
end

* ritest sampling program for permute(probabilities): one independent draw per
* cluster from the design's per-cluster marginal probabilities, copied down the
* cluster - which is what wooldid_permxp does
cap prog drop rxc_prob_sp
prog def rxc_prob_sp
  syntax , resampvar(varname) clus(varname) lev(real) plist(string) [*]
  tempvar u first
  sort `clus', stable
  by `clus': gen byte `first' = _n == 1
  gen double `u' = runiform() if `first'
  replace `resampvar' = .
  local n : word count `plist'
  forval j = 1(2)`n' {
    local c : word `j' of `plist'
    local p : word `=`j'+1' of `plist'
    replace `resampvar' = `lev' if `clus' == `c' & `first' & `u' <= `p'
  }
  by `clus': replace `resampvar' = `resampvar'[1]
end

timer clear

****************************************************************************
* X1 - EXACT: permute with enumerate on the cluster frame vs wooldid permexact.
* Tiny panel: 6 clusters, 3 at onset 4 and 3 never -> C(6,3) = 20 assignments, so
* BOTH sides enumerate the whole space and neither carries Monte-Carlo error.
* PASS (all numeric):
*   X1-1 r(n_reps) == r(n) == e(perm_npossible) == e(perm_completed) == 20
*   X1-2 sorted permute draws == sorted e(permdist_b) at reldif <= 1e-10 (bitwise
*        count reported, not asserted - saving() round-tripping is trap a/b)
*   X1-3 r(c_lower)/20 == e(perm_p_lo) and r(c_upper)/20 == e(perm_p_hi) at 1e-12
*   X1-4 r(p_twosided) == e(perm_p_att) at 1e-12 (they coincide exactly under
*        permexact, where wooldid drops the add-one and uses plain shares)
*   X1-5 reldif(r(b), e(perm_comparator_b)) <= 1e-13 (NOT bitwise: trap b)
****************************************************************************
clear
qui set obs 6
qui gen i = _n
qui gen ttre = cond(_n <= 3, 4, .)
qui expand 8
qui bysort i: gen t = _n
set seed 24601
qui gen y = 3*(t >= ttre & ttre < .) + 0.1*i + rnormal()
sort i t
tempfile TINYA
qui save `TINYA'

tempname PDX CMPX
scalar __x1ok = 0
timer on 1
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
local rc1 = _rc
timer off 1
if `rc1' != 0 {
  dis as error "X1-0 FAIL: wooldid permexact rc=" `rc1'
  global pf = $pf + 1
}
else {
  tempname CB PLO PHI
  matrix `PDX' = e(permdist_b)
  local cx = colnumb(`PDX', "main_att")
  if `cx' >= . local cx = 1
  matrix `PDX' = `PDX'[1..., `cx'..`cx']
  matrix `CB' = e(perm_comparator_b)
  local cb = colnumb(`CB', "main_att")
  if `cb' >= . local cb = 1
  scalar `CMPX' = `CB'[1, `cb']
  scalar __x1npos  = e(perm_npossible)
  scalar __x1ncomp = e(perm_completed)
  scalar __x1patt  = e(perm_p_att)
  matrix `PLO' = e(perm_p_lo)
  matrix `PHI' = e(perm_p_hi)
  local cl = colnumb(`PLO', "main_att")
  if `cl' >= . local cl = 1
  scalar __x1plo = `PLO'[1, `cl']
  scalar __x1phi = `PHI'[1, `cl']
  scalar __x1ok = 1
}
timer list 1
scalar __x1t = r(t1)

* the external side: Stata's own prefix, on the cluster-level frame
global RXCY y
use `TINYA', clear
cap frame drop rxpnl
cap frame drop rxcl
frame copy default rxpnl
frame put i ttre, into(rxcl)
frame rxcl {
  qui bysort i: keep if _n == 1
  qui gen double ons = cond(ttre >= ., 0, ttre)
  drop ttre
}
frame change rxcl
mata: rxc_reset()
tempfile ENUMF
scalar __x1eok = 0
timer on 2
cap noi permute ons _b[main_att], enumerate nodrop nodots nowarning saving(`ENUMF', replace double): wd_cl ons
local rc2 = _rc
timer off 2
if `rc2' == 0 {
  scalar __e_nreps = r(n_reps)
  scalar __e_n     = r(n)[1,1]
  scalar __e_clo   = r(c_lower)[1,1]
  scalar __e_chi   = r(c_upper)[1,1]
  scalar __e_plo   = r(p_lower)[1,1]
  scalar __e_phi   = r(p_upper)[1,1]
  scalar __e_p2    = r(p_twosided)[1,1]
  scalar __e_b     = r(b)[1,1]
  scalar __x1eok = 1
}
timer list 2
scalar __x1et = r(t2)
frame change default
cap frame drop rxpnl
cap frame drop rxcl

if `rc2' != 0 {
  dis as error "X1-0 FAIL: permute with enumerate rc=" `rc2'
  global pf = $pf + 1
}
if __x1ok == 1 & __x1eok == 1 {
  if __e_nreps != 20 | __e_n != 20 | __x1npos != 20 | __x1ncomp != 20 | rowsof(`PDX') != 20 {
    dis as error "X1-1 FAIL: counts permute n_reps=" __e_nreps " n=" __e_n " vs wooldid npossible=" __x1npos " completed=" __x1ncomp " permdist rows=" rowsof(`PDX') " (all want 20)"
    global pf = $pf + 1
  }
  else dis "X1-1 info: both engines enumerated the same 20 assignments"
  tempname EXD
  preserve
  qui use `ENUMF', clear
  mkmat _pm_1, matrix(`EXD')
  restore
  if rowsof(`EXD') != rowsof(`PDX') {
    dis as error "X1-2 FAIL: permute saved " rowsof(`EXD') " draws, e(permdist_b) has " rowsof(`PDX')
    global pf = $pf + 1
  }
  else {
    mata: rxc_setcmp(st_local("EXD"), st_local("PDX"))
    if __sc_bad != 0 {
      dis as error "X1-2 FAIL: " __sc_bad " of " __sc_n " enumerated draws differ from e(permdist_b) by more than reldif 1e-10 (max " %12.4e __sc_maxrel ")"
      global pf = $pf + 1
    }
    else dis "X1-2 info: enumerated draw sets agree, max reldif " %12.4e __sc_maxrel ", bitwise " __sc_bit "/" __sc_n
  }
  if reldif(__e_clo/20, __x1plo) > 1e-12 | reldif(__e_chi/20, __x1phi) > 1e-12 {
    dis as error "X1-3 FAIL: permute c_lower/20=" __e_clo/20 " c_upper/20=" __e_chi/20 " vs e(perm_p_lo)=" __x1plo " e(perm_p_hi)=" __x1phi
    global pf = $pf + 1
  }
  else dis "X1-3 info: one-sided shares agree: lo " __e_clo/20 " == " __x1plo ", hi " __e_chi/20 " == " __x1phi
  if reldif(__e_p2, __x1patt) > 1e-12 {
    dis as error "X1-4 FAIL: permute p_twosided=" __e_p2 " vs wooldid permexact e(perm_p_att)=" __x1patt " (reldif " %12.4e reldif(__e_p2, __x1patt) ", tol 1e-12)"
    global pf = $pf + 1
  }
  else dis "X1-4 info: two-sided p " __e_p2 " == " __x1patt " (permexact drops the add-one, so these coincide exactly)"
  if reldif(__e_b, `CMPX') > 1e-13 {
    dis as error "X1-5 FAIL: permute r(b)=" %21x __e_b " vs e(perm_comparator_b)=" %21x `CMPX' " reldif " %12.4e reldif(__e_b, `CMPX') " (tol 1e-13)"
    global pf = $pf + 1
  }
  else dis "X1-5 info: observed statistic reldif " %12.4e reldif(__e_b, `CMPX') " (decimal round trip; about 1 ulp is expected)"
  dis "X1 info: wooldid permexact " %6.2f __x1t "s, permute enumerate " %6.2f __x1et "s"
}
else dis as error "X1-1..X1-5 skipped: an X1 engine run failed (already counted)"

****************************************************************************
* X2 - permute SAMPLED on the same 20-assignment panel.
* PASS (all numeric):
*   X2-2 permute makes exactly 1 lead-in call before the replications
*   X2-1 every sampled draw within reldif 1e-10 of some member of the 20-value exact
*        set from X1 (bitwise count reported)
*   X2-3 never more than the 20 distinct cluster-level onset patterns the design
*        admits, and 0 multiset violations (every draw carries three onset-4 clusters
*        and three never-treated); ALL 20 cells hit whenever an empty cell is
*        essentially impossible at the realised count (20*(19/20)^reps < 1e-6)
*   X2-4 chi2(19) of the pattern frequencies < 43.82 (upper .001 point)
*   X2-5/6 tails against the exact-set counts recomputed under PERMUTE's OWN eps rule
*        (#{s <= s_obs + 1e-7}, #{s >= s_obs - 1e-7}), NOT wooldid's exact counts.
*        For a tail whose exact share is 0 or 1 the MC se is 0, so equality of the
*        count is asserted instead of running a degenerate 4*se z-test.
****************************************************************************
if __x1ok == 1 {
  use `TINYA', clear
  global RXCY y
  cap frame drop rxpnl
  cap frame drop rxcl
  frame copy default rxpnl
  frame put i ttre, into(rxcl)
  frame rxcl {
    qui bysort i: keep if _n == 1
    qui gen double ons = cond(ttre >= ., 0, ttre)
    drop ttre
  }
  frame change rxcl
  mata: rxc_reset()
  tempfile MCF
  timer on 3
  cap noi permute ons _b[main_att], reps(`RMC') rseed(7) nodrop nodots nowarning saving(`MCF', replace double): wd_cl ons
  local rc3 = _rc
  timer off 3
  if `rc3' == 0 {
    scalar __m_n   = r(n)[1,1]
    scalar __m_clo = r(c_lower)[1,1]
    scalar __m_chi = r(c_upper)[1,1]
  }
  timer list 3
  scalar __x2t = r(t3)
  mata: rxc_pat("X2", `RMC', (4\4\4\.\.\.), J(6,1,1), 20)
  frame change default
  cap frame drop rxpnl
  cap frame drop rxcl

  if `rc3' != 0 {
    dis as error "X2-0 FAIL: permute sampled rc=" `rc3'
    global pf = $pf + 1
  }
  else {
    if __X2_lead != 1 {
      dis as error "X2-2 FAIL: permute made " __X2_lead " lead-in call(s) before the `RMC' replications (expected exactly 1, the observed-data call); the pattern analysis below is misaligned"
      global pf = $pf + 1
    }
    tempname MCD
    preserve
    qui use `MCF', clear
    mkmat _pm_1, matrix(`MCD')
    restore
    mata: rxc_member(st_local("MCD"), st_local("PDX"), 1e-10)
    if __mem_bad != 0 {
      dis as error "X2-1 FAIL: " __mem_bad " of " __mem_n " sampled draws are not within reldif 1e-10 of the 20-value exact set (max " %12.4e __mem_maxrel ")"
      global pf = $pf + 1
    }
    else dis "X2-1 info: " __mem_n "/" __mem_n " sampled draws are members of the exact set (bitwise " __mem_bit "), max reldif " %12.4e __mem_maxrel
    * A 21st distinct pattern, or any draw whose cluster-level onset multiset is not
    * {4,4,4,.,.,.}, means the sampler left the assignment space: always a failure.
    * "ALL 20 cells were hit" is only a valid assertion when an empty cell is
    * essentially impossible - the chance is about 20*(19/20)^reps, which at the
    * shipped 500 reps is 2e-10 and at a development 60 reps is 0.64.  The shipped
    * rule is therefore unchanged; only a reduced-replication run relaxes it.
    local pempty = 20*(19/20)^`RMC'
    if __X2_nd > 20 | __X2_bad != 0 {
      dis as error "X2-3 FAIL: " __X2_nd " distinct cluster patterns (the design has 20) and " __X2_bad " multiset violations (want 0) over " __X2_R " draws"
      global pf = $pf + 1
    }
    else if `pempty' < 1e-6 & __X2_nd != 20 {
      dis as error "X2-3 FAIL: only " __X2_nd " of the 20 assignments appeared in " __X2_R " draws; the chance of an empty cell at this count is " %12.4e `pempty'
      global pf = $pf + 1
    }
    else dis "X2-3 info: " __X2_nd " distinct patterns of 20, 0 multiset violations, freq " __X2_minf ".." __X2_maxf " over " __X2_R " draws (P[some cell empty] = " %12.4e `pempty' ")"
    if __X2_chi2 >= 43.82 {
      dis as error "X2-4 FAIL: pattern-frequency chi2(19) = " %8.2f __X2_chi2 " >= 43.82 (upper .001 point)"
      global pf = $pf + 1
    }
    else dis "X2-4 info: pattern-frequency chi2(19) = " %8.2f __X2_chi2 " (.001 point 43.82)"
    mata: rxc_tails(st_local("PDX"), st_numscalar(st_local("CMPX")), 1e-7)
    local elo = __t_cloe/20
    local ehi = __t_chie/20
    forval side = 1/2 {
      if `side' == 1 {
        local nm lower
        local ex = `elo'
        local obc = __m_clo
      }
      else {
        local nm upper
        local ex = `ehi'
        local obc = __m_chi
      }
      local obp = `obc'/__m_n
      if `ex' <= 0 | `ex' >= 1 {
        if `obc' != `ex'*__m_n {
          dis as error "X2-5 FAIL: the `nm' tail share is a degenerate `ex' in the exact set, but permute counted " `obc' " of " __m_n
          global pf = $pf + 1
        }
        else dis "X2-5 info: `nm' tail share is exactly `ex' in the exact set and in all " __m_n " draws (z-test skipped: the MC se is 0)"
      }
      else {
        local se = sqrt(`ex'*(1 - `ex')/__m_n)
        if abs(`obp' - `ex') > 4*`se' {
          dis as error "X2-6 FAIL: `nm' tail " %8.5f `obp' " vs the exact set under the eps rule " %8.5f `ex' ", MC se " %8.5f `se' ", z " %6.2f (`obp' - `ex')/`se' " (tol 4 se)"
          global pf = $pf + 1
        }
        else dis "X2-6 info: `nm' tail " %8.5f `obp' " vs exact " %8.5f `ex' ", se " %8.5f `se' ", z " %6.2f (`obp' - `ex')/`se'
      }
    }
    dis "X2 info: permute sampled `RMC' reps, " %6.2f __x2t "s"
  }
}
else dis as error "X2 skipped: the X1 wooldid permexact run failed (already counted)"

****************************************************************************
* X6 - CONVENTION IDENTITIES on wooldid's own draws.  This block pins the documented
* mapping in place, so X3/X4 can compare like with like: it shows which wooldid
* p-value convention IS ritest's / permute's number and which is not.  On the
* 20-assignment tiny panel, so it costs seconds.
* PASS (every share at reldif <= 1e-12 against counts recomputed offline from
* e(permdist_b) and the exact-bit e(perm_comparator_b)):
*   exact mode  : default == min(1, 2*min(c_lo,c_hi)/M); permabsp == permrawp ==
*                 c_abs/M; permside(left) == c_lo/M; permside(right) == c_hi/M;
*                 e(perm_p_lo) == c_lo/M; e(perm_p_hi) == c_hi/M
*   sampled mode: default == min(1, 2*(1+min(c_lo,c_hi))/(n+1));
*                 permabsp == (1+c_abs)/(n+1); permrawp == c_abs/n  <- THE EXTERNAL
*                 ENGINES' NUMBER; permside(left/right) == (1+c_lo)/(n+1),
*                 (1+c_hi)/(n+1); e(perm_p_lo)/e(perm_p_hi) the same two
*   and in both modes every convention's e(permdist_b) is BITWISE identical to the
*   default run at the same rseed (a p-value convention must not move the draws)
****************************************************************************
use `TINYA', clear
tempname PDREF
local x6fail 0
timer on 5
foreach mode in exact samp {
  if "`mode'" == "exact" {
    local mopt permexact
    local mtag X6e
  }
  else {
    local mopt permreps(`RSMP') rseed(5150)
    local mtag X6s
  }
  local firstconv 1
  forval j = 1/5 {
    if `j' == 1 {
      local conv
      local cname default
    }
    if `j' == 2 {
      local conv permabsp
      local cname permabsp
    }
    if `j' == 3 {
      local conv permrawp
      local cname permrawp
    }
    if `j' == 4 {
      local conv permside(left)
      local cname permsideleft
    }
    if `j' == 5 {
      local conv permside(right)
      local cname permsideright
    }
    cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) `mopt' `conv'
    if _rc != 0 {
      dis as error "`mtag' FAIL: wooldid `mopt' `conv' rc=" _rc
      global pf = $pf + 1
      local ++x6fail
      continue
    }
    tempname PD6 CB6 L6 H6
    matrix `PD6' = e(permdist_b)
    local c6 = colnumb(`PD6', "main_att")
    if `c6' >= . local c6 = 1
    matrix `PD6' = `PD6'[1..., `c6'..`c6']
    matrix `CB6' = e(perm_comparator_b)
    local d6 = colnumb(`CB6', "main_att")
    if `d6' >= . local d6 = 1
    scalar __o6 = `CB6'[1, `d6']
    scalar __p6 = e(perm_p_att)
    matrix `L6' = e(perm_p_lo)
    matrix `H6' = e(perm_p_hi)
    local e6 = colnumb(`L6', "main_att")
    if `e6' >= . local e6 = 1
    scalar __plo6 = `L6'[1, `e6']
    scalar __phi6 = `H6'[1, `e6']
    mata: rxc_tails(st_local("PD6"), st_numscalar("__o6"), 0)
    local n6 = __t_n
    if "`mode'" == "exact" {
      if "`cname'" == "default"       local want = min(1, 2*min(__t_clo, __t_chi)/`n6')
      if "`cname'" == "permabsp"      local want = __t_cabs/`n6'
      if "`cname'" == "permrawp"      local want = __t_cabs/`n6'
      if "`cname'" == "permsideleft"  local want = __t_clo/`n6'
      if "`cname'" == "permsideright" local want = __t_chi/`n6'
      local wlo = __t_clo/`n6'
      local whi = __t_chi/`n6'
    }
    else {
      if "`cname'" == "default"       local want = min(1, 2*(1 + min(__t_clo, __t_chi))/(`n6' + 1))
      if "`cname'" == "permabsp"      local want = (1 + __t_cabs)/(`n6' + 1)
      if "`cname'" == "permrawp"      local want = __t_cabs/`n6'
      if "`cname'" == "permsideleft"  local want = (1 + __t_clo)/(`n6' + 1)
      if "`cname'" == "permsideright" local want = (1 + __t_chi)/(`n6' + 1)
      local wlo = (1 + __t_clo)/(`n6' + 1)
      local whi = (1 + __t_chi)/(`n6' + 1)
    }
    if reldif(__p6, `want') > 1e-12 {
      dis as error "`mtag' FAIL: `cname' e(perm_p_att)=" __p6 " but the documented share is " `want' " (n=" `n6' ", c_lo=" __t_clo ", c_hi=" __t_chi ", c_abs=" __t_cabs ")"
      global pf = $pf + 1
    }
    else dis "`mtag' info: `cname' p = " %9.6f __p6 " == documented share " %9.6f `want' " (n=" `n6' ", c_lo=" __t_clo ", c_hi=" __t_chi ", c_abs=" __t_cabs ")"
    if reldif(__plo6, `wlo') > 1e-12 | reldif(__phi6, `whi') > 1e-12 {
      dis as error "`mtag' FAIL: `cname' e(perm_p_lo)=" __plo6 " e(perm_p_hi)=" __phi6 " vs the documented " `wlo' " / " `whi'
      global pf = $pf + 1
    }
    if `firstconv' == 1 {
      matrix `PDREF' = `PD6'
      local firstconv 0
    }
    else {
      mata: rxc_bitcmp(st_local("PDREF"), st_local("PD6"))
      if __bc_same != __bc_n {
        dis as error "`mtag' FAIL: `cname' moved the draws - only " __bc_same " of " __bc_n " e(permdist_b) values are bitwise identical to the default run at the same rseed"
        global pf = $pf + 1
      }
      else dis "`mtag' info: `cname' draws bitwise identical to the default run (" __bc_n "/" __bc_n ")"
    }
  }
}
timer off 5
timer list 5
dis "X6 info: 10 convention runs, " %6.2f r(t5) "s"
if `x6fail' == 0 dis "X6 info: the external engines' number is permrawp (c_abs/n); wooldid's default and permabsp carry the add-one and cannot equal an external p except by coincidence"

****************************************************************************
* X5 - STRATA and PROBABILITIES against ritest.
* (i) permstrata(): 8 clusters, s=1 holds {4,4,.,.} (6 assignments), s=2 holds
*     {4,.,.,.} (4) -> 24 assignments.  wooldid permexact enumerates them; ritest
*     with cluster(i) strata(s) samples them.
*     PASS: X5-1 npossible == 24; X5-3 ritest makes exactly 2 lead-in calls;
*     X5-4 never more than 24 distinct patterns and 0 PER-STRATUM multiset
*     violations (all 24 hit when the count makes an empty cell impossible); X5-5 every
*     ritest draw within reldif 1e-10 of the exact 24-value set; X5-6 |r(p) -
*     c_abs/24| <= 4 MC se; X5-7 observed statistic at reldif 1e-13.
* (ii) permute(probabilities): 6 clusters with per-cluster onset probabilities
*     .8/.6/.4/.2/0/0; ritest runs the same table through a samplingprogram().
*     PASS: X5-10/X5-11 both engines' p within 4 MC se of the EXACT design tail
*     computed from the 15 estimable subsets and their product probabilities;
*     X5-12 per-cluster treated frequencies within 5 binomial sd of each cluster's
*     conditional marginal; X5-13 probability-0 clusters never treated.
****************************************************************************
if `HASRI' {
  clear
  qui set obs 8
  qui gen i = _n
  qui gen s = cond(_n <= 4, 1, 2)
  qui gen ttre = cond(inlist(_n, 1, 2, 5), 4, .)
  qui expand 8
  qui bysort i: gen t = _n
  set seed 777
  qui gen y = 3*(t >= ttre & ttre < .) + 0.1*i + rnormal()
  sort i t
  tempfile TINYS
  qui save `TINYS'

  tempname PDS CBS
  scalar __x5ok = 0
  cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permstrata(s)
  local rc5 = _rc
  if `rc5' != 0 {
    dis as error "X5-0 FAIL: wooldid permexact permstrata rc=" `rc5'
    global pf = $pf + 1
  }
  else {
    matrix `PDS' = e(permdist_b)
    local cs = colnumb(`PDS', "main_att")
    if `cs' >= . local cs = 1
    matrix `PDS' = `PDS'[1..., `cs'..`cs']
    matrix `CBS' = e(perm_comparator_b)
    local ds = colnumb(`CBS', "main_att")
    if `ds' >= . local ds = 1
    scalar __x5obs  = `CBS'[1, `ds']
    scalar __x5npos = e(perm_npossible)
    scalar __x5ok = 1
    if __x5npos != 24 | rowsof(`PDS') != 24 {
      dis as error "X5-1 FAIL: the stratified design has npossible=" __x5npos " and " rowsof(`PDS') " draws (want 24 / 24)"
      global pf = $pf + 1
    }
    else dis "X5-1 info: 24 stratified assignments enumerated"
  }
  if __x5ok == 1 {
    use `TINYS', clear
    mata: rxc_reset()
    tempfile STRF
    timer on 4
    cap noi ritest ttre _b[main_att], reps(`RSTR') cluster(i) strata(s) seed(303) nodots saving(`STRF', replace double): wd_rec y i t ttre, att cluster(i) fastnose
    local rc6 = _rc
    timer off 4
    if `rc6' == 0 {
      scalar __r5p = r(p)[1,1]
      scalar __r5c = r(c)[1,1]
      scalar __r5n = r(reps)[1,1]
      scalar __r5b = r(b)[1,1]
    }
    timer list 4
    scalar __x5t = r(t4)
    mata: rxc_pat("X5S", `RSTR', (4\4\.\.\4\.\.\.), (1\1\1\1\2\2\2\2), 24)
    if `rc6' != 0 {
      dis as error "X5-2 FAIL: ritest cluster(i) strata(s) rc=" `rc6'
      global pf = $pf + 1
    }
    else {
      if __X5S_lead != 2 {
        dis as error "X5-3 FAIL: ritest made " __X5S_lead " lead-in call(s) before the `RSTR' replications (expected exactly 2); the pattern analysis is misaligned"
        global pf = $pf + 1
      }
      * same reasoning as X2-3: a 25th pattern or a broken per-stratum multiset is
      * always a failure; "all 24 cells hit" is asserted only when an empty cell is
      * essentially impossible (24*(23/24)^720 = 1.2e-12 at the shipped count)
      local pempty5 = 24*(23/24)^`RSTR'
      if __X5S_nd > 24 | __X5S_bad != 0 {
        dis as error "X5-4 FAIL: " __X5S_nd " distinct stratified patterns (the design has 24), " __X5S_bad " per-stratum multiset violations (want 0) over " __X5S_R " draws"
        global pf = $pf + 1
      }
      else if `pempty5' < 1e-6 & __X5S_nd != 24 {
        dis as error "X5-4 FAIL: only " __X5S_nd " of the 24 stratified assignments appeared in " __X5S_R " draws; the chance of an empty cell at this count is " %12.4e `pempty5'
        global pf = $pf + 1
      }
      else dis "X5-4 info: " __X5S_nd " distinct patterns of 24, 0 per-stratum multiset violations, freq " __X5S_minf ".." __X5S_maxf " (P[some cell empty] = " %12.4e `pempty5' ")"
      tempname SRD
      preserve
      qui use `STRF', clear
      mkmat _pm_1, matrix(`SRD')
      restore
      mata: rxc_member(st_local("SRD"), st_local("PDS"), 1e-10)
      if __mem_bad != 0 {
        dis as error "X5-5 FAIL: " __mem_bad " of " __mem_n " ritest stratified draws are not within reldif 1e-10 of the exact 24-value set (max " %12.4e __mem_maxrel ")"
        global pf = $pf + 1
      }
      else dis "X5-5 info: " __mem_n "/" __mem_n " ritest stratified draws are members of the exact set (bitwise " __mem_bit ")"
      mata: rxc_tails(st_local("PDS"), st_numscalar("__x5obs"), 0)
      local ex5 = __t_cabs/24
      local se5 = sqrt(`ex5'*(1 - `ex5')/__r5n)
      if abs(__r5p - `ex5') > 4*`se5' {
        dis as error "X5-6 FAIL: ritest p " %8.5f __r5p " vs the exact design share " %8.5f `ex5' ", MC se " %8.5f `se5' ", z " %6.2f (__r5p - `ex5')/`se5' " (tol 4 se)"
        global pf = $pf + 1
      }
      else dis "X5-6 info: ritest p " %8.5f __r5p " vs exact " %8.5f `ex5' ", se " %8.5f `se5' ", z " %6.2f (__r5p - `ex5')/`se5'
      if reldif(__r5b, __x5obs) > 1e-13 {
        dis as error "X5-7 FAIL: ritest r(b)=" %21x __r5b " vs e(perm_comparator_b)=" %21x __x5obs " reldif " %12.4e reldif(__r5b, __x5obs)
        global pf = $pf + 1
      }
      else dis "X5-7 info: observed statistic reldif " %12.4e reldif(__r5b, __x5obs)
      dis "X5 info: ritest strata leg `RSTR' reps, " %6.2f __x5t "s"
    }
  }

  * ---- probabilities ----
  clear
  qui set obs 6
  qui gen i = _n
  qui gen ttre = cond(_n <= 3, 4, .)
  qui expand 8
  qui bysort i: gen t = _n
  set seed 2718
  qui gen y = 2*(t >= ttre & ttre < .) + 0.15*i + rnormal()
  sort i t
  tempfile TINYP
  qui save `TINYP'
  * external enumeration of the 15 estimable subsets of {1,2,3,4} (bitmask order) with
  * their design probabilities.  The all-never draw (q = .2*.4*.6*.8 = .0384) does not
  * estimate, so the design probabilities are renormalised by 1/(1-q).
  tempname QEXT QPR
  matrix `QEXT' = J(15, 1, .)
  matrix `QPR'  = J(15, 1, .)
  local pr1 = .8
  local pr2 = .6
  local pr3 = .4
  local pr4 = .2
  local qfail 0
  forval m = 1/15 {
    cap drop ttre2
    qui gen double ttre2 = .
    local pp = 1
    forval c = 1/4 {
      if mod(floor(`m'/2^(`c'-1)), 2) == 1 {
        qui replace ttre2 = 4 if i == `c'
        local pp = `pp'*`pr`c''
      }
      else local pp = `pp'*(1 - `pr`c'')
    }
    sort i t
    cap qui wooldid y i t ttre2, att cluster(i) fastnose
    if _rc != 0 {
      dis as error "X5-8 FAIL: the external probability subset `m' failed, rc=" _rc
      global pf = $pf + 1
      local ++qfail
    }
    else matrix `QEXT'[`m', 1] = _b[main_att]
    matrix `QPR'[`m', 1] = `pp'/(1 - .0384)
  }
  cap drop ttre2
  sort i t
  cap qui wooldid y i t ttre, att cluster(i) fastnose
  scalar __qobs = _b[main_att]
  mata: st_numscalar("__qtail", sum(select(st_matrix(st_local("QPR")), abs(st_matrix(st_local("QEXT"))) :>= abs(st_numscalar("__qobs")))))
  if `qfail' == 0 {
    timer on 6
    dis "X5-9 info: the exact design |s| tail is " %8.5f __qtail " over the 15 estimable subsets"
    use `TINYP', clear
    set seed 424242
    cap noi wooldid y i t ttre, att cluster(i) permute(probabilities) permreps(`RPRB') permrawp permprobs((ttre 4 .) (1 .8 .2) (2 .6 .4) (3 .4 .6) (4 .2 .8) (5 0 1) (6 0 1))
    local rc7 = _rc
    if `rc7' != 0 {
      dis as error "X5-10 FAIL: wooldid permute(probabilities) rc=" `rc7'
      global pf = $pf + 1
    }
    else {
      scalar __qwp = e(perm_p_att)
      scalar __qwn = e(perm_completed)
      local sew = sqrt(__qtail*(1 - __qtail)/__qwn)
      if abs(__qwp - __qtail) > 4*`sew' {
        dis as error "X5-10 FAIL: wooldid permrawp " %8.5f __qwp " vs the exact design tail " %8.5f __qtail ", se " %8.5f `sew' ", z " %6.2f (__qwp - __qtail)/`sew' " (tol 4 se)"
        global pf = $pf + 1
      }
      else dis "X5-10 info: wooldid permrawp " %8.5f __qwp " (n=" __qwn ") vs exact " %8.5f __qtail ", z " %6.2f (__qwp - __qtail)/`sew'
    }
    use `TINYP', clear
    mata: rxc_reset()
    tempfile PRBF
    cap noi ritest ttre _b[main_att], reps(`RPRB') seed(909) nodots samplingprogram(rxc_prob_sp) samplingprogramoptions("clus(i) lev(4) plist(1 .8 2 .6 3 .4 4 .2 5 0 6 0)") saving(`PRBF', replace double): wd_rec y i t ttre, att cluster(i) fastnose
    local rc8 = _rc
    if `rc8' == 0 {
      scalar __q5p = r(p)[1,1]
      scalar __q5n = r(reps)[1,1]
    }
    mata: rxc_pat("X5P", `RPRB', (4\4\4\.\.\.), J(6,1,1), 0)
    if `rc8' != 0 {
      dis as error "X5-11 FAIL: ritest samplingprogram rc=" `rc8'
      global pf = $pf + 1
    }
    else {
      if __X5P_lead != 2 {
        dis as error "X5-11 FAIL: ritest made " __X5P_lead " lead-in call(s) before the `RPRB' replications (expected exactly 2); the frequency analysis is misaligned"
        global pf = $pf + 1
      }
      local ser = sqrt(__qtail*(1 - __qtail)/__q5n)
      if abs(__q5p - __qtail) > 4*`ser' {
        dis as error "X5-11 FAIL: ritest p " %8.5f __q5p " vs the exact design tail " %8.5f __qtail ", se " %8.5f `ser' ", z " %6.2f (__q5p - __qtail)/`ser' " (tol 4 se)"
        global pf = $pf + 1
      }
      else dis "X5-11 info: ritest p " %8.5f __q5p " (n=" __q5n ") vs exact " %8.5f __qtail ", z " %6.2f (__q5p - __qtail)/`ser'
      mata: st_matrix("__qfrq", rowsum(editmissing(RXCREC[., (cols(RXCREC) - `RPRB' + 1)::cols(RXCREC)], 0) :> 0))
      forval c = 1/4 {
        local pc = `pr`c''
        local mu = `RPRB'*`pc'
        local sd = sqrt(`RPRB'*`pc'*(1 - `pc'))
        local ob = __qfrq[`c', 1]
        if abs(`ob' - `mu') > 5*`sd' {
          dis as error "X5-12 FAIL: ritest treated cluster `c' in " `ob' " of `RPRB' draws; expected " %8.2f `mu' " +- " %8.2f 5*`sd'
          global pf = $pf + 1
        }
        else dis "X5-12 info: ritest treated cluster `c' " `ob' "/`RPRB' (expected " %8.2f `mu' " +- " %8.2f 5*`sd' ")"
      }
      local anch 0
      forval c = 5/6 {
        if __qfrq[`c', 1] != 0 {
          dis as error "X5-13 FAIL: cluster `c' has probability 0 but ritest treated it in " __qfrq[`c', 1] " draws"
          global pf = $pf + 1
          local anch 1
        }
      }
      if `anch' == 0 dis "X5-13 info: both probability-0 clusters were never treated by ritest"
    }
    timer off 6
    timer list 6
    dis "X5 info: probabilities leg (wooldid `RPRB' + ritest `RPRB' reps), " %6.2f r(t6) "s"
  }
  else dis as error "X5-9..X5-13 skipped: the external probability enumeration failed (already counted)"
}
else dis "X5 info: SKIPPED (uncounted) - ritest is not installed"

****************************************************************************
* X3 / X4 on a real fixture: dta_binary (960 obs, 80 clusters, 12 periods) with a
* NULL outcome y0 = y - eff*(post).  The null outcome is not optional: on the
* effect-carrying outcome every engine bottoms out at its own floor (measured:
* wooldid .04 = 2/(n+1), ritest 0) and a p-agreement test has no content.
*
* X3 - EXACT (no Monte-Carlo tolerance at all).  wooldid runs permute(shuffle) and
*   writes its own draws with permsave(); those assignments are transposed into a
*   cluster-keyed file and fed back to ritest through samplingsourcefile(), so ritest
*   re-estimates THE SAME assignments with an ordinary wooldid call.  The two draw
*   sets must then agree rep by rep.
*   X3-2  the offline ddof=1 sd of e(permdist_b) == e(perm_sd) at reldif 1e-10
*   X3-3  NO draw lies within 1e-7 of |s_obs| - the precondition that makes wooldid's
*         exact tie rule and the engines' eps rule count the same draws.  If one did,
*         the p comparison would have to use #{|s| >= |s_obs| - 1e-7}/n instead.
*   X3-4  the permrawp rerun at the same rseed is BITWISE identical draw for draw
*   X3-5  its engine p == the offline c_abs/n at 1e-12
*   X3-6  the permsave file carries one onset_ column per cluster (wooldid documents
*         that a column can be silently dropped on a name collision)
*   X3-7  the transposed _rep row is 0,1,..,R
*   X3-8  rep-by-rep reldif <= 1e-10 (bitwise count reported)
*   X3-9  sd(ritest replay draws) == e(perm_sd) at reldif 1e-10
*   X3-10 reldif(r(b), e(perm_comparator_b)) <= 1e-13
*   X3-11 ritest p == the offline c_abs/n at 1e-12
*   X3-12 ritest p == wooldid's own permrawp ENGINE number at 1e-12 (engine to engine)
*
* X4 - MONTE-CARLO, independent draws.  ONE-SIDED tails only, against a binomial
*   standard error: for two independent draw sets estimating the same share,
*   se_diff = sqrt(pbar(1-pbar)(1/n1 + 1/n2)), pass at |diff| <= 4*se_diff.  SDs at
*   4*sqrt(sd1^2/(2(n1-1)) + sd2^2/(2(n2-1))).
*   X4-1/2 ritest |s| tail and null sd; X4-3 permute rc and lead-in; X4-4 every
*   permute draw holds the {25 x 2005, 25 x 2008, 30 never} onset multiset fixed;
*   X4-5 both permute tails; X4-6 permute null sd; X4-7 observed statistic.
****************************************************************************
local BIN `DATA'/dta_binary.dta
cap confirm file `"`BIN'"'
if _rc != 0 {
  dis as error "X3-0 FAIL: fixture not found: `BIN'"
  global pf = $pf + 1
}
else {
  use `"`BIN'"', clear
  qui gen double y0 = y - eff*(t >= ttre & ttre < .)
  qui keep i t ttre y0
  sort i t
  tempfile BINN
  qui save `BINN'
  qui levelsof i, local(ilist)
  local NCL : word count `ilist'
  dis "X3 info: dta_binary null outcome, `NCL' clusters, reps `RB'"

  tempname PDW CBW SDW
  scalar __w1ok = 0
  scalar __w2ok = 0
  capture erase rxc_wd.dta
  timer on 7
  cap noi wooldid y0 i t ttre, att cluster(i) permute(shuffle) permreps(`RB') rseed(11) permsave(rxc_wd) replace
  local rcw = _rc
  timer off 7
  if `rcw' != 0 {
    dis as error "X3-1 FAIL: wooldid permute(shuffle) permreps(`RB') rc=" `rcw'
    global pf = $pf + 1
  }
  else {
    matrix `PDW' = e(permdist_b)
    local cw = colnumb(`PDW', "main_att")
    if `cw' >= . local cw = 1
    matrix `PDW' = `PDW'[1..., `cw'..`cw']
    matrix `CBW' = e(perm_comparator_b)
    local dw = colnumb(`CBW', "main_att")
    if `dw' >= . local dw = 1
    scalar __wobs   = `CBW'[1, `dw']
    scalar __wpatt  = e(perm_p_att)
    scalar __wncomp = e(perm_completed)
    matrix `SDW' = e(perm_sd)
    local sw = colnumb(`SDW', "main_att")
    if `sw' >= . local sw = 1
    scalar __wsd = `SDW'[1, `sw']
    scalar __w1ok = 1
  }
  timer list 7
  scalar __x3t1 = r(t7)

  if __w1ok == 1 {
    mata: rxc_tails(st_local("PDW"), st_numscalar("__wobs"), 1e-7)
    scalar __wn    = __t_n
    scalar __wclo  = __t_clo
    scalar __wchi  = __t_chi
    scalar __wcabs = __t_cabs
    scalar __wsdo  = __t_sd
    dis "X3-2 info: wooldid draws n=" __wn ", c_lo=" __wclo ", c_hi=" __wchi ", c_abs=" __wcabs ", e(perm_p_att)=" __wpatt
    if reldif(__wsdo, __wsd) > 1e-10 {
      dis as error "X3-2 FAIL: the offline ddof=1 sd of e(permdist_b) " %18.15f __wsdo " != e(perm_sd) " %18.15f __wsd " (tol reldif 1e-10)"
      global pf = $pf + 1
    }
    else dis "X3-2 info: offline sd " %18.15f __wsdo " == e(perm_sd) " %18.15f __wsd
    if __t_mintie <= 1e-7 {
      dis as error "X3-3 FAIL: a draw lies " %12.4e __t_mintie " from |s_obs| (<= 1e-7), so wooldid's exact tie rule and the engines' eps rule no longer count the same draws and the p comparisons below are not like-for-like"
      global pf = $pf + 1
    }
    else dis "X3-3 info: the closest draw to |s_obs| is " %12.4e __t_mintie " away (> 1e-7); c_abs " __wcabs " == c_abs under eps " __t_cabse

    use `BINN', clear
    timer on 8
    cap noi wooldid y0 i t ttre, att cluster(i) permute(shuffle) permreps(`RB') rseed(11) permrawp
    local rcw2 = _rc
    timer off 8
    if `rcw2' != 0 {
      dis as error "X3-4 FAIL: the wooldid permrawp rerun failed, rc=" `rcw2'
      global pf = $pf + 1
    }
    else {
      tempname PDW2
      matrix `PDW2' = e(permdist_b)
      local cw2 = colnumb(`PDW2', "main_att")
      if `cw2' >= . local cw2 = 1
      matrix `PDW2' = `PDW2'[1..., `cw2'..`cw2']
      scalar __wraw = e(perm_p_att)
      scalar __w2ok = 1
      mata: rxc_bitcmp(st_local("PDW"), st_local("PDW2"))
      if __bc_same != __bc_n {
        dis as error "X3-4 FAIL: permrawp at the same rseed moved the draws - only " __bc_same " of " __bc_n " are bitwise identical"
        global pf = $pf + 1
      }
      else dis "X3-4 info: the permrawp rerun at rseed 11 is bitwise identical (" __bc_n "/" __bc_n ")"
      if reldif(__wraw, __wcabs/__wn) > 1e-12 {
        dis as error "X3-5 FAIL: e(perm_p_att) under permrawp = " __wraw " but c_abs/n = " __wcabs/__wn
        global pf = $pf + 1
      }
      else dis "X3-5 info: the permrawp engine p " %9.7f __wraw " == the offline c_abs/n " %9.7f __wcabs/__wn " (this is the number ritest and permute report)"
    }
    timer list 8
    scalar __x3t2 = r(t8)
  }

  local replayok 0
  if `HASRI' & __w1ok == 1 {
    * transpose wooldid's permsave file into the cluster-keyed replay file that
    * ritest's samplingsourcefile() wants: i + ttre1..ttreR.  permsave() takes a
    * PLAIN file name (replace is wooldid's own top-level option, NOT part of
    * permsave()), so permsave(name, replace) silently writes a file whose name
    * contains a comma.
    cap confirm file rxc_wd.dta
    if _rc != 0 {
      dis as error "X3-6 FAIL: the permsave() file rxc_wd.dta was not written"
      global pf = $pf + 1
    }
    else {
      use rxc_wd.dta, clear
      qui ds onset_*
      local nons : word count `r(varlist)'
      if `nons' != `NCL' {
        dis as error "X3-6 FAIL: the permsave file carries `nons' onset_ columns for `NCL' clusters (wooldid drops a column silently on a name collision)"
        global pf = $pf + 1
      }
      else {
        dis "X3-6 info: the permsave file carries `nons' onset_ columns, one per cluster"
        local replayok 1
      }
      if `replayok' == 1 {
        qui keep _rep onset_*
        qui xpose, clear varname
        local nv = `RB' + 1
        qui gen long i = real(subinstr(_varname, "onset_", "", .))
        cap assert v1 == 0 & v2 == 1 & v`nv' == `RB' if _varname == "_rep"
        if _rc != 0 {
          dis as error "X3-7 FAIL: after transposing, the permsave _rep row is not 0,1,..,`RB'"
          global pf = $pf + 1
          local replayok 0
        }
        else dis "X3-7 info: the transposed _rep row is 0,1,..,`RB'"
        qui drop if _varname == "_rep"
        qui drop v1
        forval k = 2/`nv' {
          qui rename v`k' ttre`=`k'-1'
        }
        qui keep i ttre1-ttre`RB'
        sort i
        tempfile RSF
        qui save `RSF'
      }
    }
  }
  if `replayok' == 1 {
    use `BINN', clear
    tempfile REPF
    timer on 9
    cap noi ritest ttre _b[main_att], reps(`RB') samplingsourcefile(`RSF') samplingmatchvar(i) nodots saving(`REPF', replace double): wooldid y0 i t ttre, att cluster(i) fastnose
    local rc9 = _rc
    timer off 9
    if `rc9' == 0 {
      scalar __rrp = r(p)[1,1]
      scalar __rrc = r(c)[1,1]
      scalar __rrn = r(reps)[1,1]
      scalar __rrb = r(b)[1,1]
    }
    timer list 9
    scalar __x3t3 = r(t9)
    if `rc9' != 0 {
      dis as error "X3-8 FAIL: the ritest samplingsourcefile replay failed, rc=" `rc9'
      global pf = $pf + 1
    }
    else {
      preserve
      qui use `REPF', clear
      mata: rxc_repcmp(st_local("PDW"))
      restore
      if __rc_n < 0 {
        dis as error "X3-8 FAIL: the ritest replay saved " __rc_na " draws, e(permdist_b) has " __rc_nb
        global pf = $pf + 1
      }
      else if __rc_bad != 0 {
        dis as error "X3-8 FAIL: " __rc_bad " of " __rc_n " replayed draws differ from wooldid's own by more than reldif 1e-10 (max " %12.4e __rc_maxrel ")"
        global pf = $pf + 1
      }
      else dis "X3-8 info: ritest replayed all " __rc_n " of wooldid's own draws, max reldif " %12.4e __rc_maxrel ", bitwise " __rc_bit "/" __rc_n
      if reldif(__rc_sd, __wsd) > 1e-10 {
        dis as error "X3-9 FAIL: the sd of the replayed draws " %18.15f __rc_sd " vs e(perm_sd) " %18.15f __wsd " (tol reldif 1e-10)"
        global pf = $pf + 1
      }
      else dis "X3-9 info: the replay sd " %18.15f __rc_sd " == e(perm_sd) " %18.15f __wsd
      if reldif(__rrb, __wobs) > 1e-13 {
        dis as error "X3-10 FAIL: ritest r(b)=" %21x __rrb " vs e(perm_comparator_b)=" %21x __wobs " reldif " %12.4e reldif(__rrb, __wobs) " (tol 1e-13; a decimal round trip of about 1 ulp is expected, more is not)"
        global pf = $pf + 1
      }
      else dis "X3-10 info: observed statistic reldif " %12.4e reldif(__rrb, __wobs)
      if reldif(__rrp, __wcabs/__wn) > 1e-12 {
        dis as error "X3-11 FAIL: ritest p " %18.15f __rrp " (c=" __rrc "/" __rrn ") vs wooldid's offline c_abs/n " %18.15f __wcabs/__wn
        global pf = $pf + 1
      }
      else dis "X3-11 info: ritest p " %18.15f __rrp " == the offline c_abs/n " %18.15f __wcabs/__wn " on the same draws"
      if __w2ok == 1 {
        if reldif(__rrp, __wraw) > 1e-12 {
          dis as error "X3-12 FAIL: ritest p " %18.15f __rrp " != wooldid's permrawp engine p " %18.15f __wraw " on the same draws (tol 1e-12)"
          global pf = $pf + 1
        }
        else dis "X3-12 info: ritest p == the wooldid permrawp engine p, " %18.15f __rrp " (engine to engine, no offline step)"
      }
      dis "X3 info: wooldid " %7.1f __x3t1 "s, permrawp rerun " %7.1f __x3t2 "s, ritest replay " %7.1f __x3t3 "s"
    }
  }
  else if !`HASRI' dis "X3 info: the ritest replay leg is SKIPPED (uncounted) - ritest is not installed"

  if __w1ok == 1 {
    if `HASRI' {
      use `BINN', clear
      tempfile RIF
      timer on 10
      cap noi ritest ttre _b[main_att], reps(`RB') cluster(i) seed(22) nodots saving(`RIF', replace double): wooldid y0 i t ttre, att cluster(i) fastnose
      local rc10 = _rc
      timer off 10
      if `rc10' == 0 {
        scalar __i4p = r(p)[1,1]
        scalar __i4c = r(c)[1,1]
        scalar __i4n = r(reps)[1,1]
      }
      timer list 10
      scalar __x4t1 = r(t10)
      if `rc10' != 0 {
        dis as error "X4-1 FAIL: ritest cluster(i) rc=" `rc10'
        global pf = $pf + 1
      }
      else {
        local pb = (__i4c + __wcabs)/(__i4n + __wn)
        local sd = sqrt(`pb'*(1 - `pb')*(1/__i4n + 1/__wn))
        local df = __i4p - __wcabs/__wn
        if abs(`df') > 4*`sd' {
          dis as error "X4-1 FAIL: the ritest |s| tail " %8.5f __i4p " vs wooldid c_abs/n " %8.5f __wcabs/__wn ", diff " %8.5f `df' ", se_diff " %8.5f `sd' ", z " %6.2f `df'/`sd' " (tol 4 se)"
          global pf = $pf + 1
        }
        else dis "X4-1 info: ritest |s| tail " %8.5f __i4p " vs wooldid " %8.5f __wcabs/__wn ", se_diff " %8.5f `sd' ", z " %6.2f `df'/`sd'
        preserve
        qui use `RIF', clear
        qui summ _pm_1
        scalar __i4sd = r(sd)
        restore
        local sse = sqrt(__i4sd^2/(2*(__i4n - 1)) + __wsd^2/(2*(__wn - 1)))
        if abs(__i4sd - __wsd) > 4*`sse' {
          dis as error "X4-2 FAIL: the ritest null sd " %9.6f __i4sd " vs e(perm_sd) " %9.6f __wsd ", se " %9.6f `sse' ", z " %6.2f (__i4sd - __wsd)/`sse' " (tol 4 se)"
          global pf = $pf + 1
        }
        else dis "X4-2 info: ritest null sd " %9.6f __i4sd " vs e(perm_sd) " %9.6f __wsd ", z " %6.2f (__i4sd - __wsd)/`sse'
        dis "X4 info: ritest independent leg " %7.1f __x4t1 "s"
      }
    }
    else dis "X4-1/X4-2 info: SKIPPED (uncounted) - ritest is not installed"

    use `BINN', clear
    global RXCY y0
    cap frame drop rxpnl
    cap frame drop rxcl
    frame copy default rxpnl
    frame put i ttre, into(rxcl)
    frame rxcl {
      qui bysort i: keep if _n == 1
      qui gen double ons = cond(ttre >= ., 0, ttre)
      drop ttre
    }
    frame change rxcl
    mata: rxc_reset()
    tempfile PMF
    timer on 11
    cap noi permute ons _b[main_att], reps(`RB') rseed(33) nodrop nodots nowarning saving(`PMF', replace double): wd_cl ons
    local rc11 = _rc
    timer off 11
    if `rc11' == 0 {
      scalar __p4n   = r(n)[1,1]
      scalar __p4clo = r(c_lower)[1,1]
      scalar __p4chi = r(c_upper)[1,1]
      scalar __p4b   = r(b)[1,1]
    }
    timer list 11
    scalar __x4t2 = r(t11)
    mata: rxc_pat("X4P", `RB', (J(25,1,2005) \ J(25,1,2008) \ J(30,1,.)), J(80,1,1), 0)
    frame change default
    cap frame drop rxpnl
    cap frame drop rxcl
    if `rc11' != 0 {
      dis as error "X4-3 FAIL: permute on dta_binary rc=" `rc11'
      global pf = $pf + 1
    }
    else {
      if __X4P_lead != 1 {
        dis as error "X4-3 FAIL: permute made " __X4P_lead " lead-in call(s) before the `RB' replications (expected exactly 1)"
        global pf = $pf + 1
      }
      if __X4P_bad != 0 {
        dis as error "X4-4 FAIL: " __X4P_bad " of " __X4P_R " permute draws violate the cluster-level onset multiset {25 x 2005, 25 x 2008, 30 never}"
        global pf = $pf + 1
      }
      else dis "X4-4 info: all " __X4P_R " permute draws hold the onset multiset fixed; " __X4P_nd " distinct patterns"
      forval side = 1/2 {
        if `side' == 1 {
          local nm lower
          local oc = __p4clo
          local wc = __wclo
        }
        else {
          local nm upper
          local oc = __p4chi
          local wc = __wchi
        }
        local pb = (`oc' + `wc')/(__p4n + __wn)
        local sd = sqrt(`pb'*(1 - `pb')*(1/__p4n + 1/__wn))
        local df = `oc'/__p4n - `wc'/__wn
        if abs(`df') > 4*`sd' {
          dis as error "X4-5 FAIL: the permute `nm' tail " %8.5f `oc'/__p4n " vs wooldid " %8.5f `wc'/__wn ", diff " %8.5f `df' ", se_diff " %8.5f `sd' ", z " %6.2f `df'/`sd' " (tol 4 se)"
          global pf = $pf + 1
        }
        else dis "X4-5 info: permute `nm' tail " %8.5f `oc'/__p4n " vs wooldid " %8.5f `wc'/__wn ", se_diff " %8.5f `sd' ", z " %6.2f `df'/`sd'
      }
      preserve
      qui use `PMF', clear
      qui summ _pm_1
      scalar __p4sd = r(sd)
      restore
      local sse = sqrt(__p4sd^2/(2*(__p4n - 1)) + __wsd^2/(2*(__wn - 1)))
      if abs(__p4sd - __wsd) > 4*`sse' {
        dis as error "X4-6 FAIL: the permute null sd " %9.6f __p4sd " vs e(perm_sd) " %9.6f __wsd ", se " %9.6f `sse' ", z " %6.2f (__p4sd - __wsd)/`sse' " (tol 4 se)"
        global pf = $pf + 1
      }
      else dis "X4-6 info: permute null sd " %9.6f __p4sd " vs e(perm_sd) " %9.6f __wsd ", z " %6.2f (__p4sd - __wsd)/`sse'
      if reldif(__p4b, __wobs) > 1e-13 {
        dis as error "X4-7 FAIL: permute r(b)=" %21x __p4b " vs e(perm_comparator_b)=" %21x __wobs " reldif " %12.4e reldif(__p4b, __wobs)
        global pf = $pf + 1
      }
      else dis "X4-7 info: observed statistic reldif " %12.4e reldif(__p4b, __wobs)
      dis "X4 info: permute leg " %7.1f __x4t2 "s"
    }
  }
  capture erase rxc_wd.dta
}

if $pf == 0 dis "RITEST_CROSSCHECK_ALL_PASS"
else dis as error "RITEST_CROSSCHECK_FAILS=$pf"
exit, clear
