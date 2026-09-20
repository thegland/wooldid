# Skeptic notes -- refuting PX-01 / PX-02 / PX-03 (permute x everything else)

## Approach
Own DGP (60 firms x 12 years 1998-2009, onsets 2002/2004/2006 + 15 never-treated,
unit FE + linear year trend + AR(1)-ish noise + heterogeneous dynamic effects).
Nothing from the project test suite, nothing from the accused agent's do-files.

## R1_px01.do -- PX-01 and PX-02 (DONE)
PX-01 REPRODUCES on my own data, exactly as claimed:
 - A permanalyticse alone       : rc 0, perm_analyticse=1, vcvstate=analytic-blocks,
                                  e(V)=1 e(V_es_att)=1 e(b_es_att)=1, se_analytic populated (.337...)
 - B fastnose + permanalyticse  : rc 0, perm_analyticse=1, vcvstate=ABSENT,
                                  e(V)=0 e(V_es_att)=0 e(b_es_att)=0, se_analytic column
                                  PRESENT but ENTIRELY BLANK.
                                  Header still: "se_analytic = the analytic cluster(firm)
                                  standard error, shown for reference only".
                                  Footer still: "The standard-error column is the ANALYTIC
                                  standard error (se_analytic, permanalyticse) ... e(V) and
                                  e(V_es_*) hold the analytic covariance, so test or lincom ..."
 - Estimates and p_perm in A and B are IDENTICAL -> no number is silently wrong.
 - D fastnose+permtstat         : r(198) refused
 - E fastnose+jointtests+permute: rc 0 ACCEPTED
 - F fastnose+permanalyticse+jointtests: r(198) REFUSED  <- the asymmetry they claim is real
 CORRECTION to their OBSERVED: the se_analytic COLUMN is not "entirely missing"; it is
 present in e(wdidmainresults) colnames and in the printed table, with all cells missing.

PX-02: all reachability attempts refused.
 - H cluster + customvce                : r(198) "customvce() is not compatible with robust and cluster()"
 - I no cluster + customvce             : r(198) "permute() requires cluster() with a single variable"
 - J cluster + customvce(robust) safety(off): r(198) same customvce message (guard is NOT safety-gated)
 - K customvce() empty                  : rc 0 (parses as absent; ado:128 still not emitted)
 - L customvce( ) single space          : (see log)
 Also grepped: validate never ASSIGNS customvce internally (only subinstr reads at 1418/1539),
 so there is no internal path that makes it non-empty. ado:128 is dead.

## R2_px03.do -- PX-03 (RUNNING)
Tests the accused pair AND separates the two candidate triggers:
 (a) raw outcome non-positive  vs  (b) model PREDICTION crosses zero.
`ymix` = strictly-positive predictions (mean ~12) with a few large negative raw values.

## Results in (2026-09-17)
R3  semielasticity DENOMINATOR oracle: dydx att is level-invariant (1.394985 at y+12 and
    y+120); eydx att is .10003225 at y+12 (hand: 1.394985/14.4316 = .09666) and .01143591
    at y+120 (hand: 1.394985/122.4316 = .011394).  So margins IS dividing by the FULL
    fitted value including the absorbed effects.  No wrong-denominator finding.
R4  (i) fastnose+permanalyticse+subgroups(): still no e(V), no e(V_es_att), no e(b_es_att),
        no e(numpre_att); the honest comparator (no fastnose) has all four.  No partial fill.
    (ii) section-7 item 11: plain `fastnose' DOES post e(b_es_att) and e(numpre_att)
        (and no e(V)/e(V_es_att)).  sthlp:1325's list is ACCURATE.  No doc-vs-code gap.
R5b LEVEL LADDER, min(y) vs rc for `semielasticity':
      -4.418 -> 459 ; -3.418 -> 459 ; -2.418 -> rc 0 att .23864 (se .0783, p .003) ;
      -1.418 -> rc 0 ; -0.418 -> rc 0 ; +0.582 -> rc 0 ; ... monotone to +5.582.
    => "the outcome takes non-positive values" is NOT the trigger.  Success band with
    negative y reports a meaningless fully-inferred semi-elasticity and no caution.
    Fitting att(c)=D/(m+c): D=1.0200, m=2.2092, residuals < 1.5% at every c.
    FRAME/state after a failing run is CLEAN (default frame, all user vars, N=720).
    [My earlier R5 r(111) was MY harness bug: `cap drop yc' abbreviates to `ycross'.]
R6  level() reaches the permute shift interval CORRECTLY: colnames flip to lb_90ci/ub_90ci,
    e(perm_level)=90, and reported lb/ub equal b - q_hi(draws) / b - q_lo(draws) EXACTLY at
    both 95 and 90 (hand-computed from e(permdist_b), 120 draws).
    nofastcollapse under permute(): estimate, p_perm, lb, ub bit-identical.  "no answer
    changes" holds.
R7  rc-0 HUNT: 8 independently drawn 40x8 two-onset panels, `semielasticity fastnose' on a
    sign-crossing outcome -> r(459) EVERY time.  The rc-0-with-an-all-missing-table branch
    did not reproduce (21 attempts in total across 3 panel shapes and 11 level shifts).
R8  permonsetcohorts (section-7 T1 machinery): cluster finer than i without it -> r(459)
    as documented; with it -> rc 0, e(perm_irebuilt)=1, e(perm_npossible_str)=2.8456e+33
    (= 60!/(15!)^4, exactly right), observed att bit-identical to the plain analytic run.
    At cluster()==i it is a bit-exact no-op with e(perm_irebuilt)=0.  Behaves as documented.
R9  (running) does the permute() design population follow the user's if-restriction?
R10 DECISIVE for SK-01: fitted values from a refit of wooldid's att-only specification
    (TWFE absorbing i,t with one dummy per treated (cohort, rel>=0) cell):
      shift 0: fitmin -3.362, 160 nonpos (17 treated) -> wooldid r(459)
      shift 1: fitmin -2.362,  58 nonpos ( 5 treated) -> wooldid r(459)
      shift 2: fitmin -1.362,  14 nonpos ( 3 treated) -> wooldid rc 0, att .23864, p .003
      shift 3: fitmin -0.362,   3 nonpos ( 0 treated) -> wooldid rc 0
      shift 4: fitmin +0.638,   0 nonpos             -> wooldid rc 0
    So the SUCCESS band at shift 2 really does average (dp/dx)/p over three rows with p<0.
    SK-01 stands (confidence probable: the fitted values are from my refit, not wooldid's
    own internal model, whose tempvars are gone by the time a user can predict).
R10 PART 2: `ereturn post B V' and `ereturn post B' with a missing element of B both fail
    r(504); a BARE `ereturn post' succeeds rc 0 with e(b) ABSENT. That is the ado:4082-4085
    fallback and it is the mechanism behind PX-03's rc-0/no-e(b)/e(cmd)=wooldid state.
    Under fastnose the only caution on that path is gated off (ado:4084), so the run is
    silent; otherwise the caution printed says "Point estimates remain in e(b)", which is
    false on exactly the branch that fired. -> SK-02.
FINAL: FINDINGS.md holds PX-01 CONFIRMED / PX-02 CONFIRMED(low) / PX-03 PARTIAL(medium)
    + SK-01 + SK-02 + nine negative results + the file map.
