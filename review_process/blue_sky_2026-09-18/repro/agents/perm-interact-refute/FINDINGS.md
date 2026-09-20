# Skeptic verdicts -- "permute() x everything else" (PX-01, PX-02, PX-03)

All evidence below is from MY OWN data and MY OWN do-files, in this directory.
I read the accused agent's do-files only to learn what conditions they used; I ran none
of them as evidence.

My DGP (identical across probes): 60 firms x 12 years (1998-2009); 15 never-treated and
15 each with onset 2002 / 2004 / 2006; firm effects N(0,1.4), a linear year trend, an
AR(1)-ish firm-level shock, and heterogeneous dynamic effects 0.35 + 0.18*rel + a
firm-block shifter. A plain, unremarkable staggered panel.

---

## PX-01 -- CONFIRMED (severity medium stands). Own repro: RF01_fastnose_permanalyticse.do

Reproduced in full on my own data (R1_px01.log):

| run | rc | e(perm_analyticse) | e(wdid_vcvstate) | e(V) | e(V_es_att) | e(b_es_att) | se_analytic(att) |
|---|---|---|---|---|---|---|---|
| `permanalyticse`             | 0 | 1 | analytic-blocks | yes | yes | yes | .33714813 |
| `permanalyticse fastnose`    | 0 | 1 | **absent**      | **no** | **no** | **no** | **. (blank)** |
| fast-mode default            | 0 | (absent) | absent | no | no | no | blank |
| `permtstat fastnose`         | **198 refused** | | | | | |

The accepted pair still prints the header
  "se_analytic = the analytic cluster(firm) standard error, shown for reference only"
above a column of dots, and the footer
  "The standard-error column is the ANALYTIC standard error (se_analytic, permanalyticse),
   shown for reference only; ... e(V) and e(V_es_*) hold the analytic covariance, so test
   or lincom run by hand on them are analytic, not permutation, quantities."
with no e(V) posted at all.

The asymmetry they cite is real and measured on my data:
 - `fastnose jointtests permute(...)`               -> rc 0, ACCEPTED (the `fnperm` stand-down)
 - `fastnose permanalyticse jointtests permute(...)`-> r(198), REFUSED
So the program already computes the correct predicate and already uses it to decide that
fastnose+permanalyticse is NOT the fast-mode default -- it just never applies that
decision to permanalyticse itself.

CORRECTION to their OBSERVED (does not change the verdict): the se_analytic **column is
not missing**. It is present in e(wdidmainresults) (colnames
`estimate se_analytic t_perm p_perm lb_95ci ub_95ci relyr`) and in the printed table, with
every cell missing. "Column of dots", not "no column".

CONFIRMED that no number moves: estimates and p_perm are bit-identical between
`permanalyticse` and `permanalyticse fastnose` (att 1.0306437, p_perm .06451613 in both).
The damage is a silent no-op, a false machine-readable e(perm_analyticse)=1 that a script
would branch on, and two printed sentences that instruct the reader to use objects that
do not exist.

Permute-only: the bootstrap twin refuses `fastnose` outright at ado:6405
("bootstrap() may not be combined with fastnose"), so this surface does not exist there.

---

## PX-02 -- CONFIRMED (severity low; arguably none). Own repro: RF02_customvce_permute.do

Every route I could construct to permute()+customvce() is refused, on my own data:

| attempt | rc | message |
|---|---|---|
| `cluster(firm) customvce(cluster firm)` | 198 | "customvce() is not compatible with robust and cluster()..." (ado:4667) |
| no cluster, `customvce(cluster firm)`   | 198 | "permute() requires cluster() with a single variable..." (ado:5947) |
| `cluster(firm) customvce(robust) safety(off)` | 198 | same customvce message -- the guard is NOT safety-gated |
| `customvce()` (empty)                   | 0 | parses as absent; ado:128 is not emitted |
| `customvce( )` (a space)                | 0 | same |

I also checked the only remaining escape hatch -- an INTERNAL assignment that could make
`customvce` non-empty without the user typing it. There is none: the only occurrences of
the macro outside the syntax/forwarding lines are two `subinstr ... local customvce` READS
at ado:1418 and ado:1539. ado:128 is dead code.

So the factual claim is right. As a defect it is a cosmetic one with no user-visible
consequence today, which is the brief's lowest category; I would file it as `low` at most
and would not call it a defect at all were it not a correction to the seeded lead.

---


## PX-03 -- PARTIAL (severity medium, down from high). Own repro: RF03_semielasticity.do
Evidence: R2_px03.log, R5b_semi_trigger.log, R7_rc0_hunt.log, R10_predsign.log, RF03.log

WHAT IS TRUE, on my own data. `semielasticity` on an outcome whose fitted values reach far
enough below zero fails with a bare, undiagnosed message, and wooldid says nothing in its
own words:
 - `semielasticity` alone      -> r(459) "could not calculate numerical derivatives --
                                  discontinuous region with missing values encountered"
 - `semielasticity fastnose`   -> r(459) "inconsistent estimation sample levels 0 and 1 of
                                  factor __000007"
 - `semielasticity permute()`  -> r(459), the SAME bare margins message
 - e(cmd) is left as `reghdfe` with tempvar-named coefficients on all three.
 - shifting the outcome up makes all three rc 0.

THREE CORRECTIONS.

(1) THE STATED TRIGGER IS WRONG. It is not "the outcome takes non-positive values"; it is
    how much mass the model's FITTED VALUES put at or below zero. Level ladder on one
    outcome (R5b), with the fitted values from a refit of wooldid's own att-only
    specification -- TWFE absorbing i and t with one dummy per treated (cohort, relative
    period) cell (R10):

      shift  min(y)    min(fitted)   fitted<=0 (all / treated)   wooldid semielasticity
        0    -4.418      -3.362         160 / 17                 r(459)
        1    -3.418      -2.362          58 /  5                 r(459)
        2    -2.418      -1.362          14 /  3                 rc 0, att .23864029,
                                                                 se .0783198, p 0.003
        3    -1.418      -0.362           3 /  0                 rc 0, att .19509366
        4    -0.418      +0.638           0 /  0                 rc 0, att .16414207
        5    +0.582      +1.638           0 /  0                 rc 0, att .14149061

    At shift 2 the raw outcome is non-positive on 30 of 720 rows AND the fitted value is
    non-positive on 14 of them, three of which are inside the treated-and-post subpop the
    att is averaged over -- and the run completes and reports a fully-inferred
    semi-elasticity with no caution of any kind. A non-positive outcome is therefore
    neither necessary nor sufficient for the failure the finding is built on.

(2) THE SUCCESS BAND IS ITS OWN, UNFILED DEFECT -- see SK-01 below. The finding stops at
    the crash; the more interesting behaviour is on the other side of the threshold.

(3) TWO SPECIFIC OBSERVATIONS DID NOT REPRODUCE, AND ONE OF THEM I CAN NOW EXPLAIN.
    - "with fastnose it returns rc 0 ... and posts NO e(b) at all": on my panels
      `semielasticity fastnose` ALWAYS returned r(459) -- 21 attempts across three panel
      shapes, eight independently drawn 40x8 panels and eleven level shifts. I could not
      hit the branch. But it is not a harness artefact of theirs: the branch is
      ado:4051-4085 and I verified its premise directly (R10 PART 2) --
      `ereturn post` with a coefficient vector containing a missing value fails r(504),
      and the fallback at ado:4083 is a BARE `ereturn post` with no vector at all, which
      succeeds and leaves e() with no e(b). See SK-02 below; that is the mechanism, and the
      reason the run is silent about it under fastnose.
    - "Under permute() ... r(459) 'the observed comparison statistic for att could not be
      read from the rep-0 run', blaming the comparator": on my data the permute() run died
      with the same bare margins message as the plain run, before the permutation engine
      was reached. The wrong-object message did not reproduce.

SEVERITY. What survives is an undiagnosed bare margins error on erroneous-but-plausible
input, an e(cmd) left as reghdfe, and a help entry whose only hedge ("tends to require
large samples and a y-variable that is not sparsely distributed") never says the prediction
must be positive. That is brief categories 4 and 5, i.e. medium, not high.

---

# NEW FINDINGS (mine; neither is in the seeded leads -- SK-01 is the negation of PX-03's
# stated trigger, and SK-02 is the mechanism PX-03 observed but did not identify)

## SK-01 [medium/probable] semielasticity averages (dp/dx)/p over a subpopulation in which
## some p are negative, and reports the result with se, t, p and a 95% interval, unscreened
Repro: RF03_semielasticity.do (section A) + R10_predsign.do (the fitted-value counts).

CLAIM. `semielasticity` asks margins for eydx, which divides by the prediction. wooldid
never checks the sign of what it divides by. On a panel where 14 of 720 fitted values are
non-positive -- three of them inside the treated-and-post subpopulation the att averages
over -- the run prints

    Main Results (Full Estimation Sample): Treatment Effects Presented as Semi-elasticities
             att |  .23864029    .0783198    3.05   0.003   .08192275   .39535784

with no caution. What I measured is narrower than a claim about margins' internals: the
averaging set contains observations whose fitted value is negative, the run returns a finite
fully-inferred number anyway, and nothing in the output mentions it. A "% change in y" is not
defined at those rows, so the headline figure is contaminated by them to an extent the user
is given no way to see.

OBSERVED / EXPECTED. The counts above (R10_predsign.log). The contamination is small on
this panel -- three rows at a prediction of about -1.4 against a median of +3.3 move the
average by roughly 3-4%, which is why the ladder still fits att(c) = D/(m+c) to within 1.5%
-- but nothing bounds it: the same code path with more mass near zero produces an
arbitrarily corrupted number, and the only thing standing between the user and that is the
point where margins happens to give up with a bare error instead. wooldid screens and names
comparable conditions in its own words elsewhere (the poisson perfect-separation message,
the dose-slope screen), and the cell-collapse caution that fires on THIS VERY RUN already
says in as many words that "semielasticity ... divide[s] by the prediction". Having
identified the division, the run does not look at the sign of the divisor.

CONFIDENCE: probable rather than certain. The fitted values are from my refit of wooldid's
att-only specification, not from wooldid's own internal model (its tempvars are gone by the
time the user can predict); the specification is the same one, but I did not read the
numbers out of wooldid itself.

REALISTIC: yes for a z-scored, demeaned, differenced or profit/margin outcome. Asking for a
semi-elasticity of such an outcome is a user error; reporting an inferred number for it
without a word is the program's.

## SK-02 [medium/probable] when the estimand vector contains any missing coefficient the
## run posts NO e(b) at all, returns rc 0, and either says nothing (fastnose) or prints a
## caution asserting the opposite of what happened
Repro: R10_predsign.do PART 2 (the ereturn semantics), plus the code path below.

CLAIM. ado:4051-4085 is a three-step fallback:
    cap ereturn post `eb' `ev', ...        <- fails if ANY variance cell is missing
    ... partial-identification rescue, taken only if at least one variance is identified ...
    cap ereturn post `eb' , ...            <- fails if ANY COEFFICIENT is missing
    if _rc != 0 ereturn post   , ...       <- BARE: posts an e() with no coefficients
    if `evpostmiss' > 0 & "`fastnose'" == "" dis "... Point estimates remain in e(b); ..."
Measured directly (R10 PART 2): `ereturn post B V` and `ereturn post B` with a missing
element of B both fail r(504); a bare `ereturn post` succeeds with rc 0 and leaves e(b)
absent. So on that third step the run completes successfully having posted no coefficients,
and:
  (a) under `fastnose` the only caution on the path is gated OFF, so nothing is printed at
      all -- the run returns rc 0, prints its formatted table, and has no e(b);
  (b) otherwise the caution that does print says "so no variance-covariance matrix e(V) was
      posted. Point estimates remain in e(b)", which is false precisely on the branch that
      just fired.
This is exactly the rc-0/no-e(b)/e(cmd)=wooldid state PX-03 reported and I could not
reproduce; the branch is real, its silence under fastnose is structural, and semielasticity
is only one way to get a missing coefficient into the vector.

CONFIDENCE: probable. The code path and the `ereturn post` semantics are certain (measured
in R10 PART 2, which demonstrates the PREMISE, not the full run); what I did not do is
manufacture a missing coefficient on my own data. The end-to-end corroboration is the other
agent's own log, probes/perm-interact/PX03_repro.log, which records exactly the state this
path produces: rc 0, e(cmd)=wooldid, no e(b), a printed all-missing table. The author's own source comment at ado:5017 names "the data-dependent
r(459)/missing-with-rc-0 behaviour" in passing, so the state is known to exist.

EXPECTED. A run that posted no coefficients should not return rc 0, and must not print a
caution that tells the reader the coefficients are in e(b). At minimum the caution's text
should branch on which post actually succeeded, and the `fastnose` gate should suppress only
the standard-error half of the sentence, not the whole disclosure.

---

# NEGATIVE RESULTS (probed hard, found clean)

1. **level() under permute()** -- verified against a hand oracle, not read off the header.
   `level(90)` flips the result colnames to `lb_90ci ub_90ci`, sets e(perm_level)=90, and the
   reported interval equals b - q_{hi}(draws) / b - q_{lo}(draws) computed by hand from
   e(permdist_b) (120 draws) EXACTLY at both 95 and 90. (R6_perm_oracle.log)
2. **nofastcollapse under permute()** -- the cell-collapse caution promises "no answer
   changes". Same rseed, both ways: estimate, p_perm, lb_95ci, ub_95ci bit-identical in hex.
   (R6_perm_oracle.log)
3. **permonsetcohorts (the T1 machinery known_bydesign section 7 flags as one day old and
   unassessed)** -- on i = cohort / cluster = unit: refused r(459) without the option as
   documented; with it rc 0, e(perm_irebuilt)=1, observed att bit-identical to the plain
   analytic run, and e(perm_npossible_str) = 2.8456e+33 = 60!/(15!)^4 exactly. At
   cluster()==i it is a bit-exact no-op with e(perm_irebuilt)=0. (R8_onsetcohorts.log)
4. **The permute() design population follows the estimation sample** -- oracle on the
   multinomial assignment count. Unrestricted: N=600, 2.8456e+33 = 60!/(15!)^4. With
   `if region==1`: N=400, 4.7054e+21 = 40!/(10!)^4. With a bare 0/1/2 flag: N=400,
   4.7054e+21, e(N_dropped)=200 -- BS-37's "exactly 1" rule reaches the randomization design
   too. (R9_perm_ifpop.log)
5. **known_bydesign section 7, item 11** -- plain `fastnose` DOES post e(b_es_att) and
   e(numpre_att) (and withholds e(V)/e(V_es_att)). sthlp:1325's list is accurate; no live
   doc-vs-code divergence there. (R4_px01_sharp.log)
6. **The semielasticity denominator** -- margins on wooldid's reghdfe fit divides by the FULL
   fitted value including the absorbed effects, not by an FE-free index: dydx is
   level-invariant (1.394985 at y+12 and y+120) while eydx is .10003225 and .01143591, both
   within 3.5% of dydx / mean(y | treated, post). (R3_semi_oracle.log)
7. **The permute() forwarding allow-list** -- I audited every option in the worker's syntax
   line against the rep string at ado:125-172. Everything that shapes a point estimate is
   either forwarded or refused: `contreatwithin` (ado:5832), `contreatcontrols(untreatedZero)`,
   `impvalues`, `unconditionalse`, `oldsyntax`, `pretrendpower()` and
   `suppresssgprimaryeffects` are refused in wooldid_permctl; `semielasticity` is forwarded
   only while it is still "eydx", which is right because the poisson-xb rewrite at ado:5034
   re-fires identically inside the replication; `esprelength` is emitted RAW (so validate's
   re-firing normalization cannot widen the window per nesting level) while `espostlength` is
   emitted post-normalization, which is also right because its only mutation (ado:4856) is
   self-cancelling on re-parse. The one dead entry is `customvce` (PX-02). I found no
   unforwarded option that would silently give the replications a different design from the
   observed run.
8. **State after a failed semielasticity run** -- the caller is left in the `default` frame
   with every user variable intact and N unchanged; the only residue is e(cmd)==reghdfe with
   tempvar-named coefficients, wooldid's generic post-failure e() state, guarded for replay by
   WOOLDID_LASTRC. (R5b_semi_trigger.log)
9. **fastnose + permanalyticse + subgroups()** -- the subgroup-by-horizon margins call keeps
   df(df_r) even under fastnose (ado:1981), so I checked whether the "shown for reference
   only" column comes back PARTIALLY filled. It does not: no e(V), no e(V_es_att), no
   e(b_es_att), no e(numpre_att), column still all dots. (R4_px01_sharp.log)

# FILE MAP
  R1_px01.log         PX-01 + PX-02 measurement
  R2_px03.log         PX-03 first pass, four paths x three outcomes
  R3_semi_oracle.log  semielasticity denominator oracle
  R4_px01_sharp.log   PX-01 subgroups sharpening + section-7 item 11
  R5b_semi_trigger.log level ladder + post-failure frame/variable state
  R6_perm_oracle.log  level() hand oracle + nofastcollapse invariance
  R7_rc0_hunt.log     eight panels hunting the rc-0 branch
  R8_onsetcohorts.log permonsetcohorts
  R9_perm_ifpop.log   permute design population vs the if-restriction
  R10_predsign.log    fitted-value signs in the success band + ereturn post semantics
  RF01/RF02/RF03      the three self-contained repro files (all run clean, EXIT 0)
