# hv-noparent-refute — SKEPTIC verdicts

Every repro is self-contained, generates its own data, asserts the pinned build with
`which wooldid`. The project test suite was not used.

## hv-noparent-F1 — PARTIAL (severity revised medium -> low)
Repro: r_f1_noparent.do  (+ r_f3b_semantics.do Part 3 for the paired validation leg)

BEHAVIOUR HALF CONFIRMED. Nine paired legs, identical data: all rc 0,
mreldif(e(wdidmainresults), baseline) = 0 for every one, nothing printed beyond the
ordinary table. saveplots(p1) and `makeplots saveplots(p2)` wrote NO files; the control
`makeplots esprelength(2) espostlength(2) saveplots(p3)` wrote p3_ES.dta and p3_ES_att.gph.

DOC HALF REFUTED FOR FIVE OF THE NINE. Explicitly documented as inert and silent:
  ccnoabsorb        sthlp:1121 "the option is inert and silent on a run with no continuous treatment"
  contreatpoly()    sthlp:917  "without contreat() the option is otherwise inert: nothing reads it and nothing is said about it"
  cfxplotnose       sthlp:978  "Without makecfxplots the option is inert: nothing reads it, and no message is printed"
  makecfxplotsbysg  sthlp:976  "silently inert without either: no error and no caution"
  makeplots (no horizon) sthlp:807 "No plot is drawn unless at least one event-study horizon is requested"
The finding conceded only two. Five explicit entries by the same author in the same file
make sthlp:1181 an overstatement of its own scope, not a promise these options break; its
closing clause scopes it to runs where "the run changed your specification".
The four remaining members (saveplots() alone, suppresssgprimaryeffects, cfxplottypes(),
latticeignoreweights) are inert by implication only — a doc-internal inconsistency, no
number moves.

WHAT SURVIVES AT FULL STRENGTH IS A DIFFERENT CLAIM: the cfxplottypes() VALIDATION GAP.
The token validator sits inside `if "makecfxplots" != ""` (ado:5188). Same data, same
options, one difference:
  ... contreat(dose) contreatcontrols(untreatedZero) cfxlattice(dose0) cfxplottypes(psot) safety(off)  -> rc 0, silent
  ... same + makecfxplots -> r(198) "cfxplottypes() contains an unrecognized entry: psot"
Without contreat() at all, cfxplottypes(psot), (1.5), (nosuchvariable), (postpost 2post2)
are all rc 0 and silent. sthlp:974 states the rejection with NO parent qualifier, and the
sibling entry sthlp:917 says contreatpoly() "is validated whether or not contreat() is
given". ado:5188's own VAL-G13A comment records that garbage tokens used to be swallowed;
the fix was placed inside the parent block, so without the parent the pre-fix behaviour is
still live. Doc-vs-code, medium-low, no numeric consequence.

## hv-noparent-F2 — CONFIRMED (medium stands), three surface corrections
Repro: r_f2_cet.do

Four runs, identical binary data, contreat() never given. att hex +1.d169f113af0c8X-001 in
all four; mreldif(e(wdidmainresults)) = 0 for eyex/eydx/dyex. Header moves from
"Main Results (Full Estimation Sample):" to "... Marginal Treatment Effects Presented as
Elasticities / Semi-elasticities (eydx) / (dyex)". e(wdid_cetype) dydx -> the request.

CORRECTION 1 (against the finding): the ES GRAPH TITLES are NOT relabeled. Both saved .gph
files carry ytitle("ATT Effects on y by Horizon"). ado:13308-13310 sets elasinj2, which
feeds only ame_it_text/ame_text — the AME graphs, absent without contreat(). The ATT graph
title comes from elasinj (13305-06), driven by semielasticity alone. The ES TABLE is
relabeled only because its rows sit under the same Main Results header (ado:11368).

CORRECTION 2 (against the finding): histogramcohorteffects is NOT relabeled. Ceffecttype
(ado:1930-33) is consumed only at ado:12391/12428 and 2506-2602, each gated on
`if "contreat" != ""`.

CORRECTION 3 (the finding UNDERSTATES): the relabel also reaches the per-subgroup result
headers (ado:11528) and the pairwise subgroup comparison header (ado:11588):
  "Results for Subgroup sgv == 0: Marginal Treatment Effects Presented as Elasticities"
  "Pairwise Comparisons of Different Subgroups' Estimates: ... Elasticities"

THE FINDING'S UNTESTED POISSON PREDICTION IS CORRECT: poisson without poisexpresults, no
contreat(), contreatelasticitytype(eyex) prints "Main Poisson Regression Results (Full
Estimation Sample): Marginal Treatment Effects Presented as Elasticities" with NO BS-6
caution, att .28180154 bit-identical to the plain poisson run, e(wdid_cetyperaw) not posted.

## hv-noparent-F3 — PARTIAL (core CONFIRMED; medium; confidence certain), two corrections
Repro: r_f3e_oracle.do (margins semantics, decisive); r_f3b_semantics.do (wooldid on an
unambiguous straddle); r_f3c_pole.do (the knife-edge); r_f3d_scope.do (scope legs)

CORE CONFIRMED. Panel Q (r_f3b_semantics.log): 80 units x 10 periods, unit effects
U(-4,4), resid sd 0.2, R2 = 0.997, min y over the treated-post subpop = -2.87 — so the
fitted values in the region margins evaluates unambiguously straddle zero.
  wooldid y unit t ttre, att robust                  -> .83196634
  wooldid y unit t ttre, att robust semielasticity   -> rc 0, 1.0437756 se .0431678
                                                        t 24.18 p 0.000, printed under
                                                        "Treatment Effects Presented as
                                                        Semi-elasticities", NO domain caution
Weighted legs on the same panel ([aw=w] and [pw=w], weights written BEFORE the comma):
rc 0, 1.0842342. att_it on a straddling panel (r_f3d_scope.log): rc 0, 1.2005979.
No screen exists in the code and sthlp:243 warns only about "large samples and a
y-variable that is not sparsely distributed" — nothing about the sign of the fitted outcome.

THE FINDING'S OWN OPEN QUESTION IS ANSWERED: (a), not (b). Standalone oracle, no wooldid
(r_f3e_oracle.log), on a domain that genuinely straddles zero — min yhat -1.990, max
+4.013, 70 of 200 subpop rows with yhat < 0:
  hand-computed mean(_b[x]/yhat) over the subpop  = -.24488522
  ratio-of-means alternative                      =  1.00068
  margins, eydx(x) subpop(if x==1)                = -.24488522, rc 0
Exact match. margins' eydx is the MEAN OF PER-OBSERVATION RATIOS and drops nothing. On a
tighter straddle (min yhat -.301, max +.696) it again matched exactly (.67907475, rc 0).
Note what that means: a treatment whose true effect is +1.0, on an outcome whose fitted
values average +1.0, is reported as a semi-elasticity of -0.24. The sign is wrong, the
magnitude is meaningless, and margins raises nothing.

CORRECTION 1: the r(3598) is NOT "the weighted leg". Every r(3598) I obtained is
UNWEIGHTED (r_f3c_pole.log, shifts -8.4 through -20, no weights anywhere in the run).
The prober's attribution to an aweight does not hold.

CORRECTION 2: neither party has characterized what selects among the three behaviours,
and MY first attempt at a rule was wrong too. Depth of the straddle is NOT the
discriminator: panel Q reaches roughly -2.7 in fitted value and returns a number, while
pole shift c = -8.4 (min y over treated-post only -0.745) dies r(3598). What the two
crashing configurations have in common is that the fitted values are tightly CLUSTERED
near zero rather than spread across it, which is where a numerical derivative of
ln(yhat) is worst behaved — the traceback is Stata's own numerical-derivative machinery
(_deriv() / __marg_cr_eval()), and my plain-reg oracle never crashes because that path
differentiates analytically. That is a hypothesis consistent with both datasets, not a
measured rule. The honest statement is: rc 0 with a number, r(459) and r(3598) all occur
on ordinary panels, and which one you get is data-dependent and uncharacterized.

CORRECTION 3 (knife-edge, r_f3c_pole.log): one panel, one model, one additive shift c on
y; the level ATT is .7543601 at EVERY c, so the model never changes.
  c = 0     min y over treated-post +7.655 -> .08617444 (se .0053) rc 0
  c = -6                           +1.655 -> .27606125 (se .0171) rc 0
  c = -8                           -0.345 -> 1.1572378 (se .0864, t 13.39) rc 0
  c = -8.4                         -0.745 -> r(3598) "variable __00001O not found" /
                                             st_data(): 3598 / __marg_cr_eval() /
                                             deriv__call2user_v()
  c = -8.7, -9, -10, -20                  -> r(3598)
0.4 in the units of y moves the run from a confident complete inference row to a raw Mata
traceback naming an internal tempvar. (At c = -8, min y over treated-post is only -0.345
against a residual sd of 0.3, so I did not measure whether its fitted values cross zero;
panel Q above, not this row, is the evidence for the silent-number leg.) Adding
esprelength(2) espostlength(2) jointtests to the c = -8 run turns rc 0 into r(459)
"could not calculate numerical derivatives — discontinuous region with missing values
encountered" (r_f3d_scope.log).

NOT ESTABLISHED, by them or by me: the contreatelasticitytype(eyex|eydx) half of the
family. My contreat designs died on wooldid's multicollinearity screen r(322).

## New (unseeded) probe and its honest negative result
Repro: r_new_validgate.do

The generalizable class behind F1 is "value-domain validation gated behind a parent".
Sweeping every other valued option:
 - esrelativeto(3) — a POSITIVE value, where the help documents (-int) and the gated
   message says "must be a negative integer" — is rc 0 without esfixedbaseperiod and
   r(198) with it. BUT wooldid prints "Caution: esrelativeto() has no effect without
   esfixedbaseperiod and is being ignored", so the inertness IS announced. Low.
 - pretrendpower() out-of-domain: r(198) with and without a window. Screened.
 - customweightmultiplier() with negatives: r(411) (matches ruled VAL-G06).
RESULT: cfxplottypes() is the only valued option I could find that is BOTH unvalidated and
unannounced when its consumer is not reached. That is the honest scope of F1's live part,
and it also shows the author does announce ignored options elsewhere.
