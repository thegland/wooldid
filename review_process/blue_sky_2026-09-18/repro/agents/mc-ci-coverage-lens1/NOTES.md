# LENS 1 (DOCUMENTED-BEHAVIOUR) — mcci-03-esprelength-moves-att

**Verdict: PARTIAL, severity low.** The behaviour is real and reproduces on a third
dataset, but it is documented, was found by the project itself, was ruled, and was built.

## The documentary record the prober and the skeptic both missed

The finding says the behaviour "IS documented -- but only in the Stored Results entry for
e(escomp_refN)". That is **false on the pinned build**. It is also documented, far more
explicitly, in the main help prose:

`build/wooldid.sthlp:1493-1497`, "Comparing Output to did_imputation":

> ... requests no pre-treatment horizons (`esprelength(0)`) ... **The `esprelength(0)`
> condition is easy to miss and it matters: the pooled reference is the set of horizons
> the run asks for, so requesting pre-periods removes them from that pool and shifts every
> post horizon by a constant (measured 4.8e-02 on a common-onset panel) while both commands
> still run cleanly.** `test/suite/newfeat/bjs_identity.do` pins both directions - the match
> at `esprelength(0)` and the divergence without it.

Provenance — this is the project's own prior finding, not a new one:
- `wooldid_old_bug_report_2026-09-13.md:4282` = **C.2.2** "The did_imputation equivalence
  conditions omit esprelength(0)", with the same mechanism, measured against an independent
  oracle: wooldid `esprelength(0)` vs BJS max|diff| 4.07e-08 (EQUIVALENCE_HOLDS);
  `esprelength(4)` vs BJS 4.76e-02 (EQUIVALENCE_BROKEN); and BJS's own `pretrends(4)` moves
  its post taus by **exactly 0.0** (POST_TAUS_UNCHANGED).
- `review_process/held_for_ruling_2026-09-13.md:63,704-745` = **item 13 / T-BJS1**, ruled
  "(iii) documentation wording ... needs no ruling", fix at `sthlp:1464`.
- `known_bydesign.md` §7 names items 1, 8 and 13 as the three of thirteen that **were**
  verified in code. Item 13 is one of them, and the shipped passage above exceeds the
  ruled minimum (the ruling asked only for `esprelength(0)` to be added to the list).

It is also standard econometrics: putting leads in any TWFE/imputation event-study
regression takes those cells out of the untreated pool that identifies the counterfactual,
so the omitted category moves and precision falls. PS-34 already rules the pooled reference
"default and unlabelled", which is the other half of "nothing prints at run time".

## What survives

Exactly one crumb, and it is a cross-reference, not a defect: item 13 fixed the *BJS
comparison* passage; nobody was ever asked whether the **`esprelength()` option entry**
(`sthlp:729-732`, whose `{pmore}` covers only the contiguous-slot issue BS-36) should point
at it. A user who does not use `did_imputation` has no reason to read that section.

## What falls

- "documented only in the Stored Results entry" — false (above).
- "up to 60%" — dataset-specific. Three datasets, three numbers: prober +60%, skeptic +21%,
  mine +9.5% (common onset) and +26.6% (staggered).
- "esprelength()/espostlength()" — `espostlength()` moves nothing. Confirmed **bit-exact**
  in hex on both my panels, and the help's own sentence says "pre-periods".

## My measurement (lens_espre.do / lens_espre.log), 80 units, T=12, cluster(unit), i=cg

e(N)=960 and e(df_r)=79 in every run of both panels.

| panel | spec | b | se | e(escomp_refN) |
|---|---|---|---|---|
| hom | none | 1.265880342874 | 0.148524822029 | . |
| hom | espostlength(3) | 1.265880342874 | 0.148524822029 | 360 |
| hom | esprelength(3) espostlength(-999) | 1.351063444131 | 0.162645597762 | 180 |
| hom | esprelength(3) espostlength(3) | 1.351063444131 | 0.162645597762 | 180 |
| stag | none | 1.313287406473 | 0.116928762989 | . |
| stag | espostlength(3) | 1.313287406473 | 0.116928762989 | 360 |
| stag | esprelength(3) espostlength(-999) | 1.408820394599 | 0.147998295742 | 180 |
| stag | esprelength(3) espostlength(3) | 1.408820394599 | 0.147998295742 | 180 |

`none` and `espostlength(3)` are identical in hex on both panels
(stag: b `+1.50339a7d0cd74X+000`, se `+1.def0b1cffe61dX-004`).

## ADJACENT OBSERVATION — one shipped sentence is wrong on staggered panels

Not part of mcci-03's severity; report separately. The help says the re-basing "shifts every
post horizon **by a constant**". Per-horizon shift, `esprelength(0)` -> `esprelength(3)`,
same data, same `espostlength(3)`:

| horizon | hom shift | stag shift |
|---|---|---|
| att | 0.085183101257 | 0.095532988126 |
| post0 | 0.085183101257 | 0.084752485441 |
| post1 | 0.085183101257 | 0.122836625585 |
| post2 | 0.085183101257 | 0.061364233795 |

Constant to 12 decimals under common onset; **varies by a factor of 2 under staggered
adoption** — which is the design wooldid exists for. The project's own file already recorded
this (`held_for_ruling_2026-09-13.md:735`: costag 1.17e-2 "**not** constant", dta_binary
post0 3.71e-2 "not constant"), but the sentence that shipped states the constancy
unqualified and attaches the hedge only to the parenthetical measurement. Suggested wording:
"shifts every post horizon — by a constant under common onset (measured 4.8e-02), by a
horizon-varying amount under staggered adoption".

## For the coordinator

`known_bydesign.md` §3 lists BS-7 and BS-36 on the ES options but **omits T-BJS1**. That
omission is why two agents spent a probe and a refutation cycle on ruled ground. Worth adding.
