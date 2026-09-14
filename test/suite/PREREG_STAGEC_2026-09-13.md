# Preregistered signature — Stage C Monte Carlo campaigns, 2026-09-13

Written BEFORE any of the three campaigns is run. The three drivers (`mc_adv_bs12.do`,
`mc_bands.do`, `mc_ri_ext.do`) shipped in `4280f3b` as scaffolds, each carrying an ACCEPTANCE
block that says of itself "to be preregistered before the production run". This is that
preregistration. Analyzers: `analyze_advbs12.py`, `analyze_bands.py`, `analyze_riext.py`.

Old/new is not the axis here: these measure the CURRENT build against statistical truth, because
the features they cover do not exist in v1.3 and have no external oracle. That is the whole
reason the chunk asks for correct sizing rather than an oracle comparison.

## C0 — what is preregistered, and what is honestly not

Every band below is one of three kinds, and each is labelled in the analyzer source at its
constant:
- **DRIVER ACCEPTANCE** — lifted verbatim from the driver's own block, written when the driver
  was written and before any campaign.
- **PREREGISTERED ELSEWHERE** — exactly one: the assignment-space floor of ~3000, from
  `PREREG_NEWFEAT_2026-08-24.md` line 39 ("assignment spaces chosen >= ~3000 so tie prob <= ~3%").
- **NOT PREREGISTERED** — chosen 2026-09-13 while writing the analyzers, stated here so it cannot
  later be presented as a prediction. These are: the minimum-n gate of 30 per cell; the MARGINAL
  corridors [.90,.98] for a .95 coverage claim and [.85,.94] for a .90 one (the house size
  corridors do not transfer to a coverage claim); the arithmetic self-consistency checks; and the
  requirement that every p lie in [0,1].

## C1 — the BS-12 adversarial size ladder (`mc_adv_bs12.do`)

**Why it exists.** Every adversarial size number on record — `MC_REPORT_2026-08-24.md`,
`REPORT_ADVBOUND_2026-08-26.md`, and the help at sthlp:745-752 — was measured on the
free-intercept statistic that BS-12 replaced on 2026-09-10. The ruling's own prerequisite
simulation was never run. This is that simulation.

**P0 (falsifiable).** On every row: `rc == 0`, `advabsent == 0`, `p` and `pjoint` in [0,1],
`ncomp == B` and `nfail == 0` on the `pm` and `bt` arms, `advw_big == 0` on `PM`.

**P1 — the headline.** Rejection COUNT of `p` at nominal .05 and .10, per (config, arm), inside
the central ~95% region of Binomial(M, attained-alpha). **Attained**, not nominal: the p lattice
is discrete, so the analyzer enumerates the run's own lattice and bands against P(p <= alpha) on
it. At B = 199 an add-one p attains .05 exactly; a plain k/B p attains 9/199 = .0452 and the band
moves by a rank.

**P2 — a correction to how this driver is read, recorded in advance.** The plant is
`SLOPE*(relyr+1)`, which is ZERO at the reference period `relyr = -1`, so the true profile is a
line THROUGH the reference period — inside the BS-12 null for EVERY slope. The sloped arms are
therefore a SECOND NULL, not a power arm, and their rejection rates are size numbers to be banded
exactly like the flat arm's. The driver's header says this at its lines 31-34; a brief written
2026-09-13 called them power and was wrong. **Consequence, stated plainly: this campaign measures
size across the composite null and contains NO power measurement.** Size is what the chunk asks
for ("whether the adversarial line test p-values are correctly sized"), so this is not a gap
against the ask; it IS a gap against anyone who later wants power, and adding one needs a
deviation that is not a line through the reference period.

**P3 AMENDMENT, 2026-09-14 — a precondition neither this preregistration nor the driver
recorded, and what it does to P3's claim.** `wooldid.sthlp:745` ("Size of the adversarial test")
states a condition specific to the continuous case: *"Under a continuous treatment the line test
additionally requires cells containing many observations - on the order of fifty - because each
cell's dose slope is otherwise too imprecisely estimated for the test's covariance to be
trustworthy."* AB1C's geometry is G = 30 cohorts x U = 6 units x T = 12 periods, so a treated
cohort-by-period cell holds **six** observations against the help's ~50. **AB1C sits about an
order of magnitude below the documented precondition.**

This amendment changes no band and no prediction. What it changes is what the rung can be read as
showing. As designed, AB1C **cannot certify that the continuous line test is correctly sized**; it
can only measure the test on a design the help already declares untrustworthy, which is a
confirmation of a documented caveat rather than a test of the feature. Recorded here because P3 as
written ("the DRIVER now covers this") invites exactly the wrong reading.

Note which half of the caveat AB1C *does* meet: the cluster condition. `cluster(unitid)` puts six
clusters in every cell and the help's threshold is "six or more" - which is why the `att` line
test on the very same fits is expected to be, and is measured to be, correctly sized. The
asymmetry between the two families on one fit is the help's own prediction.

**Consequence for the chunk's question.** "Whether the adversarial line test p-values are
correctly sized ... including with continuous treatments" needs a COMPLIANT geometry. The measured
arm is therefore a dose-response in observations per cell at a fixed spec (`probe_ame_u.do`, U as
the only varying quantity, an arm only), with U = 50 as the help's own number. The prediction
stated in advance: rejection at nominal .05 falls toward nominal as U rises, and at U = 50 the ame
line test is inside the Binomial(n, .05) central region while the att twin stays inside it at
every U.

**P3 — the known gap, recorded so the campaign is not later read as covering it.** The
geometries are all BINARY. The chunk asks for sizing "including with continuous treatments". A
contreat geometry, `AB1C`, was ADDED to the driver on 2026-09-14 and its null was verified before
anything else, so the DRIVER now covers this. Every CAMPAIGN run to date still does not, and no
campaign output may be cited for the continuous case until one includes AB1C. `AB1C` is an-arm
only: `cluster(unitid)` does not nest `i`, so the resampling engines refuse it and there is still
no continuous rung at the permutation or studentized arms.

## C2 — simultaneous sup-t band coverage (`mc_bands.do`)

**Why it exists.** The band shipped 2026-09-10 after every MC had already run. `test_esband.do`
pins its structure only, never its coverage, and no MC has used a cluster strictly BELOW `i` —
the geometry the help recommends and the only one on which the analytic `e(V_es_*)` the band is
simulated from is the right object.

**P0.** `rc == 0`, `band_absent == 0`, `nh == 8`, `pw_n == nh`, `ncomp == B` and `nfail == 0` on
the bootstrap arms, and the summary columns arithmetically consistent with the per-horizon ones
(`band_hits == sum(bd*)`, `pw_hits == sum(pw*)`, `band_all == 1[band_hits == nh]`). That last
family exists because a driver that mis-writes its summary columns turns a coverage number into
fiction without failing anything.

**P1 — the headline.** Simultaneous coverage per arm inside Binomial(n, .95), nominal .95 because
`mc_bands.do` passes no `level()` and wooldid defaults to 95.

**P2.** Per-horizon pointwise coverage inside Binomial(n, .95), as the control.

**P3 — directional, and the sign IS the claim.** The `pw_all` share must sit BELOW the `band_all`
share. A simultaneous band that is not wider than the pointwise intervals is not doing its job,
and this is the one check whose direction, not magnitude, carries the meaning.

## C3 — randomization-inference extensions (`mc_ri_ext.do`)

**Why it exists.** Continuous treatment, `permside(right)`, `level(90)`, `permstrata()` and a
`permtstat` null arm have no MC coverage at HEAD.

**P0.** Per phase: `rc == 0`, `nfail == 0`, `ncomp == B`, `ndraw == B`, `namesok == 1`,
`pconv` equal to the phase's convention (`right` on E2, else equal-tailed), `level` 95 (90 on E3),
`tau == 0` on the size phases, and every p in [0,1].

**P1 — assignment space.** `npos >= 3000` per phase. This is the one number preregistered
elsewhere and it is load-bearing: a small space puts a floor under the achievable p, so a size
result from a small space is not evidence about the statistic.

**P2 — the headline.** Size per (phase, p-column) at attained .05 and .10 inside
Binomial(n, alpha)'s central region, covering `p_main`, `p_pre`, `p_post0`, `p_post1`, plus on E1
`p_ame`, `p_amepost0`, `p_amepost1`, plus `jointp` and `jointCp`. The joint p's are judged on the
ONE-SIDED lattice whatever the arm's `permside`, because `e(joint_horiz)` is the add-one one-sided p.

**P3 — lattice.** `p_main` must be an exact multiple of the lattice step. Off-lattice on `p_main`
is a FAILURE: that column lost draws, so the headline p is not exact. Off-lattice on any other
column is a printed finding and is NOT scored, because a horizon can legitimately be unidentified
under some reassignments.

**P4 — E3 interval coverage.** Coverage of tau by `[ll_main, ul_main]` inside Binomial(n, .90).
The blank-interval share is reported and is OUTSIDE the coverage denominator, so that rate is
conditional and must be read as such.

## C4 — what a failure means, decided in advance

A P0 violation is a broken run: stop, diagnose, do not report a rate from it. A P1/P2 miss is a
FINDING about the statistic and goes to the author with its n and its interval, never silently
re-banded. A cell below the minimum-n gate reports INSUFFICIENT_N and is never scored as a pass —
this matters because the failure mode these campaigns exist to prevent is a family that wrote
missing values and got read as a zero result.

## C5 — the one column-level assumption no data in hand can settle

The `bt` arm's p-value convention is assumed to be add-one `k/(B+1)` by analogy with `pm`, whose
form is observable. No `bt` row exists in any sample collected so far. `analyze_advbs12.py`
therefore MEASURES the observed lattice and prints a finding naming what to change if the
assumption is wrong, rather than trusting it. Recorded here so that if the convention turns out
to differ, the band change is visibly a correction and not a post-hoc adjustment.
