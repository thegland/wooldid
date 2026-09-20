# Three screens held for ruling — RULED 2026-09-20, see each verdict below

Everything else from the blue-sky review is built (19 ruled items in `1c09efa`; the
non-number-moving batch in the commits alongside this file). These three are the residue.

They are held for the same reason, and it is not that they move numbers — none of them
moves a reported number on any run that keeps running. It is that **each one adds a screen
whose firing is a new event in a user's log**, and the evidence behind each is the
proposer's fixtures, not a measured rate on ordinary data. A screen that cries wolf is how
L2 died (see below). So the question for each is one line: **caution / r(459) / gated /
don't ship.**

---

## 1. LAT L3 — the lattice support screen

**What it does.** On `contreat()` + `lattice()` runs only, one `gegen` pass over the treated
rows counts how many (cohort-period x lattice-bin) pieces hold no more observations than the
`1 + contreatpoly` parameters that piece is given. Such a piece is fitted *exactly* — its
residuals are identically zero — so no residual-based variance estimator (robust, cluster, or
the bootstrap's analytic leg) can see any variation there.

**Why it matters.** This is the mechanism behind the worst number the review found. On the
broken fixture the reported ame standard error is **0.017 of the estimator's own conditional
sampling sd** (frozen-design Monte Carlo, `proposals_LAT.md` s0.3 — only the noise is
redrawn, so the MC sd *is* the conditional sampling sd; the ratio is calibrated against a
control arm). Today the run prints that standard error and says nothing.

**Discrimination, as measured:**

| design | pieces | saturated | treated obs in them | screen |
|---|---|---|---|---|
| broken (5-12 units per state) | 273 | 101 | 157 | **fires** |
| healthy (40 units per state) | 234 | **0** | 0 | **silent** |
| broken, no `lattice()` | — | — | — | silent (by construction) |
| broken, `safety(off)` | — | — | — | silent (gate works) |

**False positives.** Structurally none *in the claim it makes*: a saturated piece really does
have zero residuals — that is arithmetic, not a heuristic, and there is no threshold to argue
about. The open risk is **nuisance**: a lattice run with one thin cell out of 300 would print
the caution. I have no measurement of how often ordinary lattice runs have one thin cell.

**Cost.** Upper bound 5.05 s at 5,000,000 rows / 60,000 groups; well under a second on
ordinary panels. Gated to `contreat()` + `lattice()`, skipped under `safety(off)`, and kept
out of every bootstrap/permute replication, so it runs once per run or not at all.

**The proposer recommends:** caution (the code above). Refusing was considered and rejected
in the proposal — the point estimates are unbiased, it is only the standard errors that are
worthless.

**Coupling:** the help now describes this support requirement (new paragraph in the
`lattice()` discussion). I wrote it **without** the clause "wooldid cautions when it sees
this", because that clause is L3. If L3 ships, that sentence can be added back.

> **RULED 2026-09-20: DON'T SHIP.** Author: *"I think we can skip Lat L3. Let's just put a
> comment in the help that when using a lattice, you need to be careful not to over saturate
> the specification."* The screen is not built. The help carries the warning instead: the
> `lattice()` discussion now opens its support paragraph with "take care not to over-saturate
> the specification with lattice terms", and keeps the mechanism, the consequence for the
> standard errors, and the units-per-cell rule of thumb behind it.

---

## 2. FIX-PERM-1 — the inert-covariate screen under `permonsetcohorts`

**What it does.** Under `permute(..., permonsetcohorts)` with a `cluster()` finer than `i`,
the replications rebuild `i` from the onsets they draw. A covariate that is constant within
the *original* `i` is collinear and dropped in the observed fit, but is **not** collinear
inside the rebuilt `i`, so it is live in the replications. Observed and reference are then
not the same function, and the permutation p-value is not a test of anything. The screen
detects the slot and refuses.

**Measured, nine arms** (`p3_patched.log`): the four arms that must be untouched are
**bit-identical**, the two that must be caught are caught, and `safety(off)` reproduces
today's number exactly. `permute()` without `permonsetcohorts` is untouched; so is
`cluster() == i`.

**False positives.** The proposal closed the one route it could construct
(`p7_window.log`). The screened predicate is "constant within `i`", which is the exact
condition for the collinearity, so a false positive would have to be a covariate that is
constant within `i` and yet not dropped.

**Cost.** 0.156 s once, on 200,000 rows / 2,000 clusters, and **exactly zero** on every run
that does not use `permonsetcohorts`.

**Installed base.** `permonsetcohorts` + finer `cluster()` is from T1 (2026-09-16) and its
assessment suite was deferred, so there is essentially nothing to break.

**The two tiers the proposal itself offers:**
* **r(459) under default safety, caution under `safety(off)`** — the proposer's
  recommendation, on the argument that the option provably cannot change the estimate, so
  refusing costs the user only a deleted token. Matches the `bootstrap(bayes)` and
  `contreatwithin` patterns.
* **caution only** — delete the safety wrapper. Then nothing is refused, no number moves at
  all, and the user is merely told. The proposal calls this the zero-risk choice.

> **RULED 2026-09-20: DON'T SHIP, no change at all.** Author: *"let's not ship any changes,
> we chatted about this earlier, it's too much of an unusual edge case and I think users
> should be able to tolerate it."* Nothing is built and no help text is added.

---

## 3. GTE-3 — the all-zero-VCV screen tests for exact zero too

**Held because the proposal says so**, and it is right: it changes the *firing predicate* of
a message you ruled on (AMESE R4, shipped 2026-08-24).

**What it does.** `wooldid_amescreen, mode(vzero)` currently reads `max(abs(e(V))) == 0`. A
cluster-robust sandwich that is rank-deficient in *every* direction returns **cancellation
residue**, not zeros. Measured: `max|e(V)| = 6.0e-28` on a two-cluster design where every
reported standard error was ~1e-15 and this screen stayed silent; the run printed
`se = 1.60299e-15, p = 4.41339e-16`. GTE-3 replaces the exact-zero test with a floor relative
to the fit's own largest coefficient (so it is equivariant in the units of y).

**What it costs you.** The two caution texts the screen controls become false as written
("every standard error ... is missing" — they are 1e-15, not missing), so both are amended in
the same patch. That is why it needs you: it widens a ruled message and rewrites it.

**Note on ordering.** LAT L1 (shipped in this batch) rewrote the *branch structure and remedy
list* of those same two messages — it changed what they say, not when they fire, which is why
it was separable from this hold. If you take GTE-3, its two substring amendments must be
re-anchored onto L1's new branch texts; the mechanism is untouched by that.

> **RULED 2026-09-20: SHIP, with a conservative band.** Author: *"I guess use a conservative
> band here."* Built separately from this batch, with the floor placed on measured evidence
> rather than chosen by feel - see `GTE3_BAND_2026-09-20.md`.

---

## What died on measurement, for contrast — why these are being asked rather than shipped

* **LAT L2** — a free predicate on the dropped-column names. Built, run, measured: **it cried
  wolf**, firing on the healthy fixture too. Rejected on its own evidence.
* **LAT L2c** — a second free predicate (excess drops over the structural baseline). Also
  built, also measured, also rejected.
* **A relative-sd threshold for L3** (copying the `bycoh` screen) — catches 46 of 273 pieces,
  essentially the singletons only, and **misses all 56 two-observation pieces**, which hold
  112 of the 157 zero-residual observations. It asks "is this piece well conditioned?" when
  the damaging question is "is this piece saturated?".
* **FIX-PERM-2b** — the numeric clusters-per-cell predicate. Not asked: your standing ruling
  was to leave that statement qualitative, and the qualitative form is already shipped in both
  the help and the runtime caution.
