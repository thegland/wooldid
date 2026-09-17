# Runtime cost of guard checks — where a wooldid run spends time proving things are fine

Investigation opened 2026-09-16 at the author's request, as a parallel track to the
fineclus wrap-up. **Documentation only — no code changed.** The author's instruction was
explicit: "don't make the edits right away, just document as we go."

## The question, and the filter it comes with

How much of a typical run goes on testing and verifying guard conditions, and are there
checks worth dropping as redundant with later fail-safes? The author set the bar himself,
and it is a conjunction of three things. A candidate must:

1. cost **non-trivial additional computation** (not a macro comparison);
2. test for something **rare** — under roughly 1 in 100 practical uses; and
3. be **redundant with a later fail-safe** — the same failure gets caught downstream anyway,
   just later and with a less expressive message.

With an explicit carve-out: *"core estimability checks, for example, even if expensive are
too valuable to even think about skipping."* And an explicit exchange rate: don't spend one
unit of time to avoid fifty units in a failure mode that applies in fewer than 1 in 100 runs.

Criterion 3 is the one that does the work. Most of what looks like a guard in this codebase
fails it, because the failure it prevents is **silent** — a wrong number, not a crash.

## Method

Static first, then measured. A parser walks `wooldid.ado` tracking brace depth and program
boundaries, and tags every expensive operation (`gegen`, `gstats`, `gduplicates`, `gunique`,
`glevelsof`, `gcollapse`, `gisid`, `greshape`) with its **innermost enclosing condition**
rather than its enclosing program. That distinction turned out to matter more than anything
else in the investigation: a program-level count says `wooldid_permctl` is the hottest guard
site in the package, and it is, but only on runs that use `permute()`.

Scripts and their regeneration instructions: `review_process/research_2026-09-16/`
(`guardmap.py` for the static map, `count.do`, `fire.do`, `calib.do`, `si8.do` for the probes).

## Landscape

Expensive-op sites by innermost enclosing condition (top of the distribution):

| n | program | innermost condition | on a typical run? |
|---|---|---|---|
| 9 | `wooldid_worker` | `if espre+espl > 0 & wdrepmode != "1"` | yes, if an event study is requested — **already rep-gated** |
| 8 | `wooldid_worker` | `if subgroups != "" & suppresssgprimaryeffects == ""` | only under `subgroups()` |
| 7 | `wooldid_permctl` | top of program | only under `permute()` |
| 6 | `wooldid_permctl` | `if onc != ""` | only under `permonsetcohorts` |
| 5 | `wooldid_worker` | `if safetyinput != "off"` | **yes — this is the candidate** |
| 5 | `wooldid_reltime` | top of program | yes, but load-bearing (builds relative time) |
| 4 | `wooldid_worker` | top of program | yes, but load-bearing |
| 4 | `wooldid_permctl` | `if sg != ""` | only under `permute()` + `subgroups()` |
| 2 | `wooldid_fcol` | top of program | yes, but load-bearing (the collapse itself) |
| 2 | `wooldid_bootctl` | `if clus != ""` | only under `bootstrap()` |

The 26 sites in `wooldid_permctl` and `wooldid_bootctl` are **control-program setup**, run
once per run ahead of the replication loop (verified: none sits inside a `forv` over
replications). Their cost is O(1) against the engine's O(R), so they are not candidates
however expensive they look in a raw count. The same is true of `wooldid_permrun`'s sites.

## Census of the `safety()`-gated blocks

The `safety()` family is where a reader would expect the skippable checks to live, so it was
enumerated exhaustively rather than sampled. There are **24** `safety`-gated blocks. Excluding
display lines and matching only outside quoted strings (message text mentioning `margins`
otherwise scores as a `margins` call), **7** contain a genuinely executed expensive operation:

| site | program | what it is | verdict |
|---|---|---|---|
| L762 | `wooldid_worker` | F13 onset-constancy, 2 `gegen by(i)` | **candidate — rep-gate** |
| L898 | `wooldid_worker` | pre/post-period sufficiency, 1 `gstats` | core estimability — author's carve-out |
| L1341 | `wooldid_worker` | si8 variable-count, 3 `gegen tag` | **candidate — rep-gate** |
| L2399 | `wooldid_worker` | `makecfxplots` production, 4 `margins` | not a guard; output, and only under `makecfxplots` |
| L4959 | `wooldid_validate` | float-precision screen, 1 `gstats` | cheap, once |
| L5326 | `wooldid_validate` | time-range screen, 2 `gstats` | cheap, once |
| L12234 | `wooldid_sgdom` | subgroup-constancy, 2 `gegen by(i t)` | **candidate — rep-gate**, under `subgroups()` |

**The other 17 are macro or scalar comparisons and cost nothing measurable.** That is the
short answer to the survey half of the author's question: the `safety()` surface is mostly
free, and its cost is concentrated in three blocks.

## The one real finding: the whole-sample safety blocks re-run on every replication

Two blocks in `wooldid_worker` are gated on `safetyinput != "off"` and on nothing else:

- **F13, onset-constancy** (`ado:762-779`) — two `gegen ... by(i)` passes over the estimation
  sample plus an `assert`, checking that `ttre` is constant within `i`.
- **si8, variable-count** (`ado:1341-1365`) — three `gegen ... tag()` passes plus three
  `count`s, producing a caution when `i` has over 100 levels and pre-empting an `r(907)`
  when the implied variable count would exceed `c(maxvar)`.

A third block has the identical shape one level down: **`wooldid_sgdom:12234`**,
subgroup-constancy within `i`-by-`t`, two `gegen ... by(i t)` passes and an `assert`. It runs
only under `subgroups()`, but `wooldid_worker:985` calls `wooldid_sgdom` on every worker entry,
so under `subgroups()` it too is paid once per replication. It carries the same safety argument:
under bootstrap re-id the rebuilt `i` is *finer* than the original, so constancy within
`(i, t)` is implied by the constancy rep 0 already verified; and `permonsetcohorts` refuses
`subgroups()` outright, so the rebuild case cannot arise.

`wooldid_worker` is re-entered once per replication: `wooldid_bootrun` and `wooldid_permrun`
set `WOOLDID_REPMODE 1` around each nested `wooldid` call. Measured on a 240-unit × 16-period
panel with a cohort-level `i`:

| run | `wooldid_worker` entries | F13 block | si8 block |
|---|---|---|---|
| analytic `att` | 1 | 1 | 1 |
| `bootstrap(standard) reps(20)` | 22 | **22** | **21** |

So five `gegen` passes over the estimation sample are paid on every replication. At
`reps(500)` that is ~2,500 whole-sample grouped passes.

**This is not a new principle — it is an existing one with two sites missed.** The worker's
own header comment says so (`ado:533-537`):

> *"Replication results are read for their coefficients alone, so blocks that only fill
> diagnostic e() returns can be skipped there."*

And the rule is already applied: the nine event-study sites at `ado:3816` carry
`& "\`wdrepmode'" != "1"` explicitly. The two safety blocks simply never got the same gate.

### Why these two are safe to skip in a replication specifically

Not "the check is worthless" — it is valuable, at rep 0, where it runs on the observed data
and produces its full expressive message. The claim is narrower: a replication cannot
produce a violation that rep 0 did not already catch.

- **F13**: under bootstrap re-id, `i` is rebuilt as `group(__wdbns, i)`, so onset stays
  constant within the new `i` by construction. Under `permute()` an onset is drawn per
  cluster, and where a cluster is finer than `i` (`permonsetcohorts`) the rebuild *makes*
  `i` a function of the drawn onset — which is exactly what the dedicated rep-0 assertion
  already checks.
- **si8**: resampling cannot create more distinct cohorts than the data has, so the implied
  variable count cannot grow past what rep 0 measured. And the failure it pre-empts is
  `r(907)`, which Stata raises loudly on its own.

Measured, not assumed. A build instrumented at each block's **failure branch** records every
firing together with the `wdrepmode` it fired under:

| run (40 reps each) | F13 entries | si8 entries | F13 fired | si8 fired |
|---|---|---|---|---|
| `bootstrap(standard)` | 42 | 41 | never | never |
| `permute(shuffle)` | 42 | 42 | never | never |
| `permute(shuffle) permonsetcohorts` | 42 | 42 | never | never |

**And the probe was calibrated in the other direction before that table was believed**, which
is the discipline this slate learned the hard way. Breaking one unit's onset inside a cohort
so `ttre` is genuinely non-constant within `i`:

| calibration | F13 entries | fired | rc |
|---|---|---|---|
| broken data, safety on | 1 | **yes**, `rep=` (observed run) | 459 |
| broken data, `safety(off)` | 0 | no — block never entered | 0 |

So the instrument can fire, and fires exactly where it should; the zero firings across ~126
replication-time executions above are a real reading rather than a dead probe.

The honest limit of that evidence: "never fired on healthy data" is weak on its own. The
load-bearing argument is the structural one above, plus the fact that the author already
applies precisely this gate at `ado:3816` for the same reason.

## Screened out, with reasons

| site | criterion it fails |
|---|---|
| L899 pre/post-period sufficiency | **core estimability check — the author's explicit carve-out.** Also cheap: one `gstats sum`. |
| F13 / si8 on the *observed* run | fails (3): F13's failure is **silent**, not loud. Its own comment says rows of the same `i` with different onsets "are silently re-timed row by row ... with no error." A check against a silent wrong answer is never redundant with a later fail-safe. |
| `permctl` / `bootctl` setup (26 sites) | fails (1) in relative terms: once per run against R replications. |
| `wooldid_reltime`, `wooldid_fcol` top-of-program sites | not guards — load-bearing computation. |
| `wooldid_amescreen` (10 sites, under `contreat()`) | **core estimability — the author's carve-out.** The code says so itself at the `mode(bycoh)` call site: "byCohort point-identification screen ... the pre-treatment-period sufficiency errors' sibling, one layer down: rows can exist while the baseline dose slope's identifying VARIATION does not." A point-identification screen is the paradigm case of a check too valuable to skip, and its failure is silent. |
| `wooldid_validate` double-parse | not a guard for a rare thing. It is architecture: the worker parses `0`, then `wooldid_validate` re-parses it with a byte-identical syntax line, and that happens on every worker entry and so on every replication. Noted as the largest *structural* per-replication overhead, but outside the author's stated criteria — recorded, not recommended. |

## Incidental finding — now a ticket, not a footnote

The si8 maxvar pre-emption **did not fire on a case that then failed with Stata's own
`r(907)`**. On a plain 200-unit × 20-period panel with `i = unit`, the block measured
`cells=2004 ilev=200 tlev=20 demand=2424` against `c(maxvar) = 5000`, correctly stayed silent,
and the run then died on Stata's generic "maxvar too small" rather than wooldid's expressive
message. The true requirement exceeded the estimate by more than 2×. Re-verified on a build
carrying one added `dis` and no logic change, so it is not an instrumentation artifact.
Partly mitigated: `si8ilev` was 200, so the `> 100` caution fired first and did give the
coarsen-`i` advice.

This is a correctness defect and has been **opened in `open_defects.md` under OPEN**, where a
reader looking for defects will find it. It is unrelated to anything this report recommends,
and it does not weaken the rep-gate case: under the rep-gate the block still runs at rep 0
exactly as it does today.

## Rarity estimates are the weak column

Criterion 2 asks how often each condition actually bites, and the honest answer is that
**every rarity figure here would be a guess**, so none is offered as a measurement. The
recommendation below does not depend on any rarity estimate — it does not drop a check on
the grounds that its condition is rare. It removes a repetition of a check whose answer
cannot have changed. That is a strictly weaker claim and a strictly safer change.

## The one candidate

One change, offered for the author's ruling and **not made** here; the measurement section below
prices it, and on the direct measurement it **is** worth making:


> Add `& "\`wdrepmode'" != "1"` to the `safetyinput != "off"` gates at `ado:762` and
> `ado:1341`, matching `ado:3816` and the rule stated at `ado:533-537`. The same gate applies
> to `wooldid_sgdom:12234` under `subgroups()`, with the caveat that `wooldid_sgdom` does not
> currently receive `wdrepmode` and would need it passed in — a slightly larger change, and
> worth treating as a separate decision from the two worker sites.

What it buys and what it costs is measured in the next section. What it does **not** cost:
any diagnostic on the observed run. Both checks still run at rep 0, on the user's real data,
with their full messages.

## Measured cost — directly, after a differencing attempt failed

**The first attempt used the wrong instrument and I nearly shipped the wrong conclusion from it.**
Paired alternating A/B/A/B (A pristine, B rep-gated), 10,000 observations, `reps(50)`:

| pass | analytic safety on | analytic `safety(off)` | boot A | boot B |
|---|---|---|---|---|
| 1 | 0.536 | 0.600 | 32.303 | 31.619 |
| 2 | 0.587 | 1.223 | 29.200 | 34.179 |
| mean | 0.625 | 0.775 | 30.75 | 32.90 |

Both comparisons came out with the wrong sign. On that evidence I downgraded the recommendation
to "not worth making". **That was wrong.** The same bootstrap configuration varies 29.2 to 34.2
seconds run to run - about 17% - so differencing cannot resolve an effect several times smaller
than its own noise, whichever sign it has. A null from an instrument that could not have detected
the effect is not a null.

**Replacing it with a direct instrument** - `timer on/off` wrapped around each block inside a
scratch build, so the guard time is read rather than inferred - settles it. Two independent runs:

| run | total (s) | F13 (s) | si8 (s) | guards | run 2 total / F13 / si8 |
|---|---|---|---|---|---|
| analytic `att` | 0.417 | 0.010 | 0.012 | **5.28%** | 0.438 / 0.010 / 0.012 |
| `bootstrap(standard) reps(50)` | 21.893 | 0.507 | 0.647 | **5.27%** | 21.868 / 0.490 / 0.649 |
| `bootstrap(standard) reps(150)` | 64.249 | 1.504 | 1.920 | **5.33%** | 64.079 / 1.464 / 1.907 |

The share is **5.3%, and it is the same 5.3% across totals spanning a factor of 150** - which is
what it should be, since both the guards and the engine scale with worker entries. That internal
consistency, reproduced across two runs, is far stronger evidence than the differencing ever was.

### What the rep-gate actually buys

| run | guard time | removable by the rep-gate | saving |
|---|---|---|---|
| analytic | 0.022s | 0.000s (rep 0 must run) | **none** |
| `reps(50)` | 1.154s | 1.131s (50 of 51 entries) | **5.17%** |
| `reps(150)` | 3.424s | 3.401s (150 of 151) | **5.29%** |

So the change buys about **5% of every `bootstrap()` and `permute()` run**, costs nothing on the
analytic path, and removes no diagnostic anywhere: rep 0 still runs both checks on the user's real
data with their full messages.

## Recommendation

**Make the change** (still not made - this is a document for the author's ruling):

> Add `& "\`wdrepmode'" != "1"` to the `safetyinput != "off"` gates at `ado:762` and `ado:1341`,
> matching `ado:3816` and the rule the worker's own header states at `ado:533-537`.

`wooldid_sgdom:12234` is the same shape and the same argument, but its syntax takes `si()` and no
repmode, so gating it means threading a new option - a separate, slightly larger decision.

## Answering the author's question as asked

*"How much time in a typical wooldid run is spent testing and verifying guard conditions?"*

**About 5% for the two whole-sample safety blocks**, plus a few `gstats` calls elsewhere that are
smaller still. 17 of the 24 `safety()`-gated blocks are macro comparisons and cost nothing
measurable. So guards are a real but modest slice, not a hidden tax.

*"Are there checks that can be skipped as redundant with later fail-safes?"*

**Essentially none, and the reason is worth keeping.** The author's third criterion - the failure
would be caught later anyway - is what disqualifies nearly every candidate, because the failures
these guards prevent are **silent**. F13's own comment says mis-specified onsets are "silently
re-timed row by row ... with no error". A check standing between the user and a wrong number is
never redundant with a later fail-safe, however cheap the later failure would be. The core
estimability screens are the author's explicit carve-out in any case.

The one exception is not a redundant check but a **repeated** one, and the fix keeps every check
while paying for it once instead of `reps` times.

## Honest limits

- One machine, one panel shape (10,000 observations, cohort-level `i`), one estimand family. The
  5.3% share was stable across three run lengths but was not varied across data shapes.
- The timers bracket each block including its gate test; the gate itself is a macro comparison and
  contributes nothing measurable.
- `safety(off)` in the differencing table disables more than these two blocks, so that row bounds
  all guard cost from above rather than attributing it - one reason it was the weaker instrument.
- Rarity, the author's second criterion, is not measured anywhere in this report. It is not needed:
  the recommendation does not drop a check on the grounds that its condition is rare.
