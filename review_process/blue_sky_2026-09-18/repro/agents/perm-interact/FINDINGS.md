# perm-interact — findings
Agent: adversarial prober, target **permute() x everything else / the hand-maintained
forwarding allow-list**, method **interaction matrix, hunting SILENT NO-OPS**.
Build: pinned `.../bluesky/build/wooldid.ado` (v1.4, sha c1d8979), asserted with
`which wooldid` + a `nofastcollapse` parse in `p01_base.log`.

---

## PX-01 (SELF-FOUND, CONFIRMED) — `fastnose` silently overrides `permanalyticse`, and the run then ASSERTS the analytic standard errors and `e(V)`/`e(V_es_*)` exist when they do not

**Surface:** `permute()` x `fastnose` x `permanalyticse`
**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact/PX01_repro.do`
— fully self-contained: it generates its own 40-unit panel and needs nothing but the
pinned build. (`p10_fnase.do` in the same directory is the wider screen it came from;
that one reads `base.dta`, which `mkdata.do` rebuilds.)

**Claim.** Under `permute()`, `wooldid` refuses `fastnose` together with `permtstat` and
with `permtci` — r(198), "permtstat/permtci randomize on t-statistics, which fastnose does
not compute" — but ACCEPTS `fastnose` together with `permanalyticse`, which asks for the
same object. `fastnose` then wins silently: no analytic standard errors are computed.
The run nevertheless posts `e(perm_analyticse) = 1` and prints two statements that are
false on that run.

**Observed** (`p10_fnase.log`, tag `fn_ase`, 40 units / 320 rows, `cluster(unit)`,
`esprelength(2) espostlength(2) permute(shuffle) permreps(30)`):

| run | `e(perm_analyticse)` | `e(wdid_vcvstate)` | `e(V)` | `e(V_es_att)` | `e(b_es_att)` | se_analytic col |
|---|---|---|---|---|---|---|
| fast-mode default      | . | absent          | no  | no  | no  | .        |
| `permanalyticse`       | 1 | analytic-blocks | yes | yes | yes | .1690174 |
| `fastnose`             | . | absent          | no  | no  | no  | .        |
| **`fastnose permanalyticse`** | **1** | **absent** | **no** | **no** | **no** | **.** |
| `fastnose permtstat`   | r(198) refused | | | | | |
| `fastnose permtci`     | r(198) refused | | | | | |

and the run prints, over an entirely blank `se_analy~c` column:

* header — `Wooldid Estimation for Outcome: y; Standard Errors: randomization inference
  (permute(shuffle), 30 replications); se_analytic = the analytic cluster(unit) standard
  error, shown for reference only`
* footer — `The standard-error column is the ANALYTIC standard error (se_analytic,
  permanalyticse), shown for reference only; ... e(V) and e(V_es_*) hold the analytic
  covariance, so test or lincom run by hand on them are analytic, not permutation,
  quantities.`

`e(V)`, `e(V_es_att)` and `e(b_es_att)` are all absent on that run, so the footer's
instruction cannot be followed and the header describes a column that is empty.

**Expected.** Three of the four available behaviours would be fine; the shipped one is not.
(a) Refuse the pair, as the `permtstat`/`permtci` twin is refused with reasoning that
applies verbatim — the help itself treats `fastnose` as acceptable only "under the
`permute()` fast-mode default ... (no `permtstat`, **no `permanalyticse`**), where fast
mode is already in force", i.e. it does not contemplate this pair being in fast mode;
(b) let `permanalyticse` win, since the help says it "restores" the SEs and it is the more
specific request; (c) accept but say what happened. What ships is (d): the option is a
silent no-op and the run then asserts in three places — `e(perm_analyticse)`, the table
header and the footer — that it took effect.

**Mechanism** (`ado:5448`): `if "`permute'" != "" & "`permtstat'" == "" &
"`permanalyticse'" == "" local fastnose fastnose` only ADDS the fast path in the default
case; it never clears a user-typed `fastnose`. `wooldid_permctl`'s guard (`ado:5881-5883`)
gates on `ts`/`tc` only, not on `permanalyticse`. `wooldid_permwire` then keys its header
(`ado:8214`), `e(perm_analyticse)` (`ado:8207`) and the display footer (`ado:11679`) on
the requested `permanalyticse` rather than on whether a covariance was actually posted.
The true state is computed seven lines later — `cap confirm matrix e(V)` / `local wdphasv =
(_rc == 0)` at `ado:8220-8221`, which correctly records `absent` in
`e(wdid_vcvstate)` — and is never fed back into the header or into
`e(perm_analyticse)`, which is why the same run says both things at once.

**realistic:** yes, and the help invites the first half of the pair: the `fastnose` entry
says that naming `fastnose` under `permute()` is "accepted rather than an error" because
fast mode is already in force, so a user types it for speed; adding `permanalyticse`
afterwards to see a standard error is following the `permanalyticse` entry, which promises
it "restores `e(V)`, `e(b_es_*)` and `e(V_es_*)` (analytic)".

**Honest calibration of the damage.** Nothing downstream silently returns a wrong number:
`test`/`lincom` fail loudly (r(321)) because there is no `e(V)` (measured on the fast-mode
default in `p05_down.log`), and `esttab`/`estout` should be unaffected because permwire's
`!wdphasv` branch still posts `e(p)`/`e(ci_l)`/`e(ci_u)` — that branch does fire here, since
`e(wdid_vcvstate)` came back `absent` (measured), though I did not run `esttab` itself. The damage is (i) the requested option silently does nothing, (ii)
`e(perm_analyticse) == 1` is a machine-readable claim that is false, which a script
branching on it will act on, (iii) two printed sentences are false, one of which tells the
reader to run `test`/`lincom` on a covariance that does not exist, and (iv) the
`e(b_es_*)`/`e(V_es_*)` honestdid/pretrends handoff the help promises is withdrawn with no
message.
**severity:** medium. **confidence:** certain.

### PX-01 addendum — the codebase already has the predicate it needs

`ado:5396` defines exactly the right test and uses it two lines later:

```stata
local fnperm = ("`permute'" != "" & "`permtstat'" == "" & "`permanalyticse'" == "")
if "`adversarial'" != "" & !`fnperm' { ... error 198 }
if "`jointtests'"  != "" & !`fnperm' { ... error 198 }
```

Measured (`p11_gaps.log`): `fastnose jointtests permute(...)` is ACCEPTED (rc 0) but
`fastnose permanalyticse jointtests permute(...)` is REFUSED r(198) — "fastnose is not
compatible with jointtests; fast mode reports point estimates only." So the program
already decides, correctly, that `fastnose`+`permanalyticse` is NOT the fast-mode default.
`ado:5448` then reuses the same three-term condition to decide whether to ADD `fastnose`,
but no screen anywhere refuses a USER-typed `fastnose` when `!fnperm`. The one-line fix is
the mirror of `wooldid_permctl`'s existing `ts`/`tc` guard.
Also asymmetric: `bootstrap()` refuses `fastnose` outright ("bootstrap() may not be
combined with fastnose; fastnose requests point estimates only, leaving nothing for the
bootstrap to attach inference to").

---

## PX-02 (SEEDED LEAD, CONTRADICTED) — `customvce()` can never reach a permutation replication: `ado:128` is dead code, and the "permute forwards customvce, bootstrap refuses it" asymmetry has no behavioural consequence

**Surface:** `permute()` x `customvce()`
**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact/p11_gaps.do` (first two commands)

The seeded lead says the twins are asymmetric because `customvce()` "is forwarded to
permutation replications while bootstrap REFUSES customvce() with reasoning that applies
verbatim to permutation p-values". The source asymmetry is real (`ado:128` emits
`customvce(...)` into `wdp_ro`; `ado:6414` refuses it under `bootstrap()`), but the
combination is **unreachable**, so no user can be affected:

* with `cluster()`: `ado:4667` errors first — "customvce() is not compatible with robust
  and cluster(). You must manually specify the contents of reg's vce() statement in this
  case." (measured, `p03_refuse.log` tags `cvce_cl`, `cvce_nocl`)
* without `cluster()`: `wooldid_permctl` errors — "permute() requires cluster() with a
  single variable - the unit whose treatment-onset assignment is randomized." (measured,
  `p11_gaps.log` tags `cvce_truly_nocluster`, `cvce_robust_nocluster`, both r(198))

and `permute()` requires exactly one `cluster()` variable, so there is no third case.

**severity:** low (dead code + a misleading in-source asymmetry, no user-visible effect).
**realistic:** no — nobody can reach it. **confidence:** certain.
Reported because the coordinator note asks for contradictions of the record, and because a
future relaxation of either guard would silently switch the line on.

---

## Seeded lead (1) — gteffects x permute: **CLEAN**
`p05_down.log` (tags A1/A2). A plain run posts `e(gteffects)` 144x15 with
`se`/`p`/`lb_95ci`/`ub_95ci` filled (e.g. row 1: .2107, .0016, -1.1396, -.2872). The same
run under `permute(shuffle) permreps(30)` posts the SAME 144x15 matrix with the estimates
intact and those four columns **missing on every row** — not stale analytic values.
`e(perm_gt_blank)` is 1, the display prints "The gteffects cells are keyed by cohort onset,
which the replications reassign, so they have no permutation distribution: their inference
columns are blank" (`ado:11726`), and `sthlp:1150` says the same in the `gteffects` entry.
Documented, implemented, and told to the user on every run. No finding.

## Seeded lead (2) — permute() posting an analytic covariance under a permutation header: **CLEAN**
`p05_down.log` (VCHK block).
* Fast-mode default: `e(wdid_vcvstate)` = `absent`, no `e(V)`, no `e(V_es_att)`,
  `e(vce)`/`e(vcetype)` cleared. `test main_att = 0` → r(321), `lincom main_att` → r(321),
  `estat vce` → r(322). Nothing is silently computed; the footer says
  "No analytic standard errors were computed ... and no e(V) or e(V_es_*) is posted."
* `permanalyticse` / `permtstat`: `e(V)` and `e(V_es_att)` present and ANALYTIC;
  `test` gives p = 4.898e-13 and `lincom` se = .1690174, both analytic; and the footer
  states it verbatim — "e(V) and e(V_es_*) hold the analytic covariance, so test or
  lincom run by hand on them are analytic, not permutation, quantities" (`ado:11677/11679`),
  with the same sentence in the help. No finding.
  (PX-01 is the one case where that footer is printed and is FALSE.)

---

## PX-03 (SELF-FOUND, CONFIRMED) — `semielasticity` on an outcome that is not strictly positive fails inside `margins` and is never diagnosed; under `fastnose` — which is what `permute()` selects by default — the run returns **rc 0 with a printed results table and no `e(b)` at all**

**Surface:** `permute()` (fast-mode default) x `fastnose` x `semielasticity`
**Repro:** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact/PX03_repro.do`
(reproduces all four modes and both controls; its only dependency is the review's shared
data generator, included by absolute path. See also `p16_semi.do`, `p17_semi2.do`,
`p19_bracket.do`.)

**Claim.** `semielasticity` (margins' `eydx`) divides by the prediction. When the outcome
takes non-positive values, `margins` fails. `wooldid` never says so. The failure surfaces
in four different ways depending on the panel and the options, and on one of them the run
reports **success**.

**Observed** (`PX03_repro.log`; panel A is 40 units / 8 years / onsets 2003 and 2005 with
heterogeneous effects, y in −3.70 .. 6.86; panel B is a plainer inline panel,
y in −2.63 .. 5.33):

| run | rc | `e(b)` | `e(cmd)` |
|---|---|---|---|
| A, plain (no `semielasticity`) | 0 | `main_att` | wooldid |
| **A, `semielasticity fastnose`** | **0** | **ABSENT** | wooldid |
| A, `semielasticity permute(shuffle) permreps(30)` | 459 | absent | wooldid |
| A, `semielasticity` (no fastnose) | 3598 | internal `0b.__000006#...` | **reghdfe** |
| B, `semielasticity fastnose` | 459 | internal | **reghdfe** |
| A and B, the SAME outcome + 10 (strictly positive) | 0 | `main_att` | wooldid |

so the trigger is exactly "the outcome is not strictly positive": a constant shift of the
same outcome makes every row clean. A plain run on the same outcome is fine, and `fastnose`
alone on the same outcome posts `e(b)` normally, so neither `fastnose` nor the data is the
cause on its own.

What the **rc-0** case prints (`p16_semi.log`, tag `semi_fn`) — no error, no caution, and
then `matrix e(b) not found`:

```
Wooldid Estimation for Outcome: y; Standard Errors: not computed (fast mode)
N clusters = 40; N = 320 (0 obs dropped from initial sample)
R2 = 0.8989;  R2adj = 0.8178;  R2-within = 0.6799;  R2-withinadj = 0.5063

Main Results (Full Estimation Sample): Treatment Effects Presented as Semi-elasticities
             |   estimate          se         t         p     lb_95ci
main         |
         att |          .           .         .         .           .
```

The other three messages are all Stata internals, none of them a `wooldid` sentence:
* `permute()`: `Error: the observed comparison statistic for att could not be read from the
  rep-0 run; permutation inference is impossible.` + `something that should be true of your
  data is not` — which points at the comparator run and at the data, when in fact the
  OBSERVED run produced no coefficient vector either.
* no `fastnose`: `variable __00001O not found / st_data(): 3598 Stata returned error /
  __marg_cr_eval() / deriv__call2user_v() / _deriv__compute_value() / _deriv() /
  st__marg_cr_compute()`.
* panel B: `inconsistent estimation sample levels 0 and 1 of factor __000007` — Stata's own
  `margins` error, naming a `wooldid` tempvar.

On the two error paths `e()` is left holding the internal stage-1 model: `e(cmd)` is
`reghdfe` and `e(b)`'s names are tempvar factor terms. `estimates store` and `est table`
then operate on that object with rc 0 (measured), while `_b[main_att]` gives r(111).

**Expected.** `wooldid` screens for comparable conditions and says so in its own words — on
the same data under `poisson` it prints "ppmlhdfe dropped every observation of the treated
cohort-period cell unit == 17, relative-t == 4 (perfect separation: the outcome is zero on
all of them, so the cell's coefficient does not exist). This is not multicollinearity...".
`semielasticity` on a sign-crossing outcome deserves the same treatment. At an absolute
minimum, a run that produced no coefficients must not return rc 0, and must not print a
formatted table with an R² line and the header "Treatment Effects Presented as
Semi-elasticities" as if it had estimated. The help's only hedge is "this option tends to
require large samples and a y-variable that is not sparsely distributed" (`sthlp`,
`semielasticity` entry); it does not say the outcome must be strictly positive, and it does
not warn that failure can be silent.

**realistic:** yes. Asking for a semi-elasticity of an outcome that can be negative is a
specification mistake, but an ordinary one — a net, differenced or de-meaned outcome run
through a loop of specifications — and the point is that nothing tells the user.
`permute()` makes it worse twice over: `permute()` *selects* `fastnose`, so the rc-0 branch
is the one a permutation user is on, and the error they eventually get names the wrong
thing.
**severity:** high. Not a wrong number, but a run that reports success, prints a table
headed "Semi-elasticities" with an R² beneath it, and returns nothing — in a scripted sweep
that is a silently empty column rather than a visible failure.
**confidence:** certain for the behaviour and for the trigger (pinned by the +10 shift
control on two independent panels). The rc-0 branch specifically is panel-dependent: of the
six panels tried, one produces it and the others produce r(459) or r(3598); all six fail,
none is diagnosed.

---

## What else was probed and came back CLEAN

### The comparator detector (novel instrument) — 65+ configurations, all exactly 0
`wooldid_permrun` runs a mandatory rep 0 on the TRUE assignment under the rep option
string, and `wooldid_permwire` posts that vector as `e(perm_comparator_b)` (`ado:8169`).
The code comment at `ado:7651` states that under the fast-mode default "the comparator then
equals it bitwise" — and nothing asserts it. So
`max |e(perm_comparator_b) − e(b)| / |e(b)|` over name-matched columns is a direct test of
"the replications are fitting the model the user asked for", catching both dropped and
mutated options with no RNG dependence.

Result: **exactly 0.0** on every one of 65+ configurations across
`p01_base.log`, `p02_more.log`, `p04_escmp.log`, `p12_coarse.log`, `p13_contreat.log`:
`controls()`, `fe()` (on a time-varying group), `coarsecohortcontrols(c.·)`,
`interactivecontrols()`, `timetrends`, `subgroups()`, `esprelength()/espostlength()`,
`esfixedbaseperiod`, `esrelativeto()`, `espostlength(-999)`, `anticipation()`,
`oostreatedcontrols()`, `customweightmultiplier()`, `nofastcollapse`,
`emptycellsoverride`, `regtol()`, `safety(off)`, `level()`, `robust`, `gteffects`,
`jointtests`, `adversarial`, `permstepdown`, `permstrata()`, `[aw=]`, `[pw=]`, an `if`
qualifier, weights and `if` together, weights+`if`+`subgroups()`, each of `att`/`att_i`/
`att_it`/`att_itime` and all four at once, `permtstat`, `contreat()` with
`contreatcontrols(basic)`/`(byTime)`, `contreatpoly(2)`,
`contreatelasticitytype(dyex)` — and the same sweep repeated with `cluster()` STRICTLY
COARSER than `i` (20 paired clusters, unit-level `i`), which is the geometry `permute()`
exists for and which `p01`/`p04` did not cover.

**The detector has power** (positive control): `permanalyticse` makes the observed run full
while rep 0 stays `fastnose`, and there the same statistic is 7.3e-11 on a near-zero
pre-period coefficient — exactly the F9 fastnose-vs-full contract the help documents. So
the zeros above are informative, not a broken instrument.

Secondary screen: toggling each option on/off at a fixed `rseed()` and fingerprinting
`e(permdist_b)` (sum and sum of squares, in hex). Every option that moved `e(b)` also moved
the null distribution; every option whose fingerprint did not move also left `e(b)`
unmoved (`fe(region)` absorbed by the unit effects, `oostreatedcontrols(2)` with nothing
out of sample, `robust` alongside `cluster()`, `nofastcollapse` — documented to change no
answer — `gteffects`, and `esrelativeto()` without `esfixedbaseperiod`, which prints
"Caution: esrelativeto() has no effect without esfixedbaseperiod and is being ignored").
**What would have shown up:** an option that silently never reached the replications would
have left the fingerprint bit-identical while `e(b)` moved; an option forwarded in mutated
form would have shown a nonzero comparator difference. Neither occurred.

### The refusal matrix — every refusal fired, and every message named the cause
Measured in `p03_refuse.log`, `p07_ct_onc.log`, `p11_gaps.log`, `p12_coarse.log`,
`p13_contreat.log`: `impvalues`, `unconditionalse`, `contreatwithin`,
`contreatcontrols(untreatedZero)`, `suppresssgprimaryeffects`, `pretrendpower()`,
`oldsyntax`, `simultaneousbands`, multi-way `cluster()`, `fastnose`+`permtstat`,
`fastnose`+`permtci`, `permabsp`+`permrawp`, `bootstrap()`+`permute()`,
`bootstrap()`+`fastnose`, user-supplied `[fw=]`, `permonsetcohorts`+`subgroups()`,
a `cluster()` that does not nest `i` without `permonsetcohorts`, a `permstrata()` that
straddles clusters, a `cluster()` that crosses onset groups, and a second `permsave()` to
an existing file (r602). `cleanmatrices`+`permute()` runs clean under `adopath ++` loading
(the SI-17 reversion holds), and `verbose`, `summarystats`, `histogramcohorteffects` and
`makeplots`+`saveplots()`+`subgroups()` all complete, the histogram one under its
documented analytic caution.

### `permonsetcohorts` (T1, one day old) x permute
`p07_ct_onc.log`: `i`=cohort with `cluster(unit)` and no assertion is refused r(459) with a
message that names the remedy; with `permonsetcohorts` it runs, `e(perm_irebuilt)`=1, and
the comparator is exactly 0 with `esprelength/espostlength`, with `controls()` and with
`gteffects`; with `cluster()==i` the rebuild correctly reports `e(perm_irebuilt)`=0 (a
declared no-op) and `e(perm_npossible)`=6 against 1.7e17 for the finer cluster.
`permonsetcohorts`+`subgroups()` is refused with a message that explains the r(411)
mechanism.

### `contreat()` / the ame family x permute
`p13_contreat.log` (i=cohort, 6 cohorts, cluster(cohort), unit-level dose):
`contreatcontrols(basic)`, `(byTime)`, `contreatpoly(2)` and
`contreatelasticitytype(dyex)` all run, post `main_att` and `main_ame`, and the comparator
is exactly 0 on every column including `main_ame` and the 12-column ES run.
`permtstat`+`contreat` at `cluster()==i` stops r(459) — but only after the run has printed
the cell-containment caution AND "the estimation command returned an all-zero
variance-covariance matrix, so no standard error in this run is identified", which is the
documented geometry; the failure is diagnosed, not silent.

### Independent exact-enumeration oracle (novel) — wooldid's permutation distribution equals a hand enumeration BITWISE, with a working negative control
`p06_oracle.do` / `p06_oracle.log`. A 6-cluster design whose onset multiset
{2003, 2003, 2005, ., ., .} admits exactly 60 distinct shuffle assignments.
`wooldid ..., permute(shuffle) permexact` reported `e(perm_npossible)` = 60,
`e(perm_completed)` = 60, `e(perm_failed)` = 0. I then enumerated the same 60 assignments
myself (Mata `cvpermute` over the onset multiset), refitted `wooldid ... , att
cluster(unit) fastnose <the same options>` once per assignment, and compared the sorted
hand distribution against `e(permdist_b)`:

| arm | wooldid draws | oracle draws | max relative difference |
|---|---|---|---|
| no extra options              | 60 | 60 | **0.0** |
| `controls(zbig)`              | 60 | 60 | **0.0** |
| `timetrends`                  | 60 | 60 | **0.0** |
| `gteffects` (not on the allow-list) | 60 | 60 | **0.0** |
| **negative control**: the `controls(zbig)` run against the NO-controls oracle | 60 | 60 | **0.946** |

So the replications really are fitting the model the user specified — not merely a
self-consistent one — the enumeration is the right set of 60 assignments, `controls()` and
`timetrends` genuinely reach the replications, `gteffects` is genuinely a no-op there, and
the comparison would have detected a mismatch had there been one (the negative control
separates by 95%).

### `safety()` x the replication failure set — CLEAN, measured
`p08_safety.do`, a 30-unit panel with 45% late entry and 25% early exit and 200
replications: `e(perm_failed)` = 0 under `safety(on)`, `safety(off)`, `emptycellsoverride`
and with an ES window, and the p-value and null SD are bit-identical across all of them.
So on this design the estimability screens do not shape the reference distribution. What
would have shown up: a different `e(perm_failed)` between the two safety settings (a
selection effect on the null), or the same failure count with a different p. Neither.
The help documents the behaviour anyway ("each draw is the estimate wooldid would compute
under that counterfactual assignment, estimability screens included. Replications that fail
to estimate are excluded and counted (e(perm_failed)), with a caution above 10 percent").

### Replay, state leakage, downstream tools — CLEAN
`p18_replay.do`: after a `permute()` run with subgroups and an ES window, a bare `wooldid`
replays with `e(b)`, `e(perm_p)`, `e(perm_ci)` and `e(perm_kernel)` bit-identical
(max difference exactly 0) and `e(perm_p_att)` equal in hex.
`p09_odd.do`: `honestdid` refuses r(198) after a fast-mode permute run (no `e(V)`) and runs
after `permanalyticse` on the analytic covariance, which the footer discloses; the next
plain `wooldid` run has `e(perm_method)` empty, and two consecutive plain runs after a
permutation run give bit-identical `e(b)` (no state leak).
`p09`/`p12` also cleared: `makeplots`+`saveplots()`+`subgroups()`, `timetrends`+`jointtests`
+`adversarial`, `espostlength(-999)`+`subgroups()`, `level(97.5)`+`gteffects`, `permrawp`,
`permside(right)`, `cleanmatrices`, `verbose`, `summarystats`+`histogramcohorteffects`.

### known_bydesign §7 item 11 — SETTLED, and the suspicion does not hold
The coordinator's list asks: "Check whether `fastnose` actually posts
`e(b_es_att)`/`e(numpre_att)`; if it does not, this is a live doc-vs-code divergence."
Measured (`p15_item11.do` / `p15_item11.log`, no `permute()`, `esprelength(2)
espostlength(2)`):

| run | `e(b_es_att)` | `e(V_es_att)` | `e(numpre_att)` | `e(wdid_vcvstate)` |
|---|---|---|---|---|
| plain      | matrix  | matrix | scalar 2 | analytic-blocks  |
| `fastnose` | **matrix** | absent | **scalar 2** | absent |
| `impvalues`| absent  | absent | absent   | analytic-diagonal |

So `fastnose` DOES post the coefficient vector and the pre-period count without the
covariance, under the caution "no horizon of the att event study has an identified
variance, so e(b_es_att) posts without e(V_es_att)..." — exactly what `sthlp:1325` says of
the list it heads. Keeping `fastnose` on that list is correct; removing it would have made
the help wrong. And `impvalues` is already OFF that list on the pinned build (the sentence
now reads "`fastnose`, typically `timetrends`, and any run whose estimation returns an
all-zero covariance"), which matches the measurement that `impvalues` posts none of the
three. **No divergence.** Reported because the coordinator named it as open.


---

## Severity rationale (so the calibration is visible)
* **PX-03 high.** The rc-0 branch is undetectable: no rc, no message, a formatted table.
  A crash would be better. The three error branches on their own would be medium
  ("uninformative error on realistic input"); it is the rc-0 branch that lifts it.
  The specification (`semielasticity` on an outcome that can be negative) is itself a user
  error — the defect is that `wooldid` neither refuses it nor names it, and on one branch
  reports success.
* **PX-01 medium.** No number is wrong and nothing downstream silently misreports; the
  damage is a requested option doing nothing, plus three false statements about it, one of
  them machine-readable.
* **PX-02 low.** Dead code with no user-visible effect; filed only because it contradicts
  the seeded lead and because relaxing either guard would switch the line on.

## Secondary observation (not filed as a finding)
On the two `semielasticity` error paths, `e()` is left holding `wooldid`'s internal stage-1
model: `e(cmd)` becomes `reghdfe` and `e(b)`'s column names are tempvar factor terms
(`0b.__000006#0b.__000007#...`). `estimates store` and `est table` then operate on that
object with rc 0 (measured in `PX03_repro.log`), while `_b[main_att]` gives r(111).
`wooldid`'s own replay is protected (`WOOLDID_LASTRC` refuses it), and the source already
notes at `ado:7736` that "nothing in wooldid clears e() on failure", so this is a known
property rather than a surprise — but a `cap wooldid ... ; est store m` loop will store the
internal model under the user's name. Reachable from any failure after the stage-1 fit, not
just from `semielasticity`; I only measured it there, so it is recorded rather than filed.

### p-value conventions and p-to-interval coherence — MEASURED, all five conventions exact on every column
`p14_pconv.do` / `p14_pconv.log`. For every reported column I recomputed the p-value by
hand from `e(permdist_b)` (divided by `e(permdist_se)` under `permtstat`) and
`e(perm_comparator_b)`, applying the documented rule — plain shares under a complete
enumeration, add-one `(1 + count)/(n + 1)` for sampled draws, and plain `c/n` for
`permrawp` on both branches — and compared it against `e(perm_p)`. I also checked the
help's promise that "the interval excludes zero exactly when the p-value falls below the
level" on every column, for the two conventions where it is supposed to hold (the default
and `permabsp`; `permside()` keeps a two-sided interval under a one-sided p by design).

15 runs on a 6-cluster / 60-assignment design: 8 under `permexact` (default, `permabsp`,
`permrawp`, `permside(left)`, `permside(right)`, `permtstat`, `permtstat`+`permabsp`,
`level(90)`) and 7 under sampled draws with `permreps(80)` (the same five conventions,
`permtstat`, and a `subgroups(sg2)` run with 18 columns).
**MISMATCH = 0 and COHERENCE_VIOL = 0 on every run, worst absolute difference 0.000e+00
over 108 column-checks.** `e(perm_p_lo)` and `e(perm_p_hi)` also match the hand shares
exactly (e.g. exact: .666667 / .350000; sampled add-one: .641975 / .407407; on the t scale:
.700000 / .316667). What would have shown up: a convention applying the add-one where the
help says plain shares or vice versa, a `permside()` that reported the wrong tail, or an
interval quantile that disagreed with its own p — any of which is a wrong number on an
inference surface.

Corroborating code read: `wooldid_permtail` (`ado:15945-16068`) implements exactly those
formulas, and turns the share convention on **membership** (`memb` — is the observed
statistic one of the n reference values) rather than on completeness, so an element under
`permexact` whose own observed row failed keeps the add-one form and is named in
`e(perm_addone)`.
