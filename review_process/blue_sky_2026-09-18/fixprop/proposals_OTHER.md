# Fix proposals — mechanism `OTHER` (single finding: SK-FC-A)

Workdir: `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-OTHER`
Patched private copy of the ado: `.../prop-OTHER/mybuild/wooldid.ado` (pinned build + 37 lines, nothing else touched).
Probes: `pA_arms.do`, `pB_bias.do`, `pC_cover.do`, `pE_guard.do`, `pF_ident.do`, `pG_time.do` (+ `.log` beside each).

---

## 0. FRAMING CORRECTION — read before the proposals

Two things about SK-FC-A as written in the brief are wrong or misleading, and both change what the
fix should be. I measured both.

**(a) `cluster()` is not the mechanism.** The finding's own repro (`probe_unbal2.do`) runs
`att robust` — no `cluster()` at all — and the coordinator separately measured point estimates as
bit-invariant to `cluster()` (mreldif 0 on 774 draws). The bias is **coarse `i` with no sub-`i`
fixed effect on a panel that is unbalanced inside its cohorts**. `cluster()` enters only because
`sthlp:255` is the sentence that talks the user into coarsening `i`. A guard keyed on the cluster
geometry as *the cause* would be keyed on the wrong thing; the guard below is keyed on it as
*the population the help's sentence creates*, which is a different and defensible claim.

**(b) The headline 21% is not bias — it is one draw of noise.** Under MCAR raggedness the
coarse-`i` error averages to zero and the interval covers. I ran 120 missingness draws on the
finding's own zero-noise DGP (`pC_cover.do`, noise sd 0.10 added so `e(V)` exists):

| design (i=coh, cluster(state), TRUE att = 1) | mean error | t(mean=0) | 95% coverage |
|---|---|---|---|
| MCAR 8% missing, G=40 | **-0.0252** | -1.23 | **0.933** (120) |
| MCAR 8%, + `fe(state)` | +0.0010 | 0.44 | 0.892 (120) |

So the refuter's `att 1.2147  p 0.000  CI [.953, 1.477]` is an honestly-priced realization: the CI
covers 1. All 32 single-draw arms in `pA_arms.log` cover 1. **Do not ship a fix on the 21%.**

**(c) But there IS a real, large, non-covering bias — under SELECTIVE attrition.** When the
missingness is correlated with the unit's own permanent level *and* with treatment status, the
error is systematic, it
**grows with G while the sd shrinks**, and coverage collapses. Same file, same DGP, missingness rule
"treated units with `s_g>0` lose post-onset rows w.p. 0.5":

| design (i=coh, TRUE att = 1) | mean error | t(mean=0) | 95% coverage |
|---|---|---|---|
| G=40, `robust` | -0.604 | -40.3 | **1.000** (120) |
| G=40, `cluster(state)` | -0.604 | -40.3 | **0.650** (120) |
| G=40, `cluster(state)` + `fe(state)` | **+0.0008** | 0.32 | 0.883 (120) |
| G=200, `cluster(state)` | -0.824 | -58.5 | **0.000** (60) |
| G=200, `cluster(state)` + `fe(state)` | **+0.0020** | 1.24 | 0.900 (60) |

Be precise about what that rule assumes: rows are dropped only for *treated* units with `s_g>0`.
Missingness correlated with a unit's own level is ordinary; correlated with treatment status as
well is the **worst case within the realistic range** — if controls attrit the same way the time
fixed effects absorb the common composition shift and the bias falls back toward a variance cost.
Neither proposal below depends on this being the typical case, only on its being reachable.

**(d) NEW MECHANISM, and the reason this belongs in "What level to cluster at" specifically.**
Compare rows 1 and 2 of that table: *identical point estimates, identical bias, coverage 1.000 vs
0.650.* Clustering below `i` removes from the standard error exactly the variation that coarsening
`i` left in the residual — the omitted unit effect is constant within the cluster, so the
cluster-robust meat matrix nets it out across periods, while `robust` treats each row's share of it
as independent noise and prices it. Direct SE comparison on one ragged panel (`pA_arms.log`,
RAGGED, i=coh): `robust` se **0.591**, `cluster(state)` se **0.130** — 4.5x smaller on the same
point estimate. **The two halves of the recommended geometry interact: coarsening `i` puts unit
heterogeneity into the residual, and clustering at the unit stops the standard error from seeing
it.** Neither half is wrong alone. Together, without `fe()`, they produce a wrong number with a
confident interval. I have not seen this stated anywhere in the project record.

**(e) The repair is exact and free.** `fe(state)` and `controls(i.state)` each return the point
estimate to 1.000000000000 with `e(r2) = 1.0` on all four panels (balanced / attrition / late
entry / ragged) under both `robust` and `cluster(state)` — 32 arms, `pA_arms.log`. The biased arms
all carry `e(r2) = 0.14`, which is the exact-fit check (harness hazard e) doing its job.

---

## PROPOSAL DOC-1 — three anchored help/message edits [RECOMMENDED, SHIP]

`risk: low · moves no reported number · runtime cost: zero`

The behaviour is not a code defect. `wooldid` fits the model the user specified; a model without
unit effects on a compositionally-shifting panel is biased, and that is econometrics, not a bug.
What *is* a defect is that `sthlp:255` names a remedy for the standard errors that silently costs
the user the point estimate, and lists `fe()` nowhere. Three edits, which must ship together —
shipping one leaves the file contradicting itself.

### Edit 1 — `wooldid.sthlp:255`, the "What level to cluster at" paragraph

Anchor (unique, byte-exact, mid-paragraph):

```
which build their inference from replications instead. Those two remedies now compose:
```

becomes

```
which build their inference from replications instead. A finer {cmd:cluster()} is not always free: with one row per unit per period there is nothing below the unit, so reaching this geometry means moving {it:i} up to the treatment-onset grouping and clustering on the unit, and that move takes the permanent differences between units out of the model. On a panel balanced inside its cohorts it is exactly innocuous - the {it:att} point estimate is bit-identical - but where a unit is absent from periods in which its cohort is present, the composition of the cohort-by-period cells moves with the missingness and the omitted unit differences move the point estimate with it, while clustering on the unit also takes that variation out of the standard error, which {it:robust} would have priced. Specify {opt fe()} at the unit level whenever {it:i} is coarser than the unit - or {it:controls()} with your own {it:i.} prefixes on the {cmd:regress}-based paths of {it:unconditionalse} and {it:customvce()}, which cannot absorb {opt fe()} - and the point estimate returns to the one {it:i} at the unit level gives. See the {it:i} entry under {help wooldid##syntax:Syntax} for how the two choices divide. Those two remedies now compose:
```

Register check against the surrounding paragraph: no capitalised emphasis, no bullets, four " - "
asides in ~150 words against the paragraph's existing three in ~190, `{cmd:}`/`{it:}`/`{opt}`
markup as used, one back-link of the kind `sthlp:67` already makes forward.

### Edit 2 — `wooldid.sthlp:70`, the `i` synopt entry

Anchor:

```
 specified separately in {opt fe()}, so a coarse {it:i} combines with unit-level ones - and go finer only where that heterogeneity, or a
```

becomes

```
 specified separately in {opt fe()}, so a coarse {it:i} combines with unit-level ones and, on a panel not balanced inside its cohorts, requires them, since the permanent differences between the units of a cohort are otherwise absent from the model and a shifting cell composition moves the point estimate - and go finer only where that heterogeneity, or a
```

This is the one sentence in the file that pairs the three parts of the working route. It currently
states `fe()` as a *preference* ("combines with"); the measurement above says it is a
*requirement conditional on unbalance*. That is the whole edit.

### Edit 3 — `wooldid.ado:1371`, the SI-8 caution text

Anchor — the short tail is **not** unique (`ado:1375`, the r(907) message, ends the same way), so
anchor on the whole recipe clause, which occurs exactly once:

```
Coarsening `i' is normally the single largest speedup available: specify unit-level fixed effects via fe() and swap `i' to a coarser level of aggregation, such as cohorts or treatment onset groups. (safety(off) suppresses this check.)"
```

becomes

```
Coarsening `i' is normally the single largest speedup available: specify unit-level fixed effects via fe() and swap `i' to a coarser level of aggregation, such as cohorts or treatment onset groups. Keep the fe() half of that recipe: on a panel that is unbalanced inside its cohorts it is what holds the point estimate fixed, not a modelling preference. (safety(off) suppresses this check.)"
```

Optional companion: `ado:1375-1376`, the r(907) message, carries the same recipe with the same
omission. It aborts the run, so the user must act anyway and it already names `fe()`; I would leave
it alone rather than lengthen an error message.

**MERGE WARNING.** Two findings in the `structure_and_ux` lane (`structure_and_ux.md:540` and
`:603`) propose an edit to this same message, adding that `att_i` and `att_itime` change definition
when `i` is coarsened. Both clauses belong; whoever lands second must merge rather than replace.
My clause is about the point estimate of `att`; theirs is about what `att_i` means.

**Blast radius.** Help text and one `dis as text` string. No estimation path, no `e()`, no number.
The SI-8 message already fires on any panel with >100 distinct `i`, so Edit 3 reaches the users
most likely to take the advice. **Verification:** rebuild the help, `help wooldid` renders; diff
the ado for the single changed string. No re-bless.

**Alternatives rejected.** (i) Documenting it only in the `i` entry (Edit 2 alone): the user
arrives at this hazard from `sthlp:255`, which is where they were sent for the standard errors, so
the caution has to be there. (ii) Removing "a finer `cluster()`" from the remedy list: it is a
genuinely effective remedy for the standard errors — the coordinator measured SE/SD 0.970 and
coverage 0.935 in that arm — and removing it would trade a measured fix for a measured defect.

---

## PROPOSAL GUARD-1 — a runtime caution on the precondition [BUILT AND MEASURED; AUTHOR MUST RULE]

`risk: low · moves no reported number (verified, 32 arms byte-identical) · cost 0.24 s at 950k rows`

### The edit

37 lines inserted in `wooldid_worker` at `wooldid.ado:1373`, immediately after the SI-8
`si8ilev > 100` block closes and before `if \`si8demand' > c(maxvar) {`. It sits **inside** the
existing `if "\`safetyinput'" != "off" & "\`wdrepmode'" != "1"` block, so it inherits the SI-8
safety gate and the replication skip for free, and it reuses `\`si8ilev'`, already computed three
lines above. Full text: `.../prop-OTHER/patch_block.txt`; applied copy:
`.../prop-OTHER/mybuild/wooldid.ado`. The predicate:

```stata
if `: word count `cluster'' == 1 & "`fe'`ccc_absorb'" == "" & strpos(" `controls' ", " `cluster' ") == 0 & strpos("`controls'", ".`cluster'") == 0 {
  tempvar skctag skcitag
  qui gegen byte `skctag' = tag(`cluster') if `ifstatement'
  qui count if `skctag' == 1
  local skcn = r(N)
  qui gegen byte `skcitag' = tag(`cluster' `i') if `ifstatement'
  qui count if `skcitag' == 1
  local skcin = r(N)
  cap drop `skcitag'
  if `skcn' > `si8ilev' & `skcn' == `skcin' {          // strictly finer, and nested
    tempvar skntc skntI
    qui gegen double `skntc' = nunique(`t') if `ifstatement', by(`cluster')
    qui gegen double `skntI' = nunique(`t') if `ifstatement', by(`i')
    qui count if `skctag' == 1 & `skntc' != `skntI'    // ragged at the cluster level
    local skragged = r(N)
    cap drop `skntc'
    cap drop `skntI'
    if `skragged' > 0 {
      dis as text "Caution: cluster(`cluster') is nested strictly inside `i' - `skcn' clusters across `si8ilev' cohorts - and no fixed effect at or below `cluster' is specified, while `skragged' of those clusters are not observed in every period their cohort is. Clustering below `i' is the intended geometry for the standard errors, but because `i' is coarser than `cluster' the permanent differences between clusters are not in the model, and on a panel ragged at that level they move the point estimate as well as the standard error. Specify fe(`cluster') - or controls(i.`cluster') on the regress-based paths - to hold the point estimate fixed. (safety(off) suppresses this check.)"
    }
  }
  cap drop `skctag'
}
```

`skcn == skcin` is nesting (each cluster lies in exactly one `i`); `skcn > si8ilev` is strictness.
The `\`ifstatement'` restriction is the screened sample the SI-8 block already uses — deliberately
**not** `smp`, which four separate findings in the `SMP` lane show is the wrong sample.

### Measured behaviour — `pE_guard.do` / `.log`, seven arms

| arm | expected | observed |
|---|---|---|
| ragged, `i=coh`, `cluster(state)`, no fe | FIRE | **fires**: "40 clusters across 3 cohorts … 23 of those clusters are not observed in every period their cohort is" |
| balanced, `i=coh`, `cluster(state)` | silent | silent |
| ragged + `fe(state)` | silent | silent |
| ragged + `controls(i.state)` | silent | silent |
| ragged, `i=state`, `cluster(state)` (cluster == i) | silent | silent |
| ragged, `i=coh`, `robust` (no cluster) | silent | silent (**deliberate gap**, see below) |
| ragged + `safety(off)` | silent | silent |

Two deliberate false negatives, both stated here so nobody files them later. The `robust` arm is
biased by the same mechanism and the guard says nothing, because with no `cluster()` there is no
sub-`i` variable to key on (rejected alternative 4). And the predicate skips on **any** `fe()`,
so `fe(region)` sitting *above* `state` silences it even though it absorbs nothing at the cluster
level; the message's "no fixed effect at or below `cluster`" is true whenever it fires, but the
screen is coarser than the sentence. Erring toward silence is the right direction for a caution,
and testing containment of `fe()` inside `cluster()` would cost a further pass for a case in which
the user has already thought about fixed effects.

### Measured cost — `pG_time.do` / `.log`

The added sequence, timed standalone on panels with 5% MCAR missingness:

| rows | clusters | seconds |
|---|---|---|
| 237,383 | 5,000 | 0.100 |
| 949,979 | 20,000 | 0.242 |
| 3,799,827 | 40,000 | 0.923 |

Linear in N. A full `wooldid y coh t ttre, att cluster(unit)` on the 949,979-row panel took
**19.53 s**, so the guard is **1.2 % of runtime** — and it is skipped entirely, at the cost of one
`word count` and two `strpos` on macros, whenever `cluster()` is absent or multi-way, or any
`fe()`/`ccc_absorb()` is given, or the cluster name appears in `controls()`. The two `nunique`
passes run only after the two cheap tag passes have already established nesting and strictness.

### Does it move a number — no, verified

`pF_ident.do` re-runs the entire 32-arm matrix of `pA_arms.do` against the patched build. Every
`ARM` line — att to 12 decimals, se to 8, `e(r2)` to 10 — is **byte-for-byte identical** to
`pA_arms.log`. `diff` is empty. No re-bless needed.

### Honest cost: this caution will fire often

On the 949,979-row 5%-MCAR panel, 18,435 of 20,000 clusters are ragged, so it fires — and by §0(b)
that panel's bias is indistinguishable from zero. **The guard detects the precondition, not the
magnitude.** I looked for a cheap statistic that separates benign MCAR raggedness from dangerous
selective attrition and did not find one: the separating quantity is the covariance between a
unit's presence pattern and its own unobserved permanent level, which cannot be computed without
fitting the model both ways. Two candidate sharpenings were tested against the designs and both
fail — "cluster absent from the whole pre- or whole post-window" does not fire on the selective
design (units keep their onset period), and "pre/post distinct-cluster counts differ within cohort"
does not either, for the same reason.

What makes me recommend it anyway: the caution can only fire on a user who has *already* chosen a
cluster strictly inside `i` with no `fe()` — the exact population `sthlp:255` creates — and the
remedy it names (`fe(cluster)`) is weakly correct for every one of them, costs nothing on the
reghdfe path, and is standard practice. That is a better false-positive profile than SI-8's, which
tells a user their run is slow when it may not be. The author declined both a refusal and a caution
on T4 (2026-09-15), so this is a ruling, not a default — but T4's caution had no free remedy.

### Alternatives rejected

1. **Change the default — absorb sub-`i` effects automatically when `cluster()` is finer than
   `i`.** Rejected: moves the reported `att` on every such run (by exactly the amounts in §0(c)),
   overrides a modelling choice the user made, and `i`/`fe()` are documented as separate axes.
2. **Refuse the specification.** Rejected: the specification is legitimate on a balanced panel,
   where the estimate is bit-identical, and a refusal there would be a false positive with a
   correctness claim attached.
3. **Post `e(wdid_subiragged)` silently and print nothing.** Same four passes, no display noise,
   and it has the author's own precedent — `e(perm_irebuilt)` was ruled the record in place of a
   display line on 2026-09-16, "Not to be reopened". Rejected as *primary* because a scalar nobody
   reads does not reach a user who is about to publish a number, but it is the right fallback if
   the author finds the caution too noisy, and the same 37 lines produce it by replacing the
   `dis as text` with an `ereturn scalar`.
4. **Key the guard on unbalance alone, with no cluster variable, so it also covers the `robust`
   arm.** Rejected on measurement of the population, not of cost: without a sub-`i` variable the
   only detectable signal is "rows per (i,t) cell vary over t", which fires on every repeated
   cross-section — where no unit effect exists to add, so the caution would be unactionable. This
   is why arm 6 above is a deliberate gap: `i=coh` with `robust` and no `cluster()` is biased by
   the same mechanism and the guard says nothing. DOC-1 Edit 2 is what covers that user.

### Test that would prove it correct

Already run, and it needs no oracle: the exact-fit design (`pA_arms.do`) in which every term lies
in the model's span when `i` is the unit or `fe(unit)` is given, so `e(r2) == 1` and truth is
exactly 1.000000000000. The invariant to assert in a permanent test is **not** a number but an
identity: on any panel, `wooldid y coh t ttre, att cluster(u) fe(u)` and `wooldid y u t ttre, att
cluster(u)` must return the same `att` to `mreldif < 1e-10`, and the guard must be silent in the
first and fire in neither when the panel is balanced. Both legs hold in `pA_arms.log` /
`pE_guard.log`. `att_i` and `att_itime` are deliberately excluded from that identity — they average
over whatever `i` is and change definition by construction when `i` is coarsened.

---

## Not pursued, recorded for whoever wants it

On the zero-noise design, `i=coh` with `cluster(state)` posts **no `e(V)` at all** —
`_se not allowed when e(V) is not present`, r(111) — on the balanced, attrition and late-entry
panels, while the ragged panel returns `se = 0.1296` (`pA_arms.log`, the four
`vce=cluster(state) i=coh` arms print `se = .`). `known_bydesign` §5 says an all-zero covariance under `robust` or
`cluster(unit)` on realistic data is a finding; mine is on a zero-residual-variance design and is
almost certainly degenerate, so I did not chase it. Flagging it because §7 lists a neighbouring
undetermined item (`lattice(all 3)` all-zero VCV).

Also observed and **explained, not open**: the *repaired* arms (`fe(state)`, cluster strictly
finer than `i`) cover 0.883 / 0.900 / 0.892 across my three MC cells. That is below nominal, but it
is not a finer-than-`i` failure and should not be filed as one — `pC_cover.do`'s `mkpanel` assigns
onsets by `mod(state,10)==1` and `==2`, so at G=40 each treated cohort holds only **four** clusters
(fewer after attrition). Coverage of 0.88–0.90 on four clusters is the ordinary few-cluster result
and a different design from the coordinator's 0.935. Recorded so a verifier does not chase it.
