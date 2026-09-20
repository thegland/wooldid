# Fix proposals — mechanism `SG` (subgroups())

Eight critical/high findings. They fall on **three fix sites**, and every finding on a site
shares one edit.

| # | proposal | findings | moves a reported number? | risk |
|---|---|---|---|---|
| SG-1 | fast mode's continuous subgroup-by-horizon call keeps its VCE | FASTNOSE-03 | only in `fastnose` runs, toward the full run | low |
| SG-2 | refuse a partially-labelled unit (F13-PM) | fpsg-01, skp-01, skp-02 | no (refuses instead) | low |
| SG-3 | rescale the bootstrap covariance to the printed SEs | hv-bootrank-03, skep-sgexport-01, skep-threese-02, SK-03 | **yes — e(V), e(V_boot), e(V_es_\*), pretrendpower()** | medium |

Everything below was **applied to my own copy** of the pinned build
(`probes/prop-SG/ado/wooldid.ado`) and run. Logs and do-files are in
`probes/prop-SG/`. The repo was not touched.

The three patches were also verified **together**: with all of them applied, SG-1's repro comes
back byte-identical to the run in which only SG-1 was applied, and SG-2's repro comes back
numerically identical across the `wdrepmode` refactor. They do not interact.

---

## SG-1 — `fastnose` + `contreat()` + `subgroups()`: give the continuous subgroup-by-horizon call the same VCE its binary twin already gets

**Covers:** FASTNOSE-03.

### The edit

Two hunks, both in `wooldid_worker`'s continuous-treatment subgroup block.

**(a) `wooldid.ado:3285`** — one token.

```stata
* NOW
qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstextC' `marginsdof' `fcnoes'
* BECOMES  (only the penultimate macro changes)
qui margins [`weight'=`w_`horizitem''] , `contreatelasticitytype'(`contreat') subpop(if `ifstatement' & `targetflag'==1  & `_t' >= `horizonstarter' & `_t' <= `horizonend' & `subgroups' == `sgcase' ) noestimcheck   over(`_t') `poismarginstextC' `sgmarginsdof' `fcnoes'
```

This is the exact twin of the **binary** subgroup-by-horizon call at `ado:2912`, which already
threads `sgmarginsdof`. `sgmarginsdof` is built at `ado:1980-1981` and differs from
`marginsdof` **only** under `fastnose` with `subgroups()` and an ES window, where it is forced
to `df(`df_r')` so margins' *not-estimable* omission marker survives for the SI-7 filter in
`wooldid_sglbl` to read. The comment at `ado:1975-1979` states the rationale and calls it
"this ONE call"; the continuous twin was simply missed.

**(b) after `wooldid.ado:3317`** — the fast-mode blank the binary twin has at `ado:2975-2983`.
Without it, hunk (a) alone starts *publishing* delta-method standard errors in fast mode
(measured: `fn03_patch_v1.log` printed `ES_ame post0 se .21955975 t 5.08 p 0.000`), which
contradicts `sthlp:290`'s "`e(V)` is not posted ... fast mode reports no inference".

```stata
                matrix `___`horizitem'_Choriz_s`sgcase'' = ( `f9tab'', `relyrcolumn_Cs`sgcase'')
+               * BS-26 twin: sgmarginsdof keeps THIS call's variance in fast mode so the SI-7
+               * estimability markers exist for wooldid_sglbl to filter on; the variance itself
+               * is not part of fast mode's output, so discard it here and leave the
+               * se/t/p/interval cells missing, exactly as the binary twin does above
+               if "`fastnose'" != "" {
+                 forval fnr = 1/`=rowsof(`___`horizitem'_Choriz_s`sgcase'')' {
+                   forval fnj = 2/6 {
+                     matrix `___`horizitem'_Choriz_s`sgcase''[`fnr', `fnj'] = .
+                   }
+                 }
+               }
                * AMESE R1 scrub: null se/t/p/ci columns of the per-sg ame ES rows
```

(The blank is applied to the assembled matrix rather than to `f9tab`, matching the AMESE R1
scrub idiom two lines below; the binary site blanks `f9tab` before the transpose. Equivalent.)

**Comment upkeep:** `ado:1974-1979` and `ado:2935-2937` say "this ONE call" — now two.
**Help upkeep:** `sthlp:292` "the one exception is the subgroup-by-horizon margins call ... at a
small cost per subgroup" → "calls ... at a small cost per subgroup per estimand family".

Nothing else moves: the F2 b/V capture at `ado:3305-3331` stays keyed on `marginsdof`
(still `nose`), so fast mode still posts **no** `e(V_es_ame_s#)` — BS-26's other half is intact.

### Measured before / after

`probes/prop-SG/fn03_base.do` vs `fn03_patch.do` (copies of the finding's own
`repro_FN03.do`, adopath retargeted, nothing else changed).

| | analytic | fastnose BEFORE | fastnose AFTER |
|---|---|---|---|
| `e(b_es_ame_s1)` colnames, `untreatedZero` arm | `post0 post1 post2 post3` | `pre4 pre3 pre2 pre1 post0..post3` | **`post0 post1 post2 post3`** |
| `e(numpre_ame_s1)` | 0 | 4 | **0** |
| `e(b_es_ame_s2)` colnames, `basic` arm | `post0` | `pre4 pre3 pre2 pre1 post0 post1` | **`post0`** |
| `e(numpre_ame_s2)` | 0 | 4 | **0** |

* The surviving ame values are **bit-identical** before and after
  (`1.005702736  .9271277637  1.13950981  1.239007685`) — the patch removes the structural-zero
  columns and perturbs nothing.
* The **ANALYTIC** arms are **byte-identical** before and after (`ANALYTIC_IDENTICAL` in the
  run log) — `sgmarginsdof == marginsdof` off `fastnose`, so the patch is inert there.
* The **no-`subgroups()` control arm** is byte-identical before and after (`ARM2_IDENTICAL`).
* After hunk (b), the fastnose `ES_ame` rows print `se t p ci` as missing, as fast mode must.

### Cost

`fn03_cost.do`: 60,000 rows, 30 `i` values, 5 subgroups, `contreat()`, `esprelength(6)
espostlength(6)`, `set processors 2`.

| | pinned | patched |
|---|---|---|
| `fastnose` | 103.195 s | **122.129 s**  (+18.9 s, +18.3%) |
| analytic (reference) | 145.223 s | 144.708 s  (-0.4%, noise — the patch is inert here) |

The added work is **one delta-method VCE per subgroup per estimand family** — exactly the cost
BS-26 already accepted on the binary side, and exactly what `sthlp:292` already advertises
("at a small cost per subgroup"). It is bounded by the number of subgroups, **not** by N: the
margins call's subpop and weight pass are already paid; only the VCE leg is added. It does not
fire outside `fastnose` + `contreat()` + `subgroups()` + an ES window.

### Moves a reported number?

**Yes, but only inside `fastnose` runs, and only toward the full run.** Specifically:
`e(b_es_ame_s#)` loses its structural-zero columns, `e(numpre_ame_s#)` goes 4 → 0,
`e(wdidsgresults_#)` loses its `ES_ame:pre*` and `pretreatment:ame` rows, and the printed
subgroup table loses those rows. Every one of those is a row the **analytic run already
drops**, so the change *removes* a documented-vs-actual divergence rather than creating one.
No analytic number moves anywhere. I judge this **not a re-bless-worthy number move**, but it
does change `fastnose` baselines if any exist for this option combination.

### Blast radius

`fastnose` × `contreat()` × `subgroups()` × ES window only. Downstream: `permute()`'s fast-mode
default borrows `fastnose` from the shell (`ado:~5274`), so a `permute()+contreat()+subgroups()`
run now ranks the same column set the observed run reports. BS-26's note says the binary side's
column-set divergence "abort[ed] bootstrap()/permute() permanalyticse with r(459) at rep 0"
before it was fixed; the continuous side plausibly still has that failure and this cures it.
**I did not test that** — worth one run before shipping.

### Alternatives rejected

* **Leave it and fix `sthlp:292` instead** ("fast mode reports the same rows as the full run,
  except under `contreat()`"). Rejected: the row set is the *output contract* for
  `honestdid`/`pretrends`, `e(numpre_ame_s#)` feeds Roth's `numpre`, and a post-treatment
  horizon reported as an average marginal effect of exactly 0 is a wrong number, not a doc gap.
* **Strip the omission marker instead of reading it.** Explicitly rejected by SI-7's own comment
  at `ado:2954-2960`: stripping the `o` leaves a token that parses as a plain number and the
  non-estimable value is then *reported as an estimate*.
* **Force `df(df_r)` on every margins call under `fastnose`.** That is the whole cost of fast
  mode; a non-starter.

### Test that proves it

The oracle-free invariant, which is also the help's promise: **for every family `f` and
subgroup `s`, `colnames(e(b_es_f_s))` and `e(numpre_f_s)` must be identical between a
`fastnose` run and the otherwise-identical full run.** `fn03_patch.do` asserts it on three arms
(`untreatedZero`, `basic`, no-subgroups control) and it now holds on all three. Add
`contreatcontrols(lattice)` and a `poisson` arm before shipping.

---

## SG-2 — F13-PM: a unit labelled on *some* of its rows is a split unit, and must be treated as one

**Covers:** fpsg-01, skp-01, skp-02. One root cause, one fix site.

### Why one edit covers all three

The two F13 screens in `wooldid_sgdom` both fold missing through `cond(sg >= ., ., sg)`, and
`gegen min/max` then **ignore** it. So a unit labelled `1` on some in-sample rows and missing on
the rest has `min == max == 1` and passes both. But `sg` is only the *margins subpop*: the
estimand weights `w_att_i`, `w_att_it`, `w_att_itime` carry denominators (`w_per_i`,
`w_per_it`, `t_per_i`, `ado:1619-1621`) built over the **whole estimation sample**, and
`margins` renormalises within the subpop. A unit whose labelled rows carry share `s` of its
post weight therefore enters its own subgroup with total weight `s`, not 1.

* **fpsg-01** is that at the unit level (`att_i`, `att_itime`).
* **skp-01** is that at the cell level (`att_it`, and `att` — a cell labelled on 1 of 5 rows
  counts as 0.2 cells). A partially-labelled *cell* is necessarily inside a partially-labelled
  *unit*, so the unit-level screen subsumes it.
* **skp-02** is the same thing in the *pre-treatment* weight bucket: `w_per_i` is grouped
  `by(posttreatment, i)` (`ado:1620`), so a label missing only on pre rows leaves `main:att_i`
  exact while share-weighting `pretreatment:att_i`, `pretreatment:att_itime` and the pre-trend
  ES rows — one unit at weight 1.00 and 0.25 in two rows of one printed table. The F13 domain
  `ifstatement & ttre != .` already spans pre **and** post rows of treated units, so the same
  screen catches it.

wooldid **already refuses** the structurally identical case — a unit whose `sg` value *changes*
over time — with `r(459)`, and its error text at `ado:12274` names precisely this distortion
("...weighted by its share of sample weight there, rather than with the equal unit weighting
those estimands normally imply"). The project's own standard for this consequence is a refusal.

### The edit

`wooldid.ado`, in `wooldid_sgdom`, **inserted between line 12278** (the `}` closing the
split-unit `if/else`) **and line 12279** (`qui drop `f13sgimin' `f13sgimax'`). It **reuses the
existing `f13sgimin`**, so it adds **no `gegen`**.

*(Anchor re-verified on the pinned build after the fact: `sed -n '12278,12279p'` gives `  }` then
`  qui drop \`f13sgimin' \`f13sgimax''`. An earlier draft of this document said 12279/12280, one
line late; the anchor text, not the number, is what to match.)*

```stata
  * F13-PM: PARTIAL LABELLING. Both screens above fold missing through cond(sg >= ., ., sg) and
  * gegen's min/max then IGNORE it, so a unit labelled on SOME of its in-sample rows and missing
  * on the rest passes both. It must not. sg selects the margins subpop, but w_att_i / w_att_it /
  * w_att_itime carry denominators (w_per_i, w_per_it, t_per_i) built at the weights block over
  * the WHOLE estimation sample, so margins' within-subpop renormalisation enters such a unit with
  * total weight equal to its LABELLED SHARE rather than 1 - the same distortion the split-unit
  * screen above refuses, reached by a different route. w_per_i is grouped by (posttreatment, i),
  * so a label missing only on pre-treatment rows leaves main:att_i exact while share-weighting
  * pretreatment:att_i / att_itime and the pre-trend ES rows, and one unit then carries two
  * different weights in one printed table. Reuses f13sgimin, which is missing only when the unit
  * has NO labelled in-domain row, so a WHOLLY unlabelled unit still passes and sthlp:887's
  * "observations with missing subgroup values belong to no subgroup" keeps the meaning it was
  * written for. No extra gegen: one assert pass, no sort, no temporary variable.
  * gated on replication mode exactly as the onset F13 screen is (ado:762, 772): a
  * replication that REBUILDS i - permonsetcohorts - can pool a labelled unit with a wholly
  * unlabelled one into one rebuilt i, which is a resampling artefact and not something the
  * user should be told about. The observed run's screen is the one that matters, and this
  * way the assert costs NOTHING per replication under reps(999).
  if "`rm'" != "1" {
  cap assert `subgroups' < . if `ifstatement' & `ttre' != . & `f13sgimin' < .
  if _rc != 0 {
    qui count if `ifstatement' & `ttre' != . & `f13sgimin' < . & `subgroups' >= .
    local f13pmn = r(N)
    if "`safetyinput'" != "off" {
      dis as error "Error: the subgroup variable `subgroups' is missing on `f13pmn' in-sample row(s) of treated units of `i' that carry a subgroup value on their other rows."
      dis as error "A partially labelled unit is not simply outside the subgroup: att_i, att_it and att_itime divide by denominators built on ALL of the unit's rows, so the unit enters its own subgroup weighted by the share of its sample weight that carries the label rather than with the equal unit (att_i, att_itime) or equal cell (att_it) weighting those estimands name. A label missing only on pre-treatment rows distorts the pretreatment and pre-trend rows while the main rows stay exact."
      dis as error "Blank the subgroup variable on ALL of a unit's rows to leave that unit out of every subgroup, fill the label in within the unit, define `i' at the finer level over which the label is observed (e.g., egen newi = group(`i' `subgroups'), missing), or restrict the sample with if."
      dis as error "To proceed anyway, specify safety(off): the affected rows are then share-weighted, exactly as a split unit is."
      error 459
    }
    else dis as text "Caution: the subgroup variable `subgroups' is missing on `f13pmn' in-sample row(s) of treated units of `i' that carry a subgroup value on their other rows; proceeding under safety(off) using share weighting. Such a unit enters its own subgroup weighted by its labelled share of sample weight rather than with the equal unit (att_i, att_itime) or equal cell (att_it) weighting those estimands normally imply, and a label missing only on pre-treatment rows distorts the pretreatment and pre-trend rows while the main rows stay exact."
  }
  }
```

**Late correction to the remedy text — and to `sthlp:887`'s existing one.** The remedy sentence
already in the help and in the F13 cell-guard error is `egen newi = group(i sg)`. **Bare
`egen ... group()` returns missing whenever any grouping variable is missing**, so on exactly
the data this screen fires on it would set `newi = .` on every unlabelled row, and wooldid's
`markout` would then drop those rows from the estimation sample — silently changing the pooled
`att`, the never-treated control mass and `e(N)`. The `, missing` option (used above) keeps them
as their own group, which is the behaviour the sentence intends. This is a **pre-existing** flaw
in the help's remedy, independent of the patch: `sthlp:887`'s parenthetical
"(e.g. `egen newi = group(i sg)`)" should gain `, missing` whether or not SG-2 ships. I did not
run this; it follows from `egen group()`'s documented semantics, and it is cheap for the author
to confirm.

The gate needs `wdrepmode` in `wooldid_sgdom`, which is a separate program, so two more tokens:
`syntax , ... [si(string) RM(string) MAKECfxplotsbysg]` at `ado:12202`, and
`... si(`safetyinput') rm(`wdrepmode') ...` at the one call site, `ado:995`.

**Help edit (required, whatever the author rules).** `sthlp:887` currently ends
"Observations with missing subgroup values remain allowed and belong to no subgroup." — which
reads as benign and is exactly what misleads. Replace with:

> Observations with missing subgroup values belong to no subgroup and are allowed **when the
> whole unit is unlabelled**. A unit labelled on some of its in-sample rows and missing on the
> rest is refused under default safety for the same reason as a unit that changes subgroup:
> `att_i` and `att_itime` would enter it weighted by its labelled share of sample weight, and
> `att_it` and `att` would count a partly-labelled cell as a fraction of a cell.

### Measured before / after

`A_r06_attit_base.do` / `A_r06_attit_patch.do` — the finding's own noiseless design
(`e(r2) = 1.0000`, so every target is exact), plus three arms I added.

| arm | pinned | patched |
|---|---|---|
| `sgA` clean subgroup | runs; `att 6.7058823529 / att_it 6.0 / att_i 6.0 / att_itime 6.0` | **numerically IDENTICAL** |
| `sgB` partially-labelled cells | **runs silently**; `att 8.1538461538 / att_it 7.0 / att_i 7.6 / att_itime 7.0` (documented targets 6.0) | **`r(459)`**, "missing on **8** in-sample row(s)" |
| `sgB, safety(off)` | same numbers | runs with the caution; numbers **bit-identical** to the pinned default run |
| `sgC` — a **whole** unit unlabelled (the case `sthlp:887` promises) | runs | **runs, numbers IDENTICAL** |
| `sgD` — never-treated units partly unlabelled | runs | **runs** (domain is `ttre != .`) |

`A_p10_partialmissing_patch.log` (fpsg-01's repro): refused, names **5** rows.
`A_r02_mechanism_patch.log` (skp-02's repro): refused on both parts, names **6** and **7** rows.

### Cost

**One `assert` pass, O(N), no sort, no temporary variable, no `gegen`.** The block it joins
already runs four `gegen ... by()` calls over the same rows; a linear `assert` is a fraction of
one of them, and the whole block is a fraction of one `reghdfe` fit. The `qui count` that
builds the row count for the message runs **only on the failure branch**, where the run is
about to stop or caution anyway.

**Nothing is paid per replication.** `wooldid_sgdom` is called from `ado:995` on every nested
`wooldid` call, replications included, so the gate on `wdrepmode` is not cosmetic: without it a
`reps(999)` bootstrap would pay 999 extra linear passes for a screen whose answer cannot change
(resampling duplicates whole clusters and cannot manufacture partial labelling). With it, the
screen runs exactly once, on the observed data. That also removes a false-alarm route:
`permonsetcohorts` rebuilds `i` inside each replication, and a rebuilt `i` can pool a labelled
unit with a wholly unlabelled one — which would have aborted the run mid-replication.

If even one extra pass is unwanted on very large data, the success path can be held at exactly
today's cost by folding the two tests into one assert and re-testing separately on failure:

```stata
  cap assert (`f13sgimin' == `f13sgimax') & (`subgroups' < . | `f13sgimin' >= .) if `ifstatement' & `ttre' != .
```

I did **not** test that variant.

### Moves a reported number?

**No.** Every run that still passes reports bit-identical numbers (measured on `sgA`, `sgC`,
`sgD`). What changes is that a class of runs that used to print wrong numbers now stops. That
is a **behaviour change users will hit**, and `safety(off)` reproduces the old numbers exactly.

### Blast radius

`subgroups()` only, and only on data where a treated unit is partly labelled. Never-treated
units, wholly-unlabelled units and clean data are untouched. `wooldid_sgdom` is called once per
run from `ado:995`, before estimation, so nothing downstream sees a changed object.

**One over-strictness corner, stated so the author can price it.** The screen's domain is
`ifstatement & ttre != .` — every in-sample row of a treated unit, *including the omitted
reference-period rows*. Those rows sit in their own `posttreatment` bucket (`ado:1595-1604`
assigns `0` to pre-onset rows of treated units that are *not being reported* — the reference
period, and any pre period outside the requested window), so they enter **no** reported
estimand's denominator, and a unit whose label is missing *only* on rows with
`posttreatment == 0` suffers **no** distortion
yet would still be refused. Narrowing the domain is not free: `wooldid_sgdom`'s `syntax` at
`ado:12202` receives `tvar` (calendar `t`) but **not** `_t`, and `posttreatment` does not exist
yet — it is built at `ado:1599`, hundreds of lines after the `ado:995` call site — so the
refinement means threading the relative-time variable into a program that was split out
precisely for compiled-size headroom, and under the default pooled cohort-varying reference
there is no single omitted period to test against anyway. **I recommend leaving it.** The corner
is narrow (a label missing on exactly the reference rows and nowhere else), the failure mode is
a refusal with an explanatory message rather than a wrong number, and the remedy the message
already gives — blank the label on all of the unit's rows, or fill it in — costs the user one
line. It should be named in the ruling, not engineered around.

### Alternatives rejected

* **Fix the weights instead** — rebuild `w_per_i` / `w_per_it` / `t_per_i` restricted to the
  labelled rows so the subgroup estimands are equal-unit again. One extra `gegen ... by(post i
  sg)` would serve every subgroup at once, so the *cost* is not the objection. Three things
  are: (i) the pooled `att_i` must keep the all-rows denominator, so it needs a **second**
  weight variable per estimand and a rewrite of ~20 `margins` call sites to choose between
  them; (ii) it would silently change the **documented** `safety(off)` split-unit behaviour,
  which the author ruled in 2026-07-28; (iii) it presupposes a semantic answer the author has
  not given — a unit labelled on 1 of 6 post periods would then enter the subgroup with weight
  1 carrying **only that one cell's** effect, which is not obviously the intended "unit effect"
  either (fpsg-01's noisy run D measured exactly this gap: 0.9545 observed vs 0.9486
  share-weighted, the residual being precisely this).
* **Caution rather than error.** Defensible, and it is the author's call (see Open questions).
  I propose error-by-default because the split-unit case with the *identical* consequence is an
  error, and because the distortion is silent and large (measured +16.7% on `att_it`, +26.7% on
  `att_i`, 33% on the pre-treatment row, 45% on noisy data).
* **Screen at the cell level only.** Would miss skp-02 entirely (a unit can be cleanly labelled
  within every cell and still be unlabelled for whole pre-treatment periods).
* **Documentation only.** The number under a printed row labelled `att_i` is not the estimand
  that label names; `sthlp:887` promising that missing values are benign is part of the defect,
  not the whole of it.

### Test that proves it

No oracle needed — two invariants.

1. **Refusal symmetry.** For any dataset, blanking `sg` on a subset of a treated unit's rows
   and *changing* `sg` on the same rows must produce the same class of outcome (both `r(459)`,
   or both a `safety(off)` caution). Today the first is silent and the second is `r(459)`.
2. **Equal-unit identity on admitted runs.** On any run the screen admits,
   `att_i` for subgroup `s` must equal the unweighted mean over the units of `s` of their own
   `e(histogramestimates)` effects. That holds on clean data and fails on partially-labelled
   data — which is what the screen now excludes.

Plus the regression that matters: **a wholly-unlabelled unit must still run and give identical
numbers** (`sgC` arm above), because that is what `sthlp:887` promises and what users actually
do to drop a unit from the subgroup comparison.

---

## SG-3 — put the printed standard errors on the diagonal of the bootstrap `e(V)`

**Covers:** SK-03, hv-bootrank-03, skep-sgexport-01, skep-threese-02. One fix site
(`wooldid_bootfin`), plus its twin `wooldid_bootves`, a caution and a footer line.

### The problem in one sentence

`wooldid_bootfin` runs **two moment passes on different row sets** — the printed `se[k]` over
the replications where column `k` alone is available, `V[k,k]` over the replications complete
across **every** reported column — and ships both under names that imply they agree. The code's
own comment at `ado:14500-14506` says exactly this and cites A11's ruling that it is "a
documented consequence rather than a defect", on the wording that a hand test "does not EXACTLY
reproduce the table's Wald numbers". **Measured, it is not an exactness tail**: the four
findings put it at 4%, 9.7%, 19.8%, 30%, 36% and 50.4% on individual coefficients, and
`e(boot_N_V)` collapses to 14 of 60 with **zero cautions printed**. Per the coordinator note,
contradicting a ruled entry by measurement is the finding.

### The edit — `V = D R D`

**Keep the listwise pass. Use it for the CORRELATIONS only. Take the SCALE from the per-column
pass, which is exactly the `se` column the table prints.** Because `D = diag(se)` is a positive
diagonal, `D R D` is a congruence of the listwise covariance and is **positive semidefinite
whenever the listwise matrix is** — which is the whole reason listwise was chosen over
pairwise-complete (`ado:9640-9643`: "A11: PSD by construction; pairwise-complete is measurably
indefinite"). This fix therefore does **not** reopen that objection.

**(a) `wooldid.ado:14540-14559` (Mata, `wooldid_bootfin`)** — inserted after the listwise
block, before the sentinel loop:

```mata
  // A11-R: rescale the listwise covariance to the PRINTED standard errors. The listwise
  // pass supplies the CORRELATIONS (over the rows complete across every non-excluded
  // column); the per-column pass above supplies the SCALE, and that scale is exactly the
  // se column the table prints. V = D R D with D = diag(se) is a congruence by a positive
  // diagonal, so it preserves the positive-semidefiniteness that made listwise the choice
  // over pairwise-complete in the first place, while making sqrt(diag(e(V))) equal the
  // printed se bit for bit - which is what sthlp:554 promises under bootdisplay(wald),
  // and what esttab/coefplot/lincom/test/honestdid silently export. sdlw is captured
  // BEFORE the rescale because the sentinel test below asks whether the LISTWISE pass
  // found a variance, not whether the per-column pass did.
  real colvector sdlw
  real rowvector scl
  sdlw = sqrt(diagonal(V))
  scl = J(1, K, 0)
  for (k = 1; k <= K; k++) {
    if (sdlw[k] > 0 & se[k] < . & se[k] > 0) scl[k] = se[k] / sdlw[k]
  }
  if (K > 0) V = V :* (scl' * scl)
```

and the sentinel loop's predicate moves from `V[k,k] == 0` to `sdlw[k] == 0` (same question,
asked of the pre-rescale matrix), with the diagonal written from `se` directly so the equality
is **exact and not merely to rounding**:

```mata
-   if (excl[k] == 1 | cnt[k] < 2 | V[k,k] == 0) {
+   if (excl[k] == 1 | cnt[k] < 2 | sdlw[k] == 0) {
      V[k, .] = J(1, K, 0)
      V[., k] = J(K, 1, 0)
      V[k, k] = -1e100
    }
+   else V[k, k] = se[k] * se[k]
```

**(b) `wooldid_bootves` (`ado:14582`)** — the same rule, so `e(V_es_*)` (the honestdid entry
point at `ado:10068`, and what `pretrendpower()` reads under bootstrap) agrees with the table
too. New `sename` parameter; the family's own listwise pass still supplies the correlations.

```mata
- void wooldid_bootves(string scalar accname, ..., string scalar newname)
+ void wooldid_bootves(string scalar accname, ..., string scalar newname, string scalar sename)
  ...
  V = (dev' * dev) / (rows(CC) - 1)
+ sefull = st_matrix(sename)
+ sdlw = sqrt(diagonal(V))
+ scl = J(1, m, 0)
+ for (k = 1; k <= m; k++) {
+   if (sdlw[k] > 0 & sefull[1, idx[k]] < . & sefull[1, idx[k]] > 0) scl[k] = sefull[1, idx[k]] / sdlw[k]
+ }
+ V = V :* (scl' * scl)
+ for (k = 1; k <= m; k++) {
+   if (sdlw[k] > 0 & sefull[1, idx[k]] < . & sefull[1, idx[k]] > 0) V[k, k] = sefull[1, idx[k]]^2
+ }
```

**(c) `ado:10080`** — the caller passes the se vector (`bse` is already a tempname in scope in
`wooldid_bootwire`, declared at `ado:9635`):

```stata
-   mata: wooldid_bootves("`acc'", "`excl'", "`fk'", "`evold'", "`evnew'")
+   mata: wooldid_bootves("`acc'", "`excl'", "`fk'", "`evold'", "`evnew'", "`bse'")
```

**(d) `ado:10296` — a caution that actually fires.** hv-bootrank-03's sharpest point is that at
`boot_N_V = 14 of 60` the run printed **zero** cautions, because the existing `wdb_lowav`
caution is **per column** and every column was above 50%. The listwise count collapses
*multiplicatively* across columns, so it needs its own test:

```stata
  if `wdb_nlw' < . & `done' > 0 & `wdb_nlw' < 0.5*`done' dis as text "Caution: only `wdb_nlw' of the `done' completed replications are complete across every reported coefficient at once, so the CORRELATIONS in e(V) rest on that many draws. The standard errors - in the table and on the diagonal of e(V) - use each coefficient's own available replications (e(boot_navail)) and are unaffected. Adding reported coefficients, through subgroups() or a longer event-study window, shrinks this count multiplicatively."
```

**(e) `ado:11797` — the footer stops mis-describing `e(V)`:**

```stata
- "  e(V) holds the bootstrap covariance (complete-replication listwise, e(boot_N_V) replications; ...)"
+ "  e(V) holds the bootstrap covariance: its DIAGONAL is the printed standard errors squared (each coefficient's own available replications, e(boot_navail)); its correlations are complete-replication listwise over e(boot_N_V) replications. The analytic covariance, when one existed, remains in the hidden e(V_analytic)."
```

**(f) footer wording — skep-threese-02's secondary point.** The bootstrap footer groups
"esttab, coefplot and honestdid" as "the tools that read `e(V)`", but the help's own honestdid
recipe (`sthlp:1560-1567`) is `matrix V = e(V_es_att)`. After (b) the two matrices agree on the
diagonal, so the sentence is no longer *wrong*; it is still imprecise. Suggest: "...esttab and
coefplot read `e(V)`; `honestdid` should be fed `e(b_es_*)`/`e(V_es_*)` as the recipe under
`e(V_es_)` shows."

### What each finding gets

| finding | before | after |
|---|---|---|
| **SK-03** — `lincom` ≠ printed row under `bootdisplay(wald)` | `WORST_PCT_V_VS_PRINTED = 19.84125` | **`WORST_PCT_V_VS_PRINTED = 0.00000`** |
| **hv-bootrank-03** — `esttab` exports 50.4% of the printed `se` | `sqrt(diag(e(V)))` is the listwise SD | `sqrt(diag(e(V)))` **is** the printed se (exactly on the columns the rescale touches; to the last ulp on the rest — see the inertness guard) |

A note on hv-bootrank-03's magnitude, because its own refuter corrected it. The prober's "a
factor of two" is that run's `e(V)` diagonal; the refuter measured the **systematic** listwise
selection bias at roughly -6% to -15% in SE and attributed the rest of any single run's gap to
14-draw noise. **The fix does not depend on which reading is right.** It makes the diagonal the
per-column standard error whatever the split between bias and small-`N_V` noise, because the
defect it repairs is that two numbers shipped under names implying they agree do not — not that
the listwise one is biased by any particular amount. Read the 50.4% as the size of the gap a
user can meet on one run, not as the characterisation of the bias.
| **skep-sgexport-01** — adding `subgroups()` shrinks exported main-block SEs by up to 36% | diagonal depends on which other columns exist | diagonal is the per-column se, **invariant** to added columns |
| **skep-threese-02** — one coefficient, three SEs | printed / `e(V)` / `e(V_es_att)` disagree by 56% | all three agree on the diagonal |

### Measured before / after (SG-3)

`C_sk03_base.do` / `C_sk03_patch.do` — SK-03's own repro (24 states, 2-4 units each, 9 years,
`subgroups(sg3)`, `cluster(state)`, `reps(399)`, `bootdisplay(wald)`), adopath retargeted, with
a PSD check and a `V_es` dump appended. `e(boot_completed) = 399`, `e(boot_N_V) = 368`,
`K = 24` — identical in both runs, so nothing about the replications changed.

| | pinned | patched |
|---|---|---|
| `WORST_PCT_V_VS_PRINTED` over all 24 coefficients | **19.84125** | **0.00000** |
| `pretreatment_SG0_att`: printed se / `sqrt(V_kk)` | 0.269231486 / **0.215812605** | 0.269231486 / **0.269231486** |
| `ES_att_SG0_pre1`: printed se / `sqrt(V_kk)` | 0.324726026 / **0.293179186** | 0.324726026 / **0.324726026** |
| `main_att`: printed se / `sqrt(V_kk)` | 0.150555484 / **0.151332945** (0.52% off) | 0.150555484 / **0.150555484** (0.00000%) |
| `ES_att_post1`: printed se / `sqrt(V_kk)` | 0.161083438 / **0.163316021** (1.39% off) | 0.161083438 / **0.161083438** (0.00000%) |
| `lincom _b[pretreatment_SG0_att]` | se **.2158126**, t 0.81, p **0.427**, CI [-.2720803, .6208045] | se **.2692315**, t 0.65, p **0.524**, CI [-.3825857, .7313099] |
| `min(symeigenvalues)` of the non-sentinel block (24 of 24) | **+3.296e-05** | **+3.371e-05** |
| `max(symeigenvalues)` | .30782853 | .32824899 |

Two things to read off that table.

* **The promise is met exactly, not approximately.** `WORST_PCT_V_VS_PRINTED` is the worst
  relative gap over every coefficient, and it is `0.00000`. `lincom` now returns the printed
  row — and note the p-value moved across 0.5 in the process: this is the size of thing a user
  publishes.
* **PSD survives.** The minimum eigenvalue is positive before and after, and slightly *larger*
  after. That is the point of `D R D` over pairwise-complete.

`sqrt(diag(e(V_es_att)))` is **identical** in the two runs (0.141479752, 0.158863641,
0.176620916, 0.161083438) — on this design the `att` event-study family's own listwise set
already equalled each column's available set, so the inertness guard correctly left that matrix
alone. That is the guard working, not the fix failing: the coefficients that moved are exactly
the subgroup ones, where the two passes really did disagree.

### The inertness guard (important for re-blessing)

Without care, the rescale would move `e(V)` by a last-ulp amount on **every** bootstrap run,
including runs where nothing is missing, because `variance(A[use,k])` and
`diagonal(dev'dev)/(n-1)` are different floating-point paths to the same number. A
bit-identity baseline would light up everywhere. So the patch carries a guard:

```mata
  scl = J(1, K, 1)
  for (k = 1; k <= K; k++) {
    if (excl[k] == 0 & cnt[k] != nlw & sdlw[k] > 0 & se[k] < . & se[k] > 0) scl[k] = se[k] / sdlw[k]
  }
  if (K > 0 & anyof(scl :!= 1, 1)) V = V :* (scl' * scl)
```

`cnt[k] == nlw` means column `k` was available in exactly the listwise rows (those rows are a
subset of `k`'s available rows, so equal counts make the sets equal), and the two passes then
agree to the last ulp. Such a column is left **strictly alone**, and `V[k,k] = se[k]^2` is
written only where `scl[k] != 1`. **Consequence: on any bootstrap run in which every
replication delivered every reported coefficient, this patch is byte-for-byte inert.** It moves
`e(V)` only where the two passes genuinely disagree — which is exactly the case the four
findings are about. `wooldid_bootves` carries the same guard.

**MEASURED, not asserted (added on a later pass — this was the one load-bearing claim in this
document that had only been reasoned about).** `D_inert2_base.do` / `D_inert2_patch.do`:
60 units in 5 cohorts of 12, 10 years, 4 staggered onsets + a never-treated cohort,
`att att_it att_i att_itime esprelength(2) espostlength(2) cluster(unit) fe(unit)
bootstrap(standard) reps(30) rseed(777)`. `i` is the **cohort** and the cluster is the **unit**,
12 per cohort, so `P(a cohort absent from a replication) = (1-12/60)^60 ~ 1.5e-6` and no column
can blank. Result, identical in both runs: `boot_completed 30`, `boot_N_V 30`,
`min(e(boot_navail)) = max(e(boot_navail)) = 30`, i.e. `ALL_COLS_COMPLETE=1`. The **full 16x16
`e(V)`, off-diagonals included, dumped in `%21x`, is BYTE-IDENTICAL between the pinned build and
the patched build** (`D_inert2_base_V.txt` vs `D_inert2_patch_V.txt`, `diff` empty, 89 lines
each; sentinel `INERTNESS_16x16_BYTE_IDENTICAL`). The 1x1 warm-up (`D_inert_base/patch.do`,
`att` alone) likewise gives `main_att +1.2c42dd95ddf53X-008` on both sides.

**So: a bootstrap run in which every replication delivers every reported coefficient does not
need a re-bless.** Only runs where the two moment passes genuinely disagree move, which is the
class the four findings are about.

**The trade-off the guard buys, stated plainly, because it is an author choice.** With the
guard, `sqrt(diag(e(V)))` equals the printed `se` **exactly** on every column the rescale
touches, and **to the last ulp** on the columns it deliberately does not (where
`variance(A[use,k])` and `diagonal(dev'dev)/(n-1)` are two floating-point routes to the same
number). Dropping the guard — always writing `V[k,k] = se[k]^2` — would make the equality
bit-exact everywhere, at the price of moving `e(V)` by an ulp on **every** bootstrap run and
lighting up any bit-identity baseline. I chose inertness. One-line change either way.

### Cost

Free. The rescale is `K` divisions plus one `K x K` elementwise multiply, against a `K x K`
matrix multiply (`dev' * dev`) that the function already performs. `wooldid_bootves` adds one
`colsum` per family row. Nothing scales with `N` or with `reps()`. Nothing new is allocated
per replication — this all runs **once**, after the replication loop. On large data the cost is
indistinguishable from zero.

### Moves a reported number?

**Yes — this is the re-bless item, and the author must rule on it.**

* **Printed tables: unchanged.** Column 2 of `e(wdidmainresults)` / `e(wdidsgresults_*)` is the
  per-column SD and is not touched.
* **`e(V)`, `e(V_boot)`, `e(V_es_*)`: changed** on any bootstrap run with a missing replication
  cell. Magnitudes from the findings: the headline `att` SE goes from 50.4% of the printed one
  to 100% of it; `pretreatment_SG0_att` from 19.84% too small to exact; ES rows from 30-36% too
  small to exact.
* **Downstream consumers whose *printed* numbers therefore move:** `pretrendpower()` reads
  `e(V_es_*)` under bootstrap (`ado:4497`, `10381-10445`); `simultaneousbands` reads
  `e(V_es_*)` (`ado:10938-10943`, `11051-11063`). Both will print different numbers on affected
  runs. `wooldid_bootjwr`'s per-replication joint tests are built inside each replication from
  that replication's own analytic `e(V_es_*)` and are **not** affected.
* **User-side:** `esttab`, `coefplot`, `test`, `lincom`, `honestdid` all change — in the
  direction of finally matching the printed table.
* Two **new** output lines: the `wdb_nlw` caution (d) and the reworded footer (e).

### Blast radius

Every `bootstrap()` run, but only through `e(V)`-family matrices, and only where the two passes
disagree. The `-1e100` sentinel convention (PS-62) is preserved exactly: the predicate moves
from `V[k,k] == 0` to `sdlw[k] == 0`, which is the *same question asked of the same matrix*,
just before the rescale instead of after. `ereturn repost` still sees a matrix that may be
indefinite because of the sentinels — unchanged, and the comment at `ado:10052-10057` still
holds.

One behaviour worth the author's eye: today a column whose **listwise** variance is zero gets
the sentinel even though the table shows a finite `se` (the code comment at `ado:14554-14559`
calls this out). The patch keeps that. It could instead now be given `V[k,k] = se[k]^2` with
zero off-diagonals, which would be strictly better for `lincom`. I left it alone to keep the
patch to one idea.

### Alternatives rejected

* **Pairwise-complete covariance.** The natural estimator, and it would also make
  `skep-sgexport-01`'s *off-diagonals* invariant. Rejected because `ado:9640-9643` records that
  it was tried and is "measurably indefinite", and because a `D R D` rescale gets the standard
  errors exactly right without giving that up. If the author wants the off-diagonals on more
  replications too, the right follow-up is pairwise **correlations** (which satisfy `|r| <= 1`
  pairwise, so every 2x2 submatrix stays PSD and `lincom a - b` is always safe) with the same
  `D` scaling — a strictly larger change that I have not tested.
* **Disclosure only** — keep listwise, print the ratio, fix the footer. This is the honest
  fallback if the author declines the number move, and it is (d)+(e)+(f) above without (a)-(c).
  It leaves `sthlp:554`'s promise ("a `test` or `lincom` run by hand afterwards reproduces the
  reported p-values ... **this holds under** `bootdisplay(wald)`") false on its own terms, so I
  do not recommend it alone.
* **Make the printed `se` the listwise one instead** (move the table to `e(V)` rather than
  `e(V)` to the table). Rejected: that throws away replications for no gain — the per-column SD
  uses strictly more information — and it would move every printed bootstrap table in the
  package.
* **Post `diag(se^2)` as `e(V)`** (the analytic convention, PS-62). Rejected: it discards the
  genuine event-study off-diagonals that PS-62 deliberately splices in, and `honestdid` needs
  them.
* **Extract `e(V_es_*)` as a subblock of the new `e(V)`** instead of rescaling inside
  `wooldid_bootves`. Attractive — it makes `sthlp:1318`'s "the same blocks are spliced into
  `e(V)`" literally true under bootstrap, and there would be exactly one construction. Rejected
  as the default because the family's own listwise set is larger than the global one, so the
  family's correlations would get *worse*. Worth putting to the author.

### Test that proves it

1. **The identity the help promises, no oracle needed.** For every coefficient `k` on any
   bootstrap run: `sqrt(e(V)[k,k])` **==** column 2 of `e(wdidmainresults)`/`e(wdidsgresults_*)`
   for that row, exactly. `C_sk03_patch.do` computes `WORST_PCT_V_VS_PRINTED` over all 24
   coefficients; it must be 0.
2. **`lincom` reproduces the printed row.** Under `bootdisplay(wald)` the printed interval is
   `b +/- crit*se`, so `lincom _b[<name>]` must return the printed `se`, `t`, `p` and interval.
3. **Cross-matrix identity.** `sqrt(diag(e(V_es_att)))` == the printed `ES_att:*` standard
   errors == `sqrt(diag(e(V)))` on those rows.
4. **The matched-pair invariant (skep-sgexport-01), which needs no truth at all.** Two runs on
   the same data with the same `rseed()`, differing only in `subgroups()`. An option that adds
   subgroup rows is not a statement about the main event study, so `sqrt(diag(e(V)))` on the
   **shared** coefficient names must be identical. `C_pair_base.do` / `C_pair_patch.do`.
5. **PSD.** `min(symeigenvalues(S))` on the non-sentinel submatrix must be `>= -1e-12 * max`.
   Computed in `C_sk03_*.do` as `PSD_MINEV` / `PSD_RATIO`.
6. **Inertness.** A bootstrap run with no missing replication cell must produce a byte-identical
   `e(V)` before and after.

---

## Things I do NOT propose changing

* **The split-unit `safety(off)` share-weighting rule itself.** It is an author ruling
  (2026-07-28) with a caution that names the consequence. SG-2 extends its *domain*; it does not
  reopen the rule.
* **`e(V)`'s cross-family zero off-diagonals (PS-62).** Untouched by SG-3.
* **The `-1e100` sentinel convention.** Untouched; SG-3 asks the same question of the same
  matrix, one step earlier.
* **The per-column standard error itself.** It uses strictly more replications than the listwise
  one. It is the right number; the problem is that `e(V)` does not carry it.
* **`e(boot_N_V)`'s name.** After SG-3 it still means what it says — the listwise count — it just
  governs the correlations rather than the whole matrix. The footer (e) now says so.

## Open questions for the author

1. **SG-2: error or caution by default?** I propose `r(459)`, because the split-unit case with
   the identical consequence is `r(459)` and because the distortion is silent and large
   (+16.7% `att_it`, +26.7% `att_i`, 33% on the pre row, 45% on noisy data). A caution-only
   ruling is defensible and is a one-word change to the patch. There is a real class of users
   who will hit this — an industry code that starts partway through the panel, item nonresponse
   on a demographic — and for them the remedy is a one-line `egen`, but the run stops today
   rather than tomorrow.
2. **SG-3: accept the number move?** This contradicts A11 by measurement, which is exactly what
   the coordinator note says is a first-class finding, but A11 is a ruling and only the author
   can reverse it. If declined, ship (d)+(e)+(f) — the disclosure half — and amend `sthlp:554`,
   which currently promises the opposite of what the code does.
3. **SG-3(b): rescale inside `wooldid_bootves`, or make `e(V_es_*)` a literal subblock of the
   new `e(V)`?** The subblock makes `sthlp:1318` true under bootstrap and leaves one
   construction in the package; the rescale keeps the family's larger listwise row set for the
   correlations. I implemented the rescale.
4. **Should the listwise-zero-variance column (`sdlw[k] == 0`, finite printed `se`) now get
   `V[k,k] = se[k]^2` instead of the sentinel?** Strictly better for `lincom`; a separate idea,
   so I left it out.
5. **SG-1: is `permute() + contreat() + subgroups() + ES` currently aborting `r(459)` at rep 0**
   the way BS-26 says the binary side did before its fix? If so this patch cures that too and
   the finding is larger than reported. Not tested.

6. **A stale comment found while checking SG-3's premise, unrelated to any of the eight
   findings.** `ado:9654-9658` says "the displayed standard error under [the percentile] package
   is the interval inverted, `(ub-lb)/(2*zc)`". That was overruled on 2026-09-04
   (`ado:9747-9754`), and `ado:9707` writes `bse` into column 2 with no branch on
   `bootdisplay()`. The comment survives and would mislead the next reader of this code — it
   very nearly misled me into scoping SG-3 to `bootdisplay(wald)` only. Worth deleting while
   the file is open.

## What I could not test

* SG-1 under `contreatcontrols(lattice)`, `poisson`, or inside `permute()`/`bootstrap()`.
* SG-2's one-pass variant (the folded assert).
* ~~SG-3's inertness guard on a real no-missing-cell bootstrap run.~~ **Now measured** — see
  "The inertness guard" above: a 16x16 `e(V)` with `boot_N_V == boot_completed == 30` is
  byte-identical in `%21x` before and after (`D_inert2_*`). This was the last asserted-not-
  measured claim in the document.
* Anything at the scale the author actually cares about for cost (millions of rows). My cost
  claims for SG-2 and SG-3 are structural (one linear pass; one `K x K` elementwise multiply
  after the replication loop); only SG-1's is measured, and that one is the one that is not free.
---

## Artifacts

All under `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-SG/`:

| file | what |
|---|---|
| `wooldid_patched_SG.ado` | the pinned build with all three proposals applied |
| `prop-SG.patch` | `diff` of the pinned build against it (14 hunks) |
| `fn03_base.do` / `.log`, `fn03_patch.do` / `.log` | SG-1 before/after (the finding's own repro, adopath retargeted) |
| `fn03_patch_v1.log` | SG-1 with hunk (a) only — shows the SE leak that made hunk (b) necessary |
| `fn03_cost.do` / `.log`, `fn03_cost_patch.do` / `.log` | SG-1 cost, 60,000 rows |
| `A_r06_attit_base.do` / `.log`, `A_r06_attit_patch.do` / `.log` | SG-2 before/after + `safety(off)` + the two benign-case regressions |
| `A_p10_partialmissing_patch.log`, `A_r02_mechanism_patch.log` | SG-2 against fpsg-01's and skp-02's own repros |
| `C_sk03_base.do` / `.log`, `C_sk03_patch.do` / `.log` | SG-3 before/after, plus PSD and `V_es` checks |
| `C_pair_base.do` / `.log`, `C_pair_patch.do` / `.log` | SG-3 matched-pair invariant (with/without `subgroups()`) |
| `fn03_patch_stage1.log` | SG-1 verified again on the FINAL combined ado (all three fixes applied): byte-identical to the stage-1 SG-1-only run, so the three patches do not interact |
| `A_r06_attit_patch_pre_gatemove.log` | SG-2 before the `wdrepmode` gate was moved to wrap the assert: numerically identical, so the move is behaviour-preserving |
| `D_inert_base.do` / `.log`, `D_inert_patch.do` / `.log` | SG-3 inertness, 1x1 warm-up: `e(V)` `%21x` identical |
| `D_inert2_base.do` / `.log`, `D_inert2_patch.do` / `.log`, `D_inert2_base_V.txt`, `D_inert2_patch_V.txt` | SG-3 inertness, **16x16 `e(V)` byte-identical** on a run with `boot_N_V == boot_completed == 30` |
| `cmp_pair.py` | comparison script for the matched-pair logs |
| `NOTES.md` | running notes |

Reproduce any of them with
`"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do <file>.do` from that directory.
---

## Comment and help upkeep the patches make necessary

A patch that leaves the code contradicting its own comments is not finished. These all become
false and must be rewritten with the change:

| site | what it now says | why it is stale after the patch |
|---|---|---|
| `ado:1974-1979`, `ado:2935-2937` | "the subgroup-by-horizon margins call keeps its VCE ... this ONE call" | SG-1 makes it two calls |
| `sthlp:292` | "the one exception is the subgroup-by-horizon margins call ... at a small cost per subgroup" | same; and the measured cost is +18% of a `fastnose` run on my geometry, not obviously "small" |
| `sthlp:887` | "Observations with missing subgroup values remain allowed and belong to no subgroup" | SG-2 restricts that to *wholly* unlabelled units (verbatim replacement given above) |
| `ado:14500-14506` (`wooldid_bootfin` header) | "the two disagree, so `sqrt(diag(e(V_boot)))` is not the printed se column and a hand-run `test` or `lincom` does not exactly reproduce the table's Wald numbers. A11 ruled this a documented consequence" | SG-3 makes them agree; the A11 citation needs replacing with the new rule |
| `ado:9645-9651` (caller) | "CONSEQUENCE of the two different denominators ... `sqrt(diag(e(V_boot)))` does NOT equal the printed se column" | same |
| `ado:14554-14559` (sentinel) | "the branch also fires on `V[k,k] == 0`, i.e. when the LISTWISE pass found no variance ... a column can show a standard error in the table and carry the sentinel in `e(V)`" | the predicate is now `sdlw[k] == 0`; the behaviour is kept but the text must name the new variable |
| `ado:9654-9658` | "the displayed standard error under [the percentile] package is the interval inverted, `(ub-lb)/(2*zc)`" | **already stale on the pinned build**, independently of my patch: `ado:9747-9754` records the 2026-09-04 ruling that overruled it, and `ado:9707` writes `bse` into column 2 unconditionally. Worth fixing while in the file — it is the one comment that would make a reader doubt SG-3's premise |

That last row matters for SG-3's headline claim, so I checked it: **the printed standard error
is `bse` (the per-column bootstrap SD) under every non-studentized display**, written at
`ado:9707` with no branch on `bootdisplay()`. `sthlp:554`/`556` are correct and the `ado:9654`
comment is not. So "the diagonal of `e(V)` becomes the printed standard error squared" holds
under `percentile` and `bc` as well as `wald`, not only in the display I tested.

## Two blast-radius notes the measurements do not show

* **SG-1 inside `permute()`.** `permute()`'s fast-mode default borrows `fastnose` from the
  shell, and replications are nested `fastnose` calls. The added VCE is therefore paid **per
  replication**, not once — so the +18.3% I measured is a per-replication figure under
  `permute() + contreat() + subgroups() + ES`. The binary side already pays it, so this is
  parity rather than a new class of cost, but it should be priced on a real `permreps()` before
  shipping.
* **SG-3 under `boottstat`.** The studentized path carries its own `seacc` accumulator and its
  own `wooldid_boottp` package. Column 2 of the results matrices is still written from `bse` at
  `ado:9707`, so I expect the fix to behave identically there, but I did not run a `boottstat`
  arm.
