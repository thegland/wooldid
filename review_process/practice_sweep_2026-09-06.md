# Practice sweep: where wooldid is likely to look unusual to an applied microeconomist

2026-09-06. Branch `bugfixes`. Package files read at `wooldidinstall/wooldid.ado` (13,579 lines) and
`wooldidinstall/wooldid.sthlp` (1,375 lines).

## 0. Purpose and method

The author asked a single question: where is wooldid's implementation likely to strike applied
microeconomists as unusual and unlike what they see in other commonly used software or in published
papers — what is displayed, what is offered, what is missing, and which defaults were chosen. This
report is the answer as a catalogue, not as a verdict: every item states what wooldid does, what the
neighbouring package or paper does, why a reader might stumble, and what harmonizing would cost, and
several items exist only to record that a defensible wooldid choice has no precedent a reader will
recognise. Method: seven lenses were swept independently — display, bootstrap, permute, estimator,
features, syntax, terminology — each reading the ado and the help, running wooldid against the shipped
fixtures where display mattered, and comparing against the installed copies of csdid, jwdid,
did_imputation, xtevent, ritest, boottest, honestdid and Stata's own `bootstrap`, `permute`,
`didregress`, `hdidregress` and `bstat` help and source. Each finding was then handed to a separate
adversarial verifier with instructions to check both halves (the wooldid claim and the common-practice
claim) against the tree and against the comparator's own source rather than from memory. The sweep
produced **99 catalogued items** (32 display, 12 default, 17 feature-missing, 14 feature-unusual,
7 syntax, 8 terminology, 9 documentation), consolidated from 130 raw findings — two or more lenses raised
the same point independently in about thirty cases, and each such merge is noted at the head of the item,
so nothing was dropped in the consolidation. Of these, **48 carry verification evidence** — an
adversarial verifier pass, a live run of the specific behaviour, or both — and **51 were not reached before
the shared Stata process cap saturated** and are reported as the finding's author wrote them, marked "Not
verified (cap)". Among the 48: **11 confirmed on both halves**, **17 confirmed on the wooldid half with a
correction to the practice half**, **13 partly confirmed** (a correction to the wooldid half or to both),
**7 confirmed by direct observation** of the behaviour at issue, and **0 refuted outright**. Every
correction the verifiers made is carried into the text below and flagged, and the refuted sub-claims are
listed in section 5 so nothing is silently dropped. No item in this report
re-opens an author ruling of 2026-09-06; where a ruling itself is unusual relative to practice that is
noted in one line and the item is listed in section 4.

## 1. The fifteen most likely to be noticed

| id | area | what wooldid does | common practice | harmonization option | priority |
|---|---|---|---|---|---|
| PS-33 | standard errors | a bare `wooldid y i t ttre` returns iid SEs; the header prints "Standard Errors: " with nothing after the colon | csdid, jwdid, did_imputation and xtdidregress all default to clustering on the panel/cohort id | label the default ("unadjusted (iid)"), or default to `cluster(i)` | high |
| PS-34 | event study | the default reference is a pooled, cohort-varying block of pre-window periods; the plot's zero sits at a `*` tick, not at -1 | csdid `long2`, eventstudyinteract, xthdidregress normalise at a single omitted period, conventionally -1 (did_imputation and jwdid pool, as wooldid does) | print the regime under the ES block; relabel the `*` tick | high |
| PS-62 | stored results | analytic `e(V)` is `diag(se^2)` — every off-diagonal exactly 0 — so `test` and `lincom` answer rather than refuse; under `bootstrap()` it is the full bootstrap covariance | every comparator posts the sampling covariance that produced its own table | one line under the analytic table mirroring the bootstrap footer's own disclosure | high |
| PS-2 | export | under the default `bootdisplay(percentile)` the printed p and CI are rank-based while `e(V)`/`e(df_r)` carry the Wald family that esttab, coefplot and honestdid will pick up | Stata's prefix keeps the printed and stored families the same and puts percentile/BC in `estat bootstrap` | name the exporters in the footer beside `test` and `lincom` | high |
| PS-63 | sample formation | refuses (r(459)) when every in-sample unit is eventually treated; only `safety(off)` proceeds, and then the pooled regime runs fine | csdid falls back to `notyet` automatically; jwdid's default is not-yet-treated; hdidregress has `controlgroup(notyet)` | downgrade to a caution under the pooled regime; keep the refusal under `esfixedbaseperiod` | high |
| PS-64 | pre-trend tests | `jointtests` is a hard error with analytic clustered SEs when i is nested in the cluster variable — i.e. under `cluster(i)`, which the help twice recommends | did_imputation's `pretrends()` and `estat ptrends` run at panel-level clustering without a gate | cross-reference the restriction from the i and runtime entries; offer the bootstrap/permutation joint test as the automatic substitute | high |
| PS-45 | event-study inference | every interval is pointwise; the words simultaneous/uniform/sup-t appear nowhere; the only FWER device is permute-only and post-only | csdid's `wboot` defaults to uniform bands (`pointwise` opts out); R `did` plots simultaneous bands by default | say "pointwise" in the table and plot note; a sup-t band from the stored draws | high |
| PS-83 | vocabulary | "cohort" (223 uses) is wooldid's word for the panel id `i`; the timing group is `ttre`, reachable as `subgroups(ttre)` | in CS/Sun-Abraham/Wooldridge and in xthdidregress and csdid, "cohort"/"group g" IS the first-treatment period | a "note on words" box mapping i/ttre onto ivar/gvar | high |
| PS-47 | covariates | covariates enter only the outcome regression; "doubly robust", "propensity" and "IPW" appear 0 times in the help | csdid's default with covariates is `dripw` (Sant'Anna-Zhao 2020) | a paragraph stating the position: cell-saturation makes the outcome model near-nonparametric, and naming where DR would still differ | high |
| PS-19 | pre-trend display | the `adversarial` p prints in the pre-trend block; a small p there means the path is NOT a single line — i.e. usually good news — the opposite sign to the line above it | xteventtest's `linpretrend` runs the same-shaped test on pre-period coefficients only, where rejection has the conventional sign | carry the reading in the printed string; add a "how to read it" paragraph | high |
| PS-48 | post-estimation | estimands are fixed at estimation time; `estat event` returns r(321); no calendar-time aggregate; `e(histogramestimates)` rows are `1.__00000B:1.__000005` | csdid and jwdid re-aggregate off one stored fit (`estat simple/group/calendar/event`); hdidregress has `estat aggregation` | label the histogram rows with cohort values; add a calendar aggregate; an estat that re-aggregates `e(gteffects)` | high |
| PS-1 | stored results | the cluster count is `e(Nclus)`; `e(N_clust)`, `e(vce)`, `e(vcetype)` and `e(clustvar)` are never posted | csdid posts `clustvar` and `N_clust`; didregress documents `e(N_clust)`; estout keys on exactly these names | post the conventional names alongside `e(Nclus)` | high |
| PS-3 | table rendering | matlist's default format: `p 1.12e-35`, `t 13.33597`, header `se_anal~c`; `ereturn display` on the same e() prints `13.34  0.000` | every comparator routes through `_coef_table` and looks identical to `regress` | per-column `cspec()` formats, or build the table with `_coef_table` | high |
| PS-20 | banner | `impvalues` changes both the inference and the estimand while the header still reads "Standard Errors: cluster(i)" | wooldid's own bootstrap and permute banners name their method; Stata names a non-default VCE in the header | extend the existing banner mechanism to the IM path | high |
| PS-46 | permute export | under the fast-mode default no `e(V)` is posted, so esttab/coefplot/test/lincom have nothing; with `permanalyticse` the posted `e(V)` is the ANALYTIC one | ritest restores the observed run's analytic e(b)/e(V) and users `estadd` the RI p — the same hazard, quietly | document the estadd/esttab recipe; caution under `permanalyticse` | high |

## 2. The catalogue

Grouped by kind. Within each item the "wooldid" paragraph already carries the verifier's corrections;
the verification line records what changed.

### 2.1 Display

#### PS-1 — the cluster count is `e(Nclus)`; `e(N_clust)`, `e(vce)`, `e(vcetype)` and `e(clustvar)` are never posted

**wooldid.** reghdfe's `e(N_clust)` is harvested into a local (ado:1887) and reposted under the name
`Nclus` (ado:4437, `foreach scalar in N Nclus ll ll_0 ...`). A live `wooldid y i t ttre, att cluster(i)`
on `dta_binary` shows `e(Nclus) = 80`, `di e(N_clust)` missing, and `e(vce)`, `e(vcetype)`, `e(clustvar)`
all empty. A grep of the whole ado for those four names returns exactly one hit — ado:1887, the read.
The help documents `e(Nclus)` once and never mentions `N_clust`. The run also prints no
"(Std. err. adjusted for 80 clusters in i)" line, so there is no on-screen confirmation either.

**Practice.** csdid posts `ereturn local clustvar` and `ereturn scalar N_clust` unconditionally whenever
a cluster variable is in play (csdid.ado:918-921); didregress documents `e(N_clust)` (didregress.sthlp:478,
529); reghdfe posts `N_clust`, `vce`, `vcetype`, `clustvar`. `esttab`/`estout` `scalars(N_clust)`,
`estimates table, stats(N_clust)` and most lab table templates key on those names.

**Why it matters.** A table template asking for `N_clust` prints a blank cell with no error; `di e(clustvar)`
after a clustered wooldid run returns nothing. wooldid's own code has already been bitten once by the
absence (a footer keyed on `e(clustvar)`).

**Harmonization.** Post `e(N_clust)` alongside `e(Nclus)` (one line at ado:4437) and set
`e(vce)`/`e(vcetype)`/`e(clustvar)` from the `vcespecifier` already assembled at ado:1532-1534. Keeping
`e(Nclus)` preserves backward compatibility. Cost: a handful of lines plus a help-table edit.

**Verification.** Confirmed on the wooldid side (reproduced live, including `e(Nclus) = 80`). Practice
partly: csdid's printed "(Std. err. adjusted for ...)" header fires only on its wild-bootstrap path, so
the printed-header comparator is better cited to didregress/xtdidregress and reghdfe. The stored-results
contrast — the point of the item — stands.

#### PS-2 — the printed inference family and the exported one differ, and after a default `permute()` run there is nothing to export

**wooldid.** Under `bootdisplay(percentile)` (the default) the table prints rank-based p-values and
asymmetric percentile intervals while `e(V)`, `e(df_r)` and `e(V_es_*)` carry the standard-deviation
(Wald) family referred to t(G-1). Verified on a live run: `pretreatment_att` b = -.04572, sqrt(V) = .07489
so the Wald p an exporter would star is .543, while `e(boot_p_percentile)` for the same row is .44. The
footer at ado:10113 discloses this for `test` and `lincom` only; the help says the same at sthlp:503, 504,
1061, always naming `test`/`lincom`, while sthlp:1030 and 1173 advertise esttab/estimates-store
compatibility unqualified and sthlp:1044 documents the honestdid positional handoff without saying which
family `e(V_es_*)` then carries. Under the `permute()` fast-mode default `e(properties)` is `b` and no
`e(V)` is posted at all, so the permutation p, interval and null band live only in `e(perm_p)`,
`e(perm_ci)`, `e(perm_nullband)`, which no standard exporter reads.

**Practice.** Stata's `bootstrap` prefix posts `e(V)` = the bootstrap covariance and keeps
percentile/BC/BCa in a separate `estat bootstrap` table, so the printed and exported families never
diverge; didregress's `wildbootstrap()` likewise prints its bootstrap p beside an `e(V)` that still holds
the analytic vce. did_imputation posts a joint cluster-robust V and displays it. **Counterexample worth
recording:** csdid does exactly what wooldid does, with less disclosure — with `wboot` the printed
intervals are uniform bands by default while "asymptotic SE are stored in e()" (csdid.sthlp:238), and
there is no footer.

**Why it matters.** The stumble is one step downstream of a correct screen. An author who runs the
default bootstrap and then builds the paper table with esttab, the figure with coefplot, or the
robustness appendix with honestdid silently publishes the Wald/t(G-1) family while the log they checked
showed percentile ranks; honestdid's "Original" row will not match the wooldid event-study table it was
handed.

**Harmonization.** (a) One footer clause naming the tools, not just the commands: "esttab, coefplot,
honestdid and anything else reading e(b)/e(V) will report the bootdisplay(wald) family, not the columns
above." Trivial. (b) A documented one-line repost recipe that reproduces the displayed family (also the
cheap fix for the permute side — see PS-46). (c) Making `wald` the default display is a recorded ruling
and is listed only for completeness.

**Verification.** Partly. The bootstrap divergence reproduced live (the sweep's original .76/.548 pair
came from a different spec/seed and did not); the footer really does name only `test` and `lincom`. Two
corrections: the permute-side gap is already documented in the help (sthlp:263-267 names esttab
explicitly), and the "everyone else agrees" framing is undermined by csdid, which is the closest analogue
and diverges more.

#### PS-3 — numbers are rendered by matlist's general format, and one column header truncates with a tilde

**wooldid.** The tables are `matlist e(wdidmainresults)` with no `format()` (ado:9766, and likewise at
9756 summarystats, 9915 sgresults, 9960 pwsgs), so matlist's default general format governs. A plain
`att` run prints `att | .8445031  .0633252  13.33597  1.12e-35  .7201345`, while `ereturn display` on the
very same e() prints `.8445031  .0633252  13.34  0.000  .7201345  .9688717`. Under `permute()` the
standard-error header renders as `se_anal~c`. matlist also wraps the seven columns into two stacked
blocks, and `relyr` prints as `.` on the att row.

**Practice.** Every comparator routes its table through `_coef_table` (csdid.ado:214 with
`_coef_table_header` at 199; did_imputation.ado:711-712), which prints coefficients and SEs at %9.0g, the
statistic at two decimals, p at three with a leading zero, and a spanning `[95% conf. interval]` header —
so csdid, reghdfe, xtdidregress and did_imputation all look identical to `regress`.

**Why it matters.** `1.12e-35` and `t = 13.33597` in a coefficient column read as raw debug output, and a
truncated header with a tilde reads as a bug. Readers put wooldid tables next to `regress` and `csdid`
tables on the same screen.

**Harmonization.** matlist's default already IS `format(%9.0g)`, so passing that is a no-op (measured:
byte-identical output), and widening to `%12.0g` un-truncates `se_analytic` but prints ten decimals — a
single global `format()` cannot fix both, because the truncation is width-driven. The routes that work
are per-column `cspec()` (matlist.sthlp:322-393) — one string per matlist site, and the column set varies
by estimand, ES/subgroup mode and inference engine, so several context-dependent strings — or building
the table through `_coef_table`.

**Verification.** Partly. Reproduced verbatim. Three corrections carried in above: "leading zeros
dropped" is not a deviation (only the p column differs); the bootstrap `p = 0` case is a ruled item (F038)
and is already annotated in-line and in the footer, so it is excluded here; ado:9995 (gteffects) does pass
an explicit `format(%9.0g)`. The finder's proposed fix (`format(%9.0g)`) was demonstrated not to work.

#### PS-4 — the main table is ~91 characters wide, so it prints as two stacked panels, and one of its seven columns is empty on most rows

**wooldid.** The printed columns are `estimate se t p lb_95ci ub_95ci relyr` (ado:3849-3850; `z` instead
of `t` under poisson). At the batch/console default linesize 79 matlist breaks after `lb_95ci` and prints
`ub_95ci` and `relyr` in a second panel that repeats every row stub; at `set linesize 160` the whole table
fits in one 91-character panel. `relyr` is `.` on the main and pretreatment rows, and on the ES rows it is
redundant with the row stub, which already encodes the same integer.

**Practice.** Stata's `_coef_table` is sized to fit the default width with six columns and never wraps;
csdid, reghdfe, did_imputation and xtdidregress all fit. Horizon labels live in the row names, not in a
numeric column. (did_multiplegt_dyn does print descriptive `N`/`Switchers` columns, so the narrower claim
— horizons belong in the stub — is the safe one.)

**Why it matters.** The confidence interval ends up on a different screen from the point estimate, and on
event-study runs the horizon integer is in the panel that wrapped away. Copying a wooldid table out of a
log takes manual reassembly.

**Harmonization.** Drop `relyr` from the printed table (keep it in `e(wdidmainresults)`/`e(escomp)`) —
that removes ~11 characters, leaving an 80-character header line, right at the linesize-79 boundary, so
shortening `lb_95ci`/`ub_95ci` to a spanning `[95% conf. interval]` or to `lb`/`ub` may be required rather
than optional. Cost small, but it changes the printed column set the assessment baselines are blessed
against.

**Verification.** Partly. Width, split and the dot column all confirmed. Corrections: the pairwise table
(`e(pwsgs)`) has no `relyr` column at all, so the dots appear only in the main table; and the finder's
"roughly 78 characters" arithmetic was optimistic (80, not 78).

#### PS-5 — the header is a bespoke prose block rather than Stata's estimation header, and the SE label can print empty

**wooldid.** ado:9718-9732 prints three left-aligned `dis as result` lines: "Wooldid Estimation for
Outcome: y; Standard Errors: `vcespecifier'", then "N clusters = 80; N = 960 (0 obs dropped from initial
sample)", then the R2 line. `vcespecifier` is set only by `robust`, `cluster()` or `customvce()`
(ado:1532-1534, 1566), so with none of them the line ends after the colon with nothing.

**Practice.** csdid (csdid.ado:199) and did_imputation (did_imputation.ado:711) call
`_coef_table_header`, Stata's right-aligned Number of obs / F / Prob > F / R-squared block; csdid
additionally prints "(Std. err. adjusted for G clusters in var)" on its wild-bootstrap path, and Stata's
own commands surface the vcetype in the SE column header plus a cluster line.

**Why it matters.** A label with nothing after it reads as a rendering bug rather than as "classical".
See PS-33 for the substantive half.

**Harmonization.** Print an explicit default token when no VCE option was given, and post `e(vcetype)`
(PS-1). One line, independent of any behaviour change.

**Verification.** Confirmed on the wooldid side. Practice partly: csdid's explicit cluster line is
gated on the WBoot path.

#### PS-6 — the only fit statistics shown are four R2s belonging to the near-saturated stage-1 regression

**wooldid.** The header prints `R2 / R2adj / R2-within / R2-withinadj` (ado:9729-9730; "pseudo-R2" under
poisson) from the stage-1 reghdfe fit, harvested at ado:1888-1901 and re-read at 9701-9705. `e(F)`,
`e(ll)`, `e(rmse)` and `e(rss)` come from the same fit, are posted at ado:4434-4438, and are never
displayed. The figures move with the event-study request: R2 = 0.8822 plain, 0.8989 at
`esprelength(2) espostlength(2)`, 0.9104 at 3/3, 0.9353 with `esfixedbaseperiod` — same data, N = 960
throughout, more stage-1 dummies. R2-within moves further (0.4644 / 0.7058 / 0.5405 / 0.5928). The help
lists the scalars (sthlp:1130) without saying whose regression they describe.

**Practice.** csdid displays no R2 (its only `r2` tokens are Mata locals); did_imputation has no `r2` at
all; jwdid — the closest cousin, same estimator family, same saturated auxiliary regression — displays a
bare `ereturn display` with no header block and no R2. reghdfe does print R2/R2-within, but there the
regression is the object of interest.

**Why it matters.** Readers copy "R2" into tables. Here it is a property of the fitting specification, not
of the reported treatment effect, and it rises mechanically with the number of horizons requested.

**Harmonization.** Drop the R2 line from the default header (keep the scalars in e()), relabel it as the
stage-1 fit's, or show it only under `verbose`/`summarystats`. One line; the trade-off is that reghdfe
users expect to see it.

**Verification.** Partly. All numbers reproduced. One correction that softens the punchline: the att
itself also moves with the window (.8445 plain, .8555 at 3/3, .8235 at 2/2), so this is not "same
estimand, wandering R2" — it is a fit statistic that belongs to the auxiliary specification. jwdid is a
stronger comparator than the two uninstalled packages originally cited.

#### PS-7 — event-study rows are labelled `pre#`/`post#`, the relative year is in the column that wraps, and no row names the reference period

**wooldid.** Horizons print as `pre4, pre3, pre2, post0, ...` under a row-equation stub (`ES_att`,
`ES_att_it`), built at ado:3802-3817; the signed relative year appears only in the `relyr` column, which
is the last of the horizon block. Under `esfixedbaseperiod` with the default `esrelativeto(-1)` the -1
horizon is omitted and has no row; under the pooled default a `pre1` row (relyr -1) IS printed with an
estimate and a standard error. The plot, by contrast, inserts an explicit reference point in both regimes
(`rownamesvar = "reference"`, ado:10581-10617), so the picture and the table disagree about whether a
reference exists. At the default the only on-screen statement of the regime is none; the opt-in
`summarystats` prints three reference-pool lines under the pooled regime only, and `pretrendpower()`
prints "reference relative period = ..." under the fixed regime only.

**Practice.** csdid_estat labels `Tm3/Tm2/Tm1/Tp0` and did_imputation labels `pre1..preK` / `tau0..tauH`;
neither prints an explicit zero row either, so the missing base row is normal. jwdid relabels to signed
relative periods. Published figures plant a zero at -1.

**Why it matters.** `pre1` is an estimated coefficient in one regime and the omitted base in the other, so
two wooldid tables of the same command can misalign, and nothing on screen says which regime is in force.
Anyone building a figure from `e(b_es_*)` must plant the zero themselves.

**Harmonization.** (a) One line under the ES block naming the regime and the reference ("reference:
relative period -1, omitted" / "reference: pooled pre-treatment periods"). (b) Optionally print an
explicit base row as the plot already does. (c) Renaming rows to signed relative time is a breaking
change the help, plots and tests all reference — and would move wooldid away from csdid and
did_imputation, which use the same unsigned vocabulary.

**Verification.** Partly. Structure confirmed from source. Corrections: the `pre#` vocabulary is NOT the
deviation (did_imputation uses `pre#` too), so the item is really about regime ambiguity and the
table/plot asymmetry; "nothing on screen" holds at the default but not under `summarystats`; and the wrap
claim was not re-run.

#### PS-8 — joint pre-trend, across-subgroup equality and adversarial tests print a p-value with no statistic, degrees of freedom or reference distribution

**wooldid.** The line is "p-value for joint test that all pre-treatment att horizon estimates are 0:
0.6690" (ado:9790, 9793, 9796, 9799; subgroup twins at 9929-9941). At the point of computation only the p
is kept — `local joint_horiz = r(p)[1,2]` (ado:2515) from a `margins ... contrast(overjoint)` call — so the
statistic and its degrees of freedom, which margins returns, are discarded and never reach e(). The
across-subgroup chain reads `r(df)` only as a validity gate. The adversarial helper likewise returns only
`p`. wooldid's own `pretrendpower()` block, by contrast, prints "joint pre-trend test: F(k, df) = ...,
p = ..." (ado:9541, 9873).

**Practice.** `estat ptrends` after xtdidregress stores r(F), r(df_m), r(df_r), r(p)
(didregress_postestimation.sthlp:554-557); did_imputation posts `e(pre_F)`, `e(pre_p)`, `e(pre_df)`
(did_imputation.ado:703-706) and documents all three; csdid_estat prints "chi2(df) = ..." and returns the
chi2 (csdid_estat.ado:56-59).

**Why it matters.** Papers report "F(3, 79) = 0.52, p = 0.67". A bare p cannot be written into a table,
cannot be referred to a different reference distribution, and cannot be checked against a hand `test`. It
is also internally inconsistent with wooldid's own Roth block.

**Harmonization.** Harvest `r(chi2)`/`r(F)` and `r(df)` at the margins sites, post them as
`e(joint_horiz_stat)`/`e(joint_horiz_df)` and analogues, and print them — the pretrendpower wording is the
template. Under `bootstrap()`/`permute()` the analogue is the observed W, already stored as
`e(boot_jointWobs)`/`e(perm_jointWobs)`. Cost: a few lines per test site plus e() and help entries.

**Verification.** Confirmed on both halves, against the installed comparator source.

#### PS-9 — the inference footers are longer than the tables they annotate

**wooldid.** A `permute(shuffle) permreps(50)` run on a six-row table prints a 41-line, ~2,841-character
footer at linesize 79; the bootstrap footer is 25 lines / 295 words under the default percentile display
and four blocks under `bootdisplay(wald)`, with an additional ~141-word `boottstat` caution above the
table. The ado carries ~100 user-facing "Caution:" strings, the longest 992 characters. There is no
brevity or suppression switch anywhere in the file. Several permute paragraphs are the author's own
2026-09-06 rulings, so the observation is cumulative volume, not any one paragraph.

**Practice.** Stata's `bootstrap` prefix prints a 4-6 line header and puts the family in a column label
("Normal-based [95% conf. interval]"); `estat bootstrap` adds a compact table and a short key; ritest
prints a table and two notes; boottest about six lines; csdid two lines naming the outcome and treatment
models (csdid.ado:200-201). Standing explanation of what a percentile interval is lives in the help and
the manual, not in the estimation output.

**Why it matters.** Volume defeats its own purpose: the two disclosures a reader needs on a bootstrap run
— which family the p column belongs to, and that the SE column is not that family's dispersion — sit among
six other paragraphs of standing methodology and get skipped with them. A loop over subgroups or outcomes
buries the numbers under repeated identical prose, and the "not a wooldid peculiarity" sentence reads as
defensive in a log file.

**Harmonization.** (a) Split run-specific from standing: keep method, completed/failed, family, level,
resolution and any conditional disclosure; move the percentile-asymmetry explanation and the e()-matrix
inventory to the help behind one pointer line. (b) Print the full block once per session. (c) A
`bootnotes(brief|full)` switch defaulting to brief. Cost: text plus a session flag; nothing numeric moves,
but it changes the log text the assessment baselines are blessed against — and note that several
2026-09-06 rulings add clauses to this same block, so the ruled direction is currently toward more text.

**Verification.** Partly (counts corrected: 41 lines not 45; ~100 distinct caution strings not 105; the
reps-resolution line and the runtime projection belong to the bootstrap footer and the progress line, not
the permute footer). Practice confirmed for csdid from source, the rest from the tools' help.

#### PS-10 — `e(b)` carries no equation names: the row-equation colons are flattened to underscores before the post

**wooldid.** The results matrix is built with row equations (`main:att`, `ES_att:post0`,
`ES_att_it:pre2` — ado:3788, 3805, 3817, 3830, 3842, applied at 3848), and ado:4230 then does
`subinstr("`rownamespool'", ":", "_", .)` before the post, so `e(b)`/`e(V)` colnames are flat:
`main_att`, `pretreatment_att`, `ES_att_pre2`, `ES_att_it_post0`. Subgroup twins become
`ES_att_SG2_post0` (ado:4239). The printed table keeps the colons, so the table shows equation stripes
that `e(b)` does not have. The help tells event_plot users to match the flat stubs
(`stub_lag(ES_att_post#)`, sthlp:1168-1170).

**Practice.** csdid posts a multi-equation `e(b)` (`matrix coleq b1 = ...`, csdid.ado:877-886) and displays
it with `neq()`; coefplot's `keep()/rename()/eqlabels()`, esttab's `eqlabels()/drop()` and `test [eq]` are
all written around equation names.

**Why it matters.** Selecting an event-study block for a plot or a table means prefix-wildcarding rather
than naming an equation; the flattening pushes names toward Stata's 32-character limit sooner (a class of
overflow already fixed once for subgroup twins); and a user's first `lincom _b[ES_att:post0]` fails with
"equation ES_att not found", r(111), because the printed stripes are not the stored names.

**Harmonization.** Keep the colons — Stata accepts equation-named `e(b)`/`e(V)`, and `ereturn display`
would then print equation headers matching the table. Cost: every documented name changes, so it is a
one-time breaking change worth doing before wider release or not at all. Note the `test [ES_att]` motive
is weak while `e(V)` stays diagonal (PS-62): equations would buy selection convenience, not joint testing.

**Verification.** Confirmed on both halves.

#### PS-11 — a studentized interval endpoint can print as -3.77e+13 with no caution

**wooldid.** `boottstat` forms lb = b - q_hi*se_analytic and ub = b - q_lo*se_analytic from the signed
quantiles of t_r = (b_r - b)/se_r (`wooldid_boottp`, ado:11506ff; interval built at ado:8899-8900 and
9024-9025). The only per-replication guard is `if (A[r,k] >= . | S[r,k] >= . | S[r,k] <= 0) continue`, so a
tiny positive se_r is retained and yields an enormous |t_r| that becomes the quantile when reps is small.
On `dta_tiny` (8 clusters, 3 treated) with `reps(40) rseed(9)` the att row printed estimate 1.123,
se_analytic .166, t_pt 6.78, p .05, lb_95ci **-3.77e+13**, ub_95ci 1.687 — implying q_hi ≈ 2.3e14 — and no
caution fired: the two existing boottstat cautions key on availability (fewer than 2 usable replications,
ado:9163; fewer than half the completed ones, ado:9168), not on the magnitude of the quantile.

**Practice.** boottest, the Stata home of bootstrap-t inference, inverts the test over a grid and reports
the confidence set explicitly, including when it is unbounded or disconnected. Stata's bootstrap prefix
offers no percentile-t family at all. Davison and Hinkley (1997, ch. 5) warn that the studentized
bootstrap is unstable when replicate standard errors can be near zero and recommend variance
stabilization; MacKinnon and Webb (2017, 2018) recommend the null-imposed wild cluster bootstrap-t for the
few-cluster setting in which this configuration arises.

**Why it matters.** `[-3.77e+13, 1.69]` beside a p of .05 looks like the command broke. It is the
documented construction behaving as specified on a degenerate resample, but nothing on screen says so —
and few clusters with few treated is exactly the configuration a user reaches for `boottstat` in.

**Harmonization.** (a) Print the studentized quantiles the interval rests on; they are already in
`e(boot_tstud)` rows `qlo_pctl`/`qhi_pctl`, and a footer line makes the blow-up self-explaining.
(b) A caution when |q| exceeds some multiple of the t(G-1) critical value, or when
min(se_r)/se_obs over the retained replications falls below ~0.1 — `e(bootdist_se)` already holds every
se_r. (c) Winsorizing or dropping small-se_r draws trims the tail in the anti-conservative direction the
code's own comment worries about and should not be done without measurement.

**Verification.** Confirmed on the wooldid side (code path and log both). Practice partly: the
Davison-Hinkley / MacKinnon-Webb attribution was split as above.

#### PS-12 — bootstrap bias is computed and stored on every run, never shown, and never used as a diagnostic

**wooldid.** `e(bias_boot)` is posted on every bootstrap run as replicate mean minus observed, with the
comment "bias = replicate mean - observed (official convention)" (ado:9193-9197), alongside `e(b_boot)`.
Neither is printed, neither is mentioned in the footer, and nothing keys a caution on the ratio of bias to
standard error. The help lists `e(bias_boot)` by name with no gloss. The related centring diagnostic
`e(boot_z0)` (ado:9135) is likewise stored and never printed.

**Practice.** `estat bootstrap` prints Observed coefficient / Bias / Bootstrap std. err. / CI as its
standard layout; the bias column is the routine diagnostic for a badly centred bootstrap, and
Efron-Tibshirani's (1993, §10.3) rule of thumb is to worry once |bias| exceeds about a quarter of the
standard error. **Fairness note:** Stata's `bootstrap` prefix's own main table does not print bias either
— it surfaces only in the `estat bootstrap` redisplay.

**Why it matters.** This is a diagnostic the package already pays for and then withholds. It is also
exactly the diagnostic that would have surfaced the bayes `/W_c` centring defect (F001) on the first run:
a replication mean of 0.80 against an observed 1.69 is a bias of -0.89 against a standard error of 0.345.

**Harmonization.** (b) A Bias column in any `bootdisplay(all)` / re-display table (see PS-49) — this is
the practice-matching move. (a) A footer line naming the largest |bias|/se when it exceeds a threshold
(0.25 conventional), computed from `e(bias_boot)` and `e(boot_se)`, would go beyond what Stata does but is
a handful of lines.

**Verification.** Confirmed on the wooldid side. Practice partly: the finder's Stata help citation does
not support the claim (right practice, wrong lines) and Stata's main bootstrap table omits bias too, so
the deviation is the missing redisplay rather than the main table. The F001 example is now a fixed defect,
which weakens the rhetoric but shows the diagnostic class has been needed here.

#### PS-13 — permutation shift intervals print in columns still named `lb_95ci` / `ub_95ci`

**wooldid.** Under `permute()` the p column is relabelled `p_perm`, the SE column `se_analytic` and the
statistic column `t_perm` (ado:7503-7506, and 7548/7998 for the subgroup and pairwise blocks), but the
interval columns keep the analytic names — the only code that touches them is the `level != 95` branch
that rewrites the level suffix. The numbers behind them are shift intervals (the estimate minus the
equal-tailed quantiles of the null draws, the reflection of the footer's null band). The saved
event-study plot .dta carries the same column names (ado:10579, 10657-10660).

**Practice.** ritest reports no interval for the estimate: its CI column is a Clopper-Pearson interval for
the p-value itself (`cii reps c`, ritest.ado:869) with the printed note "Confidence interval is with
respect to p=c/n". Stata's `permute` prints a similarly-referenced interval. Published RI results report a
p-value, not an RI confidence interval. So the object in these columns has no counterpart in any tool the
reader knows, but carries the label of the object they know best.

**Why it matters.** A co-author reading an exported table, or a script keyed on `lb_95ci`, sees a column
named exactly like a conventional 95% CI. It is not: it is calibrated only near the null and its width
grows with the true effect (the F035 ruling documents this in text), and the column name is the one thing
that survives into a spreadsheet or a plot with the text stripped off.

**Harmonization.** (a) Rename to `lb_95ri`/`ub_95ri` under `permute()` and mark the saved plot .dta — cost:
the level-suffix convention is package-wide and scripts keyed on `lb_95ci` break under permute().
(b) Move the label to the column level; the footer already carries the full explanation, so this is a
placement fix, not a missing warning. (c) Blanking the columns behind an option is closest to ritest but
discards work and empties the ES plot bands.

**Verification.** Partly. Naming confirmed in code (line cite corrected to 7503-7511). Two corrections:
the finder's numeric example did not reproduce on the stated fixture and is dropped; and the footer
already states the shift-interval caveat in full, so the live gap is the column name, not the console text.

#### PS-14 — the step-down footer's "raw p" uses a different convention from the `p_perm` column above it

**wooldid.** `wooldid_permstep` (ado:13298-13301) computes p_raw as the share of replications with
|t_r| >= |t_obs| (add-one under sampling, plain share under `permexact`), while the table's `p_perm`
follows `e(perm_pconv)` — equal-tailed 2*min by default, or absolute/raw/one-sided under
`permabsp`/`permrawp`/`permside`. The footer prints them side by side as "horizon: raw p -> step-down p"
(ado:10044) with no convention label.

**Practice.** In rwolf/rwolf2 and in Westfall-Young / Romano-Wolf tables in applied papers the two columns
are both two-sided p-values on the same nominal scale and the adjusted one is weakly the larger, so the
pair reads monotonically as "before -> after adjustment".

**Why it matters.** The natural reading of "raw p -> step-down p" is "the p in the table above -> what it
becomes after adjusting". Here the two can disagree — at the floor by a factor of two (.02 in the table
vs .010 in the footer) — so a step-down that tightened can display as unchanged or looser. Under
`permrawp` the mismatch flips sign; under `permside` the table is one-sided while raw p is absolute.

**Harmonization.** (a) Label the convention in the footer clause; the max-|t| procedure is intrinsically
absolute-valued, so re-expressing raw p in the equal-tailed convention would create a different mismatch.
(b) Print only the adjusted p and let the table supply the unadjusted one — the only one of the two that
removes the non-monotone reading, since a labelled raw p can still print below the table's p.

**Verification.** Partly. Mechanism confirmed. Two corrections: the add-one is conditional (not applied
under `permexact`), and the "factor of two" holds only at the floor — the two conventions generally differ
by about 1/(n+1) additively. Restate as an unlabelled-convention / non-monotone-pair issue.

#### PS-15 — under `permute()` the joint pre-trend line prints the analytic wording verbatim, unlike the adversarial line beside it

**wooldid.** ado:9790, 9793, 9796, 9799 (and the subgroup twins at 9929-9941) print "p-value for joint test
that all pre-treatment att horizon estimates are 0: <p>" with no engine branch, while the adversarial line
ten lines later does branch and prints "permutation p-value (sharp null of no effect at any horizon) of
the adversarial line-fit statistic ...". The number under the unbranched line is a permutation quantity:
ado:7734 overwrites `joint_horiz*` with `wooldid_permqf`, the fixed-metric quadratic form whose null is no
effect at any horizon. The same unbranched line also prints under `bootstrap()` without `boottstat`, where
the F034 block (ado:6526) likewise replaces it with a fixed-metric rank bootstrap p.

**Practice.** csdid's `estat pretrend` names its test and prints chi2 with df; xtdidregress's `estat
ptrends` names the test and its null; did_multiplegt_dyn prints a labelled joint-placebo line. And in the
RI literature the distinction being elided matters: this is not a test of parallel pre-trends but a
permutation test of a quadratic form under the sharp null of no effect at any horizon.

**Why it matters.** A reader skimming takes 0.39 for the ordinary pre-trend F-test and concludes
"pre-trends fine". The package has already shown, on the line immediately below, that it knows how to say
otherwise.

**Harmonization.** Give the four display lines the same two-branch treatment the adversarial line has, and
branch on the bootstrap fixed-metric path too, not only `e(perm_method)`. Text only.

**Verification.** Confirmed on the wooldid side. Practice partly (csdid's is a postestimation command, not
its default display; did_multiplegt not installed). One softening: the permute footer at ado:10039 already
states the sharp-null nature of these tests, so this is a line-level labelling asymmetry rather than a
missing disclosure.

#### PS-16 — a p sitting exactly on its resolution floor prints as an ordinary number

**wooldid.** Under the default equal-tailed add-one convention the smallest attainable p is
2/(completed + 1) — exactly .02 at 99 completed draws. Nothing marks a p that has hit it: `wooldid_pfloor`
(ado:9632-9646), which prints "p = 0 above means p < ..." for the bootstrap surface, exits immediately
when `e(boot_completed)` is empty, so it has no permutation analogue, and the F036 caution fires only when
the floor exceeds .05 (permreps ≤ 38). A run can therefore print `t_perm 13.08, p_perm .02` in one row
with no marker.

**Practice.** Stata's `permute` prints c and n as columns alongside p = c/n plus the Monte Carlo error;
ritest defines p as c/N and prints the counts. Practice makes the resolution visible.

**Why it matters.** "t = 13.1, p = .02" reads like a defect, and .02 with no annotation reads as a
measured quantity rather than as the smallest number this many replications can produce. Because the
add-one convention makes p = 0 impossible, the permutation surface never triggers the F038 annotation the
bootstrap surface got — so the case that most needs a note is the one case with none.

**Harmonization.** A permutation analogue of `wooldid_pfloor`, keyed on `e(perm_completed)` and
`e(perm_pconv)`, roughly the same ten lines. It must key on `e(perm_exact) != 1`: under enumeration
p = 2/M is the true exact p, not a resolution artefact, so "raise permreps()" would be wrong there.

**Verification.** Confirmed on both halves. One nuance: the floor is 2/(completed+1), not
2/(permreps+1), and the ingredients (completed count, the convention, the resolution ladder) are on screen
in the footer and in the help — what is absent is any marker at the point of use.

#### PS-17 — `permute()` replaces the whole inference surface rather than adding an RI p beside the analytic one

**wooldid.** By ruling (ado:111-113, "only permutation p-values are ever shown under permute()"), every
p-value and interval in every table becomes a permutation quantity, the SE and statistic columns print as
`.` unless `permanalyticse` or `permtstat` is given, and a quantity with no usable permutation
distribution is blanked rather than filled from the analytic run. The footer says so.

**Practice.** ritest and Stata's `permute` are prefixes: the estimation output the researcher already
knows is printed by the command itself and the RI number arrives as a separate small table. In published
tables RI enters as one extra number per coefficient beside the clustered standard error.

**Why it matters.** The two things a referee usually wants side by side — the clustered SE and the RI p —
cannot be obtained from one wooldid call, and the second call is a full re-estimation rather than a cheap
add-on. A table with a blank SE column is not a form applied readers have seen.

**Harmonization.** (a) Extend `permanalyticse` so it keeps the analytic p and interval as extra columns
beside `p_perm` — a display and column-width question rather than new computation. (b) Document the
two-run + estadd recipe. (c) No change: mixing conventions in one table is exactly the confusion the
current design avoids. Recorded because the ruling itself is unusual, not as an open defect (F042).

**Verification.** Not verified (Stata cap).

#### PS-18 — the default output reports N, N clusters and R2 but no treated-unit, cohort or never-treated counts

**wooldid.** The header prints the outcome, the SE specifier, "N clusters = ...; N = ... (# obs dropped)"
and the four R2s (ado:9717-9736). The counts an applied reader wants exist — `e(escomp_trtN)`,
`e(escomp_refN)`, `e(escomp_refnunits)`, `e(escomp_refncoh)` — but are computed only when an event study is
requested and displayed only under `summarystats`. On the shipped fixture the ATT averages over 600
treated observations out of 960, a fact the default output never states.

**Practice.** did_imputation posts and reports the treated observation count; did_multiplegt_dyn prints
the number of switchers per horizon; csdid's aggregation tables carry group sizes; xtdidregress prints a
treated/control breakdown. The identifying sample is part of the header, not an optional diagnostic.

**Why it matters.** "How many treated units identify this estimate?" is the first question asked of a
staggered-DiD number. Meanwhile R2-within is printed prominently (PS-6).

**Harmonization.** Promote three numbers to the default header — treated units, onset cohorts,
never-treated (or out-of-window) units — from the escomp scalars already computed; or turn on the escomp
composition table whenever an event study is requested. Display only.

**Verification.** Not verified (Stata cap).

#### PS-19 — the adversarial line test rejects when the treatment effect is real, and nothing says which direction is good news

**wooldid.** With `adversarial` and `esfixedbaseperiod`, wooldid prints, immediately beneath the pre-trend
joint test and separated only by a blank line, "p-value for joint test of null that all (non-reference)
att coefficients are on a single line: 0.0007" (ado:9807, subgroup twin 9949). The line is fit through ALL
reported horizons, pre and post (sthlp:598-601, 625). On the clean binary fixture the pre coefficients are
.054/.183/.122 (p .57/.22/.43), the post are .972/.894/.864, the pre-trend joint test gives p = 0.669 and
the adversarial test gives p = 0.00075 with fitted slope .168. Nowhere in sthlp:592-615 is there a
statement of which direction of the p favours the analyst — and sthlp:592-593 points the wrong way,
describing the null as "adversarial toward finding pre-trend violations".

**Practice.** xteventtest's `linpretrend` — "requests a specification test to see if the coefficients
before the event follow a linear trend" (xteventtest.sthlp:30, 65) — is the literal published analogue and
is confined to the PRE-event coefficients, where rejection carries the conventional sign. xtevent's
`overidpre()`/`overidpost()`/`trend()` (Freyaldenhoven, Hansen, Pérez Pérez and Shapiro) keep the two
directions separate. csdid's pretrend test, did_imputation's `pre()` and every event-study plot share one
convention: a small p in the pre-trend block is bad news for the design.

**Why it matters.** Two p-values print one under the other in the same block and their signs are opposite.
Any design with a discernible break at onset drives the second to ~0. A reader skimming output — or a
coauthor reading a log — will read it as a failed diagnostic, and a reader who knows
`xteventtest, linpretrend` has a same-named test with the opposite reading.

**Harmonization.** (a) Carry the reading in the printed string: "test that the whole event-study path is a
single linear trend (small p = the path is NOT explained by a trend alone, which favours a genuine break
at onset)". One display string, twice. (c) A two-sentence "how to read it" paragraph at sthlp:592.
**Do not** use the finder's option (b): `altmethod_adversarial` shares the same reversed sign (a genuine
effect rejects the projected-line null too), so pairing them fixes nothing; only relocating the line out
of the pre-trend block would.

**Verification.** Partly. Numbers and placement reproduce exactly. Corrections: the help does use the word
"Rejection" at sthlp:604 (the accurate claim is that it never states the direction); "any DiD with a real
effect will reject" is too strong (a linearly-growing effect from onset will not); and harmonization
option (b) as written is unsound.

#### PS-20 — `impvalues` silently changes both the estimand and the inference while the banner still reads "Standard Errors: cluster(i)"

**wooldid.** The banner at ado:9718 prints only `e(wdid_vcespec)`, which is set from `robust`/`cluster()`
at ado:1532-1534 and overwritten in exactly two places — the permute path (ado:7656) and the bootstrap
path (ado:9150). No `impvalues` branch touches it. A grep for "Ibragimov|IM-style" in the ado finds only
error messages and comments, so a clean `impvalues` run prints nothing naming Ibragimov-Müller, nothing
saying the p-values come from a t-test on cluster-level estimates, and nothing saying `att` has become an
att_i-style equal-weighted average — which the help does state (sthlp:254-255, "IM-style inference
implicitly changes wooldid's target estimand"). The option is carried into e() at ado:5459, so the
information exists and simply never reaches the display.

**Practice.** wooldid's own bootstrap and permute banners name their method. Stata's estimation commands
surface a non-default VCE through `e(vcetype)` in the SE column header plus a cluster line; csdid sets
`e(vcetype) = WBoot` so the wild-bootstrap choice appears in the table.

**Why it matters.** Two runs of the "same" command with identical banners give different point estimates,
and a reader of the log cannot tell which is which. This is the situation the bootstrap/permute banners
were designed to prevent, so the inconsistency also reads as an oversight; and it is the
highest-consequence case, because an IM run's number is a different estimand, not just a different SE.

**Harmonization.** Extend the existing banner mechanism: "Standard Errors: Ibragimov-Müller cluster-level
t-test on cluster(i) (G clusters); note this makes att an equal-weighted average across clusters." One
line in the header construction plus a footer sentence. Cheap and self-contained.

**Verification.** Confirmed on both halves by code trace and help text. The finder's specific numbers
(.7647 vs .8445, 80 clusters) come from an unnamed fixture and were not reproducible this session; treat
them as illustrative.

#### PS-21 — stage-1 regression scalars are posted but never displayed, and the likelihood ones invite `estat ic`

**wooldid.** `e(ll)`, `e(ll_0)`, `e(r2)`, `e(r2adj)`, `e(r2_within)`, `e(r2_withinadj)`, `e(rmse)`,
`e(rss)` and `e(F)` are posted from the reghdfe fit (ado:4437 and the r2list above it); only the four R2s
reach the screen. `e(rank)` and `e(df_m)` are not posted. `e(F)` is present on an unadjusted run and
absent on the clustered twin — reghdfe's own behaviour, but the consequence is that wooldid's e() surface
changes shape with the VCE.

**Practice.** Stata posts `e(ll)` only for likelihood-based estimators and always alongside
`e(rank)`/`e(df_m)` so that `estat ic` and `lrtest` can work; reghdfe posts `e(df_m)`, `e(df_a)`, `e(rank)`.

**Why it matters.** A log likelihood on a two-step margins estimator invites an information criterion for
the wrong object, and a user diffing two runs sees a scalar appear and disappear with the VCE choice for
no stated reason.

**Harmonization.** Either complete the set (post `e(rank)`/`e(df_m)` and document that all of these
describe the stage-1 fit) or drop `e(ll)`/`e(ll_0)`. One line either way plus a sentence in the help.

**Verification.** Not verified (Stata cap).

#### PS-22 — `predict` and `margins` after wooldid fail with a confusing r(111)

**wooldid.** No `e(predict)` or `e(marginsok)` is posted, so Stata's generic `predict` tries to build xb
from `e(b)`'s coefficient names and errors with "variable main_att not found", r(111). `margins` gives the
identical message.

**Practice.** reghdfe posts `e(predict)`; csdid ships `csdid_p`; Stata's convention is a `_p` companion
that exits with "predict not allowed after ...".

**Why it matters.** The error points the user at a variable name they never typed and that does not exist
in their data, which reads as a wooldid bug rather than an unsupported operation.

**Harmonization.** Ship a two-line `wooldid_p` that exits with "predict is not supported after wooldid"
and post `e(predict) wooldid_p`. Self-contained.

**Verification.** Confirmed by run in the display lens (r(111) on both), not re-verified.

#### PS-23 — subgroup results print as one full titled table each, and the pairwise table carries a different column set

**wooldid.** A `subgroups(sg)` run prints "Main Results (Full Estimation Sample):", then "Results for
Subgroup sg == 1:", then "... == 2:", each a separate matlist with its own repeated seven-column header,
then "Pairwise Comparisons of Different Subgroups' Estimates:" — whose matrix has six columns (no `relyr`),
so the column set changes between adjacent tables in one run (ado:9887-9960). The across-subgroup equality
test then prints as another p-only line (PS-8).

**Practice.** `margins, over()` and the group-specific aggregations of csdid and did_multiplegt_dyn put
groups in one table as extra rows or equations. No Stata command varies its column set between two
coefficient tables of one run.

**Why it matters.** With four or more subgroups the output is several screens of repeated headers, and the
differing column sets make log diffs and scripted parsing brittle.

**Harmonization.** Stack the subgroup estimates as additional row-equations in the main table (the
machinery exists), keeping the pairwise block separate but with the same seven columns. Moderate; changes
the printed layout the baselines are blessed against.

**Verification.** Not verified (Stata cap).

#### PS-24 — the permutation p-value's Monte-Carlo precision is stored but never printed

**wooldid.** `e(perm_p_se)` holds the standard error of each p's underlying share and `e(perm_p_ci)` its
Clopper-Pearson interval (the executed F064 ruling), and the footer points at them. Neither the count c of
replications at or beyond the observed statistic nor the per-row replication count n appears in the
printed table, whose columns are `estimate se_analytic t_perm p_perm lb ub relyr`.

**Practice.** Stata's `permute` prints T(obs), c, n, p = c/n, SE(p) and a CI for p as table columns;
ritest prints the same shape.

**Why it matters.** At small `permreps` the p's Monte-Carlo noise is the dominant uncertainty and it is
invisible in the table; a reader cannot see that .039 and .078 are one draw apart without going into e().

**Harmonization.** Add SE(p) — or c and n — as a column under `permute()` (the table already varies its
column set by inference path), or append the SE to the p in the pooled-estimand summary lines. Low cost.
A presentation residual of an executed ruling, not a re-opening of it.

**Verification.** Not verified (Stata cap).

#### PS-25 — the printed bootstrap SE and the `e(V)` diagonal are different estimators, disclosed only below 50% availability

**wooldid.** `wooldid_bootfin` computes the reported SE per coefficient over that coefficient's available
replications (ado:11673) but builds `e(V)` listwise over replications complete on every non-excluded
column (ado:11690, divisor nlw-1), posting the row count as `e(boot_N_V)`. The code's comment states the
design. The only disclosure is a caution when a coefficient's availability falls below half the completed
replications (ado:9230). Between 50% and 100% availability the printed SE and sqrt(diag(e(V))) differ with
nothing said. On the sweep's fixtures the two agreed, so this is code-proven and latent rather than
reproduced.

**Practice.** In Stata output the Std. err. column is by construction sqrt of the diagonal of `e(V)`; no
Stata estimator prints a standard error its own `e(V)` does not contain. Stata's prefix records
incompleteness as `e(N_misreps)` and keeps one variance object.

**Why it matters.** The available-case rule is a deliberate improvement that the help sells as one of the
things wooldid's bootstrap does better than Stata's prefix, and it fires precisely on the long
event-study windows where the advantage is real — but once it fires, a user who checks the table against
`e(V)`, or exports with esttab, gets two different standard errors for the same coefficient.

**Harmonization.** (a) Lower the disclosure trigger: one line whenever any reported column's availability
is below the completed count, naming `e(boot_navail)` — the counts already exist. (b) Or state the rule in
the footer beside the existing `e(boot_N_V)` sentence. The code comment records pairwise-complete `e(V)` as
measurably indefinite.

**Verification.** Not verified (Stata cap).

#### PS-26 — the statistic and standard-error column headers (`t_wald`, `z_wald`, `t_pt`, `se_analytic`) exist nowhere else in Stata output

**wooldid.** Under the default percentile display the statistic column is headed `t_wald` (`z_wald` under
`bootdf(normal)`) and carries b divided by the bootstrap standard deviation — a Wald object whose tail
probability is not the rank-based p printed beside it. Under `boottstat` the columns are `se_analytic` and
`t_pt`. Observed: `estimate | se | t_wald | p | lb_95ci` on the default run, `... | t | p` under
`bootdisplay(wald)`, `... | se_anal~c | t_pt | p` under boottstat.

**Practice.** Stata's answer to the same problem is to print no statistic column at all: `estat
bootstrap`'s percentile/BC/BCa table is Observed coefficient / Bias / Std. err. / CI, with no z and no p.
No Stata command prints a column called `t_wald` or `t_pt`.

**Why it matters.** Recorded because the F053 ruling chose relabel over blank, and blank is what Stata
does. A reader exporting with esttab sees a header no documentation outside wooldid explains, and the
underlying oddity — a statistic beside a p that is not its tail probability — remains, now named. On the
other hand `se_analytic` and `t_pt` are genuinely more informative than an unlabelled column.

**Harmonization.** Nothing required; the ruling is executed. If revisited: blanking the statistic column
under percentile/bc removes the mismatch rather than annotating it; a middle option keeps `t_wald` and
adds the matching Wald p as a second labelled column.

**Verification.** Not verified (Stata cap). Listed also in section 4 as ruled.

#### PS-27 — the bootstrap coefficient p uses plain order-statistic counts, so it can print below the smallest value the tools it will be compared with can report

**wooldid.** The percentile and BC p-values are plain counts with no add-one adjustment, chosen so each is
the exact dual of its own interval (the help records 6 duality disagreements in 35,646 constructed cases
against 36 for the add-one alternative). The floor is therefore 0 at resolution 2/n rather than 2/(n+1);
the printed notes read "p = 0 above means p < 0.010 (2/200 completed replications)". The joint, equality
and adversarial p-values on the same table are add-one shares, so one printed surface carries two
conventions.

**Practice.** Davison-Hinkley (1997) and MacKinnon recommend the add-one bootstrap p; boottest and Stata's
`permute` report it, so at 999 replications the smallest p either can print is .001, never 0.

**Why it matters.** Comparability: a wooldid bootstrap p and a boottest p on the same data and the same
number of draws are not the same statistic, and a p of exactly 0.000 in a published table attracts a
question the add-one convention never gets. A reader also meets both conventions on one screen.

**Harmonization.** Nothing new to build (F038 and F033 are ruled/under research). Worth one help paragraph
placing wooldid's convention beside the add-one one — the formulas, the two floors, and the duality reason
— as the permutation engine's help already got under F052.

**Verification.** Not verified (Stata cap). Listed also in section 4.

#### PS-28 — the randomization-inference header reports clusters and possible assignments but never the onset composition

**wooldid.** The header prints "Randomization Inference (permute(shuffle)): 99 of 99 replications
completed; 1.1214e+36 possible assignments." above "N clusters = 80". The number of clusters at each onset
level — in particular the number of treated clusters — is computed internally (it is what the
multiset-permutation count is built from, ado:10949-10960) and never displayed. The few-treated-cluster
caution (F050) lives in `wooldid_bootctl` and fires only under `bootstrap()`.

**Practice.** MacKinnon-Webb organise the few-treated-clusters literature around G1; RI sections in papers
say "we permute treatment among the 50 states, 12 of which adopted". csdid and did_multiplegt print
cohort/group counts.

**Why it matters.** With few treated clusters the floors, the exchangeability risk the help now warns
about, and the power of the test all follow from the composition. "80 clusters, 1.12e+36 assignments" gives
no cue that the design might be thin.

**Harmonization.** One line under the RI header: "randomization design: k clusters per onset level (list),
n never-treated" — the counts are already in hand where M is computed. Optionally extend the F050 caution
to the permutation engine, worded for RI.

**Verification.** Not verified (Stata cap).

#### PS-29 — under `permtstat` the pooled summary line prints the estimate in outcome units beside a null band and SD on the t-ratio scale

**wooldid.** The footer prints "Pooled estimand summaries (p, 95% null band and SD of the permutation
distribution, on the t-ratio scale):" and then "estimate <b> p <p> null band [lo, hi] SD <sd>"
(ado:10067-10071). The parenthetical attaches to the band and SD; the estimate on the same line is in
outcome units. Observed: "att: estimate .82352847 p 0.0200 null band [-2.0230352, 2.310082] SD 1.136349".

**Practice.** Stata output keeps one scale per printed row; where a studentized quantity sits beside a
level the columns are labelled individually.

**Why it matters.** The line invites the reader to check whether the estimate falls inside the band — and
under `permtstat` 0.82 does fall inside [-2.02, 2.31], suggesting a null result, while the p on the same
line is .02.

**Harmonization.** Print the observed t-ratio rather than (or beside) the estimate on that line under
`permtstat`, or move the scale note into the row. One-line change; no numbers move.

**Verification.** Not verified (Stata cap).

#### PS-30 — the event-study plot is a connected line with a shaded ribbon and a zero-width reference point

**wooldid.** `twoway (connected estimate relyr, color(navy)) (rarea lb_95ci ub_95ci relyr, color(navy%50))`
with legend off, a dashed `xline(-0.5)` when any pre-period is reported, and the reference row planted as
estimate = lb = ub = 0 so the ribbon pinches to a point there (ado:10504, 10517-10521, 10586-10590). Under
the pooled default the reference tick is labelled `*`; under `esfixedbaseperiod` an asterisk is appended
to the tick's number.

**Practice.** The dominant Stata convention, and the default of Borusyak's widely used `event_plot`, is
markers plus capped vertical CI bars (`rcap`) at each event time with the omitted period shown at -1.
Shaded ribbons are for continuously indexed estimates; for a handful of discrete horizons a connecting
line implies interpolation between event times that does not exist.

**Why it matters.** It makes wooldid's own figures look unlike the figures in the papers the user is
writing, and the pinched-to-zero reference reads as "the effect is known exactly there", which is a
normalisation, not knowledge.

**Harmonization.** (a) Switch the CI layer to `rcap` and the point layer to markers, keeping the ribbon as
an option. (b) Or leave the style and note that `event_plot` reproduces the field-standard rendering from
the stored coefficients. (c) Either way, label the reference tick "ref" rather than `*`.

**Verification.** Not verified (Stata cap).

#### PS-31 — the column and row-equation vocabulary is wooldid-specific rather than Stata's coefficient-table words

**wooldid.** Every table prints columns headed `estimate se t p lb_95ci ub_95ci relyr` with equation rows
`main`, `pretreatment`, `ES_att` (and `ES_ame` under contreat). Subgroup blocks are headed "Results for
Subgroup ttre == 2005", contrasts "Pairwise Comparisons of Different Subgroups' Estimates" with row stub
`2008vs2005.ttre`, and the equality line reads "p-value for joint test that all SGs have the same att".
The interval bounds are hard-named `lb_95ci`/`ub_95ci` in the header even though `level()` is now accepted
on the analytic path.

**Practice.** Stata's coefficient tables print `Coefficient  Std. err.  t  P>|t|  [95% conf. interval]`
with the CI as a bracketed pair under one spanning header; csdid, did_imputation, boottest and reghdfe all
render that layout. "SG" as an abbreviation for subgroup appears in no other package's output.

**Why it matters.** The labels are clear on their own terms but mark the output as not-a-Stata-table,
which invites the question of whether the numbers are standard either. Two are actively confusing:
`lb_95ci`/`ub_95ci` will be wrong-looking under `level(90)`, and "SGs" reads as an acronym the user is
expected to know.

**Harmonization.** Rename the header stubs to Stata's words and make the CI header carry the actual level;
expand "SGs" to "subgroups". Keeping the wooldid names in `e(wdidmainresults)` column names while changing
only the printed header is a viable middle path.

**Verification.** Not verified (Stata cap).

#### PS-32 — error messages are paragraph-length by Stata standards

**wooldid.** The ado contains 218 `as error` messages with a median length of 117 characters and a mean of
143; 44 exceed 200 characters and 18 exceed 300, the longest being 786 (the `bootdisplay()` validity
message at ado:6148, which enumerates all four families with their rationale). The `jointtests` refusal is
570 characters over two lines. Many are two-part: an Error line plus a line of guidance.

**Practice.** Stata's own errors are one short clause ("option vce() not allowed") with the explanation
left to the help and the r(#) code; csdid, did_imputation and reghdfe follow that.

**Why it matters.** Two-sided. The long messages usually say what to do instead, which is more than most
packages manage. But in the Results window a 700-character error wraps to eight or nine lines and, inside
a loop, buries the failing iteration; users conditioned on one-line Stata errors may read the length as a
crash report.

**Harmonization.** Keep the content, change the shape: first line a short imperative clause naming the
offending option and the accepted values, subsequent lines as `as text` guidance (several messages already
do this). A mechanical pass over the 44 long messages; no behaviour change.

**Verification.** Not verified (Stata cap).

### 2.2 Defaults

#### PS-33 — no default clustering: a bare call reports homoskedastic SEs, and the SE label prints blank

**wooldid.** `local vcespecifier` is initialised empty and set only by `robust` or `cluster()`
(ado:1532-1534), so a bare `wooldid y i t ttre` passes no vce to reghdfe and returns unadjusted (iid)
standard errors, which margins then propagates by delta method. Measured on `dta_binary`: att .8445031
with se .0633252 unclustered against se .0526333 under `cluster(i)`. The header prints "Wooldid Estimation
for Outcome: y; Standard Errors: " with nothing after the colon (ado:9718), and because `e(Nclus)` is
missing the "N clusters =" line is suppressed. The string "unadjusted" appears nowhere in the ado or the
help; the `cluster()` and `robust` entries (sthlp:237, 240) never state what happens when neither is
given. (pweights imply robust SEs via reghdfe, so a pweighted bare call is not iid.)

**Practice.** did_imputation.ado:101 sets `if ("`cluster'"=="") local cluster = "`i'"`; jwdid.ado:195-197
sets `cvar = ivar` when `cluster()` is omitted; csdid's influence-function variance is unit-clustered by
construction; didregress/xtdidregress state that `vce(cluster groupvar1)` is the default
(didregress.sthlp:353-355). Bertrand-Duflo-Mullainathan (2004) and Cameron-Miller (2015) made unit/group
clustering the default expectation in this literature.

**Why it matters.** A user who copies the help's own Example 1 (`wooldid y i t ttre, att att_it
histogramcohorteffects`, sthlp:1233) gets iid standard errors on serially correlated panel data — the BDM
error — and nothing in the output says so: the label is blank rather than reading "unadjusted". Anyone
arriving from csdid, jwdid, did_imputation or xtdidregress will assume the command clustered for them.

**Harmonization.** (a) Default to `cluster(i)` when no vce option is given and say so in the header —
matches all four comparators; cost is that every previously unclustered run changes numbers, needing a
version note and a battery re-bless. (b) Keep the default but print "Standard Errors: unadjusted (iid) —
neither cluster() nor robust was specified" plus a caution. (c) Document the no-option behaviour in the
`cluster()`/`robust` entries. Verified in two independent lenses.

**Verification.** Confirmed on both halves, reproduced live. Fairness note: help Example 0 immediately
above does use `cluster(i)`, so the help demonstrates clustering — it just never says what happens without
it, and the ATT/ATT_it example a reader is likelier to copy is the unclustered one.

#### PS-34 — the pooled, cohort-varying, multi-period reference is the default, and three of the four pre-trend surfaces refuse to run under it

**wooldid.** `esprelength()`/`espostlength()` default to 0, so a plain call returns only the pooled ATT
with no horizons. When horizons are requested the default regime is `pooled_reference_period`: each treated
cohort's reference is the whole block of its pre-window pre-treatment observations, which can span
different numbers of periods for different cohorts (sthlp:154, 663-671), and the plot plants the zero at
min(relyr)-1 with the tick labelled `*` (ado:10581-10598). Under this default `pretrendpower()` requires
`esfixedbaseperiod` (sthlp:632); the adversarial p-value requires it (sthlp:592); and the direct honestdid
handoff requires "esfixedbaseperiod with the default esrelativeto(-1)", with the help warning that other
geometries hand honestdid the wrong event-time geometry **silently** (sthlp:1043-1050). `jointtests` does
run under the pooled regime (but see PS-64). A run at `espre(3) espost(3)` prints an estimated, nonzero
`pre1` row and no on-screen statement of the regime; `e(estimate_type)` is the only place it appears.

**Practice.** csdid's `long2` (and its default "short" varying base), eventstudyinteract/Sun-Abraham and
xt(h)didregress normalise at a single named period, conventionally -1, and honestdid/pretrends assume that
geometry with `numpre()` purely positional. **However:** did_imputation's pretrend reference group is "all
periods more than k periods prior to the event date" (did_imputation.sthlp:193-195), and jwdid's default
ETWFE specification pools every pre-treatment observation into the implicit baseline — so wooldid's default
is shared by the two imputation/ETWFE-family packages it is closest to, and differs from the CS/SA family.

**Why it matters.** A reader of a wooldid event-study plot sees the zero at the left edge under a bare `*`
tick, at a normalisation that differs across cohorts in calendar terms. The help itself warns this can
"obscur[e] even severe parallel trends assumption violations from detection". And the three tools an
applied reader would reach for to defend the design — Roth power, the line test's p, HonestDiD — are all
unavailable in the default configuration, so the default is the least-diagnosable regime, with the
honestdid path the one that fails silently rather than with an error.

**Harmonization.** (a) Flip the default to `esfixedbaseperiod` with `esrelativeto(-1)` — changes every
default run's numbers and moves wooldid away from the estimator family it belongs to. (b) Print a one-line
note whenever espre/espost > 0 under the pooled regime saying the zero is the pooled reference, is not
relative period -1, and that `pretrendpower()`/adversarial p/honestdid need `esfixedbaseperiod`. (c) Label
the reference tick "pooled ref" rather than `*`. (d) Refuse or hard-caution the honestdid-relevant
`e(numpre_*)` posting under the pooled regime — a small, consistent extension of an existing pattern.

**Verification.** Confirmed on the wooldid side (behaviour reproduced live). Practice partly: the pooled
multi-period reference is not unique — did_imputation and jwdid do the same — so the framing should be
"wooldid sits in the imputation/ETWFE camp, not the CS/SA camp, and the plot never says so". Separately,
the help's claim at sthlp:154 that pooled-reference event studies "always have a normalized effect of 0 in
the earliest period reported" is contradicted by the run (pre3 = .0891, se .1065): the normalized-zero
period is the unreported pool, a small help-text inaccuracy worth fixing alongside (b).

#### PS-35 — the few-cluster caution thresholds fire only well inside the region the literature already calls hopeless

**wooldid.** Two thresholds, both inside `if "`clus'" != ""` (ado:6250-6282), so an unclustered bootstrap
run gets neither: a caution below 10 total clusters (ado:6260) and, since 2026-09-06, a caution below 6
treated resampling clusters naming `permute()` and `boottest` (ado:6279, further gated on `ttrevar`).
Above those — say 25 clusters with 8 treated — a `bootstrap()` run prints confident percentile intervals
with no caution of any kind. The help does steer few-treated-cluster users to `permute()` in prose
(sthlp:38, 392) with no numeric threshold.

**Practice.** Cameron and Miller (2015) decline a sharp cutoff and cite rules of thumb around 50 clusters
below which cluster-robust asymptotics are suspect. MacKinnon and Webb (2017, 2018) show severe size
distortion with as many as 50 total clusters when only a handful are treated — over-rejection for the
pairs-type schemes wooldid implements — with simulations at treated counts in the 5-10 range, i.e. the
range wooldid's threshold stops just short of. Stata's `didregress` ships `wildbootstrap` (reps(1000),
rademacher) as the answer for this whole region, though it is an option rather than the default.

**Why it matters.** The thresholds are calibrated to "obviously broken" rather than to "the literature says
be careful here", so the common applied case — 20 to 40 state-level clusters, a handful treated — passes
silently. A reader who knows MacKinnon-Webb will notice that the tool warned on 8 clusters and said
nothing on 30.

**Harmonization.** Raise and grade: a soft note below ~50 total clusters (Cameron-Miller's own rule of
thumb) or ~10 treated, pointing at `permute()` and `boottest`, and keep the current stronger wording below
10 / 6; the message should say the pairs cluster bootstrap over-rejects there. Cost: two numbers and one
message string. Against: caution fatigue, which argues for softer wording at the higher tier.

**Verification.** Confirmed on the wooldid side (plus the unclustered-run gap the finder did not note).
Practice partly: `wildbootstrap` is an option not a default, MW2018's direction is variant-dependent, and
Cameron-Miller name no 20-50 "danger zone".

#### PS-36 — `reps()`/`permreps()` are required with no default, and `bootstrap()`/`permute()` require a mode argument

**wooldid.** `reps(integer 0)` and `permreps(integer 0)` both default to 0 and are refused:
`bootstrap()` without `reps()` exits with "Error: bootstrap() requires reps(#) with a positive number of
replications." (ado:6055-6058), the code comment recording the ruling "reps() is mandatory (author ruling:
no floor - low-rep runs are legitimate)". `permute(shuffle)` requires `permreps()` or `permexact`
(ado:5847). `bootstrap()` takes a required mode (`standard|bayes`) and `permute()` a required method
(`shuffle|probabilities`), so there is no bare `bootstrap` or `permute`. The help gives the guidance
("hundreds for standard errors, and a thousand or more before the percentile and bias-corrected intervals
are worth reading", sthlp:487) but the command proposes no number.

**Practice.** Stata's `bootstrap` prefix defaults reps(50); `permute` reps(100) with the honest note that
"100 permutations is chosen for convenience"; `didregress`'s wildbootstrap reps(1000); boottest 999; csdid
wboot 999; ritest 100. Every one runs on a bare invocation, and `vce(bootstrap)` needs no mode argument at
all.

**Why it matters.** Small friction, but it is the first thing a user meets: someone with bootstrap-prefix
habits types `bootstrap(standard)` and gets an error. More substantively, because wooldid declines to name
a number there is no canonical wooldid bootstrap the way there is a canonical 999-replication csdid or
boottest run, so two papers using wooldid are not running comparable inference unless both read the help
paragraph.

**Harmonization.** (a) Keep mandatory but make the error teach: name a starting point and the accepted
mode values in the first sentence. One string, zero behaviour change. (b) Default `reps()` — but note the
comparator arithmetic: boottest/csdid/didregress's 999-1000 are wild-bootstrap defaults where no refit
happens, while wooldid re-estimates per replication, so Stata's refit-bootstrap default of 50 is the apt
number, not 999. (c) Default `bootstrap()` to `standard` and `permute()` to `shuffle` when the parentheses
are empty.

**Verification.** Confirmed on both halves (all six comparator defaults read off the installed help files).

#### PS-37 — the default analytic SE is margins' conditional (delta) SE; matching Stata's own `hdidregress twfe` requires the non-default `unconditionalse`

**wooldid.** The help states that the default "estimates standard errors in the underlying regression and
then obtains standard errors on the treatment effects ... using the delta method. These errors are
conditional standard errors in the sense that they are conditional on the data" (sthlp:~88-92).
`unconditionalse` switches margins to `vce(unconditional)` but also forces the run off reghdfe/ppmlhdfe
onto regress/poisson and "should be slower to use on average" (sthlp:243-246). wooldid's own comparison
section says agreement with `hdidregress twfe` requires "unconditional standard errors" (sthlp:1199-1202).

**Practice.** Stata's `hdidregress`/`xthdidregress` twfe uses margins with unconditional standard errors;
csdid, drdid and did_imputation report influence-function standard errors that account for the sampling of
the estimation sample. Conditional margins SEs are unusual as a headline inference default in the DiD
literature.

**Why it matters.** A user replicating a wooldid result in `hdidregress twfe` (or vice versa) finds the
point estimates agree and the standard errors do not, with nothing in the default output naming the
difference — the header says only "Standard Errors: cluster(i)". The one place the distinction is
explained is a paragraph 1000 lines into the help.

**Harmonization.** (a) Name the concept in the header ("cluster(i), conditional (delta)" vs
"... unconditional") — one line, no numeric change. (b) A sentence in the `cluster()`/`robust` entries
pointing at the distinction and the hdidregress comparison. (c) Flipping the default is not advisable: it
would force every run off the HDFE engine.

**Verification.** Not verified (Stata cap).

#### PS-38 — `controls()` is a plain additive TWFE control; the Wooldridge-conforming specification is the separate `interactivecontrols()`

**wooldid.** "controls(): The contents of this option are passed directly to reghdfe as controls, without
any modification" (sthlp:873) — an additive, homogeneous-coefficient control. The specification that "most
closely mimics the approach to including covariates that appears in Wooldridge (2025)" is
`interactivecontrols()`, which enters the variables interacted with time dummies and, de-meaned by cohort,
interacted with every treated-cell dummy (sthlp:891-902). `coarsecohortcontrols()` is a third variant that
partials the control effect out of the reported treatment effects. The steer toward
`interactivecontrols()` appears only inside that option's own entry.

**Practice.** jwdid builds the covariate interactions by construction (the Mundlak/ETWFE regression);
csdid's doubly-robust estimators condition on covariates in the propensity/outcome models; did_imputation
distinguishes `controls()`, `unitcontrols()` and `timecontrols()` explicitly and estimates them on
untreated observations only. In all three, the natural-looking covariate slot does the method's
recommended thing. Wooldridge (2021, 2023) argues that plain additive controls in a TWFE DiD reintroduce
the contamination ETWFE is designed to avoid.

**Why it matters.** `controls()` is the option name a user will reach for first, and it gives the
specification the estimator's own author argues against, with no note in its entry pointing elsewhere. The
three covariate options also have three different semantics (in the treatment effect, partialled out of
it, additive) with no side-by-side summary.

**Harmonization.** (a) Cross-reference from `controls()` to `interactivecontrols()` and
`coarsecohortcontrols()` with one sentence saying which reproduces Wooldridge's covariate treatment.
(b) A three-row table under Controls and Fixed Effects. (c) Renaming or re-defaulting would break existing
code and move numbers.

**Verification.** Not verified (Stata cap).

#### PS-39 — under poisson the reported "ATT" is a difference in linear predictions, not the level effect

**wooldid.** "Treatment effects are presented in terms of differences in linear predictions (the poisson
predict(xb) option)" (sthlp:555-557); the level-scale version is opt-in via `poisexpresults`. The heading
merely appends " Poisson Regression" to "Main Results" and the row is still labelled `att`. `gteffects`
and `pretrendpower()` are likewise on the xb scale under poisson.

**Practice.** jwdid's postestimation calls plain `margins, subpop(...) at(__tr__=(0 1))` with no
`predict()` option, so under `method(poisson)` margins uses its default response scale and the reported
ATT is on the level (count) scale. Wooldridge (2023), "Simple approaches to nonlinear
difference-in-differences with panel data", defines the ATT for the exponential model on the levels scale.

**Why it matters.** Two commands both described as "the Wooldridge estimator with poisson" report numbers
on different scales, and wooldid's default is not the scale Wooldridge's own paper defines the estimand
on. A reader seeing a row labelled `att` under a poisson header will read it as a level or a percentage
effect; it is neither exactly — it is a log-index difference.

**Harmonization.** (a) Label the scale in the header or the row equation under poisson — display only.
(b) Make `poisexpresults` the default so `att` always names a level ATT (numeric change to every poisson
run). (c) At minimum say in the poisson entry which scale each comparator uses, so a replication mismatch
is diagnosable.

**Verification.** Not verified (Stata cap).

#### PS-40 — the ATT(g,t) table is opt-in, print-capped at 200 rows, and unavailable with `contreat()`

**wooldid.** `gteffects` returns the cohort-by-relative-period building blocks and the help says the cells
"are ATT(g,t)-style group-time effects when i indexes cohorts" and the table "is directly comparable to
the group-time tables of csdid or xthdidregress" (sthlp:953-965). It is off by default; it prints only for
tables of up to 200 rows (on the shipped fixture: "550 rows; table not printed - see e(gteffects)"); and it
is "Not available with contreat()".

**Practice.** csdid's DEFAULT output is `agg(attgt)` — the group-time table — with aggregations opt-in.
jwdid prints the aggregate but keeps `estat group/calendar/event` one keystroke away. The framing most
readers are trained on is "estimate ATT(g,t), then aggregate"; wooldid inverts it.

**Why it matters.** Readers who learned staggered DiD through Callaway-Sant'Anna look for the group-time
table first and find only aggregates; and the setting where the building blocks are most useful for
reasoning about heterogeneity (continuous treatment) is the setting where they are unavailable.

**Harmonization.** (a) Mention `gteffects` in the Description of Estimation Approach, framed as "the
ATT(g,t) table if you want it". (b) Make the 200-row print cap user-settable, or print a truncated head.
(c) The contreat() exclusion is principled (the cell coefficients are dose-zero intercepts) — say so where
a csdid user will look.

**Verification.** Not verified (Stata cap).

#### PS-41 — the honestdid handoff works only off the non-default event-study regime, and the geometry guard is positional and silent

**wooldid.** `e(b_es_type)`, `e(V_es_type)` and `e(numpre_type)` are returned and a worked honestdid
example is given (sthlp:1277-1290), but a direct positional handoff requires ALL of: `esfixedbaseperiod`
with the default `esrelativeto(-1)`; no `timetrends`; at least one pre and one post coefficient;
consecutive horizons apart from the intended -1 gap; and a non-degenerate covariance. Other geometries
"hand honestdid the wrong event-time geometry silently and need a user-supplied transformation"
(sthlp:1044-1052). The default regime is the pooled reference, which fails the condition, and
`e(numpre_att)` is posted anyway. Under `permute()` without `permanalyticse` the objects are not posted at
all.

**Practice.** honestdid ships a post-estimation path for csdid (run csdid, then honestdid with
`pre()`/`post()`), and Rambachan-Roth sensitivity is now a routine referee request. Users expect the
sensitivity analysis to follow the estimator without changing the specification.

**Why it matters.** The default specification is the one that cannot be handed to honestdid, and the
failure mode when the user does it anyway is a wrong answer rather than an error. The conditions are
documented in a Stored Results entry a user reaching for honestdid may never open.

**Harmonization.** (a) A run-time caution when `adversarial`/`pretrendpower()`/ES horizons are requested
under the pooled regime, saying honestdid input will not be available. (b) A helper (or documented
one-liner) that checks the colname signature and errors. (c) Note in the Description that the
honestdid/pretrends route requires the fixed reference regime. Neighbouring cases already refuse or
caution (PS-34), so this is the one remaining documented-but-unguarded path.

**Verification.** Not verified (Stata cap) as a standalone item; the underlying geometry facts are
confirmed under PS-34.

#### PS-42 — `bootdf()`'s default is state-dependent: t(G-1) when clustered, normal when not (more careful than practice)

**wooldid.** `bootdf()` defaults to t(G-1) when `cluster()` was given and to normal when it was not, with
G the number of resampling units; an explicit setting overrides. Observed as t(79) on an 80-cluster run and
t(7) on an 8-cluster run. But the option changes only the Wald family, and the default display is
`bootdisplay(percentile)`, whose intervals and p-values are quantile-based and unaffected — so on a default
run the entire effect of `bootdf()` is on `e(df_r)` and on stored tables nobody sees.

**Practice.** Stata's `bootstrap` prefix reports normal-based z regardless of whether `cluster()` was used;
csdid's wboot uses normal quantiles; boottest applies its own small-sample corrections. Cameron and Miller
(2015) recommend t(G-1) critical values with cluster-robust inference, so wooldid's default is the one the
literature endorses — it is simply the only package that varies the default by whether the user clustered.

**Why it matters.** Two mild stumbles: a user who reads the entry, likes the argument and leaves the
defaults alone sees nothing change; and a user comparing a wooldid bootstrap to a Stata-prefix bootstrap
finds the p-values differ for a reason invisible unless `bootdisplay(wald)` is on — which is also the
family esttab and coefplot pick up (PS-2).

**Harmonization.** Nothing to fix; this is a better-than-practice default and the help already says the
option does not move the percentile family. Optional: one footer clause under the percentile/BC displays
saying `bootdf()` affects only the stored Wald family and `e(df_r)`. **Keep as-is.**

**Verification.** Not verified (Stata cap).

#### PS-43 — `permreps()` has no default and a floor of 20 (more careful than practice)

**wooldid.** `permute(shuffle)` requires `permreps()` or `permexact`; `permreps()` below 20 is refused
outright ("the engine refuses inference on fewer than 20 completed replications", ado:5806-5808) and a
caution fires when 2/(permreps+1) exceeds .05, i.e. at permreps ≤ 38 (ado:5811-5817). The help gives the
resolution ladder: 50 → .039, 100 → .020, 199 → .010, 999 → .002.

**Practice.** ritest and Stata's `permute` both default to reps(100) and both run happily at reps(5).
wooldid's own `bootstrap()` imposes no floor on `reps()` by an explicit ruling.

**Why it matters.** This is a deviation in the safe direction: a silent 100-replication default gives a
p-value floor of .02 under wooldid's equal-tailed convention, which is exactly the trap the author has
otherwise worked to avoid, and a user who never chose the number would not know. It is recorded only
because it is asymmetric with the `reps()` no-floor ruling and because a ritest user will be surprised the
command will not run without the option.

**Harmonization.** Keep. If anything, name a sensible number in the refusal ("999 for a .002 floor; see the
resolution ladder"). Dropping the floor of 20 for consistency with `reps()` would be a step backwards.
**Keep as-is.**

**Verification.** Not verified (Stata cap).

#### PS-44 — singleton observations are dropped on the reghdfe path and kept on the regress path, so N moves with the SE option

**wooldid.** The reghdfe and ppmlhdfe calls (ado:1547, 1559) do not pass `keepsingletons` — the token
appears nowhere in the ado — so reghdfe's default singleton dropping and its degrees-of-freedom adjustment
apply. When `unconditionalse` or `customvce()` is specified, estimation switches to regress/poisson
(ado:1564-1570), which keep singleton groups. The reported N comes from the fit, and the "# obs dropped
from initial sample" figure in the header folds the engine's drops in with wooldid's own screens.

**Practice.** jwdid passes `keepsingletons` explicitly on every reghdfe call (jwdid.ado:203, 217, 243),
pinning the sample so it does not move with the engine. reghdfe's documentation treats singleton dropping
as a deliberate choice with inference consequences.

**Why it matters.** Switching from the default to `unconditionalse` — which the help itself recommends for
matching `hdidregress twfe` — can change the estimation sample and the residual degrees of freedom without
any message, so the two SE concepts are not being compared on the same data. The header's single "obs
dropped" number does not separate wooldid's screens from the engine's singleton drops.

**Harmonization.** (a) Pass `keepsingletons` so the sample is the same across SE options (numeric change on
affected designs; needs measuring first). (b) Report the engine's singleton drops separately in the
header. (c) At minimum, note in the `unconditionalse` and `customvce()` entries that the estimation sample
can differ. Practical exposure is designs with fine i or heavy imbalance; on a balanced panel with coarse i
there may be no singletons at all.

**Verification.** Not verified (Stata cap); source-only (exhaustive grep).

### 2.3 Missing features

#### PS-45 — no simultaneous (sup-t / uniform) confidence bands, on any path

*Three lenses (bootstrap, features, estimator) converged on this item independently.*

**wooldid.** Every reported and plotted interval is pointwise: the tables carry `lb_#ci`/`ub_#ci` per row
and the event-study plot is `(connected estimate relyr) (rarea lb_95ci ub_95ci relyr)` (ado:10651) over
the same per-row columns. A grep of the ado for sup-t / supt / uniform / simultaneous / pointwise returns
nothing relevant, and the same grep on the help returns nothing at all — so there is no uniform-band option
and the help never says the bands are pointwise. The material is present: `e(bootdist)` holds the full
reps × horizons matrix up to 10,000 completed replications, `e(bootdist_se)` every replicate SE, and
`e(V_es_type)` the full event-study covariance with genuine off-diagonals on both the analytic and
bootstrap paths. The only multiplicity device anywhere is the max-|studentized| step-down landed
2026-09-06 (ado:7782-7800, 13263-13322), which is `permute()`-only and filtered to post-treatment horizons
(`wdprl >= 0`); there is no bootstrap or analytic counterpart and it covers no pre-period.

**Practice.** csdid's `wboot` produces uniform confidence intervals BY DEFAULT and requires `pointwise` to
opt out (csdid.sthlp:133); R's `did` package plots simultaneous bands by default (multiplier bootstrap,
`cband=TRUE`); Callaway and Sant'Anna (2021) build the simultaneous band into their aggregation results.
Scope caveats: csdid's own non-wboot default is pointwise, and Stata's `bootstrap` prefix and other DiD
bootstrap paths are pointwise too — so the accurate statement is that the Callaway-Sant'Anna lineage makes
simultaneous bands the bootstrap default and a reader coming from there will expect them.

**Why it matters.** This is the single most visible difference between a wooldid event-study figure and a
csdid one, and it runs anti-conservative: a reader scanning six pointwise 95% bands for a pre-trend is
running six tests at 5%. A referee who knows Callaway-Sant'Anna will ask whether the bands are uniform,
and the help does not say they are pointwise. The author has already accepted the principle — the
Romano-Wolf step-down shipped for `permute()` the same day — so the bootstrap and analytic engines are now
the odd ones out inside wooldid itself.

**Harmonization.** (a) Cheapest: say "pointwise" in the event-study table header, the plot note and the
bootstrap footer. (b) A sup-t band from the stored draws — c = the (1-alpha) quantile over replications of
max_h |b_rh - b_h| / se_h, band = b_h ± c·se_h — one pass over the accumulator that already feeds the joint
test, roughly the cost of the step-down just built; expose as `bootdisplay(supt)` or `bootband()`,
defaulted off. (c) The same band simulated from `e(V_es_type)` (Montiel Olea & Plagborg-Møller 2019) would
serve the analytic path too, ~20 lines of Mata with no refit. (d) Extend the step-down to pre-treatment
horizons and to `bootstrap()`. Guard note: wooldid already flags structurally degenerate ES rows
(`e(es_ndegen_type)`) and withholds `e(V_es_ame*)` under the cell-containment caution, so a band needs the
same guards and would be unavailable exactly there.

**Verification.** Confirmed on the wooldid side by two independent verifiers (grep-exhaustive; the plot
line corrected to ado:10651; `bootdisplay(bca)` is now committed, not in-flight). Practice partly: csdid's
uniform default is conditional on `wboot`, and "nothing in the mainstream DiD toolkit prints only pointwise
bands under a bootstrap" is too broad.

#### PS-46 — no `e(V)` is posted under the `permute()` default, and `permanalyticse` restores an ANALYTIC one

**wooldid.** `permute()` without `permtstat` puts the observed run in fast (no-SE) mode (ado:111-114, 126),
so the SE column of the results matrix is missing, `diag(hadamard(se,se))` is all-missing, `ereturn post`
rejects it, and the coefficient-only fallback runs (ado:~4274-4308). Verified: `ereturn list` shows
`e(b) : 1 x 6` and no `e(V)`; the permutation quantities live in `e(perm_p)`, `e(perm_ci)`, `e(perm_sd)`,
`e(perm_se_shift)`, `e(perm_p_se)`, `e(perm_p_ci)`. The footer says "no e(V) or e(V_es_*) is posted" and
the help adds that the honestdid/pretrends handoff is unavailable (sthlp:417) and that "test, lincom,
esttab standard errors, and anything else needing a variance-covariance matrix will refuse" (sthlp:263-267).
Under `permanalyticse` or `permtstat` `e(V)` comes back but holds the ANALYTIC covariance — the footer
itself warns that `test`/`lincom` on it are analytic, not permutation, quantities.

**Practice.** Stata's `permute` is an r-class prefix that clears `e()`, so the estimation is run separately
and its own `e(b)/e(V)` exist. `ritest` is a prefix that runs the estimation itself and then RESTORES the
observed run's analytic `e(b)/e(V)` (ritest.ado:290-296, 652-653), so the standard `estadd scalar
ri_p = r(p)` / `esttab` workflow prints analytic standard errors and stars beside a randomization p-value.
That is the published-paper norm.

**Why it matters.** wooldid's `permute()` replaces the estimation run, so there is no companion `e()` to
export: a user has an estimation result they cannot post-estimate or export by any standard route, and the
analytic twin costs a second full estimation. **Correction to the finder's framing:** the
`permanalyticse` hazard the finding called "worse" is exactly ritest's own post-run state, so it is not a
deviation. The genuine deviation is the opposite one — the fast-mode default posts nothing, which is
unusual but errs conservative: it refuses to export rather than exporting analytic numbers under a
permutation banner.

**Harmonization.** (a) Post `r(table)` with the permutation p and interval in its p/ll/ul rows so tools
reading `r(table)` pick up RI numbers — cheap, and it creates no fake covariance. (b) Post
`e(V) = diag(e(perm_se_shift)^2)` with a loud note — this is effectively what ritest already gives users,
but it invites a hand `test` that is not the permutation test and the shift-implied SE is not a sampling
SE away from the null; it also touches the interval-implied-SE ruling of 09-03/09-04. (c) Cheapest and
safest: a help section with the estadd/esttab recipe, plus one caution under `permanalyticse` that `e(V)`
is analytic and export tools will show analytic inference. (c) is literally the ritest recipe.

**Verification.** Confirmed on the wooldid side (code trace plus a live `ereturn list`). Practice partly:
ritest's mechanism was corrected as above, which reverses which half of the finding is the real deviation.
The esttab/coefplot behaviour after a permute run was never executed by either lens.

#### PS-47 — no doubly-robust or IPW variant; controls enter only through the outcome regression

**wooldid.** Covariates enter the stage-1 regression via `controls()`, `coarsecohortcontrols()`,
`ccc_absorb()`, `interactivecontrols()` or `fe()`; the estimand is a margins average over the fitted
regression (sthlp:100-130, 868-895). "doubly robust", "propensity", "IPW" and "inverse probability" are 0
hits in both the ado and the help. There is no weighting-based or DR alternative and no propensity-side
diagnostic (overlap, balance).

**Practice.** csdid's default estimation method with covariates is `dripw`, the Sant'Anna & Zhao (2020)
doubly-robust IPW estimator, with `drimp`/`stdipw`/`ipw`/`reg` as alternatives (csdid.sthlp:80-99); it
falls back to `reg` only when no covariates are specified. drdid is the underlying 2×2 workhorse. Callaway
& Sant'Anna (2021) present DR as the recommended default with covariates, and Sant'Anna-Zhao is now the
modal citation for conditional parallel trends.

**Why it matters.** With conditional parallel trends a referee routinely asks whether the result survives a
doubly-robust specification, and a user comparing wooldid to csdid on the same data with covariates gets
different numbers for a reason the help never names. wooldid's help gives the user no language for
answering the question.

**Harmonization.** (a) Cheapest and probably sufficient: a paragraph in Controls and Fixed Effects stating
the position — the ETWFE design saturates cohort-by-period cells, so within a cell the outcome model is
essentially nonparametric and DR buys little when covariates are cell-constant — and naming the cases where
it does matter (covariates varying within cell, coarse i). It should be honest that
`interactivecontrols()` lets treatment effects vary linearly in the demeaned covariate, a functional-form
assumption on the outcome side precisely where DR would differ. (b) A genuine IPW/DR leg would need
cell-level propensity weights and a reweighted margins pass — a new engine.

**Verification.** Confirmed on both halves.

#### PS-48 — no `estat` aggregation surface, no calendar-time aggregate, and the cohort-effect matrix is not keyed by cohort

**wooldid.** Estimands are fixed at estimation time (`att`, `att_it`, `att_i`, `att_itime`); there is no
post-estimation re-aggregation and no `wooldid_estat.ado` anywhere, so `estat event` fails with "estat
event not valid", r(321). There is no calendar-period aggregate at all — `att_itime` is a cohort-level
average, and the only "calendar" object is the `caltime` column inside `e(gteffects)`. Replay is a bare
`wooldid` on a strictly empty command line and "accepts no options", so it is not a redisplay surface
either. The per-cohort effects from `histogramcohorteffects` come back in `e(histogramestimates)` with row
names such as `1.__00000B:1.__000005` (margins tempvar stripes), and the help concedes "Row names carry
internal temporary-variable identifiers, not cohort values; the histogram itself is the intended output"
(sthlp:1131-1133).

**Practice.** csdid offers `agg(simple|group|calendar|event)` at estimation and `csdid_estat`
(pretrend/simple/group/calendar/event/attgt) with `estore()`/`esave()` afterwards, and `saverif()` exists
precisely so other aggregations can be requested after the fit. jwdid_postestimation offers
simple/group/calendar/event with plot/esave/estore/other(). hdidregress has `estat aggregation`,
`estat atetplot`, `estat ptrends`. Reporting ATT-by-cohort and ATT-by-calendar-period from a single fit is
the standard Callaway-Sant'Anna workflow, and simple/dynamic/group/calendar is its shared vocabulary.

**Why it matters.** wooldid runs are expensive (runtime grows at about the 1.3rd-1.7th power of the
cohort-by-period cell count), so re-running the whole estimation to see a different aggregation is a real
cost. The one table that answers "which cohorts drive this?" cannot be joined back to cohort identifiers
without going to `gteffects` and redoing the aggregation by hand. And a reader who knows
simple/dynamic/group/calendar has no translation table: `att_i` is an equal-weight average over cohorts
whereas csdid's group aggregation weights cohorts by size — same name-space, different weights.

**Harmonization.** (a) Label `e(histogramestimates)` rows with the cohort values — the information is
already in `e(gteffects)`'s i column. Small, high value. (b) Add a calendar-time aggregate (an att-style
margins pass grouped on t, or a weighted read of `gteffects`). Moderate. (c) An `estat` menu that
re-aggregates `e(gteffects)` using its stored `w_att`/`w_att_it`/`w_att_i`/`w_att_itime` columns — the help
already promises those weights reproduce the aggregates exactly (sthlp:958-961), so the arithmetic is
settled and the work is inference on the re-aggregation. (d) A mapping table in the Comparison section
(att ↔ `estat simple`, att_i vs `group` with the weighting difference spelled out, ES rows ↔ `estat event`,
calendar not offered). Escape hatch that exists today: `verbose` stores the underlying regression as
`___wooldidfullmodel`, restorable with `estimates restore`, so a determined user can hand-roll a calendar
margins pass with no wooldid-side weighting or inference support.

**Verification.** Confirmed on both halves. One scope correction: jwdid's `over(varname)` is documented for
`simple` only; `other(varname)` subsets calendar/group/event.

#### PS-49 — switching the displayed bootstrap family costs a full re-run of the replications

**wooldid.** All families are computed and stored on every run (`e(boot_wald_mainresults)`,
`e(boot_pctl_mainresults)`, `e(boot_bc_mainresults)`, posted as hidden matrices at ado:8862-8865), but
`bootdisplay()` is an estimation-time option and replay is gated on a strictly empty command line
(ado:20-26; sthlp:1175, "accepts no options"). Typing `wooldid, bootdisplay(wald)` after a bootstrap fit
returns "Error: command line input should feature all of y i t ttre; too few variables specified." (r(198))
— the ordinary too-few-variables parse error at ado:4699, not a replay-specific one. To see the Wald or BC
package as a formatted table the user re-runs every replication; the numbers are otherwise reachable only
by matlisting the hidden matrices, whose names the footer and the Stored Results section do at least
disclose. The plots and the `saveplots()` datasets carry only the family displayed at estimation time
(sthlp:485) and have no such fallback.

**Practice.** Stata's whole postestimation model is that re-display is free: `estat bootstrap` prints
alternative families side by side after the fit at no cost, and every estimation command re-displays with
`level()` and formatting options on a bare replay. boottest re-uses a stored model.

**Why it matters.** A 1,000-replication run on a real panel is minutes to hours; "let me also look at the
BC interval" should not cost that again, and a user who tries the obvious thing is told their command line
is missing variables, which reads as a syntax error rather than "replay takes no options". The plots are
the genuinely irrecoverable part.

**Harmonization.** (a) Trivial: make the replay refusal say what it means when `e(cmd)==wooldid` and
options were supplied. (b) Moderate: a replay branch accepting `bootdisplay()` (and `level()` for the
quantile families) that re-renders from the stored family matrices, optionally `bootdisplay(all)` printing
the three side by side. The matrices exist; the work is the render path and deciding what `e()` looks like
afterwards. Note BCa could not be free either way — it needs the jackknife requested at fit time, exactly
as Stata's `estat bootstrap, bca` requires `bca` on the prefix.

**Verification.** Confirmed on the wooldid side. Practice partly: `estat bootstrap`'s default is `bc`, it
prints intervals not p-values, and its BCa carries the same fit-time restriction.

#### PS-50 — no `strata()` for the bootstrap resample: failed-because-no-treated-cluster replications are excluded after the fact rather than prevented

**wooldid.** The only resampling control is `cluster()`. Clusters are drawn with replacement from the whole
set, so a draw with no treated cluster fails and is excluded; the help states the consequence
("replications that draw no treated cluster fail and are excluded, so with few treated clusters the
retained distribution is conditional on drawing at least one", sthlp:485) and a caution fires below 6
treated resampling clusters (ado:6272-6279). `strata()` is not an option: "option strata() not allowed",
r(198). The permutation engine, by contrast, has `permstrata()`. Under `bootstrap(bayes)` no cluster ever
leaves the sample, so this failure mode does not arise there.

**Practice.** Stata's `bootstrap` prefix has `strata(varlist)` and `size(#)` precisely for this
(bootstrap.sthlp:163-177), so the standard fix — draw G1 treated and G0 control clusters separately — is
one option away; `ritest` also has `strata()`.

**Why it matters.** The conditioning-on-success problem the 2026-09-06 ruling chose to disclose with a
caution is, in the treated/untreated case, partly preventable rather than merely disclosable. A user who
knows `bootstrap, strata()` will look for the option, find nothing, and read the caution as "this tool
cannot do the standard thing."

**Harmonization.** Add `bootstrata(varname)` drawing G_s clusters with replacement within each stratum
(default off, so no baseline moves); the natural stratifier is treated-vs-never-treated, which could be
offered as `bootstrata(treated)` computed internally. Cost: the stratified draw plus the nesting check
`cluster()` already gets; accumulator, families and display untouched. Minimum alternative: one help
sentence saying stratified resampling is not offered and why. **Do not overclaim the benefit:** under
staggered adoption with per-horizon estimands a binary stratified draw does not guarantee every cohort
appears, so cohort-horizon cells can still be lost, and removing the conditioning event does not fix the
few-treated-cluster sizing problem.

**Verification.** Confirmed on the wooldid side. Practice partly: the Cameron-Miller attribution for
"stratified pairs resampling as the textbook remedy" could not be substantiated — what that literature
recommends for few treated clusters is the wild cluster bootstrap-t. The evidence cited for a
no-treated-cluster failure (1 of 100 excluded) came from a thin-panel long-horizon run and does not
isolate the mechanism.

#### PS-51 — only the G-1 t reference; no CR2/HC2-type correction and no few-cluster warning on the analytic path

**wooldid.** Every margins call is passed `df(e(df_r))` from the stage-1 fit (ado:2239), which under
reghdfe with `cluster()` is G-1; an archived footer confirms "referred to the t(11) distribution" on a
12-cluster fixture. There is no bias-corrected cluster-robust option: `customvce()` can pass an arbitrary
vce string to regress, but the default reghdfe path offers only cluster/robust. The few-cluster cautions
live inside the bootstrap controller only (ado:6259, 6278), so an analytic run with 6 clusters prints no
warning at all.

**Practice.** Stata's `didregress` family ships `vce(hc2 [clustvar])` and `vce(hc3 [clustvar])` implementing
the Bell-McCaffrey (2002) degrees-of-freedom adjustment with the Hansen (2025) scale option
(didregress.sthlp:340-395), precisely because G-1 t with few clusters is undersized. csdid offers `wboot`
as a first-class alternative. Cameron-Miller (2015) and MacKinnon-Webb (2017, 2018) treat a cluster count
below roughly 30-50 as a warning condition.

**Why it matters.** The risk the bootstrap cautions describe is identical on the analytic path, which is
the default path: a user with 8 clusters and `cluster(i)` gets a t(7) p-value and complete silence, while
the same user asking for `bootstrap()` gets two explicit cautions. And a reader used to Stata 18/19's
didregress will look for `vce(hc2)` and not find it.

**Harmonization.** (a) Move the <10-clusters and <6-treated-clusters cautions up to fire on any clustered
run, not only bootstrap ones (cheap). (b) Expose HC2/CR2 through the existing `customvce()` route by
documenting `customvce(hc2 clustvar)`, and measure whether the cell collapse stays exact there (the
collapse currently refuses `customvce()` outright). (c) Full parity — a native CR2 on the reghdfe path — is
a much larger job and probably not worth it given `permute()` and `bootstrap()` exist.

**Verification.** Not verified (Stata cap).

#### PS-52 — the only reference distribution is a label permutation of the observed onset vector; the placebo designs applied papers run are unreachable

**wooldid.** `permute(shuffle)` permutes the observed onset labels across clusters, holding the number of
clusters at each onset fixed (a per-stratum jumble of the observed onset column, ado:10986-10995; the exact
branch enumerates the same multiset). `permute(probabilities)` draws each cluster's onset from
user-supplied probabilities — but the `permprobs()` parser requires the listed level set to be a bijection
with the data's in-sample onset level set (ado:11196-11207), so an onset date no unit actually has is a
hard r(198).

**Practice.** wooldid's shuffle is exactly the MacKinnon-Webb (2020) / Athey-Imbens random-adoption-date
design, and also the Roth & Sant'Anna (2023) `staggered` randomization test (with `permtstat` the
studentized version) — a good default and a matched one. What is unreachable: Bertrand, Duflo and
Mullainathan (2004) generate placebo laws by drawing an adoption date uniformly from the sample period for
randomly chosen states — dates no state actually had, and a single shared date rather than a permutation of
the observed multiset. Conley and Taber (2011) build the reference from placebo assignments among the
never-treated units only — the construction wooldid's own help explicitly contrasts itself against ("This
is not the Conley-Taber construction").

**Why it matters.** A referee who asks "what happens if you assign placebo adoption dates at random over
the sample period, as in BDM?" cannot be answered inside wooldid; the user must build the design outside
the package and re-run by hand. And a reader who has just been told what wooldid is not doing
(Conley-Taber) has no way to ask for it.

**Harmonization.** (a) Smallest and highest-value: relax the `permprobs()` level check to allow levels
absent from the data (keeping the grid and never-treated checks), which buys the arbitrary-date half of the
BDM design; supporting evidence that this is safe — `wooldid_permxp` writes the level straight into the
assignment column, the zero-probability caution already anticipates a design under which the observed
assignment has probability zero, and `lvlstr` appears nowhere outside the parser. Unchecked: whether the
replication rewiring copes with a cohort whose onset date no unit has. It does not buy BDM's shared-date
structure. (b) Conley-Taber control-only reference was tabled as F035(b) and not chosen. (c) A
not-yet-treated conditional design is a new assignment generator, comparable in size to the deferred
`permute(dose)`.

**Verification.** Confirmed on the wooldid side. Practice partly: Roth-Sant'Anna (2023) was miscited as an
unreachable design when it is in fact the design wooldid implements, and the "not-yet-treated conditional
randomization" practice claim had no citable source — so the tally is one implemented, one genuinely
unreachable (BDM), one deliberately declined (Conley-Taber).

#### PS-53 — `permprobs()` must be typed one group per cluster, and there is no way to hand the engine a design from a file, matrix or program

**wooldid.** `permprobs()` is `(onsetvar level level ...)` followed by exactly one `(clustervalue p p ...)`
group per in-sample cluster value; the parser refuses any other count and each group must carry exactly K
probabilities (ado:11190-11235). `permstrata()` takes a single numeric variable. There is no option to read
assignments from a .dta, a matrix, or a user-written program. For a 50-state panel with 6 onset levels the
option string is 51 groups and about 350 tokens typed by hand.

**Practice.** ritest offers three routes: automatic resampling with `strata(varlist)`/`cluster(varlist)`;
`samplingsourcefile()`/`samplingmatchvar()`, where the user supplies a dataset of pre-computed
permutations; and `samplingprogram()`/`samplingprogramoptions()`, where the user writes the randomization
as a Stata program — its own worked example is a matched-pair design done that way. That is how people
implement rerandomization, matched-pair, covariate-constrained and two-stage designs.

**Why it matters.** The design most users want ("each cluster equally likely to get any onset, within
region") is already covered by shuffle + `permstrata()`, so the gap bites only for genuinely
non-exchangeable designs — but for those, `permprobs()` is unusable at realistic panel sizes and there is
no escape hatch. Blocking on two variables also requires an `egen group()` first, where ritest takes a
varlist.

**Harmonization.** (a) Accept `permprobs(matrix(name))` or `permprobs(using file)` with clusters in rows
and levels in columns — a small parser branch reusing the same validation. (b) Accept a varlist in
`permstrata()` (egen internally). (c) A `samplingsourcefile()` analogue accepting a .dta of assignment
columns keyed on the cluster value is the general escape hatch and would also make wooldid's RI
reproducible from an archived design file; cost is a new read path plus validation that every column is a
legal assignment.

**Verification.** Not verified (Stata cap).

#### PS-54 — no plot of the permutation distribution

**wooldid.** wooldid draws event-study plots, continuous-treatment plots and a histogram of cohort effects,
and under `permute()` redraws the ES and cfx plots with the shift-interval bands. Nothing plots the
permutation distribution itself; there is no kdensity or histogram call anywhere in the permutation path.
The draws are available — `e(permdist_b)` when completed replications are at most 10,000 (ado:224), and
`permsave()` writes them — but the user must draw the figure.

**Practice.** ritest has `kdensityplot`, producing a density of the replicated statistic with the observed
value as a vertical line. BDM (2004), Buchmueller-DiNardo-Valletta, Chetty et al. and essentially every RI
appendix show exactly this figure. It is how a reader judges whether the reference distribution is sane
(skew, bimodality, an observed value inside the bulk).

**Why it matters.** The footer already reports the null band and the SD of the null distribution — numbers
whose whole purpose is to describe a shape — while the shape itself, which practice reports as a figure, is
left to the user. It is also the cheapest available defence against the exchangeability problem the help
now warns about: a reference distribution distorted by an atypical treated cluster is usually visible in
the picture and invisible in the SD.

**Harmonization.** Add `makepermplot` (or fold it into `makeplots` under `permute()`): one histogram or
kdensity per pooled estimand from `e(permdist_b)` with an `xline` at the observed statistic, saved by
`saveplots()` like the others. Small — the matrix is already in e() and the plot paths exist. Should be
opt-in.

**Verification.** Not verified (Stata cap).

#### PS-55 — only the sharp null of zero can be tested; there is no non-zero null and so no exact inversion

**wooldid.** Every permutation p-value ranks the observed statistic against draws generated by reassigning
onset on the untouched outcome, i.e. tests tau = 0 (`wooldid_permtail`, ado:13030-13060). The interval
offered instead is the shift interval, which the help now says is "calibrated only near the null" and whose
width grows with the true effect. There is no null-value option anywhere in the permutation syntax.

**Practice.** ritest's `null(outcome value)` subtracts a hypothesized effect from the outcome before
permuting, so H0: tau = c can be tested directly; repeating over a grid of c is the standard exact-inversion
route to a Fisher randomization interval, and it is what Imbens-Rubin and Rosenbaum describe as the way to
get a confidence interval from a randomization test. Conley-Taber (2011) likewise deliver an interval by
inverting a placebo-based test rather than by shifting a null distribution.

**Why it matters.** A referee asking "can you reject an effect of 0.5?" or "what does randomization
inference give as an interval, not a shifted null band?" has no answer inside the package, and the interval
the package does print is the one it has just told the reader is only calibrated near the null. This is
the one place where the RI surface is weaker than a two-line ritest workflow.

**Harmonization.** (a) `permnull(#)` — subtract # × (post-treatment indicator) from the outcome before the
replications, exactly ritest's `null()`; one pre-loop transform reusing the whole engine, and it makes an
equivalence-style claim reportable ("we can reject an effect as large as x"). Care needed: defining the
shift for staggered designs and for the ame family (a per-horizon constant on the treated cells), and the
joint/adversarial statistics would need recentring or a refusal under `permnull()`. (b) The full grid
inversion on top of (a) was costed on 2026-09-05 and ruled to text (F035c).

**Verification.** Not verified (Stata cap).

#### PS-56 — treatment must be absorbing: no reversal, no switching on and off

**wooldid.** "The program does not currently accommodate cases where treatment is reversible — once
treated, cohorts must remain treated" (sthlp:104-105), and `ttre` must take a single value within each i
over the estimation sample, enforced with a hard error (sthlp:66). The restriction appears in the
Description prose and in the `ttre` entry, not in the option list or the Syntax block.

**Practice.** de Chaisemartin & D'Haultfoeuille's `did_multiplegt_dyn` is built for exactly this case
(treatments that switch on and off, non-binary, non-absorbing) and is one of the two most-cited modern DiD
packages. csdid, did_imputation and jwdid are also absorbing-only, so wooldid is in normal company — but
the case is common enough (licensing regimes, minimum-wage bindingness, policy repeals) that users check.

**Why it matters.** A user with a repealed policy discovers the limitation only after reading the
Description paragraph, or after hitting the onset-constancy error with a message about data integrity
rather than about design coverage.

**Harmonization.** Text only: state the restriction in the Syntax/`ttre` synopsis where users look first,
and name `did_multiplegt_dyn` as the tool for reversible or non-absorbing treatments (the Comparison
section already names hdidregress and did_imputation, so the precedent exists). Two sentences. Implementing
reversal is out of scope for this estimator family.

**Verification.** Not verified (Stata cap).

#### PS-57 — composition change across horizons is diagnosed but cannot be fixed: no balanced-cohort event study

**wooldid.** `e(escomp)` reports, per horizon, ncohorts / nunits / nclusters / nobs / sumw_att /
wshare_att / wshare_att_it, and the help says "Composition shifts across horizons are a common cause of
jumpy event-study profiles" (sthlp:1078-1090). The table is displayed only when `summarystats` is
specified. No option restricts the sample to cohorts observed over the whole reported window.

**Practice.** R's `did` exposes `balance_e` in `aggte()` for exactly this: restrict dynamic effects to
cohorts observed for at least e periods post-treatment. Callaway & Sant'Anna (2021) discuss it explicitly;
applied papers routinely report the unbalanced and balanced-composition event studies as a robustness pair.

**Why it matters.** wooldid gives the user the diagnostic and then no remedy. Once escomp shows ncohorts
falling from 2 to 1 across horizons, the natural next step is the balanced version — which the user must
build by hand with an if-restriction, silently changing the estimand and the sample-formation screens as
they go.

**Harmonization.** (a) A `balancede(#)` option keeping only cohorts observed across the reported window,
implemented as an internal sample restriction plus a note on the estimand change; moderate but
self-contained, since escomp already identifies the cohorts. (b) Cheaper: document the hand recipe beside
the escomp entry and print escomp by default on any event-study run rather than only under `summarystats`.

**Verification.** Not verified (Stata cap).

#### PS-58 — subgroup contrasts are reported pooled only, never per horizon

**wooldid.** `subgroups()` takes a single non-negative-integer variable whose levels must specify
collections of cohorts i; a cohort split across subgroups is a hard error, as is a subgroup varying within
an i-by-t cell or over time within a unit (sthlp:695-712). The pairwise comparison table reports contrasts
for the pooled `att` and `pre_att` only — no per-horizon rows and no `relyr` column — even when an
event-study window is requested. Per-subgroup event-study objects `e(b_es_*_s#)`/`e(V_es_*_s#)` ARE
returned, so the per-horizon contrast is computable but not reported.

**Practice.** jwdid_estat's `over()`/`other()` produce aggregations across subgroups off one fit. In applied
papers heterogeneity is normally shown as two event-study paths with the difference plotted or tested
horizon by horizon; a triple-difference is by construction a per-horizon contrast.

**Why it matters.** The pooled contrast answers "do the groups differ on average", which is rarely the
question after an event study has been shown. Users must build the per-horizon difference by hand from the
returned per-subgroup blocks, and the constant-within-i requirement forces a redefinition of i whenever
heterogeneity is at a finer level than the cohort.

**Harmonization.** (a) Report per-horizon subgroup differences: with the cross-subgroup covariance from the
same margins pass this is a lincom on a stacked block; moderate, and the display slot exists. (b) Cheaper:
document the hand recipe beside the `subgroups()` entry, and put the `egen newi = group(i sg)` remedy in
the option text rather than only in the error message.

**Verification.** Not verified (Stata cap).

#### PS-59 — triple-difference is offered as a recipe rather than a feature, with no DDD-specific inference or caveat

**wooldid.** The help gives the recipe `subgroups(sgvarname) fe(i.sgvarname##i.t)` and states "The
difference in SG level estimates returned will be the triple diff estimate" (sthlp:718-722). There is no
`ddd()` option, no DDD-specific diagnostics, and (per PS-58) no per-horizon DDD path — only the pooled
subgroup contrast. No example demonstrates the recipe.

**Practice.** DDD is ubiquitous in applied micro. The recent methodological literature
(Ortiz-Villavicencio & Sant'Anna 2025) shows that the naive difference of two DiDs is not generally valid
with covariates and staggered timing, and dedicated Stata implementations have appeared. Stata's
didregress/hdidregress do not offer DDD either, so wooldid is not alone — but a recipe stated without
caveat reads as an endorsement.

**Why it matters.** A user following the recipe gets a number and a p-value with no signal that the
estimand is the difference of two DiD estimands, that its validity with covariates and staggered onset is
an open methodological question, or that the reported contrast is pooled over horizons.

**Harmonization.** (a) Text: state that the recipe delivers the difference-of-DiDs estimand, name the
conditions under which that is the DDD parameter, and cite the recent DDD literature so a referee query has
an answer. (b) A first-class `ddd()` option with its own inference is a larger project and should wait on
the (a) framing.

**Verification.** Not verified (Stata cap).

#### PS-60 — no leave-one-out, influence or leverage diagnostics, and no jackknife

**wooldid.** "jackknife" appears nowhere in the help. `customvce()` forwards a `vce()` string to
regress/poisson but is ineligible for the cell collapse and incompatible with `bootstrap()`. Nothing
exposes per-cluster influence, leverage, or delete-one-cluster estimates. (A delete-one-cluster jackknife
leg now exists as the acceleration input for BCa, so the machinery has arrived for another purpose.)

**Practice.** Young (2019, QJE) made leverage and influence diagnostics for cluster-robust DiD a standard
referee ask. csdid can write out the influence functions for every ATT(g,t) and for the aggregation weights
with `saverif()` (csdid.sthlp:141-147), which is exactly what people use for leave-one-out checks.
`vce(jackknife)` is a one-word option on most Stata estimators, and MacKinnon, Nielsen & Webb (2023)
advocate the jackknife (CV3) CRVE with few clusters.

**Why it matters.** "Is the result driven by one state?" is a question every DiD paper answers, and wooldid
gives the user no mechanised way to answer it — despite investing heavily in inference machinery
elsewhere. It also leaves the few-cluster user without CV3.

**Harmonization.** (a) Expose the BCa jackknife leg's G estimates in `e()` (and optionally print the
largest-influence clusters) — the marginal cost over BCa is near zero and it delivers the Young-style
diagnostic. (b) Optionally a jackknife SE from the same G fits. (c) Cheapest interim: document that
`e(sample)` plus a loop over clusters is the supported route, and note `customvce()`'s limits.

**Verification.** Not verified (Stata cap).

#### PS-61 — requested horizons are silently truncated to what the data support, and there is no "all available horizons" shortcut

**wooldid.** The post window is set as `local espl = min(max_t, espostlength)` (ado:882, 1009) with no
message when the requested length exceeds what the data support; `esprelength` is treated analogously. Both
options take a plain integer with no `all` keyword, and both default to 0, so the default run reports no
horizons and hence no pre-treatment evidence at all.

**Practice.** did_imputation offers `allhorizons`, which "picks all non-negative horizons available in the
sample" (did_imputation.sthlp:148), alongside explicit `horizons()`; csdid estimates every ATT(g,t) and
`estat event` reports the full dynamic profile without the user naming a length; did_multiplegt_dyn takes
`effects(#)` but reports how many were estimable. did_imputation additionally refuses to report a
coefficient resting on too little treated mass by default (`minn(30)`).

**Why it matters.** A user who asks for 8 horizons and gets 5 has no signal in the output that three were
dropped — the table just has fewer rows, and the analysis described in the paper is not the analysis that
ran. Without an `allhorizons` equivalent, getting the full profile means guessing a number and checking
whether it was honoured.

**Harmonization.** (a) Print a caution when `min()` bites: "espostlength(8) requested; 5 post-treatment
horizons are estimable and are reported" — display only. (b) Accept a sentinel (`espostlength(all)` or a
negative value) meaning "all available". (c) Optionally a `minn()`-style screen flagging horizons resting
on very few treated units; the count is already in `e(escomp)`'s nunits column.

**Verification.** Not verified (Stata cap).

### 2.4 Unusual (but often deliberate) features

#### PS-62 — `e(V)` is a diagonal matrix on analytic runs and a full covariance under `bootstrap()`; `test` and `lincom` accept both without complaint

*Raised independently by the display, features and syntax lenses.*

**wooldid.** On the analytic path `e(V)` is built as `diag(hadamard(se,se))` — the reported standard errors
squared, every off-diagonal exactly zero (ado:4264/4266) — and posted at ado:4275. The help's
"Incompatible Post-Estimation Commands" paragraph says so plainly ("these are for display purposes only.
The variance-covariance matrix has all off-diagonal elements set to 0", sthlp:1186-1189) and directs users
to `e(b_es_type)`/`e(V_es_type)` "carrying genuine off-diagonal covariances, unlike e(V)" with an explicit
`ereturn post` recipe. The author's own comment above the post records the intent: results go in "without
giving any information about covariance (which we did not estimate anyway) ... you won't be able to test
estimates against each other". Because the matrix is diagonal, `test` and `lincom` necessarily run rc = 0
and return an SE built only from the two variances: after an analytic `att att_it cluster(i)` run,
`test main_att = main_att_it` returned F(1,79) = 0.00, p = 1.0000 and `lincom main_att - main_att_it`
reported Std. err. .0744 for a contrast identically zero on that fixture. There are **four** regimes, not
three: analytic-diagonal; bootstrap-full ("e(V) holds the bootstrap covariance", ado:10109);
permute-with-`permanalyticse`-or-`permtstat`, where `e(V)` is posted and holds the ANALYTIC covariance
(ado:10028-10029, with a footer saying `test`/`lincom` on it are analytic, not permutation, quantities);
and permute fast-mode, where no `e(V)` is posted at all. Nothing at the point of use distinguishes them —
though the bootstrap and permute paths each print a run footer saying what `e(V)` holds, so the analytic
path is the one with no line at the point of use.

**Practice.** Every comparator posts the sampling covariance that produced its own printed table: csdid
posts a full influence-function VCE (`ereturn post b1 s1`, csdid.ado:896) displayed through
`_coef_table`; did_imputation does `ereturn post b V` with a joint cluster-robust V and documents
`lincom tau_1-tau_0`; reghdfe and xtdidregress/didregress likewise. Stata's estimation-command contract
treats `e(V)` as the covariance, and `test`, `lincom`, `nlcom`, `margins`, `suest` and honestdid all assume
it. **Precedent worth recording:** jwdid and eventstudyinteract do not put their headline estimates in
`e(b)` at all (jwdid leaves reghdfe's regression `e(b)`/`e(V)` in place and produces every aggregate
through `margins` inside `jwdid_estat`), so keeping the aggregates outside `e(b)` has company — what has no
company is an `e(V)` that is knowingly false for the contents of the `e(b)` beside it.

**Why it matters.** Comparing two event-study horizons, averaging post-treatment coefficients, or testing
an `att` against an `att_it` is routine, and `test`/`lincom` are the reflex. They return plausible-looking
numbers rather than an error, and the same commands are Wald-correct after `bootstrap()`, analytic-but-
mislabelled under `permanalyticse`, and impossible under the permute default — so whether the answer is
valid depends on an option the table does not mention. A diagonal `e(V)` is silently wrong rather than
loudly absent, and the documentation is the only guard.

**Harmonization.** (a) Cheapest: print one line under the analytic table mirroring the bootstrap footer's
own disclosure ("e(V) is diagonal on this run; test and lincom across estimates are not valid — use
e(b_es_*)/e(V_es_*)"). The bootstrap and permute paths already print exactly this kind of footer, so the
analytic path is the odd one out within wooldid itself. (b) Post the full covariance for the main-table
rows — the margins calls that produce them already return one — which removes the deviation but is real work
and changes what esttab exports. (c) Make `test`/`lincom` refuse instead of answering: note that they gate
on `e(V)` existing, not on `e(properties)`, so the lever that works is posting `e(b)` with no V, which the
code already does as a degenerate fallback — at the cost of esttab's SE column, which is why the zeros
exist. Whether `ereturn post` tolerates missing OFF-diagonals is untested; the code only establishes that
it rejects a missing diagonal cell.

**Verification.** Confirmed on both halves by three verifiers, against the tree and the installed
comparator source. Corrections carried above: the four-regime count (the finding said permute posts no
`e(V)`, true only of fast mode); jwdid and eventstudyinteract do not post aggregates in `e(b)` either; on
the partial-identification rescue path `e(V)` can carry a negative sentinel or be omitted, so "every
off-diagonal is zero" holds wherever `e(V)` exists but `e(V)` is not guaranteed to exist.

#### PS-63 — refuses to run when every unit is eventually treated in sample

**wooldid.** Under default safety wooldid counts in-sample rows whose internal `ttre` is missing and,
finding none, prints "Error: you have no never-treated control group and no not-yet-treated control group
whose treatment falls outside the sample window." and exits r(459) (ado:1029-1039). The only bypass is
`safety(off)`. The help states it as a flat requirement (sthlp:66). Measured: recoding the fixture's
never-treated units to the last observed onset produces exactly that error; with `safety(off)` the pooled
default runs cleanly (att = .869, se .066, N = 880 with 80 obs dropped) and a full event study estimates
(pre2 -.013, pre1 -.079, post0 .842, post1 1.057). Under `esfixedbaseperiod` the same data fails with
"insufficient observations", r(2001) — so the refusal is an identification fact in the fixed-base regime
and a policy in the pooled one, exactly as the finding argued.

**Practice.** csdid's default is never-treated only, and "If there are no never treated observations,
notyet is used automatically" (csdid.sthlp:55-56) — an automatic fallback, not merely an available option.
jwdid's default IS not-yet-treated (jwdid.sthlp:27, 65). Stata's hdidregress/xthdidregress take
`controlgroup(notyet)`. did_imputation and did_multiplegt need no never-treated group. A panel in which
every unit is eventually treated within the window is the leading case in Callaway-Sant'Anna,
Borusyak-Jaravel-Spiess and de Chaisemartin-D'Haultfoeuille.

**Why it matters.** "All states adopt eventually" is one of the most common staggered designs in applied
micro. A user coming from csdid or did_imputation meets a hard error where those packages simply run, and
the suggested remedy — `safety(off)` — simultaneously disables every other screen in the package, including
the dropped-coefficient identification check that is wooldid's most valuable one (PS-66). The error text's
own vocabulary is a second stumble: wooldid uses "not-yet-treated" to mean treated-after-the-window, while
in the csdid/jwdid/CS sense every pre-onset row in this panel IS not-yet-treated — and the pooled regime
demonstrably uses those rows, which is why `safety(off)` works. A csdid user reads the message as asserting
something about their design that is false in their own terminology.

**Harmonization.** (a) Keep the refusal under `esfixedbaseperiod`; under the pooled default, downgrade it
to a caution whenever at least one treated cohort-period still has a not-yet-treated comparison after the
period screen (cost: one reachability check plus a battery config; the trimmed-window estimates are already
computed correctly). (b) Text-only fallback: have the error name which reference regime is binding, say
that the pooled regime may be estimable on a shortened window, and define its use of "not-yet-treated".

**Verification.** Confirmed on the wooldid side, including the regime split, by live runs. Practice partly:
csdid is stronger than claimed (automatic fallback); Stata's option is
`hdidregress ... controlgroup(notyet)`, not a bare `notyet`; did_multiplegt estimates switcher effects
rather than ATT(g,t).

#### PS-64 — `jointtests` is refused in exactly the configuration the help elsewhere recommends

**wooldid.** Under default safety, `jointtests` combined with analytic clustered SEs where i is nested
within the cluster variable is a hard error: "jointtests is not compatible with wooldid's analytic
clustered SEs when clustering at the level of the cohort variable i or above ... The joint test p-values
are severely oversized in this configuration" (ado:1520-1524, gate at 1498, nesting test at 1506-1518;
sthlp:582-590), downgraded by `safety(off)`. Meanwhile sthlp:64 advises "Simulations suggest the estimator
is most statistically efficient when i is set at the level at which you cluster" and sthlp:1008 advises
"Coarsening i is normally the single largest speedup available", with no cross-reference to the
restriction. The adversarial test is documented as degrading in the same design ("by an order of magnitude
when each cell lay entirely within a single cluster", sthlp:606) but is deliberately not gated; with
`contreat()` the ame-family SEs are set to missing there (sthlp:907-918). `bootstrap()` and `permute()` are
both exempt from the screen.

**Practice.** A joint Wald pre-test on pre-treatment coefficients with cluster-robust SEs is standard and
unrestricted elsewhere: did_imputation's `pretrends(k)` reports an F "corresponding to the cluster-robust
Wald test" (did_imputation.sthlp:188) with its default `cluster(i)`; `estat ptrends` after hdidregress
tests that all pretreatment periods are zero under its default `vce(cluster group)`; csdid reports
pre-treatment ATT(g,t) clustered at the panel level. **Qualification:** neither comparator is structurally
in the configuration wooldid screens — did_imputation's pretrend test is a separate lean unit/time-FE
regression on nontreated observations only, and `estat ptrends` runs off a non-saturated specification, so
the mechanism wooldid documents (cell dummies force each cell's residuals to sum to zero within its own
cluster) does not arise there. They demonstrate an unrestricted pre-trend test at the same user-facing
clustering choice, not in the same design.

**Why it matters.** The single most common wooldid configuration a user will land on — coarse cohort i,
`cluster(i)`, because both the syntax help and the runtime help point there — is the one where the package
refuses the pre-trend test, warns the line test is oversized by an order of magnitude, and (with contreat)
blanks the ame standard errors. Composing the help's own Example 0 (`cluster(i)`) with its Example 4
(`jointtests`) produces a hard error on the first try. A reader comparing against did_imputation or
`estat ptrends`, which print a pre-trend p happily at the same clustering, will read the refusal as a
wooldid limitation.

**Harmonization.** (a) Cross-reference from the i and runtime entries to the jointtests/adversarial
restriction, so the efficiency and speed advice carries its inference cost; optionally change Example 4 to
cluster finer than i so the two examples compose. (b) Offer the fixed-metric bootstrap or permutation joint
test as the automatic substitute when the analytic one is refused, rather than erroring — both engines are
already exempt, so this is "make an existing substitute automatic", not "build one". (c) State in the help
that the oversizing is a property of wooldid's saturated cell-dummy design under cluster-at-i — not a
general property of cluster-at-i Wald pre-tests, which is what the comparators' behaviour would otherwise
imply.

**Verification.** Confirmed on the wooldid side, verbatim. Practice partly: the "same design" claim was
corrected as above, which also weakens the finder's harmonization option (c) as originally worded.

#### PS-65 — `e(b)` holds only the reported aggregates, and the stage-1 coefficient vector is not `e(b)`

**wooldid.** After `wooldid y i t ttre` `e(b)` is 1×1 (`main_att`, .8445031) and grows only with requested
estimands and horizons; `e(properties)` is `b V`; `e(V)` is `.00401008` = se². The stage-1 vector is not
`e(b)` (it is retrievable only under `verbose`, as `___wooldidfullmodel`). Genuine covariances live in the
separate `e(b_es_type)`/`e(V_es_type)` matrices the user must `ereturn post` by hand.

**Practice.** did_imputation posts `e(b)` = estimates + pre-trend + control coefficients with the
corresponding V and documents `lincom tau_1-tau_0`; csdid posts a full influence-function VCV; xtdidregress
posts a full `e(V)`. But jwdid and eventstudyinteract keep their headline estimates out of `e(b)` too, so
wooldid's choice has precedent.

**Why it matters.** Posting `e(b)`/`e(V)` is a promise. The distinctive part is not that the aggregates sit
in a small `e(b)` — it is that `e(V)` beside them is deliberately false (PS-62), so the promise is kept in
form and broken in substance.

**Harmonization.** See PS-62; there is nothing to do here beyond that item, except to note in the help that
`e(b)` holds the reported aggregates only and where the stage-1 fit can be recovered.

**Verification.** Confirmed on the wooldid side (live run). Practice partly (jwdid/eventstudyinteract
miscited as full-e(b) posters).

#### PS-66 — `safety()` is a single global on/off switch, and several analyst choices are hard errors

**wooldid.** `safety(on|off)` governs every screen at once: the dropped-coefficient identification check,
the ttre-varies-within-i check, the fractional-offset check, the no-control-group refusal (PS-63), the
subgroup-crossing-cohorts check, the jointtests-under-cluster-at-i refusal (PS-64), the >100-distinct-i
caution and maxvar halt, and the negative `customweightmultiplier` check. There is no per-check bypass.
Separately, an explicit `fastcollapse` with an ineligible option is refused and "safety(off) does NOT
override that refusal", and several option pairs are hard incompatibilities (`bootstrap()` with
`impvalues`/`unconditionalse`/`customvce()`/multi-way `cluster()`; `permute()` with
`impvalues`/`contreatwithin`/`pretrendpower()`/`suppresssgprimaryeffects`/`bootstrap()`). Some refusals are
worded as judgements about the user's reasoning — `contreat()` without `contreatcontrols()` exits with
"you must specify this - substantive judgement is required." A note on volume: ordinary runs print no
cautions at all (att+cluster, four estimands, event study, impvalues, fastnose each printed zero), so this
is a register and granularity concern, not a noise concern.

**Practice.** Stata and the DiD packages overwhelmingly warn and proceed: reghdfe notes dropped singletons,
csdid prints notes, boottest warns on few clusters, did_imputation offers `autosample` to drop
un-imputable observations while recommending review, and `minn(30)` suppresses coefficients resting on too
few treated observations rather than refusing the run. csdid reports unidentifiable ATT(g,t) cells as
missing and continues. A mid-script hard stop on a documented-but-defensible specification is unusual; so
is an escape hatch that disables all data-integrity screens at once to bypass one statistical screen.

**Why it matters.** Each refusal is documented as a demonstrated wrong answer, which is arguably better
practice than the field norm. The unusual part is granularity: a user who wants the oversized joint-test
p-value for comparison must turn off safety entirely, which also disables the onset-constancy and
subgroup-integrity screens that have nothing to do with their choice — trading one known risk for several
unknown ones.

**Harmonization.** (a) Let `safety()` accept a list of check names (`safety(off: jointtests)`) so the
escape hatch is scoped; the checks are already named internally by their error text. (b) Cheaper: in each
hard-refusal message, state which other screens `safety(off)` also disables. (c) List the complete set of
safety-gated screens in one place in the help, and document the design philosophy once so the refusals read
as a stance rather than as brittleness. (d) For error voice: lead with the mechanical fact and put the
reasoning second, so the first line matches what Stata users scan for.

**Verification.** Not verified (Stata cap) for the enumeration; the two individual refusals it rests on
(PS-63, PS-64) are confirmed.

#### PS-67 — a Bayesian bootstrap as a first-class mode is unique among DiD packages

**wooldid.** `bootstrap(bayes)` draws Dirichlet weights per cluster and reweights rather than resampling,
with `bootalpha(equal)` — Rubin's Dirichlet(1,...,1) across clusters — as the default and `bootalpha(size)`
available for `bayesboot` compatibility. It is refused with `boottstat` and poisson (ado:6142-6144), and
refused for cell-weighted estimands under default safety because reweighting cannot resample cells
(ado:6213-6219). The help carries measured coverage (85.1% at nominal 90% on a 24-cluster balanced
simulation, against 89.4% for `bootstrap(standard)`) and says plainly it is not a remedy for few clusters.
It is one of exactly two modes, and Example 14 features it.

**Practice.** The Bayesian bootstrap is Rubin (1981) and well established, but in Stata its only homes are
the `bayesboot` community command and Stata 18's `rwgen bayes` used with the `bootstrap` prefix. No DiD
package offers it: csdid, did_imputation, jwdid, didregress and did_multiplegt all have exactly one
bootstrap flavour or none. Applied DiD papers essentially never report Bayesian-bootstrap standard errors.
The options applied papers do use are the pairs cluster bootstrap, the wild cluster bootstrap-t and
randomization inference.

**Why it matters.** Purely a "this will look unusual" item. A referee seeing "standard errors from a
Bayesian bootstrap" on a DiD paper will ask why, and the honest answer — no cluster ever leaves the sample,
so replication estimability is stable and the comparison group is never lost to an unlucky resample — is
stated in the example rather than where the user chooses a mode. There is also a mild risk that a user
reaches for it because it is the exotic-sounding option in a few-cluster setting, which is the one thing
the help says it will not fix.

**Harmonization.** Keep it; it costs nothing to users who ignore it and has a documented rwgen/bayesboot
cross-check. Present it explicitly as the estimability-stability fallback rather than as a co-equal
alternative, and repeat the measured coverage comparison at the point of choice, immediately before the
"not a remedy for few clusters" line so the two are read together. One or two sentences.

**Verification.** Not verified (Stata cap).

#### PS-68 — automatic re-identification, per-estimate availability accounting, and a posted draws matrix (more careful than practice)

**wooldid.** Three places where wooldid does more than the ecosystem. (1) It chooses the resampling
mechanics from the estimand: frequency reweighting for att/ame (documented as matching Stata's prefix
replication for replication under a common seed), and a physical resample with re-identification of
duplicated clusters for the cell-weighted estimands and for `boottstat` with `cluster()`, refusing with
r(459) when `cluster()` does not nest i. (2) It cautions on any failed replication at all — "1 of 100
bootstrap replications failed to estimate and were excluded (official-bootstrap convention: complete
replications only)" — and posts `e(boot_failed)` plus per-estimate availability in `e(boot_navail)`. (3) It
posts the replication matrix itself in `e(bootdist)` up to 10,000 replications.

**Practice.** Stata makes re-identification the user's job: `cluster()` plus `idcluster(newvar)`, with the
help noting that omitting `idcluster()` gives wrong results for panel commands — one of the commonest
cluster-bootstrap errors in applied work. Stata records incompleteness only as `e(N_misreps)` with no
per-parameter breakdown, and posts no replication matrix at all (draws are available only through
`saving()`).

**Why it matters.** These are the reasons to prefer wooldid's engine to a hand-rolled `bootstrap` prefix
call, and the automatic leg selection in particular removes a footgun that has produced published wrong
standard errors elsewhere. **Keep.**

**Harmonization.** Keep all three. One small residual: the leg actually taken is a shell local
(`wdb_reid`, ado:354) and is never recorded in `e()` or announced, even though the help says
re-identification "roughly doubles the cost of each replication". Posting `e(boot_leg)`
(reweight | reid | physical) and naming it on the dots line would let a user understand a run that suddenly
takes twice as long, and make the choice auditable. A few lines.

**Verification.** Not verified (Stata cap).

#### PS-69 — a studentized bootstrap of the aggregated DiD functional exists at all, with a five-condition refusal set

**wooldid.** `boottstat` bootstraps t(r) = (b(r) - b)/se(r) for every coefficient in `e(b)`, for the
gteffects cells, for pairwise subgroup contrasts and for the joint, equality and adversarial statistics,
with each replication computing its own analytic standard error. It is refused with poisson, with
`bootstrap(bayes)`, with `bootdisplay(wald)`, with `bootdisplay(bc)`, and with a `cluster()` variable that
does not nest i — five separate refusals, each with its own paragraph-length message and stated reason
(ado:6138-6144, 6160-6170, 6242). The help closes with an unusually candid disclaimer: "What the help does
not claim: ... whether the resulting intervals attain nominal coverage has not been measured."

**Practice.** No Stata DiD package offers a studentized bootstrap of the aggregated ATT. boottest gives a
wild bootstrap-t of a linear restriction on a stored regression, a different object (and the assessed route
through boottest on the stage-1 model was ruled infeasible). Cameron-Miller (2015) and Davison-Hinkley
(1997) both single out the bootstrap-t as the best-performing pairs-bootstrap variant, so the feature is
well motivated even though nothing else implements it here.

**Why it matters.** Two opposite reactions are both plausible. A reader who knows the bootstrap-t
literature will regard this as the best thing in the option set. A reader who simply tries `boottstat` may
hit two or three refusals in a row — poisson and bayes users are turned away outright, and anyone who then
tries `bootdisplay(bc)` or `(wald)` is turned away twice more — and conclude the feature is fragile rather
than carefully scoped. See PS-11 for its one uncautioned failure mode.

**Harmonization.** Nothing to change in behaviour. One "when boottstat is available" block in the help,
listing the five conditions with one line each, would let a user check compatibility before running rather
than discovering it one error at a time. The candid uncoverage disclaimer should stay; the F058 calibration
simulations are what would let it be replaced with a number.

**Verification.** Not verified (Stata cap).

#### PS-70 — exact enumeration is capped at 10,000 assignments

**wooldid.** `wooldid_permx` sets capexact = 10000 and capwarn = 1000 (ado:10972-10981): above 1,000
assignments a caution prints, above 10,000 `permexact` is refused with "permexact requires enumerating
every assignment, but there are <M> (ceiling: 10000). Use permreps()." Run-verified on the 80-cluster
fixture: refused with "there are 1.1214e+36". Where it does run, the exact p has a hard floor of 2/M under
the default convention (1/M under `permabsp`) — with one treated cluster among twelve, the exact default p
cannot fall below .167.

**Practice.** Exact randomization tests are usually described as feasible whenever the assignment count is
in the thousands to millions, because the per-assignment cost is a cheap statistic. Here each assignment is
a full wooldid refit through margins, so 10,000 is already a very long run — on the sweep's own timing, 99
sampled replications took about six minutes, so 10,000 would be roughly ten hours. The cap is honest about
the real constraint rather than about the combinatorics.

**Why it matters.** A user who has read that randomization inference is "exact" will reach for `permexact`
and, on any realistic panel, be told to use sampling instead — correct, but the opposite of the impression
the literature gives. And in the small designs where it does run, the exact floors are coarse enough that a
reader may see a minimum p of .167 and think the command has failed.

**Harmonization.** Keep the cap; name the reason in the error ("each assignment is a full re-estimation; at
the observed replication cost this design would take about N hours"), which the new runtime-projection
machinery (F066) has the ingredients for. Raising the cap without a cost estimate would trade a clean
refusal for a run the user abandons.

**Verification.** Not verified (Stata cap) beyond the refusal itself, which the permute lens ran.

#### PS-71 — randomization-design diagnostics and refusals no comparable RI tool offers (more careful than practice)

**wooldid.** The possible-assignment count M is always computed and reported, with a caution when
`permreps()` exceeds 10 percent of it and a refusal when the design admits only one assignment
(ado:10949-10981). Replications that fail to estimate are excluded and counted in `e(perm_failed)`, with a
caution above 10 percent and a refusal below max(20, half) completed. A degenerate permutation distribution
(every draw equal) draws a caution; a distribution whose central band excludes zero has its shift interval
blanked and named rather than printed. The sharp-null pre-period caveat states plainly that pre-period
permutation p-values rise toward 1 as the true effect grows. Every coefficient, horizon, subgroup level and
pairwise contrast gets its permutation p from the same replications at no extra cost, and `permsave()` plus
`e(perm_rngstate)` make the surface reconstructable and replayable.

**Practice.** ritest reports none of these: no assignment count, no single-assignment refusal, no
completion floor (a run with two successful replications reports a p), no degeneracy or off-centre
detection, and its RI p applies only to the expressions the user listed. Stata's `permute` is similar. The
MacKinnon-Webb and Young literature discusses several of these hazards but no Stata implementation checks
them.

**Why it matters.** These catch errors that would otherwise be invisible — a stratification that leaves
nothing to randomize, a design where 30 percent of draws are inestimable, a quantity whose null
distribution is not centred at zero so the shift interval would not even contain the estimate. Any move
toward "looking more like ritest" should not cost them. **Keep.**

**Harmonization.** Keep, and surface them: the assignment count and completion counts belong in the compact
footer if the footer is ever trimmed (PS-9), and the single-assignment and completion-floor refusals are
worth naming in the help's opening paragraph so users meet them as design guidance rather than as errors.

**Verification.** Not verified (Stata cap) beyond the header text and the `permexact` refusal, both of
which the permute lens ran.

#### PS-72 — the cohort-by-period cells carry no inference at all under `permute()`

**wooldid.** `gteffects` cells are keyed by cohort onset, which is exactly what a replication reassigns, so
they have no permutation distribution: their standard error, p and interval columns are blank and a note
says "under permute() these cells have no permutation distribution (their keys are cohort onsets, which the
replications reassign), so their inference columns are blank". The point estimates remain.

**Practice.** csdid reports ATT(g,t) with standard errors and (with `wboot`) bootstrapped simultaneous
bands; did_imputation and eventstudyinteract report inference on every reported coefficient. A table of
estimates with empty inference columns is not a form these users have seen.

**Why it matters.** The blanking is correct — a (g,t) cell is not a well-defined quantity across
reassignments, so there is no honest permutation distribution to report, and filling the columns from the
analytic run would violate the engine's own coverage-by-subtraction rule. But a user whose paper reports
the ATT(g,t) grid and who turns on `permute()` to answer a few-treated-clusters referee gets a grid with no
inference, arriving as blank columns rather than as guidance.

**Harmonization.** Keep the blanking; make the note prescriptive rather than descriptive — say that the
relative-time event-study rows are the permutable version of the same information and that a user wanting
RI on the (g,t) grid should report the event-study horizons instead. Text only. **Keep the behaviour.**

**Verification.** Not verified (Stata cap).

#### PS-73 — a self-updating `update` option that pulls from GitHub, and an `oldsyntax` option that changes no syntax

**wooldid.** `wooldid, update` runs `net install wooldid, from(https://raw.githubusercontent.com/...)
replace` and exits (ado:509-513) — a network write triggered by an option on an estimation command, placed
before validation, with no version comparison and no confirmation, so it reinstalls unconditionally. The
option is exempted from the replay path. There is no `wooldid, version`. `oldsyntax` does not change syntax
at all: it reverts the event-study coefficient suffix for relative period 0 from `post0` back to `contemp`
(ado:3800/3812/3825/3837, 4092/4178) — the entirety of what it does.

**Practice.** Stata's package-update mechanism is `adoupdate`/`net install` typed by the user; none of
csdid, did_imputation, jwdid, reghdfe, boottest or ritest ships a self-update option. Where a package
exposes its version it is via a `version` subcommand or option (csdid declares one).

**Why it matters.** Applied-micro replication packages must state which version produced the results. A
command that can silently overwrite itself mid-project, outside adoupdate's bookkeeping and with no version
recorded in `e()`, makes that harder — and the help elsewhere goes to considerable lengths about
bit-identical reproduction (`fullprecisionslow`, rngstate), so the two stances sit oddly together.
`oldsyntax` is a plain misnomer: a reader of the option list will expect it to restore an older command
syntax, not a coefficient label.

**Harmonization.** (a) Record the installed version (and, if available, the commit) in `e()` and in the
printed header, so a log identifies its own vintage; add `wooldid, version`. (b) Have `update` check the
installed version against the remote first and report "already current" instead of reinstalling, and point
at `net install ..., replace` as the pinnable route. (c) Add an accurately named alias for `oldsyntax`
(e.g. `oldcoefnames`), keep the old name working, and reword its help entry's first clause. Keeping the
`update` option itself is fine.

**Verification.** Not verified (Stata cap).

#### PS-74 — estimation options mutate global session state: `cleanmatrices` wipes every Stata matrix and the whole Mata workspace

**wooldid.** `cleanmatrices` executes `clear mata` and `matrix drop _all` (ado:542-545), destroying every
matrix and Mata object the user has in memory, not just wooldid's own; the help describes it as clearing
"all matrices storied in STATA's main memory and in mata". Separately, `level()` is implemented by
`set level` for the duration of the call with the caller's value restored through a global (ado:516-521
plus the shell epilogue); `___wooldidfullmodel`/`___wooldidrepmodel` are reserved stored-estimate names
overwritten on every run; and `WOOLDID_*` global macro names are reserved.

**Practice.** No comparator touches user matrices or Mata: reghdfe, csdid, did_imputation, boottest and
ritest all confine themselves to tempnames. Stata's convention is that `set` values changed inside a
program are program-scoped and restored automatically, and that estimation commands never drop user
objects. Where a package needs a reset it ships a separate command (e.g. `reghdfe, cache(clear)`), not an
option on the estimation call.

**Why it matters.** `cleanmatrices` is documented and is a genuine escape hatch for reghdfe's Mata wedging,
but its name sounds like housekeeping and its blast radius is the user's whole workspace — a user with a
hand-built matrix or a Mata session loses it with no confirmation. The `set level` route is a smaller
instance: a break or a hard crash mid-call can leave the session's level changed.

**Harmonization.** For `cleanmatrices`: gate it behind a confirmation or narrow it to wooldid's own
leavings (drop matrices matching wooldid's name patterns; `mata drop` the reghdfe objects), and rename it
to signal scope (e.g. `resetworkspace`). For `level()`: pass the level into the display routines instead of
`set level`, or keep `set level` and note in the help that a break can leave it altered. Of the
three, the `cleanmatrices` narrowing is the narrowest in scope.

**Verification.** Not verified (Stata cap).

#### PS-75 — every event-study run prints a pooled "pretreatment att" row beside the main ATT

**wooldid.** Whenever pre-treatment horizons are requested, the main table carries an equation block
`pretreatment` with a single `att` row — e.g. "pretreatment | att | .0143048 (se .0531374, p .7884753)"
directly under "main | att | .8554682". It is a pooled average of the pre-onset placebo effects, presented
in the same table with the same column layout as the treatment effect itself. The help mentions
"pretreatment" once, inside the `boottstat` paragraph; nowhere does it say what the row averages over,
whether it respects `esprelength`, or how it relates to the joint test.

**Practice.** csdid reports per-horizon pre-treatment ATTs under `agg(event)` (labelled Pre_avg/Post_avg in
its aggregation table) and a separate pre-trend test; did_imputation reports pre1..preK and a joint Wald
test; Sun-Abraham and `event_plot` workflows show leads individually. None collapses the pre-period into a
single pooled placebo estimate presented alongside the ATT with no label distinguishing them.

**Why it matters.** A reader scanning the table sees two rows both labelled `att` with very different
numbers and no header explaining that one is a placebo. It risks being read as a pre-trend test in its own
right — a small pooled pre-effect can hide offsetting leads, precisely the failure mode the per-horizon
display and `jointtests` exist to prevent. Under `permute()` the row additionally carries the sharp-null
caveat that its p rises toward 1 with the true effect (sthlp:421), easy to miss in the main table.

**Harmonization.** Keep computing it; relabel so it cannot be misread — head the block "pre-treatment
(placebo) average" and name the row `att_pre` rather than a second `att`. Add one sentence to the Event
Study section saying what the row is, that it is a summary and not a substitute for the per-horizon rows or
`jointtests`, and (under permute) that its p is a sharp-null quantity. See also PS-96.

**Verification.** Not verified (Stata cap).

### 2.5 Syntax

#### PS-76 — the positional varlist matches did_imputation, but `i` means cohort (not unit) and `ttre` is the onset date (not time-to-treatment)

*Raised by the syntax, estimator and features lenses.*

**wooldid.** `wooldid y i t ttre` (sthlp:59). The second variable is the "Cohort identifier; must consist of
integers greater than 0. Cohorts can be fine (i = individual units) or coarse (i = collections of
individual units), but should be no coarser than to consist of treatment onset cohorts" (sthlp:62), and the
help advises setting it "at the level at which you cluster"; unit fixed effects go separately into `fe()`.
The fourth is the calendar date of treatment onset (sthlp:64). Both must be numeric. There are no named
aliases (`ivar()`, `gvar()`, `onset()`) anywhere. The same variable defines the absorbed cohort FE, the
(i,t) treatment cells that carry the coefficients, and the cell dimension the estimands average over — and
the choice of i changes what is refused (PS-64) and what the gteffects table means (unit-by-period cells
rather than ATT(g,t)). Default safety cautions once i takes more than 100 distinct values, advising the
user to move unit variation into `fe()`, and runtime grows at roughly the 1.3rd-1.7th power of the
cohort-by-period cell count built from i.

**Practice.** did_imputation uses the same four-variable positional form `Y i t Ei` but `i` is
unambiguously the unit id and `Ei` the onset (did_imputation.sthlp:25-28). csdid and jwdid name them:
`ivar()` (unit) and `gvar()` (onset cohort) are separate options, so the distinction is explicit in the
call; jwdid even offers `group` to switch the FE from unit to cohort as an explicit performance choice.
xtevent takes `panelvar()`. "Time to treatment"/"event time" in the literature means t - g, not g.

**Why it matters.** Two stumbles at once. Positional slot 2 looks like did_imputation's unit id, and
passing a unit id there is legal but changes the estimator's granularity, the cluster advice, the runtime
and which features are refused. And `ttre` reads as "time relative to treatment" to anyone fluent in
event-study vocabulary; a user who supplies t - g gets a silently different design. The guidance that
resolves it lives under Runtime and under `safety()`, not next to the argument.

**Harmonization.** (a) Text: lead the `i` entry with "cohort/group identifier — this is NOT the panel id;
put unit fixed effects in fe()", name the csdid `ivar`/`gvar` correspondence, and list the consequences of
each choice (cells, runtime, jointtests availability, escomp interpretation). (b) Accept optional named
aliases (`ivar()`, `gvar()`/`onset()`) alongside the positional form; rename `ttre` in the docs to "onset"
while keeping the token. (c) Optionally an `unit()` argument passed to `fe()` automatically when i is
coarse. Renaming the positional slots would break every existing script.

**Verification.** Partly. The core — colliding positional slots with different semantics, no aliases, i's
role in what is refused — is confirmed. Two corrections to the finder's supporting claims: `egen group()`
produces 1-based ids so the "rejects a 0-based id" example is not a real scenario; and wooldid does not in
fact enforce the "integers greater than 0" rule at parse time (no sign check exists), while internally
`_i = 0` is the sentinel for "outside the treated window" (ado:878-880) — so a 0-valued cohort id is an
undiagnosed collision hazard rather than a rejection, which arguably strengthens the case for an explicit
parse-time error or an internal recode.

#### PS-77 — `[in]` is not accepted

**wooldid.** The syntax line declares `[if/]` only (ado:506), so `wooldid y i t ttre in 1/50, att` exits
with r(101) "in range not allowed". The help's syntax diagram likewise shows `[if]` with no `[in]`.

**Practice.** `[if] [in]` together is the [U] 11 contract for essentially every Stata estimation command,
and all three DiD comparators take both: csdid.ado:38 `varlist(fv) [if] [in] [iw]`, jwdid.ado:53
`varlist(fv ts) [if] [in] [pw]`, did_imputation.ado:8 `varlist(min=4 max=4) [if] [in] [aw iw]`.

**Why it matters.** `in` is how people take a quick slice while debugging or timing a specification;
hitting r(101) reads as a bug, and the error text is Stata's generic one, so nothing points at wooldid.

**Harmonization.** Adding `[in]` means threading a row range into the `ifvar` construction (ado:603-616)
and into the `ifr()`/`ifv()` plumbing carried into `wooldid_permctl`, `wooldid_bootctl` and the cfx paths
(ado:5475, 5493, 5724, 6012), so it is more than a one-word change. A cheap alternative is a targeted
error: catch `in` at parse time and say "wooldid does not support in; use if _n <= 50".

**Verification.** Confirmed by run in the features lens (r(101)); comparator syntax lines verified.

#### PS-78 — no `vce()`; the inference interface is split across five options

**wooldid.** `wooldid y i t ttre, att vce(cluster i)` exits r(198) "option vce() not allowed". The VCE is
chosen among `cluster()`, `robust`, `customvce(string)` (which passes text straight into regress/poisson's
`vce()` and is refused alongside `cluster()`/`robust`, ado:4754-4755), and the resampling engines
`bootstrap()` and `permute()` (ado:506). Which underlying estimator runs also depends on the choice:
`cluster()`/`robust` route through reghdfe/ppmlhdfe, `customvce()` and `unconditionalse` force
regress/poisson.

**Practice.** Stata's standard is one `vce()` option covering ols/robust/cluster/bootstrap/jackknife, and
hdidregress/xtdidregress follow it; reghdfe accepts both `vce()` and the legacy `cluster()`. did_imputation
and csdid use a bare `cluster()` like wooldid, so the split is not unique — but wooldid is alone in having
five entry points to the same decision.

**Why it matters.** `vce(cluster x)` is muscle memory. The r(198) is Stata's generic "not allowed" with no
hint that `cluster()` is the name to use, so the user must open the help. The five-way split also means the
mutual-exclusion rules have to be learned rather than inferred.

**Harmonization.** Accept `vce()` as a thin alias that dispatches: `vce(robust)` → `robust`,
`vce(cluster x [y])` → `cluster(x y)`, `vce(bootstrap[, reps()])` → `bootstrap(standard) reps()`, anything
else → `customvce()`. About 15 lines in `wooldid_validate` plus the byte-identical copy at ado:4667. Even a
bespoke error naming `cluster()` would remove most of the friction.

**Verification.** Confirmed by run in the syntax lens (r(198)).

#### PS-79 — 83 options in one command, 25 of them with no minimum abbreviation, in three naming styles

**wooldid.** The syntax line declares 83 options. Twenty-five carry no capitalized prefix and so must be
typed in full, including `fullprecisionslow`, `emptycellsoverride`, `altmethod_adversarial`,
`makecfxplotsbysg`, `cfxplotnose`, `cfxlattice`, `fastcollapse`, `customvce`, `oldsyntax`, `safety`,
`regtol`, `rseed`, `reps`, `dots`, `ties`, `att_it`, `att_i`, `att_itime`, `att`, `fe`. Names mix
run-together (`contreatelasticitytype`, `suppresssgprimaryeffects`), underscored (`att_it`, `ccc_absorb`,
`altmethod_adversarial`) and prefixed-family (`perm*`, `boot*`, `cfx*`, `contreat*`) styles. Longest:
`SUPPRESSsgprimaryeffects` (24), `HISTogramcohorteffects` (22), `altmethod_adversarial` (21, no
abbreviation).

**Practice.** did_imputation declares 28 options; csdid's live syntax about 21; jwdid about 13; boottest
and reghdfe are the large ones in this ecosystem at roughly 40-50. Stata's style is to give every option of
more than a few characters a documented minimum abbreviation.

**Why it matters.** Size alone is defensible — wooldid does more than its comparators — but the abbreviation
gaps are pure friction: `fullprecisionslow` and `emptycellsoverride` are long enough that a typo is likely
and the resulting message is Stata's bare "option ... not allowed". Mixed styles make guessing a name
harder than it needs to be.

**Harmonization.** Adding minimum abbreviations to the 25 bare options is a one-line edit to the syntax at
ado:506 and its byte-identical twin at 4667 (guarded by `test/check_validate_sync.py`), plus the help's
synopt lines; it is purely additive since existing full spellings keep working. Style unification would
break scripts and is not worth it.

**Verification.** Not verified (Stata cap); counts measured from ado:506.

#### PS-80 — the `perm*`/`boot*` namespace renames settings that ritest, boottest and Stata's own prefixes already have names for

**wooldid.** Inference is configured through 11 `perm*` options (`permreps`, `permexact`, `permstrata`,
`permprobs`, `permtstat`, `permtci`, `permsave`, `permanalyticse`, `permabsp`, `permside`, `permrawp`) and
8 `boot*`-family options (`reps`, `ties`, `bootalpha`, `bootdf`, `boottstat`, `bootdisplay`, `bootsave`,
`rseed`). `reps()` serves `bootstrap()` only and is refused as a stray under `permute()` with a message
naming `bootstrap()` but not `permreps()`; `seed()` is not an option at all (`rseed()` is), so it fails
with Stata's generic r(101). `nodots` is not accepted — dot suppression is `dots(0)`. The p-value convention
is selected by `permabsp`/`permrawp`; the reference distribution by `bootdf(normal|clusters|#)`; the
interval method by `bootdisplay(percentile|bc|bca|wald)`.

**Practice.** Stata's `bootstrap`/`permute`/`jackknife`/`simulate` prefixes and ritest use `reps()`,
`seed()`, `saving()`, `strata()`, `cluster()`+`idcluster()`, `nodots`, `left`/`right`, `bca`. boottest uses
`reps()`, `seed()`, `ptype(symmetric|equaltail|lower|upper)`, `weighttype()`, `boottype()`. Stata's
bootstrap postestimation selects intervals with `estat bootstrap, percentile bc bca normal`. `rseed()`
itself is not unusual (csdid's wboot and hdidregress's `wildbootstrap(rseed())` use it). So the concepts
wooldid renames — equal-tailed vs symmetric p, interval method, stratification, replication count — all have
prior names in the two packages a wooldid user is most likely to have installed.

**Why it matters.** A user who knows boottest looks for `ptype()` and finds `permabsp`/`permrawp`; a user
who knows `estat bootstrap` looks for `bca` and finds `bootdisplay()`; the first thing a ritest user types
is `reps()` and `seed()`, and both fail. The `reps()` failure is the worst because the error tells them to
add `bootstrap()`, the opposite of what they want, while listing the companions `permute()` does accept.
Because the engines are mutually exclusive, `reps()` could have served both. The help does the translation
once, well, for the p-convention (naming boottest's `ptype(equaltail)` and ritest's c/n rule) and nowhere
else.

**Harmonization.** (a) One clause in the stray-companion message: "reps() belongs to bootstrap(); the
permutation count is permreps()." Free, and it captures most of the value. (b) A cross-reference table in
the Standard Errors section — [wooldid option | ritest/permute equivalent | boottest/estat bootstrap
equivalent] — covering the p-convention, stratification, interval method, reference distribution, seed and
save options. (c) Accept the familiar names as synonyms where they do not collide: `nodots` for `dots(0)`,
`strata()` for `permstrata()`, `saving()` for `bootsave()`/`permsave()`, `seed()` for `rseed()`, `reps()`
for `permreps()` under `permute()`. Purely additive, ~10 lines in validate. Keep the `perm*`/`boot*` names
as the documented primaries.

**Verification.** Confirmed for the two refusals the bootstrap lens ran (`nodots` → r(198) generic;
`permute` + `reps()` → a message naming `bootstrap()`); the rest source-verified, comparator defaults read
off the installed help files.

#### PS-81 — the rank counts compare doubles exactly, with no tie tolerance

**wooldid.** `wooldid_permtail` counts m0 = sum(v :<= sobs), m1 = sum(v :>= sobs) and
ma = sum(abs(v) :>= abs(sobs)) on raw doubles (ado:13039-13046). There is no tolerance option and no
strict/non-strict choice. Under `permexact` the enumerated assignment identical to the observed one has to
reproduce the rep-0 comparator's statistic bitwise for the observed value to count as tied with itself, and
the comparator is a separately computed fast-mode refit.

**Practice.** ritest defaults to `eps(1e-7)` and defines its count as c = #{|T| >= |T(obs)| - eps}
explicitly "to account for numeric imprecisions", and offers `strict` (with `eps(0)`) for the
anti-conservative variant. The conservative-versus-strict choice under ties is a standard part of the
Fisher/Rosenbaum apparatus.

**Why it matters.** Discrete outcomes, small samples and symmetric designs all produce statistics that are
equal in exact arithmetic but not in floating point after two separate margins evaluations. A near-tie
counted the wrong way moves a small-M exact p by a whole 1/M step, and in exact mode the self-match is the
difference between the 2/M floor and something below it. Small in the cases wooldid mostly runs, but it is
the one part of the p-value machinery where the package is less defensive than the tool it cites.

**Harmonization.** (a) A tolerance scaled to the statistic (|v - sobs| <= 1e-9 · max(1, |sobs|) counts as a
tie, resolved conservatively) — a few characters in `wooldid_permtail`, and the same in `wooldid_permstep`
and `wooldid_permqf`. (b) Expose it as `permeps(#)` with a ritest-compatible default and note the
convention in the help beside the p-convention paragraph. Cost: it changes exact-mode numbers in tied
designs, so it needs a cert pin and a battery re-bless.

**Verification.** Not verified (Stata cap); source-verified, and structural rather than measured — no live
tie was constructed.

#### PS-82 — assorted parser deviations: weights, a silent option override, and operators in `controls()`

**wooldid.** (a) The ado accepts `[aw fw pw]` but refuses fweights with a custom message explaining the
reason (ado:580); the help's syntax diagram advertises only `[aw/pw=weights]`, and iweights are rejected by
Stata's generic r(101). (b) Under poisson, aweights are refused ("poisson only accepts pweights, not
aweights", ado:4707). (c) `robust` is silently overridden when `cluster()` is also given — the vcespecifier
is built by two unguarded sequential assignments, cluster last (ado:1532-1534) — and the help's `robust`
entry says only that it is "relevant only if cluster() is not specified". (d) `controls()` and `fe()` take
`varlist fv` so factor notation works, but time-series operators are refused (r(101)); meanwhile
`interactivecontrols()`/`ccc_absorb()` take `varlist numeric` (no fv) while `coarsecohortcontrols()` takes
a raw string the help says must carry factor notation — three control options, three declarations.

**Practice.** csdid and jwdid take `fv` on their main varlist (jwdid also `ts`). Stata refuses
contradictory VCE requests loudly because `vce()` can be given only once, rather than letting the last one
win silently. Weight types accepted are conventionally listed accurately in the syntax diagram (csdid `iw`;
did_imputation `aw iw`; jwdid `pw`).

**Why it matters.** None is serious alone. (a) the diagram under-reports what the parser accepts, so a
reader cannot tell `fw` is deliberately caught rather than accidentally allowed; (c) `robust cluster(i)`
runs and reports clustered SEs with no note, where a user may believe they asked for something else;
(d) `L.y` in `controls()` is a natural thing to try on a panel and the refusal comes from Stata with no
wooldid context.

**Harmonization.** (a) list `fw` in the diagram as "refused, with an explanation"; (b) no change needed —
a Stata constraint; (c) refuse `robust` with `cluster()` outright, or print a one-line note that
`cluster()` takes precedence (three lines); (d) add `ts` to the `controls()`/`fe()` declarations if the
downstream markout can carry lagged variables, otherwise catch the r(101) and say "generate the lag as a
variable first".

**Verification.** Confirmed by run in the syntax lens for the weight and operator legs (fweight → r(198)
custom, iweight → r(101), `controls(i.gvar)` runs, `controls(L.y)` → r(101)); the `robust`/`cluster()`
precedence is source-verified.

### 2.6 Terminology

#### PS-83 — "cohort" names the `i` variable, not the treatment-timing group the literature calls a cohort

**wooldid.** sthlp:62 defines the second positional argument as "i: Cohort identifier; must consist of
integers greater than 0. Cohorts can be fine (i = individual units) or coarse (i = collections of
individual units), but should be no coarser than to consist of treatment onset cohorts." So "cohort" is
wooldid's word for the panel unit-of-analysis id, and the timing group is a separate variable, `ttre`,
referred to variously as "treatment onset cohort" (5 uses), "onset group" (3) and "treatment onset group"
(6). "cohort" appears 223 times in the help. To obtain effects by treatment-timing group the user writes
`subgroups(ttre)`, whose output header reads "Results for Subgroup ttre == 2005" (ado:9913). wooldid itself
also uses "cohort" in the timing sense in places (sthlp:66 "never-treated cohorts"; sthlp:1256's example 6
is headed "separate treatment effects calculated by treatment onset cohort"), confirming the collision.

**Practice.** In the staggered-DiD literature and in every neighbouring Stata command, "cohort"/"group g"
IS the treatment-timing group: Stata's xthdidregress says "The cohort variable is a categorical variable
indicating the period when the unit is first treated" and its `cohort()` option takes exactly what wooldid
calls `ttre`; csdid's `gvar()` is the "Variable identifying treatment groups or cohorts" with the unit id
carried separately in `ivar()`; Callaway-Sant'Anna index ATT(g,t) by g = first treatment period;
Sun-Abraham call it the cohort; Wooldridge (2021/2025) — the source paper for this estimator — uses cohort
in that timing sense throughout.

**Why it matters.** This is the single largest reading hazard in the package. A reader who knows Wooldridge
(2025) will read `cluster(i)`/`subgroups(i)` as "by timing cohort" and get something else, and will not
think to look for ATT(g) under an option called `subgroups()`. It also makes the help's own sentences hard
to parse — "the cohort variable i or above", "cohorts treated after sample end", "onset cohorts" — because
two different objects share the noun.

**Harmonization.** (a) Cheapest: a one-paragraph "A note on words" box at the top of the Description saying
explicitly that wooldid's "cohort" = the panel unit-of-analysis id `i` (CS's i or a coarsening of it), and
that what CS/Sun-Abraham/xthdidregress call the cohort or group g is wooldid's `ttre`, reachable via
`subgroups(ttre)` — or, when i is coarsened to onset groups, via `att_i` and `histogramcohorteffects`, which
is the second path to ATT(g) and likely the origin of the naming. Add the mapping to the Comparison
section. (b) Middle: rename the prose noun to "unit group" or "i-group" and reserve "cohort" for onset
groups (large help edit, no code change). (c) Most invasive: an alias option `cohorts(ttre)` for
`subgroups(ttre)`.

**Verification.** Confirmed on both halves (the 223 count reproduces exactly; the comparator definitions
verified from the installed help files). Fairness note carried into (a): setting i = `ttre` makes wooldid's
"cohort" and the literature's cohort coincide, which softens the "will not think to look" claim slightly.

#### PS-84 — the not-yet-treated vs never-treated control-group choice is real but is never named as such

**wooldid.** Under the default pooled reference regime, treated cohorts' pre-window periods carry no own
dummies and therefore serve as comparison observations for other cohorts: the help's own comparison section
says the pooled default corresponds to `hdidregress twfe` with `notyet` (sthlp:1198-1199). Under
`esfixedbaseperiod`, "Z is modified to include all (i,t) pairs for treated cohorts, except the pair
corresponding with the period immediately prior to treatment ... the regression is saturated with i-t
indicator variables" (sthlp:154-156) — the never-treated-only design. The help discusses the switch
exclusively as a reference-period "transparency-precision tradeoff" (sthlp:667-690) and never uses the
phrase "control group" for it (that phrase occurs only in the summarystats and continuous-treatment
sections). Measured: the two regimes move the headline number on the shipped fixture — att = .8445031
pooled vs .92238285 under `esfixedbaseperiod` — and the fixed-base figure is identical whether or not an
event-study window is requested, so this is a global comparison-design switch that moves the ATT even when
no event study is asked for, not an event-study normalisation detail.

**Practice.** Every peer makes this a named, first-class option: jwdid's `never` — "Request using Never
treated as control group. Default is using not yet treated ... will exclude g-1 period from the
interactions" (jwdid.sthlp:27), the same construction wooldid reaches through `esfixedbaseperiod`, and
jwdid adds that the results are then "equivalent to Callaway and Sant'Anna (2021) and csdid, and the
proposition by Sun and Abraham (2021)" (jwdid.sthlp:65-66); csdid's `notyet` (default never-treated only);
hdidregress twfe's `controlgroup(notyet)`. Marcus & Sant'Anna (2021) and Callaway & Sant'Anna (2021) treat
the choice as a distinct identifying assumption.

**Why it matters.** "Does it hold if you use only never-treated units as controls?" is among the most
common referee robustness requests in staggered DiD, and a wooldid user searching the help for
"never-treated control" or "not-yet-treated" finds the phrase only in the hdidregress comparison paragraph.
Two users can end up running two different comparison designs while believing they only changed how the
baseline is normalised.

**Harmonization.** (a) Text only: one paragraph in Two Approaches to Estimation saying explicitly which
comparison group each regime uses, that the switch moves the pooled ATT and not only the event-study
geometry, and cross-referencing jwdid's `never` and csdid's `notyet`. (b) Add documented aliases `notyet` /
`nevertreated` that set the corresponding regime, so the option a migrating user types exists. (c) If the
mapping is not exact (out-of-sample-treated cohorts, `oostreatedcontrols` interactions), say where it
departs.

**Verification.** Confirmed on both halves, including the numeric pair, by live runs.

#### PS-85 — the estimand names carry no map onto the literature's aggregation names, or onto ATT(g,t) and ACRT

**wooldid.** The four estimands are defined purely by their internal weighting recipe (sthlp:206-229):
`att_it` is "the simple, unweighted mean of all of the cohort-period-specific estimates", implemented by
reweighting by w/W_it; `att_i` by w/W_ipre or w/W_ipost; `att_itime` by w/(W_it · T_iprepost). The
continuous-treatment derivative is called `ame` throughout (88 whole-word uses in the ado) and prints as a
row labelled `ame` beside `att`. The building-block table option is `gteffects`, headed
"Cohort-by-Relative-Period Effects"; the phrase ATT(g,t) appears exactly once in the whole help and
"group-time" twice, both inside that one option entry (sthlp:960-961), and the option name is never
expanded. "ACRT" appears nowhere.

**Practice.** csdid names its aggregations `agg(simple)`, `agg(group)`, `agg(calendar)`, `agg(event)` and
prints those words; did_imputation speaks of horizons and tau_h; Callaway-Sant'Anna define ATT(g,t) and the
simple/group/dynamic/calendar aggregations by name. For continuous treatment, Callaway, Goodman-Bacon and
Sant'Anna (2024) name the dose derivative ACRT (average causal response on the treated) and the level effect
ATT(d); de Chaisemartin-D'Haultfoeuille speak of dose-response slopes. "AME" is standard Stata-margins
vocabulary but not standard DiD vocabulary.

**Why it matters.** A reader who wants "the csdid simple aggregation" or "the group aggregation" cannot find
it and must reverse-engineer the weighting formulas; a reader looking for an ACRT will not search for
"ame"; a reader who wants ATT(g,t) will not guess the option is called `gteffects`. And `att_it`/`att_i`
read like ATT(g,t) vocabulary while denoting weighting schemes — a csdid user reading `att_i` expects "the
ATT for cohort i".

**Harmonization.** A short table in the Description: one row per estimand, columns [wooldid name |
weighting | nearest literature name | nearest csdid/did_imputation equivalent] — att ↔ observation-weighted
overall ATT ≈ csdid `agg(simple)`; att_i ↔ cohort-equal-weighted ≈ csdid `group` **with the weighting
difference spelled out**; ES rows ≈ `estat event`; ame ↔ ACRT; `e(gteffects)` ↔ ATT(g,t). Spell out
"gteffects (group-time effects)" at the option name. Text only; a supported alias `acrt` is optional and
probably unnecessary. Renaming the estimands is expensive (they are e(b) equation names, matrix row names
and option keys), so the display stub or a table note is the cheaper route.

**Verification.** Partly. All quotes verified. Corrections: the gteffects citation is sthlp:960-961; the
"190 uses of ame" figure was an inflated substring count; and the help is not wholly silent on cross-package
mapping — it carries "Comparing Output to did_imputation" and an hdidregress comparison (sthlp:1194-1225),
and the csdid group-time reference at 960. The real gap is narrower but genuine: no csdid aggregation-name
mapping, ACRT/ATT(d) absent entirely, and the one ATT(g,t) mention buried in an option entry.

#### PS-86 — `bootstrap()` takes {standard, bayes}; `wild` and `pairs` are rejected with a message that never says the wild bootstrap is unavailable

**wooldid.** The token screen is `if !inlist("`bootstrap'", "standard", "bayes")` (ado:6050), so both
`bootstrap(wild)` and `bootstrap(pairs)` draw the identical message: "Error: bootstrap() must be either
bootstrap(standard) — cluster resampling with replacement — or bootstrap(bayes) — Bayesian cluster
reweighting. Received: bootstrap(wild)." Neither mentions that the wild cluster bootstrap is not
implemented, nor names boottest. Elsewhere wooldid's own prose calls the standard mode "the pairs cluster
bootstrap" (ado:6263-6266, 6279), so the token and the prose use different names for the same thing. The
help's `bootstrap()` entry (sthlp:485) does state outright "wooldid does not offer a wild cluster
bootstrap", explains why, and names both boottest on the stage-1 regression and `permute()` — so the reader
who reads the help gets the answer.

**Practice.** "Pairs cluster bootstrap" is the name in Cameron-Miller (2015) and MacKinnon-Webb (2017);
"wild cluster bootstrap-t" is what boottest (`weighttype(rademacher)`, reps 999) and Stata's `didregress
wildbootstrap` (`errorweight(rademacher)`, reps 1000) implement, and what Cameron-Miller (2015) and
MacKinnon-Webb (2017) recommend over the pairs bootstrap generally. A user typing `bootstrap(wild)` on an
ordinary design is following the mainstream recommendation, not guessing.

**Why it matters.** The absence of the wild cluster bootstrap is a ruled and documented decision, so this
is not a re-raise of that. The residual is that the one moment when the user is definitely asking for it —
typing `bootstrap(wild)` — is the one moment the tool does not answer the question, and a user with 40
clusters never sees the caution that would have told them (PS-35).

**Harmonization.** (a) Extend the token-refusal message: name the wild cluster bootstrap explicitly as not
implemented, and point at boottest on the stage-1 regression and at `permute()` — the same pointer the help
and the caution already carry. One string. (b) Accept `bootstrap(pairs)` as a documented synonym for
`standard`, so the token matches the name the help and the literature use. Two lines and a help clause.

**Verification.** Partly. Both refusals reproduced. Corrections: the wild-bootstrap pointer is not confined
to the <6-treated-cluster caution — the help's bootstrap() entry documents the absence and names the
alternatives; only the runtime message omits it. And MacKinnon-Webb (2018) should not be cited as
recommending WCB-t for few TREATED clusters — that paper documents the wild cluster bootstrap's own failure
there (WCR under-rejects, WCU over-rejects).

#### PS-87 — `ame` names the continuous-treatment estimand in Stata margins jargon, not in the continuous-DiD literature's terms

**wooldid.** The dose effects are `ame`/`ame_i`/`ame_it`/`ame_itime`, printed under stripes `ES_ame` and
`main_ame` (ado:3827, 3839), with `contreatelasticitytype()` offering dydx/eydx/dyex/eyex; the cohort
histogram is titled "Cohort-Specific i-t Cell Average Marginal Treatment Effect" (ado:2217). The help calls
them "average marginal effects" throughout (sthlp:726-760) and does spell the acronym out — what is absent
is the DiD-literature name.

**Practice.** The continuous-DiD literature names these objects differently: Callaway, Goodman-Bacon &
Sant'Anna (2024) use ATT(d) for the level effect at dose d and ACRT/ACR — average causal response on the
treated — for the derivative; did_multiplegt_dyn speaks of effects per unit of dose. "Average marginal
effect" does not appear in that literature.

**Why it matters.** A reader who knows the continuous-DiD literature cannot tell from the option list
whether wooldid's `ame` is the ACRT, and cannot tell whether wooldid's `att` under `contreat()` is ATT(d) at
the observed dose or an average of ATT(d) over the treated dose distribution. Referees now ask specifically
about ACRT vs ATT(d) and about the strong-parallel-trends assumption the continuous case requires.

**Harmonization.** Text only: a short glossary at the head of Continuous Treatments mapping `ame` → ACRT/ACR
and att-under-`contreat()` → ATT(d) averaged over the observed treated dose distribution, and stating which
parallel-trends assumption each needs. Keep the option names. A paragraph plus one citation.

**Verification.** Not verified (Stata cap).

#### PS-88 — standard words used with wooldid-internal meanings: family, kernel, package, comparator, twin

**wooldid.** "family" carries at least three distinct senses in the same document: a bootstrap
interval/inference type ("the Wald family", "the percentile and bias-corrected families",
`e(boot_jt_family)`); an estimand block ("for each family and subgroup", "the continuous-treatment
family"); and the multiple-testing sense ("controls the family-wise error rate across the horizons") — all
three inside the step-down and bootstrap paragraphs. "kernel" means a stored per-element table of
ingredients (`e(perm_kernel)`, "the studentized kernel is e(boot_tstud)"), 4 uses, none of them a kernel
function. "package" is used for a set of table columns ("Use bootdisplay(wald) for the normal/t package in
which they do agree") inside a Stata help file where "package" otherwise means an installable ado.
"comparator" means the observed (rep-0) statistic the draws are ranked against (8 uses), not a comparison
group. "twin" means a parallel hidden e() matrix (7 uses).

**Practice.** In this literature "family" is reserved for the multiple-testing sense (FWER, Romano-Wolf,
Westfall-Young); Efron-Tibshirani and Davison-Hinkley call percentile/BC/BCa/percentile-t interval
"methods" or "types", and Stata's bootstrap postestimation calls them normal/percentile/bias-corrected
"confidence intervals", never families. "Kernel" means a smoothing kernel or a matrix null space. In DiD,
"comparator"/"comparison group" means the untreated units.

**Why it matters.** "family" is the worst because wooldid uses two senses within a few lines of the
step-down paragraph, so a reader working out whether the step-down adjusts across estimand families or
across bootstrap families has to guess. "kernel" and "comparator" send a reader looking for a bandwidth and
a control group respectively. "package" is a live ambiguity in a Stata help file.

**Harmonization.** Pure text for most: use "interval method"/"interval type" for percentile/BC/Wald
(keeping "family" for FWER only); "estimand block" for att/att_it/ame; "set" or "column block" instead of
"package"; "observed statistic" instead of "comparator" in prose; "companion matrix" for "twin". The
shipped e() names (`e(boot_jt_family)`, `e(perm_kernel)`, `e(perm_comparator_b)`) can keep their names with
a glossed description.

**Verification.** Not verified (Stata cap).

#### PS-89 — a dozen coined terms carry central meaning with no glossary and no mapping to the literature

**wooldid.** "adversarial" (34 uses) is wooldid's name for the composite test that all event-study
coefficients lie on some line; "shift interval" (8) for the permutation interval; "null band" (3) for the
quantiles of the null draws; "off-centre" (1) for a blanked interval; "fixed metric" (3) for the
permutation/bootstrap-covariance quadratic form; "cell collapse" (4) for the margins-on-collapsed-rows
speedup; "safety" (32) for the global check switch; "lattice" for dose bins; "contreat" for the dose
variable; "ttre" for the onset variable, an abbreviation never expanded anywhere in the help. Option names
extend the pattern: `fastnose`, `fullprecisionslow`, `fastcollapse`, `permabsp`, `permrawp`, `boottstat`,
`bootdf`, `bootalpha`, `altmethod_adversarial`.

**Practice.** The nearest literature terms exist and are widely used: the line test is a linear-pre-trend
or pre-test (Roth 2022; Bilinski-Hatfield; the linear-trend restriction in Rambachan-Roth 2023) —
"adversarial" is a term of art from robustness/ML, not from DiD; permutation intervals are usually
"randomization confidence intervals" or "inverted permutation intervals"; dose bins are "dose bins"/"dose
grid" (Callaway-Goodman-Bacon-Sant'Anna 2024); csdid/xthdidregress call the onset variable gvar/cohort;
and Stata's convention for "skip my checks" is `force` or `nocheck`, not `safety(off)`.

**Why it matters.** Each coinage is individually defensible — several name objects that genuinely have no
standard name (the cell collapse, the fixed metric) — but collectively a reader must learn a private
vocabulary before they can read the option list, and they cannot search for wooldid's feature under the
name they already know. A user looking for a pre-trend robustness test will not search for "adversarial"; a
user looking for `force` will not find `safety(off)`.

**Harmonization.** One "Terminology" section (~30 lines) mapping each coined term to the nearest standard
term and, where one exists, to the package that uses it: adversarial ↔ linear pre-trend / on-a-line
pre-test; shift interval ↔ inverted permutation interval; lattice ↔ dose bins; contreat ↔ dose/treatment
intensity; ttre ↔ gvar/cohort/first-treatment period; safety(off) ↔ force. Optionally accept `force` as an
undocumented synonym for `safety(off)`.

**Verification.** Not verified (Stata cap).

#### PS-90 — the Roth pre-trend power diagnostics print coined column names, and "nis" is never expanded

**wooldid.** The block is headed "Pre-Trend Power Diagnostics (Roth 2022):" with lines "size of the
any-pre-period-individually-significant pre-test: 0.1192" and "weighted mean post-treatment relative period
(hbar) = 2.454545", and a table with columns `power  mdt_wald  mdt_nis  bias_wald  bias_nis`. The footnote
explains what each powers ("mdt_nis powers the test that any pre-period coefficient is individually
significant (the convention of Roth's pretrends package)"), but the letters "nis" are never expanded as an
acronym in the help or the ado — the string appears only as a column suffix or inside those sentences.

**Practice.** Roth's own `pretrends` package reports "slope for 50% power", a Bayes factor and a likelihood
ratio, using words rather than abbreviations; Roth (2022) speaks of the "minimum detectable slope" and of
pre-testing on individually significant pre-period coefficients. "MDT" and "NIS" are not standard
abbreviations in that literature.

**Why it matters.** A reader who copies this block into an appendix has five column names none of which is
self-describing, and cannot look up "NIS" anywhere. Since these are precisely the numbers a referee will
scrutinize (minimum detectable trend, implied bias), opaque labels cost more here than elsewhere.

**Harmonization.** Expand the acronyms once at first use in both the help and the printed footnote —
"mdt = minimum detectable trend; nis = 'any pre-period coefficient individually significant' pre-test" —
and consider renaming the printed columns to `mdt_joint`/`mdt_indiv` (or `mdt_F`/`mdt_pretrends`), which
are guessable without the footnote, keeping the e() matrix column names unchanged. Also spell "hbar" as
"mean post-treatment horizon (hbar)" at first use.

**Verification.** Not verified (Stata cap).

### 2.7 Documentation

#### PS-91 — the modern staggered-DiD literature is absent from the help, the Comparison section names two estimators, and the bibliography omits ten works the text cites

**wooldid.** Grepping the 1,375-line help: Callaway 0, Sant'Anna 0, Sun 0, Abraham 0, Chaisemartin 0,
Goodman-Bacon 0, "negative weight" 0, "forbidden comparison" 0, "clean control" 0, "dynamic" 0. `csdid`
appears once, inside the `gteffects` entry; `did_multiplegt`, `eventstudyinteract` and `jwdid` zero times.
The "Comparison to Related Estimators" section (sthlp:~1189-1224) compares wooldid to exactly two things:
`hdidregress twfe` and `did_imputation`. The Citations list (sthlp:1342-1358) has six entries: BJS,
Rambachan-Roth, Roth, Correia, Ibragimov-Müller, Wooldridge (2025). Meanwhile the body cites by name at
least ten works that never reach it: MacKinnon (sthlp:394, 455, 481, 485), Roodman (485), Davison (481),
Phipson (481), Young (481, 483), Westfall (483), Romano (483), Rosenbaum (481), Conley (413), Efron (504).

**Practice.** boottest's help names MacKinnon, Webb, Roodman, Davidson, Nielsen, Cameron, Gelbach and
Miller repeatedly; csdid cites Callaway and Sant'Anna and asks users to cite CS 2021 and Sant'Anna-Zhao
2020; jwdid credits Wooldridge and Rios-Avila and itself names csdid and Callaway/Sant'Anna. jwdid is the
other Stata implementation of the same Wooldridge ETWFE estimator and already names wooldid in its own
Remarks, so the pairing is visible from the other side. **Counterexample:** did_imputation's help mentions
only Borusyak, has a one-entry reference list, and names no rival package — so a short bibliography is not
itself unusual in this corner of the ecosystem; and csdid, while it cites its own authors, does not name
its rivals either.

**Why it matters.** Two stumbles. A reader cannot tell where wooldid sits relative to the estimators they
already use, and the absence of jwdid from a page about the Wooldridge estimator reads as an omission
rather than a choice — a stated cross-check between two independent implementations of the same estimator
is exactly the kind of evidence referees find reassuring. And the mismatch between an
inline-citation-rich body and a six-entry bibliography looks like an oversight in an otherwise obsessively
complete document: a referee checking "MacKinnon and Webb 2020" finds no reference.

**Harmonization.** (a) Extend Citations to every work named in the body (~10 additions). (b) Add three
short Comparison paragraphs: vs jwdid (same estimator, different implementation and aggregation vocabulary,
ideally with a demonstrated numeric agreement on `mpdta`), vs csdid (under what settings wooldid's att
should match an aggregated ATT(g,t) — never-treated controls, no covariates, `method(reg)` — and where it
will not), vs did_multiplegt/eventstudyinteract in one line each. (c) Add jwdid and csdid to the
`vieweralsosee` block. (d) Optionally one sentence in the Description placing the Wooldridge estimator in
the Goodman-Bacon/negative-weights problem that motivated the whole class. No code change; two verification
runs for (b).

**Verification.** Confirmed on the wooldid side (all greps and counts reproduce). Practice partly: the
"everyone else maps the landscape" premise holds for boottest's citations and only weakly for csdid;
did_imputation is a counterexample.

#### PS-92 — the help is ~4× the longest comparable help, is written in an engineering register, and has no Stored results section or quick start

**wooldid.** 29,206 words / 213,593 bytes. Comparators: boottest 7,534 words, did_imputation 4,157, csdid
2,202, ritest 1,753, honestdid 1,239, jwdid 988 — so ~3.9× the longest and 7-24× the rest. The register:
"measured" 15 uses, "roughly" 16, "refused" 17, "caution" 62, 14 quantities in 1e- notation (0 of each in
all five comparators). Examples include the continuous-treatment collapse paragraph ("no cell moved by more
than 1e-9 in relative terms ... the largest absolute move anywhere was 4.8e-10, in the tenth significant
digit of a t statistic of 39.2", sthlp:307-310); dated changelog lines inside option entries ("impvalues
... left this list on 2026-09-06", sthlp:341; "Before 2026-09-05, plain basic with a lattice inserted no
baseline dose control", sthlp:757-759; "one was offered before 2026-09-06 and has been withdrawn",
sthlp:500) — eleven dated strings in all; one internal repo path
(`review_process/bootstrap_report_2026-08-30.md`, sthlp:486) plus an unresolvable "see the design record"
(sthlp:357); and a paragraph about the help's own truthfulness ("What the help does not claim: ... whether
the resulting intervals attain nominal coverage has not been measured", sthlp:500). Structurally, the
`{title:...}` markers run Syntax / Description / Options / Runtime / Working With Output / Comparison /
Examples / Bugs / Citations / Author / Disclaimer — no `{title:Stored results}`, no Remarks, no Also see, no
quick start, and Examples begin 1,228 lines in. None of the five comparators has a changelog line inside an
option entry, a path to a file the user does not have, or a what-the-help-does-not-claim paragraph; Stata's
convention puts version history in a separate whatsnew and numerical validation in Methods and formulas or
a Stata Journal article, and all five use `{title:Stored results}` with the standard four-way table.

**Why it matters.** The content is unusually honest and much of it is load-bearing, but three specific
things trip readers: a repo path is a dead end for anyone who installed from SSC; a dated "before
2026-09-05 this was wrong" line inside an option entry makes a user wonder what else silently changed and
whether their archived results are affected; and 29k words with no glossary or quick start means the
practical instruction (which option do I use for few clusters?) is buried under measurement appendices. The
missing `{title:Stored results}` heading also costs discoverability — it is where a Stata user looks and
what `help` search anchors on.

**Harmonization.** No content need be lost. (1) Move the ULP/collapse measurement paragraphs, the coverage
simulations and the size-of-the-test paragraphs into a "Numerical accuracy and simulation evidence" section
at the end, leaving one summary sentence plus a jump link at each option. (2) Collect every dated "before
<date> this behaved differently" line into a single Version history section. (3) Replace internal repo
paths with the number itself, a citation, or a link that resolves for an SSC installer. (4) Add a 20-line
Quick start at the top: baseline ATT, event study, inference under few clusters. (5) Rename or duplicate
"Working With Output" as `{title:Stored results}` with the standard four-way table, keeping the prose as
remarks. Cost: a day of help editing, no code change.

**Verification.** Partly. Every quotation, count and structural claim verified on the current tree. One
correction carried above: the finder's "~7× longer than the longest comparable help" is wrong — against
boottest, the longest comparable, wooldid is ~3.9× by both words and bytes; 7× is the ratio to
did_imputation.

#### PS-93 — Example 13 contradicts the `bootstrap()` option text about joint tests

**wooldid.** Example 13 (sthlp:1296) reads "... and joint pre-trend tests formed from the replications' own
test statistics (boottstat; without it the joint tests are absent under bootstrap())". The `bootstrap()`
entry (sthlp:486) says the opposite: "With jointtests the joint pre-trend tests ... are bootstrap
quantities in both forms: under boottstat each is the bootstrap of its test statistic ..., and without
boottstat it is the fixed-metric rank test the permutation engine also uses." The ado confirms the entry,
not the example: ado:6522, 6570 post `e(boot_jt_family) fixedmetric` and ado:10114 prints the matching
footer, i.e. F034 is built.

**Practice.** n/a — this is an internal contradiction, not a deviation.

**Why it matters.** A user reading Example 13 will add `boottstat` solely to obtain joint tests, paying the
re-identification cost and the full-variance replications for nothing.

**Harmonization.** Fix the parenthetical to say the joint tests are the fixed-metric rank test without
`boottstat` and the studentized bootstrap of the statistic with it — a one-line edit, worth a grep for other
pre-F034 sentences making the same claim.

**Verification.** Confirmed (line cite 1296, not 1297).

*Related but distinct (medium, and an enhancement rather than a deviation): the Examples section's 15
entries include no table export (esttab/estout), no publication event-study figure (coefplot/event_plot/
`graph export`), no `permstrata()`, and no triple-difference example, though the subgroup prose describes
the DDD recipe at sthlp:718-721 and Example 10 already shows a four-line honestdid handoff. Adding a
baseline-table, event-study-figure and robustness-panel trio would mirror a paper's workflow — but csdid's
and did_imputation's own Examples blocks are also one-liners that delegate figures and export to companion
commands, so this is an enhancement beyond peer practice, not a gap relative to it.*

#### PS-94 — the `bootsave()` file is not a `bstat` dataset, and the help does not say so

**wooldid.** `bootsave(name)` writes `name.dta` with one row per replication plus a `_rep = 0` row carrying
the observed estimates, the replication coefficients under their `e(b)` names, and — with `cluster()` — one
`w_<value>` column per cluster holding that cluster's multiplicity (or Dirichlet multiplier); without
`cluster()` the resampling unit is the observation and no such columns are written (sthlp:507; built by
`wooldid_bootsavemk`, ado:12002-12070). On a clustered run over 80 units with `reps(30)` the file is 31
observations by 82 variables: `_rep`, `main_att`, and `w_1 ... w_80`. `bstat using` that file returns
"option stat() required"; the word `bstat` appears nowhere in the help. The file does carry
characteristics — `_dta[wooldid_engine]`, `_dta[wooldid_seed]`, `_dta[wooldid_rngstate]` — just not the ones
bstat looks for.

**Practice.** Stata's `bootstrap` prefix writes `saving(filename)` as a bstat-ready dataset: replicate rows
only, with the observed values riding in the file as variable characteristics under `_dta[bs_version] == 2`
(`char <var>[observed]`, plus `[acceleration]` for BCa) and `_dta[N]`/`_dta[N_cluster]`/`_dta[command]` for
the header. The standard workflow is a bare `bstat using file` — which is exactly what fails on wooldid's
file. `stat()` is the fallback for a hand-built replicate file, not the normal route.

**Why it matters.** The name and shape invite the bstat workflow and then do not support it, and if a user
forces it with `stat()` the failure is silent rather than loud: the `_rep = 0` observed row would be counted
as a draw and every `w_*` column read as a statistic. Separately, the 80-to-1 ratio of resample columns to
estimate columns on a clustered run is not what anyone expects from a file described as "the bootstrap
replication draws" — though the `w_*` columns are a genuinely useful addition Stata saves no equivalent of.

**Harmonization.** (a) One help sentence: the file is for the user's own analysis, is not a bstat dataset,
and carries the observed estimates on `_rep = 0` rather than through bstat's `stat()`. (b) Optionally a
suboption to omit the `w_*` columns, or a note that they scale with the cluster count. (c) Full bstat parity
is cheaper than it looks — writing `_dta[bs_version] = 2` plus `char <var>[observed]` per coefficient would
make the file readable — but it still requires dropping the `_rep = 0` row and would have the `w_*`/`se_*`
columns read as statistics, so the disclosure remains the right call.

**Verification.** Confirmed on the wooldid side. Practice partly: the finder's mechanism for how Stata
supplies observed values was wrong (characteristics, not `stat()`), which sharpens rather than weakens the
item — the zero-effort workflow is the one that breaks.

#### PS-95 — repeated cross-sections are supported in substance but never named, and the screens are phrased in unit terms

**wooldid.** The phrase "repeated cross-section" appears nowhere in the help. The design wooldid is
optimised for — coarse i with many rows per cohort-by-period cell — is exactly the repeated-cross-section
case, and the runtime guidance actively recommends it (sthlp:1003-1015). But the `subgroups()` and
`contreat()` screens are written in panel language ("Each treated unit of i should lie wholly within one
subgroup"; "A subgroup variable that changes over time within a unit is likewise refused"), and nothing
tells an RC user which screens apply to them.

**Practice.** csdid documents panel vs repeated-cross-section estimators and offers `rc1`
(csdid.sthlp:97-99); xtevent has `repeatedcs` (xtevent.sthlp:54); jwdid has `trgvar()` for RC data. Users
with pooled CPS/ACS-style data search a help file for the phrase before deciding whether the command
applies.

**Why it matters.** An RC user cannot tell from the help whether wooldid applies to them, and the most
natural reading of the unit-phrased screens is that it does not. This is a discoverability loss on a design
the command is actually good at.

**Harmonization.** Text only: a short paragraph saying a repeated cross-section is handled by setting i at
the cohort level (rows are then observations within a cohort-period cell), stating what "unit" means in the
screens in that case, and noting that the att vs att_it weighting choice is where it shows up.

**Verification.** Not verified (Stata cap).

#### PS-96 — the results table's row-equation structure is displayed but never documented

**wooldid.** The main results matrix prints with roweq blocks `main | att`, `pretreatment | att`,
`ES_att | post0, post1, ...` (roweq handling at ado:6714, 6860, 7705; real archived output in
`test/advg_t2.txt:14-34`). Grepping the whole help for "pretreatment" returns one hit, inside the
`boottstat` paragraph; "ES_att" one hit, inside the `event_plot` paragraph. Nowhere does the help say what
the `pretreatment` row averages over, whether it respects `esprelength`, how it relates to the joint test,
or how to read the `relyr` column (blank on the main and pretreatment rows).

**Practice.** csdid's `estat event` labels the analogous quantities Pre_avg and Post_avg and documents
them; hdidregress prints an ATET table whose rows are described in the postestimation help; every Stata
estimation command documents the layout of its own output table.

**Why it matters.** A user's first question about a wooldid table is "what is the pretreatment att row?" —
a pooled average of the pre-treatment effects, presumably, but over which horizons, under which weights,
and is it the object the joint test tests? — and the help has no answer. Anyone writing up a wooldid result
has to reverse-engineer the row from the ado. See PS-75 for the labelling half.

**Harmonization.** (a) A short "Reading the results table" subsection under Working With Output describing
the three row equations, the pretreatment row's window and weighting, the `relyr` column's blanks, and the
pre#/post# naming convention — text only, and it removes a class of user questions. (b) Optionally rename
the row to `pre_avg` for parity with csdid, at the cost of breaking anyone parsing row names.

**Verification.** Not verified (Stata cap).

#### PS-97 — the help's own canonical idiom, `cluster(i)`, is the configuration several features refuse or degrade

**wooldid.** The help recommends setting i "at the level at which you cluster" (sthlp:64) and Example 0 is
`wooldid y i t ttre, cluster(i) att`. But with `cluster()` at or coarser than i, default safety hard-errors
`jointtests` with r(198) (a 570-character message); the ame-family standard errors under `contreat()` are
set to missing with a caution; `permtstat` draws a caution saying it is unprotected there (ado:5901); and
the analytic pairwise-contrast SEs can be degenerate under `boottstat`. There is no cross-reference to any
of this from the `i` synopt or the Runtime section.

**Practice.** csdid's `estat pretrend` and did_imputation's `pretrends(#)` run at the package's default
(unit-level) clustering with no gate; the few-cluster problems those gates address are handled in that
literature by switching inference method (boottest, ritest), not by refusing the diagnostic.

**Why it matters.** Composing the two idioms the help itself teaches — Example 0's `cluster(i)` and Example
4's `jointtests` — produces a hard error on the first try. The gate is right on the statistics; the
collision is with the documentation's own house style, which does not warn the reader that the recommended
clustering is the constrained one.

**Harmonization.** Documentation is enough: state in the `i` entry and near Example 0 that clustering at i
buys efficiency but forecloses the analytic joint tests and the analytic ame SEs, and point at
`bootstrap()`/`permute()` there rather than only inside the error. Optionally change Example 4 to carry a
`cluster()` finer than i so the two examples compose. See PS-64.

**Verification.** Confirmed for the jointtests leg (PS-64); the rest source-verified.

#### PS-98 — refusals are worded as judgements about the user's reasoning, and cautions run to full paragraphs

**wooldid.** Running `contreat()` without `contreatcontrols()` exits with "Error: contreatcontrols type not
specified; you must specify this - substantive judgement is required." Cautions are multi-sentence: the
`jointtests` refusal under `safety(off)` prints eight lines ending "Consider wooldid's built-in
randomization inference (the permute() option) or bootstrap-based clustered SEs for joint-type inference
here, or cluster at a finer level."; the AME cell caution runs five lines. The ado carries ~100 distinct
"Caution:" strings, and the help uses "refused" 17 times and "caution" 62 times. **In fairness:** ordinary
runs print none of them — a plain `att cluster(i)` run, a four-estimand run, an event-study run, an
`impvalues` run and a `fastnose` run each printed exactly zero cautions.

**Practice.** Stata error messages are terse and mechanical ("option x() required", "factor variables not
allowed"); packages that want to steer users do so in the help, not the error text. csdid, did_imputation
and boottest print no paragraph-length cautions and their helps contain essentially none of this vocabulary
("caution": 0 / 2 / 0).

**Why it matters.** The register is unusual in both directions. A refusal that tells the user their
judgement is required reads as the package declining to have a default rather than as a syntax error, which
some users will find helpful and others presumptuous; and multi-paragraph cautions in the middle of
estimation output are hard to distinguish from failures.

**Harmonization.** Mostly leave the trigger logic alone and edit the wording. For refusals, lead with the
mechanical fact and put the reasoning second — "option contreatcontrols() required with contreat(); see
help wooldid##continuous for how to choose a type" — so the first line matches what Stata users scan for.
For long cautions, print a one-line headline plus a help jump, keeping the full explanation in the help
entry it points to (this composes with PS-9 and PS-32). Say in the help, once, that wooldid deliberately
prints more cautions than neighbouring packages, and what a Caution does and does not mean.

**Verification.** Not verified (Stata cap) as a whole; the individual messages quoted were observed in the
lens runs.

#### PS-99 — the continuous-treatment approach is sourced to a tweet

**wooldid.** sthlp:161-164: "Note that this approach is a bit speculative and is not covered in Wooldridge
(2025), though Wooldridge has {browse "https://twitter.com/jmwooldridge/status/1695115782739435922...":
informally suggested} handling continuous treatments using interactions between treated cell indicators and
the continuous treatment variable." The continuous-treatment surface is labelled [experimental] in the
opening paragraph.

**Practice.** Package helps cite working papers, journal articles or the package's own methods note: csdid
cites Callaway-Sant'Anna, did_imputation cites Borusyak et al., honestdid cites Rambachan-Roth. The
continuous-treatment DiD literature now has citable sources — Callaway, Goodman-Bacon and Sant'Anna (2024)
on ATT(d)/ACRT, and de Chaisemartin-D'Haultfoeuille on dose-response — neither of which the help mentions
(both are 0 hits).

**Why it matters.** The honesty is admirable and the [experimental] flag is the right call, but a tweet URL
is not a durable citation (links rot, the platform renamed) and a reader deciding whether to trust the
`ame` estimates has nothing else to check. It also leaves the impression that no literature exists on
continuous-treatment DiD, which is no longer true.

**Harmonization.** Keep the candour and the flag; add a citation to Callaway, Goodman-Bacon and Sant'Anna
(2024) (and optionally de Chaisemartin-D'Haultfoeuille) with one sentence on how wooldid's within-cell dose
slope relates to their ACRT and one on where the construction differs. Keep the tweet as a secondary note
or drop it. Text only; composes with PS-87.

**Verification.** Not verified (Stata cap); the greps were run.

## 3. Where wooldid is more careful than common practice

These came out of the sweep as deviations, but in the conservative direction. Listed so that any move
toward "looking more like the neighbours" does not trade them away.

1. **Automatic re-identification of duplicated clusters (PS-68).** Stata makes this the user's job
   (`cluster()` + `idcluster()`), and omitting it is one of the commonest cluster-bootstrap errors in
   applied work. wooldid chooses the resampling mechanics from the estimand and refuses when `cluster()`
   does not nest i.
2. **Per-estimate availability accounting and a caution on any failed replication (PS-68, PS-25).** Stata
   records only `e(N_misreps)`; wooldid posts `e(boot_failed)` and per-coefficient `e(boot_navail)` and
   names the exclusion convention on screen.
3. **The replication draws are posted in `e(bootdist)`/`e(permdist_b)`** (PS-68), where Stata's prefix
   offers them only through `saving()`.
4. **Randomization-design diagnostics no RI tool offers (PS-71):** the possible-assignment count, a refusal
   when the design admits one assignment, a completion floor, over-sampling and degeneracy cautions, and an
   off-centre detection that blanks and names a shift interval rather than printing a misleading one.
5. **A replication floor and a required `permreps()` (PS-43).** ritest and `permute` default to 100, which
   under wooldid's equal-tailed convention would silently impose a .02 p-floor on a user who never chose
   the number.
6. **`bootdf()`'s t(G-1) default when clustered (PS-42)** — the convention Cameron-Miller recommend, and the
   only package that adapts the reference distribution to whether the user clustered.
7. **Blanking rather than faking inference where none exists:** the gteffects cells under `permute()`
   (PS-72), the ame SEs under cluster-at-i, the off-centre shift intervals. Filling those columns from the
   analytic run would be the easy and wrong move.
8. **The disclosure footers themselves.** The bootstrap and permute footers state what `e(V)` holds and
   what `test`/`lincom` will reproduce — a disclosure csdid's own `wboot` path (printed uniform bands,
   asymptotic SEs in e(), no footer) does not make. The problem identified in PS-2 and PS-62 is that the
   analytic path has no such line and the footers do not name the export tools, not that the footers are
   wrong.
9. **A studentized bootstrap of the aggregated DiD functional at all (PS-69),** with a candid "whether the
   resulting intervals attain nominal coverage has not been measured" — more honest than anything
   comparable in the ecosystem.
10. **The dropped-coefficient identification screen** (the check that a treatment-cell variable was not
    dropped for collinearity). No comparator screens for this; it is the single most valuable thing
    `safety(off)` turns off, which is the argument for PS-66's scoped switch rather than against the screen.
11. **`se_analytic` and `t_pt` as column names (PS-26)** are more informative than the unlabelled or absent
    columns Stata prints for the same quantities, even though no other package uses those words.
12. **`w_<cluster>` multiplicity columns in `bootsave()` (PS-94)** — Stata saves no equivalent, and they are
    what makes a saved run auditable.

## 4. Items the sweep found that were already ruled on 2026-09-06

Listed only; none is re-opened here, and where the ruling itself is unusual relative to practice the
relevant catalogue item says so in one line.

| ruled item | how the sweep met it |
|---|---|
| F038 — bootstrap p prints 0, "< 2/n" annotation | PS-3 (excluded from the format item), PS-27 |
| F053 — `t_wald`/`z_wald` relabel under percentile/bc | PS-26 (recorded: Stata blanks the column instead) |
| F035 — shift interval calibrated near the null; Conley-Taber attribution dropped | PS-13, PS-52, PS-55 |
| F050 — caution below 6 treated resampling clusters | PS-35, PS-50, PS-28 |
| F057 — observation-level bootstrap on panels: do nothing | PS-35 (the unclustered-run gap is noted there) |
| F042 — blank `se_analytic`/`t_perm` under fast mode: no action | PS-17 (recorded as "the ruling is unusual") |
| F030 — low-reps caution: do nothing | PS-36, PS-27 |
| F034 — fixed-metric rank joint tests without `boottstat` | PS-15, PS-93 (the stale Example 13 parenthetical) |
| F036 — permreps floor caution; floor stays 20 | PS-16, PS-43 |
| F064 — one-sided p, SE(p), Clopper-Pearson CI | PS-24 |
| F065 — max-\|studentized\| step-down under `permute()` | PS-14, PS-45 |
| F067 — BCa interval | PS-45, PS-49, PS-60 |
| F052 — p-convention paragraph beside ritest | PS-27, PS-80 |
| F049/F051 — exchangeability caveat and RI-t guidance | PS-28, PS-54, PS-71 |
| F054 — `level()` on the analytic path | PS-31 (the hard-named `lb_95ci` header) |
| F062 — impvalues under the cell collapse | PS-20 (a different impvalues issue: the banner) |
| F066 — runtime projection | PS-70 (proposed as the text of the `permexact` refusal) |
| F001 — bayes `/W_c` split | PS-12 (cited only as the defect a bias column would have caught) |
| F059 — rngstate record | PS-71, PS-94 |
| F063 — bootsave resample columns | PS-94 |
| F044 — `level()` into `wooldid_ptp` | PS-90 (same block, different issue) |
| F043 — `permute(dose)`: deferred | PS-52 (listed as the comparable cost of a new assignment generator) |
| F058 — calibration simulations | PS-69 (what would replace the candid disclaimer with a number) |
| F031 — worker headroom, treat as frozen | relevant to the cost of PS-8, PS-15, PS-19: new display text belongs in the display subprogram, not the worker |

## 5. Refuted findings

**No finding was refuted outright.** Every item that reached a verifier survived on its wooldid half; the
verifiers' work went into correcting one half or the other. For completeness, the specific sub-claims that
were **refuted or materially corrected**, and are therefore not asserted anywhere above:

1. "under `permute()` no `e(V)` is posted at all" — true only of the fast-mode default; with
   `permanalyticse` or `permtstat` an ANALYTIC `e(V)` is posted (PS-62).
2. "ritest posts no `e(V)`, so wooldid's `permanalyticse` hazard is worse than practice" — ritest restores
   the observed run's analytic `e(b)/e(V)`, so that hazard *is* the practice (PS-46).
3. "csdid, did_imputation, reghdfe and xtdidregress all keep printed and exported families the same" —
   csdid's `wboot` path prints uniform bands while storing asymptotic SEs, with no footer (PS-2).
4. "csdid prints '(Std. err. adjusted for G clusters)' on every run" — only on the wboot path (PS-1, PS-5).
5. "passing `format(%9.0g)` to matlist fixes the number rendering" — that is already matlist's default; the
   proposed fix is a no-op, and no single global format fixes both the decimals and the header truncation
   (PS-3).
6. "wooldid rejects a 0-based `egen group()` id" — `egen group()` is 1-based, and wooldid enforces no sign
   rule at all; the real hazard is that `_i = 0` is an internal sentinel (PS-76).
7. "did_imputation's help discusses BJS vs CS vs Sun-Abraham and points at did_multiplegt" — its reference
   list has one entry and it names no rival (PS-91).
8. "no other package normalises to a pooled multi-period reference" — did_imputation's pretrend reference
   group and jwdid's default ETWFE baseline are the same construct (PS-34).
9. "the `pre#`/`post#` vocabulary is the deviation" — did_imputation uses `pre#` too; the deviation is the
   regime-dependent meaning and the table/plot asymmetry (PS-7).
10. "Roth-Sant'Anna (2023) randomize conditional on not-yet-treated status" — their randomization test
    permutes adoption timing, i.e. it is the design wooldid implements (PS-52).
11. "MacKinnon-Webb (2018) recommend the wild cluster bootstrap-t for few treated clusters" — that paper
    documents WCB's own failure there (PS-86, PS-11, PS-50).
12. "Cameron-Miller (2015) name a 20-50 cluster danger zone" / "recommend stratified pairs resampling" —
    neither is in the paper; their rule of thumb is around 50, and their remedy is the wild cluster
    bootstrap (PS-35, PS-50).
13. "did_imputation and `estat ptrends` print a pre-trend p in the same design wooldid refuses" — same
    *clustering choice*, but from a non-saturated specification, so the oversizing mechanism does not arise
    there and it is not a general property of cluster-at-i Wald pre-tests (PS-64).
14. "the help is ~7× the longest comparable help" — ~3.9× (boottest); 7× is the ratio to did_imputation
    (PS-92).
15. "the adversarial help never uses the words reject or evidence" — sthlp:604 does; the accurate claim is
    that it never states the direction, and sthlp:592-593 points the wrong way (PS-19).
16. "`altmethod_adversarial` has the conventional sign, so printing it instead fixes the reading" — it
    shares the same reversed sign (PS-19).
17. "the joint pre-trend line's sharp-null nature is undisclosed under `permute()`" — the footer at
    ado:10039 discloses it; the gap is line-level labelling (PS-15).
18. "the step-down raw p differs from the table's p by about a factor of two" — they differ by about
    1/(n+1) additively; the factor of two occurs only at the floor (PS-14).
19. "the permute footer already lacks an explanation of the shift interval" — it carries a full paragraph;
    the live gap is the column name (PS-13).
20. "Stata's bootstrap `saving()` file supplies observed values through `bstat, stat()`" — it carries them
    as `char <var>[observed]` under `_dta[bs_version] = 2`, so the failing workflow is the bare
    `bstat using` (PS-94).
21. "`jackknife` appears nowhere in the package" — 1 hit in the help and 17 in the ado, all unrelated to the
    DR/IPW point the grep was cited for (PS-47).
22. "wooldid's bootstrap bias column is unusual because Stata prints one" — Stata's `bootstrap` prefix main
    table omits bias too; the deviation is the missing redisplay (PS-12).
23. "the `pretreatment` block and the seven-column layout include `relyr` on the pairwise table" — the
    pairwise matrix has six columns and no `relyr` (PS-4, PS-23).
24. "csdid's examples run estimate → `estat event` → `csdid_plot`" — its Examples block is a list of
    one-line calls; `csdid_plot` appears once, as a listed companion (PS-93).
25. "jwdid's `over()` drives the whole aggregation menu" — `over()` is documented for `simple` only;
    `other()` subsets the rest (PS-48).

---

*Sweep scratch and probe logs: `/Users/hegland/.claude/jobs/be9644d0/tmp/sweep/` (subdirectories per lens
and per verifier). Fixtures used: `test/data/dta_binary.dta`, `test/data/dta_cont.dta`,
`test/data/dta_tiny.dta`. Comparator sources read at `/Applications/Stata/ado/base/` and
`~/Library/Application Support/Stata/ado/plus/`.*
