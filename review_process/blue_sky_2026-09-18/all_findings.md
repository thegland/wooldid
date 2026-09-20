# All findings — the complete corrected ledger

**223 findings, regenerated 2026-09-19 from the corrected dataset.** The published version of this
file was built from a join that collapsed distinct findings sharing a generic id (`SK-01`, `SKEP-01`,
`SK-1`), so some verdicts were attached to the wrong rows. This one is keyed on unique ids.

**Every finding here carries at least two independent reviewer opinions** (118 with 2, 94 with 3,
11 with 5). Severity is the **median** of all opinions — the author's and every reviewer's — not the
most severe, so a single alarmed voice cannot inflate a row.

**Read the verdict provenance.** `CONFIRMED`/`PARTIAL`/`REFUTED` are reviewer verdicts, but the
reviewers are not interchangeable. 118 findings were checked by a prober and an independent refuter
who wrote its own repro; 16 of those also got a third-lens reviewer; **89 were checked only in the
later adjudication rounds**, because their refuter cells were lost to the weekly usage limit. Each
entry names which. "3 reviewers" does not always mean prober+refuter+lens.

| severity | count |
|---|---|
| critical | 4 |
| high | 48 |
| medium | 87 |
| low | 72 |
| none | 12 |

---

## `SG` — subgroups() (19)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `SK-01` | critical | CONFIRMED x4 | pr+4late | The ruled subgroup-degeneracy guard misses on ordinary ragged panels, and wooldid then prints a subgroup-contrast p of 0.000 where the truth is 0.25 |
| `FASTNOSE-03` | high | CONFIRMED | pr+1ref | fastnose publishes subgroup ame event-study rows the full run drops as non-estimable — including a POST-treatment horizon reported as exactly 0 — cont… |
| `SKEP-01` | high | CONFIRMED x4 | pr+4late | permute()'s fast-mode DEFAULT publishes a whole subgroup ame event study the analytic run skips entirely, as zeros carrying permutation p-values and s… |
| `SKEP-02` | high | CONFIRMED x3 PARTIAL | pr+4late | bootstrap() and permute() permanalyticse abort r(459) on contreat() + subgroups() + an event-study window, with a 'this is unexpected - please report… |
| `fpsg-01` | high | CONFIRMED x2 | pr+1ref+1lens | A partially-missing subgroup variable silently turns subgroup att_i / att_itime into a share-weighted average |
| `hv-bootrank-03` | high | PARTIAL | pr+1ref | The listwise complete-case rule makes e(V) a factor of two NARROWER than the standard errors wooldid printed, on a realistic state-policy event study,… |
| `skep-sgexport-01` | high | CONFIRMED PARTIAL | pr+2late | subgroups() silently rescales every exported standard error in e(V) - identical estimates, identical printed table, up to -36% on the export |
| `skp-01` | high | CONFIRMED x2 | pr+2late | att_it under subgroups() is cell-SHARE-weighted, not cell-equal-weighted, when an (i,t) cell is only partially labelled |
| `BT-05` | medium | PARTIAL | pr+1ref | With cluster() at or above i, every subgroups() row's analytic standard error is numerically identical to the pooled one (and so is e(boot_ase_obs), t… |
| `EC-3` | medium | PARTIAL CONFIRMED | pr+2late | A run that reports r(198) can leave a COMPLETE, CORRECT wooldid e(), and the replay refusal then tells the user something false about it |
| `SK-03` | medium | CONFIRMED x2 | pr+2late | The printed bootstrap standard error and e(V) are built on different completeness masks, so under bootdisplay(wald) lincom does not reproduce the prin… |
| `hv-sgname-01` | medium | CONFIRMED | pr+1ref | Subgroup e(b)/e(V) column names are an infix `_SG<value>_`, not the help's suffix `_s#`; the true convention appears nowhere in wooldid.sthlp |
| `skep-threese-02` | medium | PARTIAL CONFIRMED | pr+2late | One coefficient, three standard errors: printed, e(V), and e(V_es_att) - and the footer points honestdid at the wrong matrix |
| `skp-02` | medium | CONFIRMED x2 | pr+2late | Missing subgroup values on PRE-treatment rows silently share-weight the pretreatment and pre-trend rows while main:att_i stays exact, so one unit is w… |
| `EC-5` | low | PARTIAL CONFIRMED | pr+2late | Failure-site enumeration: only 3 of 20 sites destroy the prior e(), one of them is new, and the replay-staleness message never fires for exactly those… |
| `EC-7` | low | CONFIRMED x2 | pr+2late | The help names e(b) subgroup coefficients with a suffix the code never produces |
| `SAFETY-04` | low | CONFIRMED | pr+1ref | The negative-weight screen refuses on rows the run discards anyway (raw-if domain, not the estimation domain) |
| `hv-permsave-comma-02` | low | PARTIAL | pr+1ref | Under boottstat the studentized joint/equality statistics are not ereturned, though the help lists them; boot_sgW is computed and discarded |
| `SKEP-03` | none | REFUTED x2 | pr+2late | Where a subgroup's continuous event study is 'skipped entirely', e(b_es_ame_s#) is still posted as a 1x1 matrix of MISSING under the junk colname c1,… |

## `UZ` — contreatcontrols(untreatedZero) (7)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `fp-contreat-E` | critical | CONFIRMED x2 | pr+1ref+1lens | contreatcontrols(untreatedZero) silently gains a baseline dose control on the customvce()/unconditionalse/poisson paths, and the reported ame is arbit… |
| `sk-untreatedzero-regresspath` | critical | CONFIRMED x2 | pr+2late | contreatcontrols(untreatedZero) is silently ignored on every non-reghdfe estimation path: customvce(), unconditionalse and poisson estimate the contre… |
| `FP-CTZERO-ENGINE` | high | PARTIAL x2 | pr+1ref+1lens | contreatcontrols(untreatedZero): the default and explicit-dummy engines return main_att and main_ame of opposite sign, and only the explicit arm warns |
| `SK-UZ-POISSON` | high | CONFIRMED PARTIAL | pr+2late | poisson + contreatcontrols(untreatedZero) breaks on the DEFAULT engine too: a true xb-scale ame of 0.400 is reported as 0.100 |
| `HV-EPOST-01` | medium | CONFIRMED | pr+1ref | A wooldid failure after the internal regression leaves reghdfe's (or ppmlhdfe's) saturated model in e(), destroying the user's previous results |
| `skep-poisctrl-fallback` | medium | CONFIRMED x2 | pr+2late | poisson + contreat() always trips the ppmlhdfe rc-3200 fallback, on a one-level fixed effect wooldid itself constructed, and the caution blames the us… |
| `hv-F5-verify` | low | PARTIAL | pr+1ref | VERIFIED (true, but consequence-different): contreat() ships b=0 with e(V) diag -1e100 while e(V_es_ame) carries exactly 0 on the same horizons |

## `IMP` — impvalues (5)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `mcci-01-impvalues-se` | critical | CONFIRMED x2 | pr+1ref+1lens | impvalues (Ibragimov-Muller) confidence intervals cover 5-29% at nominal 95% wherever several clusters share the same levels of i |
| `skep-01-impvalues-imunit` | medium | CONFIRMED PARTIAL | pr+2late | impvalues under-covers (0.82-0.86 at nominal 0.95) even at i = unit with cluster(unit), the geometry where Ibragimov-Muller is textbook-legitimate |
| `HV-FN1` | low | CONFIRMED | pr+1ref | REFUTED: the 2026-09-13 item-11 ruling WAS built — as option (b) for fastnose plus the (a) half for impvalues |
| `skep-02-im-degeneracy-message` | low | CONFIRMED PARTIAL | pr+2late | The IM degeneracy message blames 'deterministic or noiseless data' on a panel whose idiosyncratic sd is 0.55; the real cause is the cluster geometry |
| `SKEP-N4` | none | CONFIRMED x2 | pr+2late | NEGATIVE RESULT: impvalues is clean on this axis |

## `SMP` — replication engines' sample is not e(sample) (22)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `FCM-1` | high | CONFIRMED | pr+1ref | The bootstrap's two-regime unit rule is a single GLOBAL flag on a PER-UNIT property: one straddling unit out of 96 cuts the att_i bootstrap SE by 20-2… |
| `FCM-2` | high | CONFIRMED x2 | pr+1ref+1lens | The regime predicate is evaluated on `smp`, which is NOT the estimation sample, so cluster values on rows the regression never sees change the reporte… |
| `FCM-4` | high | CONFIRMED x2 | pr+1ref+1lens | A cluster with no rows in e(sample) is still a resampling unit and still sets bootdf(clusters), leaving e(df_r) inconsistent with e(N_clust) |
| `FCW-01` | high | PARTIAL | pr+1ref | One straddling unit silently converts the cluster bootstrap of att_i into a presence bootstrap, cutting its SE 23% |
| `FCW-02` | high | CONFIRMED x2 | pr+1ref+1lens | Phantom straddler: a row the estimator discards decides which variance the bootstrap computes |
| `HV-BOOT-01` | high | CONFIRMED | pr+1ref | The bootstrap engine's cluster universe is the raw-`if` sample, so every few-cluster safety guard is silently suppressed and bootdf(clusters) is wrong |
| `HV-BOOT-04` | high | CONFIRMED | pr+1ref | The phantom bootstrap cluster needs no missing data at all: a singleton unit creates one |
| `HV-BOOT-05` | high | CONFIRMED | pr+1ref | Worst case measured: the bootstrap reports t(11) where t(1) is correct, with every few-cluster warning silenced |
| `HV-BOOT-06` | high | CONFIRMED | pr+1ref | The one-day-old T2/T3 nesting-regime decision is made on the resample frame, so one dropped row flips the regime and collapses att_i's bootstrap varia… |
| `SKEP-N1` | high | CONFIRMED x2 | pr+2late | The bootstrap's cluster count comes from the DESIGN sample, so a cluster whose outcome is never observed still buys a degree of freedom — and the same… |
| `SKEP-N1` | high | CONFIRMED x2 | pr+2late | The phantom bootstrap cluster needs no missing data and no singleton: a cohort treated before the estimation window opens mints one, and the obvious u… |
| `BOOTX-3` | medium | PARTIAL | pr+1ref | the posted bootstrap e(V) is singular whenever completed reps <= the number of reported coefficients, nothing records it, and the help promotes a hand… |
| `BT-02` | medium | CONFIRMED | pr+1ref | bootdf(clusters) counts clusters that contribute no rows to e(sample); one printout says 'N clusters = 17' and 't(19)' |
| `SK-2` | medium | CONFIRMED x2 | pr+2late | bootstrap()'s few-cluster and few-TREATED-cluster cautions are suppressed by a cluster that contributes no estimable row |
| `SKEP-N2` | medium | CONFIRMED x2 | pr+2late | A zero survey weight mints a phantom bootstrap cluster the same way |
| `SKP-1` | medium | CONFIRMED x2 | pr+2late | A cluster with no rows in the estimation sample counts toward every cluster-count screen, so both few-cluster cautions go silent on the run that most… |
| `BR-3` | low | CONFIRMED | pr+1ref | Without cluster(), the count leg and the re-identification leg draw DIFFERENT resamples at the same rseed(), so att's reported SE moves when a second… |
| `MCBP-GUARD` | low | PARTIAL x2 | pr+2late | The few-cluster caution keys on the TREATED cluster count, so a 10-cluster panel with 6 treated and 4 never-treated is warned about nothing while its… |
| `PERMEX-PX-1` | low | CONFIRMED | pr+1ref | permstepdown's p_raw/p_singlestep/p_stepdown are documented and printed as add-one shares, but permexact computes plain exact shares |
| `SK-3` | low | PARTIAL CONFIRMED | pr+2late | The unclustered bootstrap is not reproducible across a re-sort of the same data at a fixed rseed(), while the clustered bootstrap is - and nothing say… |
| `SK-FC-B` | low | CONFIRMED x2 | pr+2late | The package has a named few-treated-cluster caution with a numeric threshold, but only on the engine whose own help says it is not a remedy for few tr… |
| `SKEP-N2` | low | PARTIAL x4 | pr+4late | The unclustered leg T2 opened returns standard errors ~36% too wide on a panel, and the revisit the ticket asked for was never done |

## `GUARD` — safety guards and preflights (16)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `SAFETY-01` | high | CONFIRMED | pr+1ref | The si8 maxvar preflight falsely refuses an ordinary wide panel, contradicting its documented zero-false-positive design |
| `SKEP-01` | high | CONFIRMED x3 PARTIAL | pr+4late | At i = unit the analytic ATT standard error degrades with the SIZE OF THE CONTROL POOL, reaching SE/SD 0.21 and 37% coverage at nominal 95% |
| `bb-01-bayes-cellwtd-se-invariant` | high | PARTIAL | pr+1ref | At cluster()==i, bootstrap(bayes)+safety(off) reports standard errors 4.4x apart for a bit-identical point estimate, and the cell-weighted SE is exact… |
| `hv-ctw-ipos` | high | CONFIRMED | pr+1ref | contreatwithin + contreat() silently skips the 'i must consist of positive, non-zero integers' refusal, turning a named r(411) into a raw r(452) on an… |
| `sk-izero-preweight` | high | CONFIRMED x2 | pr+2late | A treated cohort whose identifier is literally 0 is silently downweighted in the pretreatment att_i / att_itime (and ame_i / ame_itime) rows |
| `FW-01` | medium | CONFIRMED | pr+1ref | A weight written as a +/- expression is pasted verbatim into the estimand-weight arithmetic: silently wrong att_it/att_i/att_itime on the customvce()/… |
| `SAFETY-02` | medium | CONFIRMED | pr+1ref | contreatwithin silently disables two data-integrity screens the help says gate on the user's own safety() value |
| `SAFETY-03` | medium | PARTIAL | pr+1ref | The onset-commensurability bypass silently relabels the onset year as pre-treatment and manufactures a 14-sigma pre-trend |
| `hv-ctw-float` | medium | CONFIRMED | pr+1ref | contreatwithin + contreat() silently suppresses the float-above-2^24 silent-merge caution — a regression introduced by CS-98 (e26b3dd, 2026-09-09) |
| `sk-maxvar-fp` | medium | CONFIRMED x2 | pr+2late | The maxvar preflight is not a lower bound and DOES fire falsely: si8demand adds the ABSORBED i and t levels, so an ordinary many-unit panel is refused… |
| `state-leak-01-missing-control-outcome` | medium | PARTIAL | pr+1ref | An outcome missing on every never-treated row is refused as a DESIGN failure, and the remedy the message names substitutes a different sample with no… |
| `SAFETY-05` | low | PARTIAL | pr+1ref | The esfixedbaseperiod pre/post refusal degrades to a bare Stata r(2001) under safety(off) |
| `SKEP-02` | low | CONFIRMED x2 | pr+2late | The printed relyr column rounds distinct fractional horizons onto the same integer, showing -2 / -0 / 0 for -1.5 / -0.5 / +0.5 |
| `bbr-N3-bayes-cellwtd-gate-is-regime-blind` | low | PARTIAL x2 | pr+2late | The bootstrap(bayes) x cell-weighted refusal is regime-blind: it refuses where bayes is accurate and lets safety(off) through where bayes is worst, al… |
| `hv-safety-case-1` | low | CONFIRMED | pr+1ref | safety() is validated ~500 lines after it is first consulted, case-sensitively — but nothing is silently accepted; the cost is one misleading refusal |
| `sk-ctw-bare198` | low | PARTIAL x2 | pr+2late | contreatwithin + contreat() + esfixedbaseperiod on a panel with no usable pre-period burns a full reghdfe and then dies with a bare r(198) whose entir… |

## `ELAS` — elasticity/semielasticity denominators (9)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `FP-SEMIELAS-FE` | high | CONFIRMED x2 | pr+1ref+1lens | semielasticity divides the effect by a prediction that excludes the absorbed unit and time fixed effects, so the default engine and the explicit-dummy… |
| `FUZZ-01` | high | CONFIRMED | pr+1ref | contreatelasticitytype() relabels a binary-mode ATT table as elasticities while changing nothing |
| `PX-03` | high | PARTIAL | pr+1ref | semielasticity on a non-positive outcome fails in margins and is never diagnosed; under fastnose — what permute() selects by default — the run returns… |
| `fp-contreat-A` | high | CONFIRMED x2 | pr+1ref+1lens | contreatelasticitytype(eydx)/(eyex) divide by a prediction that omits the absorbed fixed effects on the default reghdfe path, so the reported elastici… |
| `FUZZ-02` | medium | CONFIRMED | pr+1ref | Semi-elasticity forms die with a raw Stata/Mata error (tempvar leaked) when the fitted values are not strictly positive |
| `hv-noparent-F2` | medium | PARTIAL | pr+1ref | contreatelasticitytype() without contreat() is not inert — it relabels a plain level ATT as an elasticity |
| `hv-noparent-F3` | medium | PARTIAL | pr+1ref | semielasticity / contreatelasticitytype(eyex\|eydx) has no domain screen: a silent inference row, r(459), or an opaque Mata r(3598), depending on data… |
| `SK-SEMIELAS-R3598` | low | CONFIRMED x2 | pr+2late | semielasticity on the default engine lets margins' raw Mata traceback through: r(3598) 'variable __00001O not found' with no wooldid-level diagnostic |
| `SK-01` | none | PARTIAL x2 REFUTED x2 | pr+4late | semielasticity averages (dp/dx)/p over a subpopulation containing negative predictions and reports the result with se, t, p and a 95% interval, unscre… |

## `BSAVE` — bootsave()/permsave() (15)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `BR-4` | high | CONFIRMED x2 | pr+1ref+1lens | bootstrap() resamples the `smp` universe, not e(sample): units with no usable outcome count as resampling clusters, inflating e(boot_df) and silencing… |
| `BT-04` | high | CONFIRMED | pr+1ref | bootsave()'s _rep = 0 row for jointW_* holds an event-study coefficient, not the observed joint statistic |
| `HV-PERM-01` | high | PARTIAL | pr+1ref | The permutation design assigns onsets to clusters that contribute nothing, and the design is misreported to the user by four orders of magnitude |
| `hvr-bootsave-jointwobs-01` | high | CONFIRMED x2 | pr+2late | bootsave() under boottstat writes an observed COEFFICIENT, not the observed Wald statistic, into the _rep = 0 row of jointW_* and advW |
| `HV-BOOT-03` | medium | PARTIAL | pr+1ref | bootsave()'s observed row is on the wrong scale, contradicting the help by a measured factor |
| `NF-SKEP-01` | medium | CONFIRMED PARTIAL | pr+2late | saveplots() never got the pre-flight file-name check the author added to bootsave()/permsave(), so rerunning a do-file loses a completed run -- and un… |
| `SK-2` | medium | CONFIRMED x3 REFUTED | pr+4late | permute() builds its randomization design on the user-if universe too, so units with no estimable row are shuffled treatment onsets: 18 assignment slo… |
| `fpsg-02` | medium | CONFIRMED | pr+1ref | sthlp:1307 documents the wrong e(b)/e(V) row names for subgroups (_s# suffix vs the actual _SG<value>_ infix) |
| `hv-permsave-comma-01` | medium | CONFIRMED | pr+1ref | permsave() has no comma check; permsave(f, replace) writes a file literally named "f, replace.dta" at rc 0 |
| `BR-1` | low | CONFIRMED | pr+1ref | permsave() has no comma check and silently writes a file literally named "<name>, replace.dta" |
| `FCW-05` | low | PARTIAL | pr+1ref | Nothing in e(), the table or bootsave() records which bootstrap regime ran |
| `PERMREF-01` | low | CONFIRMED x2 | pr+2late | permsave()'s failed replication rows are not "all missing", though the help says they are |
| `PERMSEM-01` | low | PARTIAL x2 | pr+1ref+1lens | e(permdist_b) and e(perm_jointW)/e(perm_advW) are not row-aligned when replications fail |
| `SK-1` | low | CONFIRMED x3 PARTIAL | pr+4late | The RF-6 comma refusal for bootsave() prints a remedy that itself writes a silently mis-named file |
| `bb-02-bootsave-size-obsrow-wrong-N` | low | CONFIRMED | pr+1ref | bootsave() under bootalpha(size) writes 1/(rows in the user-if sample) on the observed row where the help promises 1/(estimation sample) |

## `PERM` — permute() (22)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `HV-PERM-03` | high | CONFIRMED | pr+1ref | permute() hard-refuses r(459) because of single rows wooldid itself drops, while the ordinary run on the same file succeeds |
| `SK-FC-C` | high | PARTIAL x2 | pr+2late | The joint pre-trend test is severely oversized in the exact geometry the help and the runtime caution name as its REMEDY, and prints nothing there |
| `fineclus-perm-FP1` | high | PARTIAL | pr+1ref | permonsetcohorts: a covariate keyed to the OBSERVED cohort is inert for the point estimate but rewrites the permutation null distribution (p moves 2x-… |
| `AR-1` | medium | PARTIAL | pr+1ref | replay drops the cautions that disown the numbers it reprints, including the runtime half of two shipped rulings (PS-64/BS-9 and AMESE R1) |
| `HV-PERM-02` | medium | PARTIAL | pr+1ref | permexact reports 'exact shares with no Monte-Carlo error' over a reference set it pruned, and the p moves 4x with a preprocessing step that cannot ch… |
| `NFC-01` | medium | PARTIAL | pr+1ref | With a continuous covariate in controls(), the default cell collapse and nofastcollapse return different VARIANCES (up to 5.5e-6 relative), reaching t… |
| `PI-1` | medium | PARTIAL | pr+1ref | A cluster whose outcome is missing on every row stays in the permutation design, so rows the estimator discards change the p-value and the interval |
| `SK-01` | medium | PARTIAL CONFIRMED | pr+2late | Under the default cell collapse a plain wooldid run is not reproducible across a re-sort of the same data: the printed SEs, t, p and CI bounds move by… |
| `SKEP-FLOOR-1` | medium | CONFIRMED x2 | pr+2late | The help's p-value MAP gives the wrong floor for permexact's default convention: it says 1/M, the measured floor is 2/M |
| `fineclus-perm-FP2` | medium | PARTIAL | pr+1ref | The i-keyed-covariate convention is documented and cautioned for the bootstrap re-id leg, and neither documented nor cautioned for the permonsetcohort… |
| `AR-2` | low | CONFIRMED | pr+1ref | replay is not byte-for-byte: it prints an extra blank line before the permutation/bootstrap inference block |
| `FCW-06` | low | PARTIAL | pr+1ref | T1's permute() relief does not reach the configuration T1's own motivation names |
| `FUZZ-06` | low | CONFIRMED | pr+1ref | dots(-5) is refused under bootstrap() and silently accepted under permute() |
| `PERMEX-PX-2` | low | PARTIAL x2 | pr+1ref+1lens | The permexact add-one caution describes a branch the code does not take: at complete-enumeration-but-not-a-member the share is plain and the p can be… |
| `PERMSEM-02` | low | CONFIRMED | pr+1ref | e(perm_sd) silently switches denominator inside one permexact table |
| `PX-02` | low | CONFIRMED | pr+1ref | customvce() can never reach a permutation replication: the customvce forwarding line in the permute rep string (ado:128) is dead code, and the seeded… |
| `SK-1` | low | PARTIAL CONFIRMED | pr+2late | permute()'s nesting refusal r(459) is diagnosed on rows that are absent from e(sample), with a count dominated by consistent live rows and no remedy f… |
| `hv-permsave-comma-03` | low | CONFIRMED | pr+1ref | dots(-5) is r(198) under bootstrap() and silently accepted under permute(), where it behaves as dots(0) |
| `hv-reps-02` | low | CONFIRMED | pr+1ref | The permute engine still prints its resolution caution for a run that then aborts — the very defect 813c414 claimed to fix, left unfixed in the other… |
| `skep-reps-N1` | low | PARTIAL x2 | pr+2late | RF-2's own rule is still broken in BOTH engines: the low-draw caution precedes the caller's r(459), the very shape 813c414 was written to fix |
| `skep-reps-N2` | low | CONFIRMED x2 | pr+2late | wooldid.sthlp:636 says 'neither engine has a hard minimum' three lines after saying permute() requires permreps() >= 20 - and the code errors r(198) |
| `SKEP-N3` | none | REFUTED x4 | pr+4late | permstrata() strata are built on the same unscreened cluster set |

## `PTP` — pretrendpower() (14)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `BOOTX-1` | high | PARTIAL | pr+1ref | pretrendpower()'s joint pre-trend test under bootstrap() is an uncorrected k-dimensional Wald referred to F(k, e(df_r)) instead of F(k, n-k) |
| `BOOTX-2` | high | PARTIAL | pr+1ref | pretrendpower()'s rank guard stops firing at completed = k+1, where the reported pre-trend slope flips sign, its CI excludes zero on the wrong side, a… |
| `PTP-1` | high | PARTIAL | pr+1ref | pretrendpower()'s minimum detectable trend does not have the stated power on a 16-cluster state-level DiD, and no caution fires there |
| `BX-N1` | medium | CONFIRMED x2 | pr+2late | Two differently-referenced "joint pre-trend test" p-values print in the same bootstrap run, and under boottstat the only one printed is the one the he… |
| `BX-N3` | medium | PARTIAL CONFIRMED | pr+2late | The singleton-treated-cell degeneracy reaches plain robust, rank-deficits the analytic event-study covariance by one with an exact unit correlation be… |
| `EC-4` | medium | PARTIAL CONFIRMED | pr+2late | e(V_es_<fam>) is singular with a strictly positive diagonal, and e(es_ndegen_<fam>) is a diagonal-only scan, so the help's "its absence marks a clean… |
| `PTP-2` | medium | CONFIRMED | pr+1ref | e(roth_nissize_*), documented as 'the actual size' of the NIS pre-test, evaluates a normal probability at a t critical value |
| `PTP-3` | medium | CONFIRMED | pr+1ref | Under poisson, pretrendpower() takes the chi-squared shortcut its own comment calls optimistic, even with 16 clusters, and sthlp:789's description of… |
| `PTP-5` | medium | PARTIAL | pr+1ref | At cluster() == i the MDT appears to deliver its 80% power only because the joint pre-trend test rejects 28% of the time under the null |
| `SKP-NEW-1` | medium | PARTIAL x2 | pr+2late | pretrendpower() uses known-covariance references throughout, so the joint pre-trend test's size, the GLS slope CI's coverage and the mdt_wald column a… |
| `hv-adv-ame-rank` | medium | CONFIRMED | pr+1ref | The adversarial ame p-value silently drops every pre-treatment horizon: it is bitwise identical to the statistic computed from the three post-treatmen… |
| `BX-N2` | low | PARTIAL x2 | pr+2late | The Monte-Carlo component of the printed Roth slope, its CI, mdt_wald and the joint p is unrecorded, and at reps(50) it is about half the printed stan… |
| `PTP-4` | low | CONFIRMED | pr+1ref | The help's only worked example of pretrendpower() puts the user on the geometry the run then says is not a calibrated test |
| `SKP-NEW-2` | low | PARTIAL x2 | pr+2late | Under bootstrap(), pretrendpower()'s reference distribution is keyed only on bootdf() and never reflects reps(); under the default bootdf(normal) it i… |

## `OTHER` — other / single-surface (54)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `SK-FC-A` | high | PARTIAL CONFIRMED | pr+2late | The documented remedy geometry (coarsen `i` so that cluster() is finer than `i`) silently biases the POINT ESTIMATE on an unbalanced panel |
| `hv-honestdid-Minvariant` | high | PARTIAL | pr+1ref | The advertised honestdid handoff on e(b_es_ame)/e(V_es_ame) returns a robust confidence set IDENTICAL at every M - it reports the continuous-treatment… |
| `AR-3` | medium | CONFIRMED | pr+1ref | the replay-staleness guard is self-poisoning: a bare `wooldid` typed when another command owns e() errors about a missing varlist AND records itself a… |
| `BR-2` | medium | CONFIRMED | pr+1ref | The reps() resolution Caution is unreachable without cluster(), contradicting the help's unconditional promise |
| `EC-6` | medium | PARTIAL CONFIRMED | pr+2late | [aw=] and [pw=] leave no trace in e(), and e(vce) reads 'ols' under pweights although the standard errors ARE robust |
| `FCM-3` | medium | CONFIRMED | pr+1ref | Nothing anywhere records which of the two bootstrap unit regimes ran |
| `FCW-03` | medium | PARTIAL | pr+1ref | The 149/150 over-coverage at cluster()==i is an ADDED cohort-composition term, not the effect the record names |
| `FCW-04` | medium | PARTIAL | pr+1ref | Regression: the reps() p-value-resolution caution is unreachable without cluster(), on a path T2 just opened |
| `FUZZ-04` | medium | PARTIAL | pr+1ref | estimates replay after wooldid fails r(198) and blames a command line the user never typed |
| `FUZZ-05` | medium | PARTIAL | pr+1ref | The requested event-study PRE window shifts the headline ATT by 5-11% through the reference-normalization pool |
| `HV-BOOT-02` | medium | CONFIRMED | pr+1ref | bootalpha(size): the Dirichlet concentration is read off the unscreened frame, understating the bootstrap SE by a further sqrt(live fraction) |
| `HV-EPOST-02` | medium | PARTIAL | pr+1ref | The replay-staleness guard is unreachable for exactly the failure class it was written for, and the user gets an error about a command line they never… |
| `HV-VAL-01` | medium | CONFIRMED | pr+1ref | wooldid refuses runs r(411) because of rows it would itself have dropped, and the message names no remedy |
| `MC-FC-1` | medium | PARTIAL | pr+1ref | The finer-than-i remedy is incomplete: the analytic-SE understatement is a gradient in clusters-per-treated-cell, not a binary in the nesting relation |
| `MC-FC-2` | medium | PARTIAL | pr+1ref | sthlp:251's "roughly 0.66 to 0.75" is not the size of the at-or-above-i SE understatement on realistic data; it runs to 0.18, and coverage to 0.29 |
| `MCBP-BAYES` | medium | PARTIAL CONFIRMED | pr+2late | bootstrap(bayes) at G=10 rejects a true sharp null 16-20% of the time and its SE is statistically indistinguishable from the analytic SE it is offered… |
| `MCBP-BOOTSTD` | medium | PARTIAL CONFIRMED | pr+2late | bootstrap(standard) at G=10 halves the level distortion but does not remove it (.122 at nominal 5%) |
| `MCBP-ESFBP` | medium | PARTIAL x2 | pr+2late | The analytic att SE understatement is estimator-specific, and the lever that fixes it is `esfixedbaseperiod` — an event-study display option, not an i… |
| `SK-1` | medium | CONFIRMED x3 PARTIAL | pr+4late | A stale $WOOLDID_LEVEL0 replaces the true return code with 198 and disarms the shell's rc-gated cleanup epilogue on every later parse-stage failure |
| `SK-FC-D` | medium | CONFIRMED PARTIAL | pr+2late | A COARSER cluster() can lower the reported standard error, reversing the conservatism every user expects from over-clustering |
| `SKEP-N6` | medium | PARTIAL CONFIRMED | pr+2late | The worker's own onset-uniqueness screen IS markout-aware, so wooldid contradicts itself about which sample a screen should use - and the pattern for… |
| `SKEP-NEW-1` | medium | PARTIAL x2 | pr+2late | esfixedbaseperiod silently deletes any treated unit that does not OBSERVE the exact relative period named by esrelativeto() |
| `bb-03-bayes-coverage-at-cluster-eq-i` | medium | PARTIAL | pr+1ref | At cluster()==i the bayes percentile interval for att missed 7 of 30 at nominal 95% where bootstrap(standard) missed 2 on the same panels at equal wid… |
| `bbr-N2-documented-se-range-is-the-homogeneous-case` | medium | CONFIRMED x2 | pr+2late | The only quantification on record for the analytic SE understatement at cluster()>=i (0.66-0.75) is the no-unit-heterogeneity case: the reported se is… |
| `hv-cwm-01` | medium | PARTIAL | pr+1ref | customweightmultiplier() silently reweights the treatment cells, so att_it / att_i / att_itime stop being the estimands their help entries define |
| `hv-reps-01` | medium | CONFIRMED | pr+1ref | The reps() p-value-resolution Caution is unreachable without cluster() — a 2026-09-13 regression, and a doc-vs-code divergence |
| `hv-sgname-02` | medium | CONFIRMED | pr+1ref | ES_att_i_* and ES_att_itime_* do not exist; sthlp:1305-1306 contradicts sthlp:1315 in the same help file |
| `hvS-cwm-miss-01` | medium | CONFIRMED PARTIAL | pr+2late | A missing customweightmultiplier() removes the row from the REGRESSION, and a participants-only multiplier makes wooldid state something false about t… |
| `FP-CWM-DOC` | low | CONFIRMED | pr+1ref | The att_it/att_i/att_itime weight formulas in the help omit customweightmultiplier(), while still calling the result a simple unweighted mean |
| `FP-NOPRE-LATEENTRY` | low | CONFIRMED | pr+1ref | The no-pre-period drop is per-unit against that unit's own first observed period, not against the start of the sample, and the help says otherwise |
| `FUZZ-07` | low | CONFIRMED | pr+1ref | `update` on a full command line runs the package installer and silently estimates nothing |
| `FUZZ-08` | low | PARTIAL | pr+1ref | Reversing the sort order of the data moves analytic SEs in the 8th significant digit |
| `FW-03` | low | PARTIAL | pr+1ref | The fweight refusal prescribes aweights or pweights, neither of which reproduces fweight inference on collapsed data, and neither the message nor the… |
| `GTE-1` | low | PARTIAL | pr+1ref | customweightmultiplier() silently turns att_it / att_i / att_itime from the documented simple means into weighted means |
| `HV-BOOT-07` | low | CONFIRMED | pr+1ref | bootdisplay(bca): the delete-one-cluster jackknife runs one fit per FRAME cluster, so dead clusters enter the acceleration |
| `HV-FN3` | low | CONFIRMED | pr+1ref | e(numpre_ame) still posts when the void screen withholds e(b_es_ame) |
| `MCBP-BOOTTSTAT` | low | PARTIAL CONFIRMED | pr+2late | boottstat is correctly sized but pays for it with intervals ~2.7x too wide at the median and 22x at the worst draw, on a path whose coverage the help… |
| `MCBP-S7COVER` | low | PARTIAL x2 | pr+2late | Contradicts a §7 known_bydesign entry: bootstrap(standard) at cluster()==i does NOT over-cover att_i on a neighbouring geometry |
| `NF-SKEP-02` | low | CONFIRMED x2 | pr+2late | The replay guard tells the user that complete, valid results of the failing call 'predate that failed call' |
| `SKEP-N3` | low | CONFIRMED x2 | pr+2late | The explanation the 2026-09-16 prereg retracted 'so it does not propagate' is still standing in the file that records the measurement |
| `hv-bootrank-02` | low | CONFIRMED | pr+1ref | Three code comments and the PS-62 ruling record justify `ereturn repost` with an `ereturn post` behaviour that does not reproduce — and on the one rea… |
| `hv-cwm-02` | low | CONFIRMED | pr+1ref | REFUTED: the bootctl cell-weighted gate's reasoning is correct; the recon's "either cwm is a no-op or the gate is wrong" is a false dichotomy |
| `hvnpr-N1` | low | CONFIRMED x2 | pr+2late | cfxplottypes() token validation is gated behind makecfxplots, so the VAL-G13A fix does not apply without the parent and sthlp:974's unconditional reje… |
| `hvnpr-N2` | low | PARTIAL REFUTED | pr+2late | esrelativeto() accepts a positive value - which its own gated error calls impossible - whenever esfixedbaseperiod is absent |
| `hvr-dots-sentinel-01` | low | CONFIRMED x2 | pr+2late | dots(-1) collides with the internal not-given sentinel: it is silently reinterpreted as the default and escapes the stray-companion screen |
| `mcci-03-esprelength-moves-att` | low | PARTIAL x2 | pr+1ref+1lens | Requesting event-study horizons silently moves the headline att and inflates its confidence interval by up to 60% |
| `sk-amename-02` | low | PARTIAL REFUTED | pr+2late | sthlp:239 makes the same false event-study naming promise for the ame family (ES_ame_i / ES_ame_itime do not exist) -- a third help location, extendin… |
| `sk-handoff-asymmetry` | low | REFUTED PARTIAL | pr+2late | Of the two third-party handoffs the help advertises on e(b_es_*)/e(V_es_*), pretrends refuses a degenerate pre block and honestdid does not |
| `state-leak-02-verbose-leaves-variable` | low | PARTIAL | pr+1ref | `verbose` adds a variable to the user's dataset and leaves it there |
| `EC-8` | none | REFUTED x2 | pr+2late | After any failed wooldid, `estimates restore` gives back a valid result set that wooldid then refuses to replay, with a false explanation |
| `MCBP-JOINTSIZE` | none | REFUTED x2 | pr+2late | The analytic joint horizon test rejects a true sharp null 75.5% of the time at nominal 5% — far outside the 11-36% oversizing the project documents |
| `SKEP-N5` | none | PARTIAL REFUTED | pr+2late | NEGATIVE RESULT: an event-study window does not drop rows, so one of the three triggers ado:9367 names does not exist |
| `fp-contreat-D` | none | REFUTED x2 | pr+1ref+1lens | fe() with customvce()/unconditionalse: wooldid cautions, then enters the fe() variables as CONTINUOUS regressors and completes with rc 0 |
| `state-leak-03-mata-namespace-leftovers` | none | REFUTED | pr+1ref | A plain run leaves three globals in the user's Mata workspace, one of them named `k` - attribution to wooldid NOT established |

## `LAT` — lattice()/contreat() (10)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `fp-contreat-B1` | high | CONFIRMED PARTIAL | pr+1ref+1lens | lattice() + contreat(): the all-zero covariance is an artifact of the absorbing engine, not of the design -- plain regress on the identical fit return… |
| `fp-contreat-B2` | high | PARTIAL x2 | pr+1ref+1lens | On the same lattice fixture the ame point estimate is 45x the truth, and no path gives the reader a usable warning: the default path has no standard e… |
| `fp-contreat-C` | medium | CONFIRMED | pr+1ref | The two all-zero-VCV cautions blame analytic clustered standard errors that were never requested and point the reader at a caution that did not print |
| `hv-lat-01` | medium | PARTIAL | pr+1ref | The maxvar preflight omits lattice() and contreatpoly(), so its bound is short by a multiplicative factor and a textbook continuous-treatment run dies… |
| `hv-lat-02` | medium | PARTIAL | pr+1ref | The omission is a FAMILY, not a lattice special case: interactivecontrols() and coarsecohortcontrols() break the same bound with no lattice at all |
| `hv-noparent-F1` | medium | PARTIAL | pr+1ref | Nine options are accepted, do nothing and say nothing without their parent — contradicting the help's own stated policy |
| `BOOTX-4` | low | CONFIRMED | pr+1ref | bootstrap() + timetrends + cluster() at the level of i: the printed table is fully rescued with bootstrap inference but e(V) is never posted, and the… |
| `SKEP-INERTPOLICY` | low | PARTIAL x2 | pr+2late | Inert options are announced inconsistently: one prints an 'is being ignored' Caution, six are silently accepted, four are refused outright |
| `hv-lat-03` | low | CONFIRMED | pr+1ref | REFINEMENT of the seeded claim's UX half: 'no wooldid message' is literally true, but 'uninformative' overstates it -- Stata's own r(907) text names t… |
| `hv-lat-04` | none | REFUTED | pr+1ref | The granular-i caution tells the user runtime, not the variable limit, is what will bite under lattice() -- on both designs measured the variable limi… |

## `GTE` — gteffects (6)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `SKEP-NEW-2` | high | CONFIRMED x2 | pr+2late | gteffects reports numerically-zero cell variances as real standard errors: se ~ 1e-15, p down to 8e-258, confidence intervals of width ~1e-14 |
| `bbr-N1-gteffects-cell-ci-covers-45pct` | high | CONFIRMED PARTIAL | pr+2late | e(gteffects) cell standard errors, p-values and confidence intervals are understated ~3.5x at cluster() at-or-above i: the cell's printed 95% interval… |
| `GTE-3` | medium | CONFIRMED | pr+1ref | A treated i that enters an unbalanced panel at or after its own onset is dropped silently, and the help describes the rule at the sample level |
| `GTE-4` | low | CONFIRMED | pr+1ref | e(gteffects) weight columns are populated for estimands that were not requested |
| `GTE-5` | low | CONFIRMED | pr+1ref | "under timetrends ... the se column is missing" is not what happens on a cohort-level i with the default VCE |
| `FASTNOSE-02` | none | REFUTED | pr+1ref | fast mode prints "SEs and related materials are not estimated in fast mode" directly above a full gteffects table of standard errors, p-values and 95%… |

## `BAND` — simultaneousbands / interval shape (5)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `BT-01` | high | PARTIAL | pr+1ref | simultaneousbands under boottstat yields a band NARROWER than the pointwise percentile-t interval, which sthlp:525 says is impossible |
| `BT-03` | medium | PARTIAL | pr+1ref | The simultaneous band's critical value never sees bootdf(), so under bootdisplay(wald) -- where the pointwise interval IS symmetric -- the band still… |
| `SK-02` | medium | PARTIAL x2 | pr+2late | The boottstat simultaneous critical value has no guard against a near-zero replication standard error, and the printed band can be 1e14 wide with no c… |
| `mcci-02-band-z-not-t` | medium | CONFIRMED | pr+1ref | simultaneousbands: the analytic sup-t band covers 0.752 at nominal 0.90 and 0.843 at nominal 0.95 with 12 clusters |
| `BOOTX-5` | low | PARTIAL | pr+1ref | the bootstrap simultaneous band's 95% sup-t critical value is literally the maximum of the completed draws at low reps, and no caution is keyed on the… |

## `LEVEL` — level() (4)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `AR-4` | high | CONFIRMED | pr+1ref | level(): 99.9 and 99.99 (the help's own upper bound) are refused r(198); the refusal also fires with no level() option, reports a SUCCESSFUL run as a… |
| `EC-2` | medium | CONFIRMED x2 | pr+2late | level(99.9) - and a session `set level 99.9` with no level() option at all - make wooldid exit r(198) with a message that contradicts what the user ty… |
| `FUZZ-03` | medium | PARTIAL | pr+1ref | level(99.9) and level(99.99) refused with a message that contradicts both the input and wooldid's own range check |
| `SKEP-LEVELSTATE` | medium | CONFIRMED x2 | pr+2late | wooldid re-renders the confidence level at 16 significant digits and feeds it to set level, so ~13% of legal session levels make every call fail r(198… |

## `VCV` — e(V) construction / rank (9)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `skep-dfmismatch-03` | high | PARTIAL x2 | pr+2late | e(V) carries e(boot_N_V)-1 degrees of freedom but is referred to e(df_r)=G-1; on a plausible spec N_V is 4 of 60 and a nominal-95% interval covers 86.… |
| `FW-02` | medium | PARTIAL | pr+1ref | pweights produce a robust covariance but the run is labelled OLS: e(vce)="ols" and e(vcetype) unposted, contradicting wooldid's own stated contract at… |
| `sk-prefix-01` | medium | CONFIRMED x2 | pr+2late | wooldid's estimand names are prefix-nested, so Stata's default varabbrev silently redirects _b[main_att] to main_att_i -- a 2.55x wrong number, sign-f… |
| `EC-1` | low | PARTIAL REFUTED | pr+2late | The bootstrap e(V) is rank-deficient at EVERY reps(), and at low reps `test` silently drops constraints and returns a wildly anti-conservative p |
| `SKEP-01` | low | CONFIRMED x2 | pr+2late | On the unconditionalse path e(vce)/e(vcetype) name the stage-1 regression's vce, not the covariance e(V) actually is: one matrix, two names |
| `hv-bootrank-01` | low | CONFIRMED | pr+1ref | SEEDED CLAIM: facts confirmed, stated mechanism REFUTED — rank tracks the draw count, not the `ereturn repost` bypass |
| `hv-bootrank-04` | low | CONFIRMED | pr+1ref | The footer line that would disclose e(boot_N_V) prints the token, not the value — `e(boot_N_V)` sits inside the string literal |
| `hv-noVsentinel-marker` | low | CONFIRMED | pr+1ref | No e() scalar records the e(V) sentinel count, so a replayed or restored estimate cannot report it |
| `skept-roworder-01` | none | REFUTED x2 | pr+2late | The same rows in a different physical order change an event-study standard error by 3.7e-08 relative, on the reghdfe path with continuous weights only |

## `FN` — fastnose / fast-mode (6)

| id | sev | verdicts | reviewed by | title |
|---|---|---|---|---|
| `FASTNOSE-01` | medium | CONFIRMED | pr+1ref | permanalyticse is a complete no-op under fastnose, and the run prints a header and a footer asserting the analytic standard errors are in the table an… |
| `HV-FN2` | medium | PARTIAL | pr+1ref | Under fastnose the ame dose-slope void screen never runs: wooldid prints and now returns eight fabricated zeros the analytic twin refuses, under a cau… |
| `PX-01` | medium | CONFIRMED | pr+1ref | fastnose silently overrides permanalyticse under permute(), and the run then asserts in three places that the analytic SEs and e(V)/e(V_es_*) are ther… |
| `RF-N1` | medium | CONFIRMED PARTIAL | pr+2late | permute()'s fast-mode default silently disables the only diagnostic wooldid owns for "this control has no within-cohort-period variation" |
| `SK-02` | medium | PARTIAL x4 | pr+4late | A missing coefficient anywhere in the estimand vector costs the entire e(b): the run returns rc 0 having posted no coefficients, silently under fastno… |
| `GTE-2` | low | CONFIRMED | pr+1ref | fastnose DOES post e(b_es_att)/e(numpre_att) - sthlp:1325 is accurate and the 2026-09-13 ruling item 11 was not built |
