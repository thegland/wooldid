# Change audit, 2026-09-22 - the record for the author's rulings

**Status line, dated: opened 2026-09-22 late evening; every ruling in sections 2 and 2b was built by 2026-09-23 (shas in the state column).** This file
holds (1) what was built tonight from the open tickets, (2) the rulings the author owes, and
(3) the audit's confirmed findings once they land. Straight bugs the audit confirms are fixed
in their own commits and listed here with the sha; anything that moves a reported number, or
is a design choice, waits here for a ruling. No caution or error wording ships from this
review: proposals go to the review page as a collapsed "shorter wording" suggestion on the
card, and the author writes the text.

## 1. Built tonight

| what | where | verification |
|---|---|---|
| Help and guide no longer claim a caution marks the histogram outputs | 76f2f38 | text only |
| Joint-test screen judged on treated rows; bare-dot equality line; one-cluster poisson guard | 52fd319, 0e7dc0d, bcffd07 | three cert legs; battery IDENTICAL |
| reps() caution on unclustered runs; saveplots() pre-check; lattice labels on engine cfx figures; plot-data frame and update fixes | 2be7e16, 54440de, 5ab3e0b, e203b93 | cert legs G39f, G41, L7-L9; battery IDENTICAL |
| Bootstrap re-id leg: covariates that name i are re-keyed to the fresh unit code (open_ticket_bootreid_covariates) | c24ab48 | gates m10-m13; six bootstrap suites; battery IDENTICAL |
| Save-failure ticket residue: guarded test reads, forced-failure gates G38/Q37, help and guide clauses (open_ticket_savefail_caution) | (sha filled in below) | gates + permute3 |

**What the bootreid probe established before building.** Drawn copies of a cluster are exact
duplicates, so pooled and separate fixed effects give the same fit: on the percentile paths the
substitution changes the draws by at most 3.6e-13 (noise), and the reported numbers of every
existing battery config are unchanged. Under boottstat the replication standard errors move,
design-dependently: a per-unit slope entered as a regressor, controls(c.x#i.unit), moved them on
both probe designs; an absorbed unit trend, fe(unit#c.t), moved them on the 24-unit design and
not at all on the 12-unit one. So the number-moving footprint is exactly: boottstat runs where
cluster() nests i and a covariate names i. No shipped config takes that path.

## 2. Rulings the author owes

| # | question | the two readings | where it stands |
|---|---|---|---|
| R-A | Cards #29 (impvalues with a finer or cross-cutting cluster) and #214 (bootstrap(bayes) with cell-weighted estimands): you unticked "Silenced by safety(off)" in round 3 | build both (the untick is the ruling; test and help changes ride along: test_fastcollapse F15, test_bootstrap_gates G6, test_bootstrap_std:77, test_bootstrap_joint:152/229, test_bootstrap_bayes:332, help ~241/250/268/357, guide ~177; #214 also kills caution #215) / keep the escapes (safety(off) stays the documented way past both) | ruled 2026-09-23, see 2b |
| R-B | Subgroup containment guards S1 and S2 (the worker's split-unit and within-cell errors): the ticket says harden all three; your page edits hardened S3 only and kept "specify safety(off)" in S1/S2's text | harden S1/S2 too (their escape sentences go) / leave them liftable (ticket closes as S3-only) | ruled 2026-09-23, see 2b |
| R-C | Save failure, Q2: each saver re-runs the directory-part r(198) and file-exists r(602) checks after the loop; they are still fatal, so a file that appears mid-run still throws the finished run away | convert the r(602) leg to the same caution / delete the duplicates as redundant / leave | ruled 2026-09-23, see 2b |
| R-D | Save failure, Q3: the caution prints above the deferred display and scrolls off | move it below the tables / leave | ruled 2026-09-23, see 2b |
| R-E | Card #409: a caution when a bootstrap run completes more than 10,000 replications and e(bootdist) is silently not posted, the twin of #408 on the permutation side | approve the drafted wording (then it is built) / edit it / delete the card and leave the bootstrap side silent | ruled 2026-09-23, see 2b |

## 2b. Rulings given 2026-09-23 (multiple choice, in session)

| item | ruling | state |
|---|---|---|
| R-A #29 impvalues finer/crosscut cluster; #214 bootstrap(bayes) with cell-weighted estimands | **Build both**: safety(off) no longer lifts either; caution #215 goes | built 61f246e |
| R-B subgroup guards S1, S2 | **Leave them liftable**; ticket closes as S3-only | done (ticket text) |
| R-C save failure Q2, the post-loop duplicate file-exists check | **Leave it fatal** | no change |
| R-D save failure Q3, caution position | **Move below the tables** | built 61f246e; refined by 98ab2fb (save first, caution last) |
| R-E card #409, bootstrap draw-cap caution | **Approve Claude's wording** | built 61f246e; moved below the tables in 98ab2fb |
| esband-normal-law-vs-t-pointwise | **Multivariate t with df_r** (the floor at the pointwise t quantile was built in 6d8bedd, then reversed; see below); re-bless band configs | built 6d8bedd |
| ptp-nis-tcrit-normal-law (blue-sky PTP-2) | **Roth: normal throughout**; help adds "under the normal reference" | built 6d8bedd; qualifier c7bad2c |
| permtstat-cfxplotnose-blank | **Override cfxplotnose**, caution approved verbatim: "Caution: cfxplotnose was dropped; permtstat needs the plot cells' standard errors." | built a45480b |
| sgdeg-inert-under-suppresssg | **Make the guard work** under suppresssgprimaryeffects | built 637062e |
| amedetect-onerow-cell | **Exclude one-row cells** from the containment predicate | built b6f1a54 |
| ccscreen-contreatwithin-amevoid | **Run the void screen under contreatwithin** | built b6f1a54 |
| ccref-nevergroup-false-refusal | **Keep the refusal** | no change |
| perm-window-labels | **Recode late onsets in the permutation engine, and caution about it** (wording to the author first) | built 8443694 |
| val-permanalyticse-fastnose | **Refuse the pair** (widened message to the author first) | built 8443694 |
| fcol-lattice-cfx-key-liftable | **Lift the refusal**: the bin joins the collapse key | built f1a2b94 |
| perm-window-labels caution text | approved verbatim: "Caution: `n' cluster(s) with onsets after the last period are permuted as never treated." | built 8443694 |
| val-permanalyticse-fastnose text | approved: "Error: fastnose is not compatible with permtstat, permtci or permanalyticse." | built 8443694 |
| perm-save-comma | **Add bootsave's refusal**: "Error: permsave() should only contain a file name." | built 8443694 |
| permctl-partial-label-no-twin | **Add the screen, unliftable**, text approved: "Error: the subgroup variable `subgroups' is missing for some in-sample row(s) of units of `i' that carry a subgroup value on their other rows. This plausibly is an error, and is not compatible with non-att estimands." | built 8443694 |
| perm-esttab-analytic-under-ase | **Always permutation**: e(p)/e(ci_l)/e(ci_u) posted from the kernel on every permute() run | built a45480b |
| perm-sgeq-uncentred | **List flagged contrast rows in e(perm_offcentre)** (existing footer names them); the equality test's p = 1 goes to research | built a45480b; research answered (sharp null), labelled in 66d6680 |
| saveplots-no-precheck | **Keep Stata's own r(602) text** for the early refusal | built 54440de (before the ruling; it already used Stata's text) |
| fcol-customvce-holdout-liftable | **Measure first** (one settling run), bring the result back | measured; ruled below |
| customvce() collapse (after measurement) | **Lift it, fix the help**: its differences are margins' own finite-difference noise, which the uncollapsed and eligible cluster() runs share; the help's machine-precision sentence is corrected | built cd46002 |
| permute() equality test and pairwise contrasts | **Label as the sharp null now**; texts approved: equality line under permute() "permutation p-value (sharp null of no effect at any horizon) for joint test that all SGs have the same `item': "; footer "Pairwise subgroup contrasts are permuted as coefficient differences (in both modes) and, like the equality tests, test the sharp null of no effect at any horizon, so a common nonzero effect alone can reject them; their kernel is e(perm_pwsgs_kernel)." | built 66d6680 |
| permonsetcohorts after the late-onset recode | **Keep the old check at cluster() = i**, where no rebuild happens | built 66d6680 |
| r(498) in wooldid_cfxstripe on all-omitted dose slopes | **Fix the parser** | built 143e457 |
| save order | **Save first, caution last**: the draws file is written before the plots, only its failure caution prints below the tables; #409 moves below the tables to match #408 | built 98ab2fb |
| pretrendpower() note | **Add the qualifier**: "Both pre-tests are at the 5% level under the normal reference, whatever level() is in force; level() moves the slope CI only." | built c7bad2c |
| analytic band floor | **REVERSED by the author the same day: no floor on any path.** "Better users see the SB approach isn't working as expected than a more sensible but fundamentally different type of result." The floor 6d8bedd added to the t branch is removed; the multivariate t law stays | built c7bad2c (floor removed) |
| D1 fastnose-nohorizon-caution-spam (card #52) | **Silence it under fastnose**; it still prints when a real variance is degenerate | built f5cd92b |
| D2 val-sg-namelength-late-crash | First ruled **dense internal index**; re-ruled **keep the caution** once Claude corrected the size (about forty stored-result names, the ES file names and both engines' resampling carry the value) | no change |
| D3 val-i-fvcap-unscreened | First ruled **dense internal index**; re-ruled **refuse in validate** once Claude corrected the size; text approved: "Error: `i' must not exceed 2,147,483,647 on treated units; recode it first (egen newi = group(`i'))." | built 5fc3ce5 |
| D4 permute-pairwise-se-analytic-undisclosed | **Reuse the analytic path's degeneracy detection** under permute(), disclosure only; text approved: the existing caution's first sentence alone | built 189af00 |
| D5 cfxplotnose-engine-asymmetry | **Extend permute()'s carve-out to bootstrap()**: the cells keep six columns and get bootstrap inference | built 516038c |
| D5 follow-up: boottstat with cfxplotnose (the carve-out made caution #225 fire on a new path) | **Drop cfxplotnose under boottstat**, the twin of the permtstat ruling; text approved: "Caution: cfxplotnose was dropped; boottstat needs the plot cells' standard errors." | built 7849439 |
| D3 follow-up (found by Claude's advisor): the cap screen refused a treated row whose i is missing, under the cap text | straight bug: the screen reads only non-missing i, as before 5fc3ce5 | fixed 07f63ae |

## 3. The audit

**What it covered.** Sixteen reviewers each took one part of the code at 8965f20 and read our
changes since v1.3 (main, 301b706), with Stata probes. They filed 152 findings: 11 high, 39 medium,
102 low. Every high and medium finding went to a separate verifier told to refute it: 29 were
confirmed, 20 confirmed in part, 1 refuted. The low findings were not verified and nothing was done
about them; they are in `change_audit_2026-09-22/findings.json` with every reviewer's evidence.

**Numbers that moved.** Only the ruled changes moved a battery config, and 929ae6c re-blessed the
five of them: band_analytic (the band's t law), ptp_two (the NIS pre-test on the normal reference),
perm_tstat (e(p), e(ci_l), e(ci_u) now posted), and cont_lat and perm_cfxnose (the collapse now
covers lattice() and makecfxplots; one coefficient moves by 3.6e-15). Every later commit left the
battery identical, and the customvce() collapse lift moved no config.

### 3a. What happened to each verified finding

| finding | severity | verdict | what the reviewer found | outcome |
|---|---|---|---|---|
| amedetect-onerow-cell | high | confirmed | amescreen detect treats a one-row treated cell as cluster containment and withholds every ame standard error, unliftably | built b6f1a54 (ruled) |
| ccref-nevergroup-false-refusal | high | in part | ccref screen refuses byTreated/byTtre whenever the never-treated group has a zero or constant dose, though that group's baseline slope is inert | ruled: keep the refusal; no change |
| esband-normal-law-vs-t-pointwise | high | confirmed | Analytic simultaneous band is simulated from a normal law while the pointwise intervals it is drawn over use t(df_r): with few clusters the band is narrower than the pointwise interval and under-covers | t law 6d8bedd; floor removed c7bad2c (ruled) |
| fastnose-binary-sg-horizon-vce-cost | high | in part | Fast mode keeps a full delta-method VCE on every binary subgroup-by-horizon margins call (and so does every permute replication) while the contreat twin synthesizes the marker | no cheap mechanism found; FYI |
| perm-window-labels | high | in part | Onsets beyond the sample window are distinct labels in the permutation design although the estimator recodes them to never-treated: M overstated, permexact refused on a working design, permonsetcohorts refused on standard not-yet-treated coding, the few-treated caution silenced | built 8443694; permonsetcohorts check 66d6680 (ruled) |
| permctl-partial-label-no-twin | high | confirmed | permute() screens the design population for four subgroup conditions but not partial labelling; replications then share-weight a partially labelled unit silently | refused, 8443694 (ruled) |
| permtstat-cfxplotnose-blank | high | confirmed | permute() + permtstat + cfxplotnose leaves every cfx cell's inference blank with no message | built a45480b; parser 143e457 (ruled) |
| saveplots-no-precheck | high | in part | saveplots() is the one saver with no pre-check: a stale file is found after the replications and the run is lost | fixed 54440de (Stata's own r(602) text, ruled) |
| sgdeg-inert-under-suppresssg | high | in part | The BS-15 degenerate-contrast guard is silently inert under suppresssgprimaryeffects | built 637062e (ruled) |
| sgdom-safetyoff-withincell-caution-false | high | confirmed | safety(off) skips the within-cell screen and prints the over-time caution, which claims att/att_it are unaffected | help fixed 29b5952; caution text: note on card #368 |
| u16-perm-atti-r198 | high | in part | Battery pin perm_atti asserts a permreps() refusal, not the att_i-under-permute() run its commit claims | fixed aae69dd |
| adv-boot-remedy-measured-inert | medium | refuted | Help and guide name bootstrap() with boottstat as the remedy for the uncalibrated adversarial p under containment, the one geometry where the author's own MC found the bootstrapped statistic inert (0/40 power, 0/200 size), and bootjnt posts it there without a caution | refuted; nothing to do |
| bootctl-respcaution-clusteronly | medium | confirmed | reps() resolution caution is nested inside the cluster() block and never fires on an unclustered bootstrap | fixed 2be7e16 |
| ccscreen-contreatwithin-amevoid | medium | confirmed | The dose-slope (amevoid) screen never runs under contreatwithin, so the ame rows print structural zeros with no caution naming them | built b6f1a54 (ruled) |
| cfxplotnose-engine-asymmetry | medium | in part | cfxplotnose: permute() keeps the 6-column cfx shape and fills it from replications; bootstrap() truncates to 1 column and the cell rewiring skips it | guide sentence 3b5712a; carve-out extended to bootstrap(), 516038c (ruled D5) |
| cfxplots-userlattice-axis | medium | confirmed | Engine-path cfx figures lose the user lattice() variable's name and value labels | fixed 5ab3e0b |
| containment-remedy-fails-boottstat | medium | confirmed | The cell-containment remedy 'bootstrap()' gives the ame rows no inference under boottstat | guide fixed 6fe9494 (help was already right); caution text: note on card #383 |
| customvce-collapse-heldout | medium | in part | customvce() is held out of the default cell collapse 'pending measurement'; measured, it is as exact as cluster() | measured, ruled, built cd46002 |
| display-eqtest-bare-dot | medium | confirmed | A withheld across-subgroup equality test prints as 'p-value for joint test that all SGs have the same att: .' | fixed 0e7dc0d |
| docs-bootstrap-reid-covariate-stale | medium | confirmed | Guide's bootstrap() block says i-named covariates keep their codes and are cautioned; help, guide ch3 and the code re-key them | docs fixed 4a6a5b9 |
| docs-collapse-ineligible-list | medium | confirmed | Guide's nofastcollapse entry enumerates two ineligible conditions and calls itself the complete list; help, guide ch11 and the code have seven | docs fixed 4a6a5b9; customvce() left the list in cd46002 |
| docs-fastnose-bootstrap-refusal | medium | in part | Help and guide say fastnose is refused with bootstrap(); the code accepts it with a caution and drops it | docs fixed 4a6a5b9 |
| docs-fastnose-impvalues-refusal | medium | confirmed | Help and guide say fastnose is refused with impvalues; the code cautions and clears it | docs fixed 4a6a5b9 |
| docs-histogram-caution-claim | medium | confirmed | Help permute() and two guide sites promise a histogram caution that no longer exists, and call the histogram 'analytic' | docs fixed 76f2f38, 4a6a5b9 |
| docs-version-stamp-stale | medium | in part | The build is stamped version 1.4 / 2026-08-19 in the ado, the help and the pkg, but 401 commits including the bootstrap engine landed after that date | your release step; FYI |
| fastnose-nohorizon-caution-spam | medium | in part | fastnose prints 'no horizon ... has an identified variance' once per event-study family (2+2*nsg lines) on a run that computed no variance by contract | built f5cd92b (ruled D1) |
| fastnose-sg-binary-vce | medium | in part | Fast mode keeps the full delta-method VCE on the binary subgroup-by-horizon call while the continuous twin synthesizes the estimability marker | no cheap mechanism found; FYI |
| fcol-customvce-holdout-liftable | medium | confirmed | customvce() is held out of the collapse 'pending measurement'; measured, it is exact | built cd46002 |
| fcol-lattice-cfx-key-liftable | medium | in part | lattice() and makecfxplots are refused although adding the binning variable to the cell key makes them exact | built f1a2b94 (ruled) |
| help-eb-rownames-sg-suffix | medium | confirmed | The help still says subgroup e(b)/e(V) columns take the suffix _s#; the code names them with the infix _SG<value>_ | docs fixed 4a6a5b9 |
| jtscreen-domain-nevertreated | medium | confirmed | Joint-test cluster caution tests i-nesting over never-treated rows, so the guide's pooled never-treated coding silences it whenever the treated i levels nest in the cluster | fixed 52fd319 |
| perm-esttab-analytic-under-ase | medium | confirmed | esttab reports the permutation p-values under the fast-mode default but the ANALYTIC p-values under permanalyticse and permtstat, against the help's 'p-values and intervals remain permutation-based' | built a45480b (ruled) |
| perm-save-comma | medium | confirmed | permsave(name, replace) writes a file literally named 'name, replace.dta' at rc 0; bootsave() refuses the same idiom | refused, 8443694 (ruled) |
| perm-sgeq-uncentred | medium | confirmed | Under permute() without permstrata() the subgroup-equality test and the pairwise contrast of the same two subgroups print p = 1.000 and p = 0.010 on one run; the pairwise row's off-centre flag is set internally but never reported, with or without strata | rows listed a45480b; researched; labelled as the sharp null 66d6680 (ruled) |
| permute-pairwise-se-analytic-undisclosed | medium | in part | Under permute() with permtstat/permanalyticse the degenerate analytic pairwise se is printed as se_analytic with no disclosure; bootstrap boottstat discloses the same number | disclosed under permute(), 189af00 (ruled D4) |
| pfloor-resolution-format-underflow | medium | confirmed | The p = 0 resolution line and both bootstrap footer sentences print 'p < 0.000' once completed replications exceed 4000 | fixed aae69dd |
| ptp-nis-tcrit-normal-law | medium | in part | pretrendpower() NIS pre-test uses a t(df_r) critical value inside a normal-law quadrature, so the reported '5%' pre-test size is the size of no test wooldid or Roth runs (2.4% at df 9, k=1); BS-11 ('pin the pre-test at five percent') is half-implemented | ruled, built 6d8bedd; note c7bad2c |
| sgdeg-inert-under-suppress-small-fix | medium | confirmed | suppresssgprimaryeffects disables the pairwise rank-deficiency screen although the pwcompare pass already computes the level margins it needs | built 637062e |
| shell-cfxtwin-frame | medium | confirmed | The analytic cfx plot-data twin leaves the user in the default frame when the write fails, with rc 0 | fixed e203b93 |
| shell-plotfail-vs-savers | medium | in part | After a complete replication loop a plot-file collision is still fatal and refuses replay; the sibling bootsave/permsave collision became a caution with a posted rc on 2026-09-21 | stale files now refused before the loop (54440de); a late plot-write failure is still fatal: **open, D6** |
| shell-update-stale-cfxsave | medium | confirmed | wooldid, update re-runs the analytic plot-data write on the PREVIOUS run's e() | fixed e203b93 |
| testgap-partial-label-screen | medium | confirmed | No cert test pins the partial-labelling screen shipped as R16 on 2026-09-20 | pinned aae69dd |
| u16-ledger-comment-truncation | medium | confirmed | b0a783f truncated 29 ledger comment lines at their first comma, including the bound policy and the mechanism note it introduced | fixed aae69dd |
| u16-si20-lat-byt-df-rows | medium | confirmed | si20_lat_byt: five unbounded ledger rows and a 1e-2 V bound rest on a df-factor mechanism the degrees-of-freedom fix removed | re-bounded aae69dd |
| uz-cellconstant-dose-refuse-vs-amevoid | medium | in part | A cohort-constant dose is voided-and-continued under basic but hard-refused r(322) under untreatedZero, against the help's unqualified promise | refusal stands under your R5 ruling; help and guide now say so, 3b5712a |
| val-fastnose-boot-docs | medium | confirmed | Help and guide still say fastnose is refused with bootstrap() and with impvalues; the code cautions and clears it | docs fixed 4a6a5b9 |
| val-i-fvcap-unscreened | medium | in part | i above 2,147,483,647 passes validate's assert and dies late with Stata's misleading 'noninteger values' message | refused in validate, 5fc3ce5 (ruled D3) |
| val-permanalyticse-fastnose | medium | confirmed | permanalyticse with a user-typed fastnose runs, posts e(perm_analyticse)=1 and reports blank analytic SEs | refused, 8443694 (ruled) |
| val-sg-namelength-late-crash | medium | in part | 8-digit subgroup values: sgdom cautions before estimation, then the finished run dies r(198) building a result name | ruled D2: keep the caution; no change |
| vzero-double-caution | medium | in part | An all-zero covariance prints two cautions for one event, and the amescreen copy carries mechanism text | caution pass: note on cards #38 and #398 |

### 3b. Still open

**Decisions (put to the author by multiple choice; D1-D5 ruled and settled 2026-09-23, see 2b; D6-D8 are left as they stand unless he asks).**

- **D1, fastnose repeats a caution.** Under fastnose the "no horizon ... has an identified variance"
  caution (card #52) prints once per event-study family, 2 + 2 x (number of subgroups) lines, though
  fastnose never asked for a variance and the help says so. Silencing it there needs no new text.
- **D2, 8-digit subgroup values.** sgdom cautions, then some option sets (att_itime with a pre-window,
  adversarial, makecfxplotsbysg from 7 digits) die r(198) after estimation. A dense internal index
  for the subgroup values would need no new text.
- **D3, i codes above 2,147,483,647.** A treated unit's code past Stata's factor-variable cap passes
  validate and dies r(452) with a message naming a tempvar and non-integer values. A dense index for
  the internal _i would need no new text; a named refusal would need the author's.
- **D4, residue pairwise se under permute().** With permtstat or permanalyticse and cluster() at i,
  e(pwsgs) shows the cancellation-residue analytic se with no disclosure, where the analytic path
  cautions. Reusing the analytic path's detection and caution needs no new text.
- **D5, cfxplotnose under bootstrap().** The plot cells carry point estimates only, while permute()
  keeps the cells and fills them. Extending permute()'s carve-out to bootstrap() would give them
  bootstrap inference at no cost to the observed run.
- **D6, a late plot-write failure.** Stale plot files are now refused before the replications, but a
  plot write that fails after them (a read-only file, a full disk) is still fatal and refuses replay,
  as the 2026-09-05 ruling set and R-C kept for the savers. Making it a caution needs his text.
- **D7, the analytic cfx figure on an omitted bin.** marginsplot draws an omitted lattice bin at 0,
  while e(cfx*) and the _cfx.dta leave it blank (143e457). The all-slopes caution explains the blank.
- **D8, the guide's two identical "From the help file" blocks** (chapters 4 and 11).

**Text for the caution pass (on the review page).** Cards #38 and #398 (one event, two cautions),
#368 (a false "unaffected" clause under safety(off) when the subgroup varies within a cell), #383
(its bootstrap remedy fails under boottstat), and #408 (the permutation draw-cap caution, drafted by
Claude, still awaiting approval). Each card carries the audit's note.

**For the release step.** The ado, help and .pkg are stamped version 1.4 / 2026-08-19. Fast mode keeps
a full delta-method VCE on the binary subgroup-by-horizon margins call; no cheap way to skip it was
found. On the pathological Q33 fixture (one row per cell, every dose slope dropped),
cfxplotnose without permute() posts the ame cfx cells as zeros rather than blank; left alone under
the pathological-data cap.

### 3c. What the reviewers checked and found sound

- **Shell and replay (U01).** level() is restored after every engine and every failing call; live
  and replayed output are byte-identical on the analytic runs; every stored-result name in a
  245-name inventory over nine configs is documented.
- **Sample and regression (U02).** anticipation() is bit-exact against shifting ttre by hand; the
  onset, impvalues-cluster and one-cluster screens behave as documented; the regression calls match
  v1.3 apart from ledgered changes.
- **Estimands, event study and subgroups in the worker (U03).** An 11-config dual build against v1.3
  on the reviewer's own configs shows only ledgered differences; the e(V) splice reproduces the
  event-study test statistic exactly; IM degrees of freedom are right.
- **Posting, gteffects and cleanmatrices (U04).** e(gteffects) reproduces the main and event-study
  rows to 1e-16 on messy unbalanced data; cleanmatrices runs clean before and after every engine;
  the collapse keeps the realized sample.
- **Validation (U05).** Worker and validate syntax lines are identical and the hand-off closure
  holds; every refusal list in the help matches the code; cautions do not print twice under the
  engines' re-validation.
- **Cell collapse (U06).** The cell key is exact; nine messy configurations agree with nofastcollapse
  to 1e-10 in every stored cell; the degeneracy screen stays quiet on healthy designs and falls back
  where it should.
- **Permutation engine (U07).** Exact enumeration, strata and probability rows preserve what they
  should; replications rebuilt by hand under permonsetcohorts match the saved rows bitwise;
  permtstat's comparator and se columns match their sources bitwise.
- **Bootstrap engine (U08).** Every documented refusal matches the code; the count and re-id legs
  agree to 6e-16 under weights; the percentile, BC, BCa and studentized families follow their
  stated conventions; bootsave's weight columns are right under both modes.
- **Pre-trend power, adversarial test, bands (U09).** Both adversarial nulls match the help; the
  Roth pre-test shares the joint test's statistic; bands take the level in force on both engines.
- **Display (U10).** Table rendering, rules and formats are right on every table; replay is
  byte-identical; the p-floor lines print only where the convention needs them.
- **Screens (U11).** The collinearity screen is observationally equivalent to v1.3; the dose,
  containment and all-zero covariance screens fire on the right domains.
- **Plots (U12).** Plotted points and bands in the saved datasets equal the posted matrices under
  both engines; subgroup band plots draw; figure metadata pairs correctly across 12 cfx matrices.
- **Continuous treatment (U13).** The collapse is exact for dose designs; every control type leaves
  the att unchanged and moves only the ame reference; untreatedZero's poisson form is now correct
  (ame .31 against a truth of .30, where v1.3 gave .19).
- **Subgroups (U14).** Chains, the degeneracy screen, labelling and the collapse chain read sound;
  both engines permute or resample contrasts as documented; every subgroup number-mover traces to
  a ledger row.
- **Documentation (U15).** Every option on the syntax line has a help entry and a guide mention, and
  every documented e() name is posted; nine refusals and fourteen cautions were run and behave as
  written.
- **Tests and claims (U16).** Commit claims with pins assert what the messages say; every ledger row
  cites a real commit; no bounded row is exceeded.

## 4. Round 2 (2026-09-23): the low findings and the untested surface

**What it covered.** The author chose to launch the work round 1 left: five verifiers re-checked all
102 low findings at 998154a (68 confirmed, 20 in part, 13 fixed by a later commit, 1 refuted), with
the critic's read-only claims settled by runs; five sweeps ran what no reviewer had run (Stata's
prefixes, a missing gtools, failing-call hygiene, the help examples, poisson everywhere, options under
the engines, Stata's bootstrap prefix as an oracle, and band coverage at 10 clusters). The sweeps filed
23 findings; each high or medium one went to a separate refuter. Raw output:
change_audit_2026-09-22/round2_results.json; triage: round2_bins.json.

**Found sound by the sweeps.** bootstrap _b: and jackknife _b: reproduce the point estimates to
4e-16, and the bootstrap prefix's e(V) equals wooldid's own bootstrap(standard) on the same seed
(1.4e-9, the prefix's float storage); wooldid's percentile, BC and BCa endpoints equal estat bootstrap
on identical draws to 3e-15; the help's examples 0-14 all run as written; poisson+subgroups runs with
the z labels and headers 4173588 and 154d5dc promised; every generated lattice type equals a lattice
built by hand; emptycellsoverride and ccnoabsorb agree with their controls.

**Data fix, done at once.** Three of the 66 round-1 wording proposals sat on the wrong review-page card
(the line-to-card map was wrong for 1647, 9500 and 10928): #34's moved to #33, #277's first to #276,
and #304's DELETE proposal to #301's second line. 17 cards got round-2 notes.

**Rulings (multiple choice, 2026-09-23).**

| item | ruling | state |
|---|---|---|
| latticeignoreweights bins frozen on bootstrap's count leg and bayes | **Fix it**: replications re-bin on the resampled data; bootstrap SEs on those runs move | built a977d58; REPLATW cleared bfbaa81 |
| contreatpoly(2) all-zero VCE on an uncentred dose (v1.3 the same) | **Leave it** | no change |
| band covers 0.90 at 10 clusters (measurement) | **Leave the docs** | no change |
| bootstrap esttab posts the Wald SD family | **Same as permute()**: post the displayed family | built 1f11d5b (footer text approved) |
| ccref passes the 1e-5..1e-3 fragility band silently | **Caution, as bycoh does**; text approved: "Caution: under contreatcontrols(`ccrtype'), `n' group(s) have reference-period or control-group variation in `contreat' below 1/1000 of its sample standard deviation, so their baseline dose slopes are nearly unidentified and the reported effects may be wrong rather than noisy; examine or drop those group(s)." | built 2d3621e |
| bare wooldid after Stata's bootstrap/jackknife prefix replays the analytic tables | **Show the prefix's table** (hand off to Stata's prefix display, as regress does) | built 23ea9e9; docs a3e2ce0 |
| e() after a call that fails past the stage-1 fit | **Keep the promise**: hold the caller's e() and restore it on failure | built 7e150db |
| replay refusal marker after a failed call | **Replay a restored store**; a mistyped replay line no longer blocks the next | built 4130833 |
| adversarial p under permute() for a subgroup whose pre horizons were dropped | **Skip it**, as the analytic run does | built b86b08b |
| hbar passed through string() under bootstrap() | **Full precision** | built 1a66218 |
| poisson + unconditionalse + cluster() prints pseudo-R2 '.' | **Spell vce(cluster)** | built 9e22cfb |
| straddle caution (#277) on att-only boottstat runs | **Print only with a cell-weighted estimand** | built 9d841dc |
| #236 under boottstat with any cluster() | **Print only on realized pairwise degeneracy**, as D4 | built 8884e44 |
| analytic cfx plot-data write failure posts no return code | **Post e(wdid_cfxsave_rc)** | built 0d44407 |
| e(bootdist_pwse) above 10,000 replications | **Cap it** like e(bootdist) | built fb27352 |
| #179's dead safety(off) escape | **Drop the gate and the escape**; text approved: "Error: the subgroup variable varies within i-t cells. This is unlikely to be valid with permute()." | built 402b324 |
| saveplots() takes a name | **Keep names; fix the help** | built d7f3500 |
| poisson dyex/dydx titles reversed | **Relabel under poisson**: dyex "Elasticities (dyex)", dydx "Semi-elasticities (dydx)" | built 066b6d1 |
| permtci redundant with permtstat permabsp | **Keep it** | no change |
| within-cell screen a subset of the split-unit screen | **Keep both** | no change |
| gteffects table formatting | **Apply PS-3 formatting**; no floor line | built a49499b |
| p-floor line under bootdisplay(bc) | **Omit it under bc** | built 27e277d |
| poisson + contreat() basic: the first ppmlhdfe fit fails on wooldid's own one-level term and the 3200 fallback caution prints | **Investigate** specifications that do not fail, benchmarked against today's results, with a 10-15 agent workflow. If one works, implement it; if not, change nothing and leave the caution | investigated (13 agents, scratchpad ppml/REPORT.md): no absorb() spelling passes ppmlhdfe 2.3.0; v5 built 7647081, leg and docs b6301b0 |
| lattice() run with an all-zero VCE: following the bootstrap() advice posts no e(V) | **Leave it** | no change |
| adversarial test reported as a line without a p-value (card #25 and siblings) | **Eliminate the line-only cases**: e(adversarial_line_*) is not posted when no p-value can be had; texts approved: "Caution: the adversarial pre-trend test failed: it needs esfixedbaseperiod.", "Caution: the adversarial pre-trend test failed[ for subgroup `sg']: its covariance is degenerate.", "Caution: the adversarial pre-trend test failed: under bootstrap() it needs boottstat." | built bab1f93; paths D and E 7c32fdb; ledger reasons f647cc1 |
| stale WOOLDID_* globals change numbers silently | **Overwrite them**: the outermost call clears the switches when it exits, so none outlives a call; replications keep the engine's settings | built 029b96c |
| ledger X-est-b and X-est-w bounds | **Widen both to 1e-3** (the ledger's 10x rule) | built e62a949 |
| ppmlhdfe investigation result (no absorb() spelling passes ppmlhdfe 2.3.0's preconditioner; v5 builds the explicit ringer spelling where the spec is built) | **Build v5** | built 7647081, b6301b0 |
| ccc_absorb() on the plain-basic poisson path posts an all-zero V | **Keep ccc_absorb() absorbed**, then **reverted to explicit** once the mechanism was corrected | REVERSED 2026-09-24 to explicit, once the all-zero V was found to be an ereturn-post rounding accident; not built |
| new crash: poisson + contreatpoly(2) with lattice()/byX, or ccc_absorb() of 2+ variables, stops r(3201) (the same ppmlhdfe bug) | **Fix now**, keeping the same model (absorbed where it can be) | built 0126d32 (poly explicit as ccnoabsorb; ccc 2+ split per variable) |
| failed replication loop: previous e() or the loop's own observed run | **Previous e() on every failure** (the 2026-09-05 loop exception is retired) | built 7e150db |
| bootstrap footer sentence once e(p)/e(ci_*) carry the displayed family | **Name esttab separately**; text approved: "... rather than the columns printed above - and so do coefplot and honestdid, which read e(V) without being asked. esttab reports the p-values and intervals printed above, from e(p), e(ci_l) and e(ci_u)." | built 1f11d5b |
| two more line-only adversarial paths | **Approved**: "Caution: the adversarial pre-trend test failed: too few usable bootstrap replications." and "Caution: the adversarial pre-trend test failed: too few usable permutations."; no line posted | built 7c32fdb |
| p-floor line under bootdisplay(bca) | **Omit it too** | built 6adca1b |
| poisson + untreatedZero: ppmlhdfe posts an all-zero V for the correct model | **Investigate a fix now**; build only if point estimates are unchanged and SEs match Stata's poisson | investigated: no fix within the bar (Stata ereturn post zeroes a 3.3e-16-asymmetric singular V; fixture has 45 of 63 constant cells); nothing built |
| pooled dose slope silently dropped by ppmlhdfe on a large-level dose (found by the v5 probe; predates v5) | **Detect and warn**; text approved: "Caution: ppmlhdfe dropped the pooled `contreat' slope as collinear with the fixed effects, so the ame is wrong; subtract its mean from `contreat' and rerun." | built efc0671 (keyed on the linear pooled term; a second check spares a correct drop under a unit fe()) |

**Straight fixes (no number moves, no new text):** applied as 47 commits, ea11d36..a0c9bf3, with guide
chapter 14 regenerated in afe9e37. Held: the poisson 3200 fallback skip (the author asked for an
investigation instead) and the help's 'returned' wording for the adversarial line (superseded by the
line-only ruling). Not built: the lattice e(V) re-post (ruled leave). The ccref band legs came with the
ccref caution build (test_amescreen B2c and B2d, 2d3621e).

**Left as they stand (FYI).** Fast mode's full VCE on the binary subgroup-by-horizon call (no cheap
mechanism); the NIS k <= 4 cap (documented); the contreat subgroup reads without the fcdeg screen (no
divergence reproduced); untreatedZero summary rows (v1.3's design); 18 bounded ledger rows without a
cited sha (housekeeping); the reg-path maxvar preflight false positive (left on purpose in R18);
jointmarginsdof reaching only two sites (no effect); unconditionalse's small-sample factor (it is
regress + margins' own); a gtools that is present but cannot load (optional); the debug globals (a help
sentence is in the docs group).

**The v5 probe (the author asked for extra certainty).** Six angles compared v5 with Stata's official
poisson and glm: v5 is bit-identical to the old fallback everywhere tested and matches official poisson
on ordinary data (points about 1e-10, SEs about 1e-9; the larger gaps are solver tolerance). Report:
scratchpad v5probe/REPORT.md, summarised in the change-audit doc.

**New findings, filed, not built.** A failed customvce() call can leave state that moves the next call's
V by about 1e-7 in the same session (predates round 2). A near-collinear ccc_absorb() pair now stops at
the treated-cell collinearity screen, r(322), whose message does not name ccc_absorb(). fe(pid) with
contreatpoly(2) on a large-level dose drops half the cell polynomial slopes (a caution prints). The
dropped-slope caution also prints on byTime + lattice(all 2) + ccnoabsorb on a shifted dose, where the
ame moves only 6e-4 SE and zero-VCV cautions already print (left under the pathological-data cap).
Under boottstat esttab's se and t still come from e(V). ppmlhdfe posting an all-zero V on some designs
where Stata's poisson gives finite SEs is an ereturn-post accident inside ppmlhdfe (ticket only).

**Review page synced 2026-09-24 (efc0671).** Eight cards took their new built text (#25, #179, #229, #230,
#264, #381, #382, and #409's renamed count); eight took corrected triggers where a ruling changed when a
message prints (#1, #34, #46, #100, #192, #236, #277, #415); two new cards carry the approved texts (#417 the
dropped-slope caution, #418 the ccref band caution). The page itself was republished (Version 7): its
embedded card list had stopped at #408, so #409-#416 had never shown on it, and it now also adds any card
the shared record holds that the embedded list lacks.
