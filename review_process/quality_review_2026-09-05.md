# Quality review of the 2026-08-30 to 09-05 code: bootstrap(), the cell collapse, efficiency work, permute()

Date: 2026-09-05/06. Branch `bugfixes`. Reviewed tree: `310fa61` (the tree the eight review lenses were pinned to); fixes landed on top of `5aabe1e` (the SI-20 close-out) and are listed by commit below.

## 0. How this was done, and how to read it

Eight independent review agents each read the pinned tree through one lens - bootstrap code (Fable), fast-collapse code (Fable), efficiency work, econometrics of the bootstrap, econometrics of the permutation engine, applied-micro convention, missing features, and Stata idiom/robustness - each with its own Stata probes. Their 100-odd raw findings were deduplicated into 67 (F001-F067, kinds BUG / JUDGMENT / DEVIATION / MISSING_FEATURE), and every one of the 67 was then handed to two adversarial verifiers (one told to reproduce it, one told to refute it). The verdict column in the tables below records that pass; where it says "pending" the verifier had not reported when this file was written and the entry will be updated.

The rule you set: straight bugs are fixed and moved on from; judgment calls, convention deviations and missing features come to you. Section 2 is what was fixed. Sections 3-6 are yours. Section 7 is the verification record, and section 8 says what still stands between this tree and the full assessment run.

Two items sit at the top because they change reported numbers or a default behaviour, not just a message:

- **F002 (fixed, 482bf53).** Under the physical re-id bootstrap leg (cell-weighted estimands with a nesting cluster, and boottstat with cluster()) the gteffects cells were harvested from each replication by the replication's *fresh* unit codes. With original codes 1..N the wrong cells were silently matched (relabelling the same data by +100 changed every bootstrapped cell column); with any other coding no cell matched and the columns were blank. Every gteffects bootstrap column produced by that leg before this fix is wrong or empty.
- **F001 (bayes /W_c split) and F003 (BC-t direction)** are the two remaining bugs that change reported numbers, and both touch a design you ruled on. They are held for your word (section 3) rather than fixed unilaterally, with the evidence and the recommended action.

## 1. Headline findings by severity

| id | kind | sev | title | status |
|---|---|---|---|---|
| F001 | BUG | high | bootstrap(bayes) with the default bootalpha(equal) bootstraps the equal-cluster-weighted estimand, not the reported one | HELD for ruling (s3) |
| F002 | BUG | high | gteffects bootstrap under the re-id leg keys cells on unit codes | FIXED 482bf53 |
| F003 | BUG | med | BC-t (boottstat + bootdisplay(bc)) shifts the studentized interval the wrong way | HELD for ruling (s3) |
| F030 | JUDGMENT | high | Default bootdisplay(percentile) with no reps floor and no low-reps caution prints confident wrong-looking inference | s4 |
| F031 | JUDGMENT | high | wooldid_worker has ~3% headroom under Stata's per-program ceiling | s4 |
| F049 | DEVIATION | high | No exchangeability caveat: with atypical treated clusters the permutation test is badly oversized, and permtstat does not rescue it when cohort = cluster | s5 |

## 2. Straight bugs fixed (one behaviour change per commit; text items bundled)

Every fix carries a certification pin that fails on the pristine tree. Battery: 34/34 configs bit-identical to the `5aabe1e` baseline at rtol 0 after wave 1 (pre-registered: none of the wave-1 changes touches a battery config's numbers; the battery's permutation configs are sampled, not exact). Suites re-run green after every wave (section 7).

| finding | commit | defect | fix | pin |
|---|---|---|---|---|
| F005 | 6e82d56 | Footer and help said the pre-period permutation p-values *fall* as the true effect grows. Measured on a 12-unit staggered fixture at effect 30: null SD 0.62 -> 4.83, p .49 -> .69, joint p .88 -> .98; no movement with a single onset date | Caveat rewritten: in staggered designs the reassigned onsets move real post-treatment outcomes into the draws' pre-treatment horizons, the null distribution widens and the p-values rise toward 1; a large pre-period p is not evidence of no pre-trend | permute3 Q13 |
| F013 | 9d1047e | bootsave() name clash or bad name refused only after every replication had run | bootctl checks the plain-name / cwd / replace rules before the loop (permsave() already did) | gates G17 |
| F059 (part) | d81fb85 | rseed() refused as a stray bootstrap companion under permute(); no e() record of the seed | Validated before the stray screen, applied in permrun before the first assignment draw; e(perm_seed) | permute3 Q14 |
| F006 | f8671ed | cluster() finer than i under permute(): a unit received two onsets per replication (every replication failed, or a wrong design under permexact) | permctl requires cluster() to nest i (the bootstrap re-id rule); r(459) with the reason; cluster == i and coarser unaffected | permute3 Q15 |
| F011 | 4484493 | Under a complete exact enumeration the fixed metric appended the observed comparator to draws that already contained it (counted twice; exact-mode joint/equality/adversarial p's moved) | Comparator appended only when the enumeration is incomplete (permqf + permadv) | permute3 Q16 (draws-only reconstruction) |
| F008, F026, F047 | 13e873a | 32-character overflows: bootsave() pwse<j>_/se_ variable names (r(3300) after the run); hidden perm/boot subgroup twins for 9+/13+ digit subgroup values (abort after the run); permcellsx se twins silently dropped for 22-23 character cfx names | Truncation (collision-free via the leading index), skip-the-twin guards, gate on the longest prefix | gates G18, permute3 Q17 |
| F014 | 579a608 | makeplots + saveplots at level() != 95: the saved _ES.dta carried lb_95ci/ub_95ci names for 90% bounds | Alias undone before both save sites | permute3 Q8 (it had pinned the wrong name) |
| F028 | 97a9491 | Each replication stored its regression as ___wooldidfullmodel: a verbose bootstrap run's retained model was the LAST replication's reweighted fit (or deleted) | Replications store ___wooldidrepmodel and drop it at exit; shell cleanup and wooldid_advstat follow the name | gates G19 |
| F009, F010 | f7c1fd2 | Post-loop failures: the replay-staleness marker was cleared as soon as the loop returned, so a failure during the e() rewiring left a half-rewired e() that a bare `wooldid` replayed; plot / pretrendpower failures under the engines were swallowed (rc 0, stale file, frame possibly left) | Marker armed ("postloop") until the display completes and named in the replay refusal; plot rc propagates after the frame restore | gates G20 |
| F016 | 5ca3f16, 57dfbb3 | dots(#) printed a mark every replication whatever # was | One mark every # replications at all three loop sites, count after 50 marks | gates G21, permute3 Q18 |
| (idiom 6) | 290b988 | boottstat pairwise-se caution printed on unclustered runs (gated on nclus, the observation count there) | Gated on cluster() | gates G22 |
| F015, F029, F024, F007 (help), F017, F022 + MC harness | 89b59ab | gteffects note under bootstrap credited the stage-1 VCE; synopt omitted bootdisplay()/bootsave(); bootdf paragraph claimed joint tests widen to F(k,df); example 13 promised joint tests without boottstat; "close to calibrated"; help said the Wald family is not formed under boottstat (all three coefficient-draw tables ARE stored - verified on a boottstat run); stale comments (bayes default, boottp add-one); MC harness null-band coverage band Binomial(400,.90) though the band follows level() = 95; mc_bootse arms read the percentile table for a normal-CI comparison | Corrected; mc_bootse reads e(boot_wald_mainresults) | text |
| F002 | 482bf53 | gteffects cells under the re-id leg keyed on the replication's fresh unit codes (see s0) | Harvest translates each replication row's unit code back to the original i through the replication frame's (__wdbnewi, i) columns; first copy of a twice-drawn cluster kept; undrawn units missing for that replication | gates G23 (relabelling invariance, missing-aware; pristine 40 populated cells vs 0) |
| F004 | 97742e2 | Every bias-corrected endpoint took v[ceil(n*p)] where the percentile family and Stata's bstat average the adjacent order statistics at integer n*p; with z0 = 0 (about 5% of runs at reps 200) the BC interval differed from bstat's e(ci_bc) and from the percentile interval it must collapse to | One Mata helper (wooldid_bootqs) applies bstat's rule at all six sites | gates G24 (against bstat itself) |
| F023, F021 | 97742e2 | Footer said "cluster bootstrap" on observation-level runs; contrast line said "coefficient differences" under boottstat | Conditioned on e(clustvar) / e(boot_tstat) | gates G25 |
| F027 | 97742e2 | pretrendpower() under bootstrap cautioned "bootstrap event-study covariance unavailable" for families never estimated | Existence checked first | gates G26 |
| F012 (help), F019, F025 | 97742e2 | Help: subgroup-screening sentence claimed every subgroup margins call is screened (the contreat() reads are not); permutation stored results missing e(perm_*Wobs), e(perm_adv), e(perm_comparator_se), the hidden cell twins, e(perm_offcentre); permsave() implied the plot cells are written and that row 0 is the comparator under permanalyticse | Corrected | text |
| F020 | d81fb85 | Stray-companion error named only bootstrap() for level() | Message notes level(), dots(), rseed() are accepted with permute() | - |

Not fixed, by design of this review:

- **F018 (BUG low)** - the cell-surface BC z0 counts draws below the full-run observed cell value while the coefficient surface counts against the rep-0 fastnose comparator. Numerically nil for continuous statistics; a fix needs a rep-0 cell comparator plumbed through the bootstrap path (the permutation engine has one, `cell0`; the bootstrap engine does not). Listed as a consistency note for you.
- **F012 (the screen itself)** - the continuous-treatment subgroup reads run on the collapsed rows without the wooldid_fcdeg screen the binary reads get. The help now says so; adding the screen needs a flag on the Mata predicate (j0 = 1 for dydx(contreat)). PLAUSIBLE, not reproduced (needs an ame-singular / att-healthy direction). Your call whether to add it.

## 3. Bugs held for your ruling (they change reported numbers under a design you ruled on)

### F001 - bootstrap(bayes) + bootalpha(equal): the /W_c split centres the replications on the equal-cluster-weighted estimand

What the code does (wooldid_bootdw, bootrun header): each cluster c draws Dirichlet mass p_c (Gamma(1) under equal, Gamma(n_c) under size), normalised across clusters, and every row of the cluster is multiplied by p_c / W_c, W_c the cluster's base-weight total. So the cluster's total replicate weight is p_c whatever its size or the user's weights.

Consequence (econ_boot t1_bayes.do, 12 states of 1,1,2,2,3,4,6,8,12,16,24,32 units, tau_s = 0.1 n_s, cluster(state), att, reps(200)): observation-weighted truth 1.604, equal-cluster truth 0.800, observed att 1.6915. bootstrap(standard): replication mean 1.514, percentile CI [.18, 2.38]. bootstrap(bayes) default: replication mean 0.7997 (= the equal-cluster truth), se .345, z0 = 1.96, and the default percentile CI [.199, 1.678] EXCLUDES the point estimate. bayes size: mean ~1.6, se .044 (the documented understatement). A balanced twin (9 units per state) shows no shift - which is why the 08-30 MC (balanced clusters) could not see it. A hand loop with multiplier p_c and no /W_c recovers the centre (1.5108 vs standard's 1.514).

Why it is held: the /W_c split is the "locked 2026-08-30" semantics, it is what makes the bayes engine match rwgen priorpowers(-n_c) bitwise (cert test_bootstrap_bayes O2), and bootalpha(equal) is your ruled default on the MC evidence - evidence measured on balanced clusters. Rubin's Bayesian bootstrap with clusters as the exchangeable units multiplies every row of cluster c by p_c (cluster total n_c p_c), which is what the plan (s3) wrote.

Recommendation: drop the /W_c division under bootalpha(equal) (row multiplier = user weight x p_c, rescaled so the multipliers average one), keep the p_c/n_c split only under bootalpha(size) (bayesboot compatibility), re-pin O2 (its oracle becomes rwgen's weights x n_c), and re-run the 08-30 MC with unequal cluster sizes before trusting the bayes calibration numbers in the help. Until then the help should say bootalpha(equal) targets the equal-cluster-weighted estimand when cluster masses differ. Verifier verdict: pending at the time of writing.

### F003 - BC-t (boottstat + bootdisplay(bc)) shifts the studentized interval against the median bias

The studentized kernel applies Efron's z0 shift to the quantiles of the centred t-draws and then reflects them (the percentile-t construction), which moves the interval AGAINST the bias and undoes the correction percentile-t already makes; no BC-t construction exists in Efron-Tibshirani or Davison-Hinkley, and the studentized bootstrap is already second-order accurate. Measured: with t_r = z - 1 (2000 draws), b = 1, se = .5, percentile-t gives [.522, 2.478] (shifted up, correct), BC-t [-.145, 1.480] (shifted down), and the BC-t upper endpoint coincides with the naive percentile-of-b_r upper endpoint - the interval BC was invented to repair. On dta_tiny (boottstat reps 199) every BC-t interval lies below its percentile-t twin.

Recommendation: refuse bootdisplay(bc) with boottstat exactly as bootdisplay(wald) is refused there, drop the e(boot_tbc_*) twins and the "bias-corrected ranks" sentences in the help, re-pin test_bootstrap_tstat. This removes a family you designed on 2026-09-04, hence held. Verifier verdict: pending.

## 4. Judgment calls

Grouped by theme. "Recommend" is my view; the numbers are the reviewers' measurements.

### 4a. Display defaults and p-value conventions (bootstrap)

- **F030 (high) - the default percentile display with no reps floor and no low-reps caution.** At reps(50) the same 50 draws print pretreat p .04 (percentile), .091 (wald), .018 (bc), .16 (boottstat) against analytic .122; at reps(2) and reps(5) every p prints 0 and every CI excludes zero with no caution. The help's "close to calibrated" rested on the Wald interval only (now quoted precisely, 89b59ab); percentile/BC coverage is unmeasured (F058). Recommend: a caution under percentile/bc display when completed reps < 200 naming the order statistics the endpoints rest on; or make wald the default (the family whose coverage was measured; Stata's prefix, did_multiplegt and csdid default to a normal/symmetric package).
- **F033 - two p conventions on one surface**: coefficient/contrast p without add-one, joint/adversarial/equality p with add-one, inside one boottstat table. The bootjwr comment's justification (W_r is "simulated under the null" unlike the coefficient draws) does not hold - both are recentred null draws. Recommend one convention per engine, documented once.
- **F038 - bootstrap p prints as exactly 0** at resolution 2/n (every reps(200) run in the econ_boot probes). Recommend printing "< 2/n" (or adopting the (R+1) convention, which never prints 0 and keeps duality).
- **F053 (deviation) - under percentile/bc display the table prints t = b/se_boot (a Wald object) beside a rank-based p**: "t -1.85, p .04" in one row. Stata's estat bootstrap prints no t/p for the percentile family. Recommend blanking or relabelling the statistic column (t_wald) under percentile/bc.
- **F040 - bootdf(clusters) without cluster() resolves to t(N-1) silently** (e(boot_df) = 47 on a 48-row fixture). Recommend refusing or messaging.

### 4b. Joint tests under the bootstrap

- **F034 - joint and across-subgroup equality tests are blanked without boottstat although e(V) is the full listwise bootstrap covariance and the Wald test is one `test` away** (T5: e(boot_jt_absent) = 1 while a hand `test` on the same run's e(V) returns F(3,79) = .40, p = .75; under boottstat the W_r test gives p .63 vs Wald-on-V_boot .70). Stata's prefix and did_multiplegt's placebo test are Wald tests on the bootstrap covariance, so the blanking is a wooldid peculiarity that pushes users to boottstat (several times the cost). Recommend posting the Wald-on-V_boot tests when boottstat is absent, labelled as the Wald family, keeping W_r under boottstat as the refinement. Note F045: wooldid_bootjw (the retired draw-covariance Wald) is dead Mata code that this ruling would revive or delete.

### 4c. Permutation inference design

- **F035 - the shift interval's width grows with the true effect; the "Conley-Taber" label is inaccurate.** 12-cluster staggered fixture, analytic se .377 at every tau: coefficient-mode null SD .52 / 1.45 / 5.33 and 95% shift-interval width 1.87 / 4.69 / 18.3 at tau 0 / 5 / 20; a pairwise contrast on the same estimate goes from [-2.3, .2] to [-14.6, 9.3]. Under permtstat the t-draw spread is flat in tau for post rows but not pre rows. Conley-Taber take their reference from placebo assignments among CONTROL units only, which carry no tau. Recommend (a) now: say in help/footer that the interval is calibrated only near the null and drop the Conley-Taber attribution; (b) the actual Conley-Taber control-only reference as the interval's null distribution when #never-treated >= #treated; (c) the sharp-null grid inversion (the Fisher inversion costed on 2026-09-05).
- **F052 (deviation) - the default p convention (equal-tailed, add-one) differs from ritest / Stata permute / Young 2019 (|s|, c/n, no add-one); no wooldid setting reproduces their number** (permabsp is |s| but still add-one). Same draws: .0392 equal-tailed, .0196 permabsp, 0 under c/n. Floors 2/(n+1), 1/(n+1), 0. Recommend keeping the default (duality with the shift interval; Davison-Hinkley / Phipson-Smyth) but documenting the alternatives side by side, naming the commands, and optionally adding a plain c/n cross-check.
- **F036 - p floors undocumented; at permreps 20-39 a 5% rejection is impossible under the default** (permreps(20): every row prints p .0952). Recommend raising the floor to 40 or cautioning when 2/(permreps+1) > .05, plus a reps-for-resolution line in the help.
- **F051 (deviation) - RI-beta default with no guidance on RI-t.** MacKinnon-Webb 2020 recommend the studentized statistic when treated clusters may be atypical; the help describes permtstat only as "randomize on t-ratios". Recommend keeping the default for speed but stating what each statistic assumes (see F049).
- **F042 - the fast-mode default leaves se_analytic and t_perm blank in the main table.** A referee expects an SE column; permanalyticse costs one full observed run, p unchanged. Recommend either flipping the default or a one-line note under the table.
- **F041 - legacy e(perm_lb_est)/e(perm_ub_est) now hold the null band**, not an interval, under CI-looking names. Recommend new names with the old ones as documented aliases.
- **F043 - permute() with contreat() reassigns onset only; the dose is never re-randomized.** A permute(dose)-style mode would be a new assignment generator. Your call whether wanted.

### 4d. Bootstrap design

- **F039 - pretrendpower() under bootstrap re-estimates the Roth slope with the bootstrap covariance** (the GLS slope m'S^-1 b / m'S^-1 m moves between an analytic and a bootstrap run of the same data, against "point estimates are unchanged"). Recommend the fixed-weight variant (analytic weights, bootstrap se) proposed in coverage entry 4.
- **F044 - wooldid_ptp takes no level()**: the Roth slope CI prints 95% under level(90).
- **F057 (deviation) - observation-level bootstrap on panel data when cluster() is omitted, no caution, header reads as an upgrade.** Recommend a caution when i repeats over t and cluster() is absent.
- **F050 (deviation) - pairs cluster bootstrap with few TREATED clusters; no wild cluster bootstrap; no pointer to boottest.** With 2 treated among 12 clusters, 11.5% of replications drew no treated cluster and were excluded (matching (10/12)^12), the retained distribution is conditioned on that event, and the only caution keys on total clusters (< 10). Recommend a caution keyed on the number of treated resampling clusters naming permute() and boottest as the remedies, and recording the F1 infeasibility verdict (boottest needs a refit) in the help.
- **F048 - the bootsave()/permsave() note is display-only**; nothing in e() records the file.
- **F046 - bootcellsx indexes cfx accumulators by position in e(wdid_cfxmats) while bootrun numbers them by existence count** (agrees only because the worker posts every listed name; the permute engine passes names explicitly). Recommend passing names() as the permute block does.

### 4e. Fast collapse and efficiency

- **F032 - dyex excluded from the poisson analytic-Jacobian gate on a rationale the file itself contradicts** (the comment says both paths "take ln of the prediction"; dyex logs the regressor). Measured cost 5.9-31.2x on three poisson-contreat-dyex configs, with max relative movement 3.1e-7 when the gate is extended (largest on a subgroup t statistic); the knock-on cell-collapse ban for poisson + dyex rests on the retained numerical form. Recommend extending the gate under the commit's own J*b == margin acceptance test and re-examining the ban; which side of the 3.1e-7 is more accurate is not adjudicated.
- **F037 - the verbose collapse-decision line misattributes or omits the reason under fullprecisionslow.**
- **F031 (high) - wooldid_worker has ~3.6-4.0 kB (~3%) of headroom under Stata's per-program ceiling** (100 injected one-token statements load; 113 fail r(1000)). Treat the worker as frozen: route new features through subprograms; consider moving a self-contained block out; record the injection probe as the measurement method in process_checklist.md.
- **F056 (deviation, low) - the F8-skip commit message's execution count is off by one for bootstrap** (measured reps + 2 pre-skip for every bootstrap flavour; the gate itself is proven correct and bit-neutral).
- **F055 (deviation, low) - the fweight/physical studentized legs' comment claims BITWISE-identical draws; measured 4.3e-16** (the load-bearing claim, identical analytic se, holds). Soften the comment.

## 5. Deviations from common practice (the ones not already under s4)

- **F049 (high) - no exchangeability caveat; with atypical treated clusters the permutation test is badly oversized and permtstat does not rescue it when cohort = cluster.** 30 null datasets, 12 clusters, 3 treated with error sd 3 vs 1: p <= .05 in 6/30 (RI-beta) and 7/30 (RI-t); P(>= 6 of 30 at nominal .05) ~ .003; the homogeneous twin is calibrated. Mechanism: label permutation mixes clusters of different variance (MacKinnon-Webb 2020); with cohort = cluster and one observation per cell the treated cohort-period dummies saturate those cells, so the CRVE cannot see the treated clusters' variance and the studentization is incomplete (the same mechanism as the standing pairwise-se degeneracy under cluster-at-i). Recommend: (1) the help states plainly that permute() presumes exchangeable onset assignment across clusters (given permstrata), is not robust to treated clusters atypical in size/variance, and that permtstat helps only when the cluster-robust se is informative about the treated clusters; (2) a caution when cluster() == i and permtstat is requested; (3) longer term, RI-t with a se that sees the treated clusters.
- **F054 - level() is refused on analytic runs; the analytic tables are hard-wired to 95%.** Every Stata estimation command honours level(); the engines already relabel lb_<lvl>ci. About ten lines to thread level() into the Mata table builder and the display header, post e(level), and drop it from the stray list.

## 6. Missing features (for the standing-issues list or a ruling)

- **F058 (medium) - no calibration simulation** for percentile/BC intervals, the BC p, boottstat, the bootstrapped or fixed-metric joint tests, the permute shift intervals, permabsp, or per-horizon/subgroup permutation p. mc_bootse records the Wald-CI coverage for att only (its own comment says percentile/BC could be recorded at no extra cost); mc_ri reads e(perm_p_att) only. This is the itemised ticket behind F030 and F035, and behind the bayes MC re-run F001 needs.
- **F059 (medium) - reproducibility record.** rseed() under permute() is now accepted (d81fb85). Still open: neither engine posts c(rngstate) at the start of its draws (Stata's prefix posts e(rngstate)), e(boot_seed) exists only when rseed() was given, and the seed is not written into the bootsave/permsave file. About five lines per engine.
- **F060 (medium) - the collapse refuses customvce() "pending measurement"; measured exact** for robust and cluster(var) under regress (zero groups at rtol 1e-12) and cluster(var) under poisson (zero at 1e-10, residual inside the clustered-poisson noise floor). Recommend lifting the refusal for robust / cluster(varname) strings and converting test_fastcollapse F4 to an equivalence leg. Other vce strings (hc2/hc3, bootstrap, jackknife) unmeasured.
- **F062 (low) - the collapse refuses impvalues; adding the cluster variable to the cell key makes it exact even for a crossed cluster** (2,878 cells vs 60 on the admin panel, zero at 1e-12; without the key change the crossed case is BROKEN: b .649 -> .628).
- **F061 (low) - realised cell count not returned in e()** (e(wdid_fcol) = 1 even when the collapse is a no-op). One hidden scalar.
- **F063 (low) - bootsave() omits the resample itself** (cluster multiplicities / Dirichlet weights) and the per-replication joint, adversarial and cell statistics.
- **F064 (low) - no one-sided p in either engine; no MC SE/CI on p** (the one-sided tail shares are computed but not stored; ritest and Stata permute print SE(p) and a binomial CI).
- **F065 (low) - standard RI companions absent**: max-|studentized| step-down multiplicity across horizons (exactly valid from the same accumulator), the Conley-Taber control-only reference plus grid inversion, a not-yet-treated placebo loop for pre-trend tests.
- **F066 (low) - no runtime projection** for the replication loops.
- **F067 (low) - no BCa interval** (needs a delete-one-cluster jackknife: G extra fastnose fits).

## 7. Verification record

- Static: check_validate_sync OK; coverage_check COVERAGE_OK.
- Battery (test/runbattery.do, 34 configs): 34/34 IDENTICAL at rtol 0 against the `5aabe1e` baseline after wave 1, as pre-registered, and again 34/34 IDENTICAL on the final tree (97742e2). No battery config exercises the re-id leg with gteffects, and no baseline BC interval sits at z0 = 0, so the battery cannot see F002 or F004; s8 asks for a config that does.
- Certification suites run green on the final tree: bootstrap_gates (G17-G26 new), bootstrap_std, bootstrap_df, bootstrap_tstat, bootstrap_joint, bootstrap_bayes, permute (P), permute2, permute3 (Q13-Q18 new; Q8 re-pinned), ptpreplay, fastnose, gteffects, cfxlat, sgdefer, advstat, advguard, esfsupport, ptpdisp, refmargin. Every new pin was shown to fail on the pristine tree by construction or by probe (G23: 40 populated cells vs 0; G24: BC upper endpoint 39 vs bstat 39.5; Q13d: the old footer wording; Q16: the comparator-appended metric).
- Verifier pass over the 67 findings: pending completion at the time of writing (see the verdict column, to be filled).
- A lesson recorded for the process: editing a suite's do-file while Stata is running it corrupts the run (Stata reads do-files incrementally); the first permute3 re-run was lost that way and repeated.

## 8. What still stands between this tree and the full assessment run

1. **Rulings on F001 and F003** (section 3): both change reported numbers; the bayes one also needs the 08-30 MC re-run with unequal cluster sizes before the help's calibration claim for bayes can stand.
2. **A battery config that exercises the re-id leg with gteffects** (att_it + cluster(i) + gteffects + bootstrap), so the F002 class is visible to the battery in future; the current 34 configs are bit-identical before and after the fix.
3. **The calibration simulations (F058)** for the families the assessment will print: percentile/BC coverage and the shift-interval coverage at a planted effect are the two the report's own numbers now put in doubt (F030, F035).
4. **The display-default rulings (F030, F053, F034)**: whichever way they go, they change what the assessment battery prints, so they should land before the baselines are blessed.
5. **Help text for F049/F051/F052** (exchangeability, RI-t guidance, p-convention paragraph): text only, but the assessment reports will cite the help.
6. The standing items from earlier sessions that this review did not touch: the author's battery/oldnew full run, the analytic pairwise-se degeneracy ruling, and the SI-20 judgment calls in review_process/si20_ship_2026-09-05.md.
