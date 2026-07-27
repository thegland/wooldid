# Ongoing Work Game Plan

Here's what we're going to do - we are going to chug forward on the edits and issues we identified here, going through them piece by piece. We are going to use a specific process. First, let's proceed through issues one by one, covering only those I flagged below, but not below the 'Items to address in due time' line. Second, when addressing an issue, we want to use a very careful process - the worst outcome here isn't failing to fix an issue, it's breaking existing functionality for existing wooldid users. In terms of designing that process, here's what I would suggest. First, for each issue, let's spin up a set of 3 independent opus 5-class subagents to assess the problem (feel free to do this in parallel, though don't run more than a dozen and a half some subagents at once - we have only so many tokens. Also, don't let the set of all subagents run more than 8 instances of Stata at once - we have only so much compute!). Their job should be to skeptically assess whether or not the flagged issue really exists, or whether some other subtlety in how the program works actually means that the problem isn't a genuine issue. This should include code review and testing in Stata to replicate the issue with sample data (consider the simple_example file I provided if ultra simple examples will work, or construct examples as needed). If the subagents conclude the issue is a genuine issue, they should also proceed to consider the listed fix and/or propose a fix of their own. Critically, they should assess the risks associated with the fix: would it compromise / affect something else that the file works with? It is critical to confirm that the fix won't break something else - being mindful that it must not break things either in the binary or continuous treatment modes. (Note: fixes that work only w/ a binary or a continuous treatment can be implemented provided the code is gated to only trigger in those cases.) Once these subagents complete their work, you should conduct your own analysis of the issue along similar lines, without looking at the work done by the opus subagents. Once you complete your review, please then look at the opus subagent work - if you and they are all in agreement, then you should proceed. If there is disagreement but a pathway to a workable fix exists, implement that fix. If there is disagreement and no clear absolutely safe pathway to a workable fix exists, call off working on this specific item and produce a report noting that you abandoned this fix (and why) - this is a perfectly acceptable outcome, not all issues will be feasibly fixable. If you choose to proceed, carefully implement (and document) the changes made, using git appropriately and responsibly along the way. Carefully consider and test your changes, looking to see (a) that your fix worked, and (b) that it didn't compromise anything else in the program, not even when the program is forced to address rare edge case situations (much work on my end was done considering strange edge cases, degenerate datasets with few observations per gorup, etc.). Once you are confident that this is the case, both from code review first principles and based off of tests, spin up a set of 3 new opus-class subagents to check your work and confirm no adverse impacts on program functionality. Once they complete their reports, please review the 3 reports and consider their observations yourself. If they raise genuine issues that would require significant work, revert you changes addressing this issue and report on what happened. If they concur the fix is good and you agree, then move on to the next issue. This will be a difficult and time consuming process, likely, but we want to be very careful with our changes -- we very much want to approach this with a 'First, Do No Harm!' mindset. If we make no changes at all to the code because of that mindset, that would be better than making 99 good changes and 1 bad one!

One thing to note when it comes to thinking about the program's robustness under different circumstances. When it comes to probing edge cases, certain edge cases just aren't that important. For example, the edge case 'user supplies a subgroup value of 1e32' may indeed cause the program to fail, but we don't care about that failure. That only happens if the user intentionally is trying to break the program. More pernicious would be a failure where reasonable looking values break the program for some reason, perhaps as a consequence of an unusual set of options being flaggead and unusual features of the data. This can especially occur if improperly structured data is fed in - things with deeply imbalanced panels, strange numbers of observations per group, etc. Though to be clear, the 32 character name limit in stata is a pernicious thing that does normally need to be worried about, as it can easily crash things when long names show up for other reasons. 

### Process upgrades (adopted 2026-07-21, after tranches 1-2)
The B1-B23 campaign's lessons are codified in `review_process/` and `test/`; the process above now runs with these additions:
- **`review_process/agent_brief.md`** — the standing brief every subagent must read first: read-only repo, pin-to-commit sourcing (`git show <sha>:path`, never the working tree), one Stata instance per agent + a `pgrep -x stata-mp` check against the shared 6-instance cap, and a liveness protocol (findings written to disk as they go; every Stata launch gets a deadline; partial reports beat silent waiting).
- **`review_process/assessment_grid.md`** — the mandatory coverage grid (estimands x treatment modes x weights x estimator paths x inference x if-restrictions x event-study x subgroups x plots x environment x data shapes). Every assessment and verification report ends with it, untested cells marked honestly; one verifier per wave is assigned the union of untested cells. The data-shape dimension includes the continuous-treatment variants and the pathological archetypes (no proper control group, cohorts lacking pre- or post-periods, cohorts missing specific relative years, one observation per ttre-year cell, deeply imbalanced panels, poisson separation), with prebuilt datasets in `test/data/dta_path_*` (generator: `test/mkpathological.do`) and their baseline behavior recorded.
- **`review_process/harness_template.do`** — the known-good Stata test harness (dump/compare/sentinel patterns); agents copy it instead of reinventing fragile machinery.
- **`review_process/process_checklist.md`** — the per-issue lifecycle: coordinator's independent analysis locked before reading agent reports; evidence-beats-vote adjudication of design disagreements; battery signature pre-registered in writing before each run; commit-message claims as testable assertions; incremental re-verification for wave-requested one-clause follow-ups vs full waves for new scope; migration rules for solved vs abandoned items.
- **`test/`** — the permanent regression battery (see its README for usage, the baseline-update policy, and the build-out roadmap). Note: while this is the baseline, selected bespoke tests should be considered based on the issue at hand.. 







### Please work on these items:


 
# Please assess if the following items are still issues now: 



### B16 (main-block half) — the SG-crash half was RESOLVED in `4483758` (see wooldid_solved_issues.md)
The two MAIN (non-SG) `altmethod_adversarial`+`esfixedbaseperiod` blocks compute the pre-trend-fit restriction as `espre + r(min) - 1`. When the fixed reference sits below the reported window (`esrelativeto < -1` with `esprelength <= |esrelativeto|-1`, so the ref row becomes the frame minimum), this excludes the relyr -1 coefficient from the fit and tests it as if post-treatment — a possible off-by-one. Deferred per the author's directive to fix only the rare SG-crash scenario in the skiplocation item; the MAIN blocks are the common path and were left untouched. Not bit-neutral if changed (would shift `adversarial_pval_*` on affected runs) -> needs its own careful wave + a re-blessed baseline config.

**ASSESSMENT (2026-07-26, two-agent wave at 0c901c4; assess-only — nothing implemented): CONFIRMED, STILL OPEN, and materially worse than "possible off-by-one."** The restriction text is byte-unchanged since 4483758 (verified by git log -S); every refactor since only moved line numbers (MAIN binary 1623-1624, continuous 2101-2102 at HEAD).
- **Mechanism, proven with worked tables + instrumentation:** `r(min)` is taken over the coefficient frame AFTER the reference row is appended, so `espre + r(min) - 1` equals the correct "relyr < 0" cutoff (`min_t`) only when the reference falls inside the reported window. In user terms the broken regime is exactly **`esprelength() < |esrelativeto()|` as typed**. There, the reference row IS the frame minimum, the cutoff lands one short, and the relyr -1 coefficient is dropped from the pre-trend line fit and tested as if post-treatment. With `timetrends` (two reference rows) it is off by TWO: relyr -1 and -2 are both mistested and the "line" is fitted through nothing but the two zero reference rows (slope 0; the reported p comes back missing). With `esrelativeto(-2) esprelength(1)` the restriction selects a single row and the run aborts r(2001) — the crash item below is this same defect's endpoint. The SG blocks inherit the off-by-one because 4483758's empty-`skiplocation` fallback is precisely this MAIN expression, and empty skiplocation is exactly the below-window regime.
- **Hand verification from e() exports alone** (`e(b_es_att)`/`e(V_es_att)`/`e(df_r)`): on `esrelativeto(-3) esprelength(2) + adversarial altmethod_adversarial`, HEAD's `adversarial_pval_att = 1.666e-25` reproduces as F(4,79) testing relyr {-1,0,1,2} against a two-point line (slope -0.0611); the documented behavior gives the three-point line (slope -0.0913) and F(3,79) over relyr {0,1,2}, p = 1.511e-22. The direction of the error is not one-signed across configs.
- **Fix (recommended for a future implementation wave):** replace the cutoff with `t <= \`min_t'` — since relyr = t - min_t - 1, this is exactly "all pre coefficients plus the reference row(s)" in every regime, with or without timetrends, independent of which horizons are estimable. Apply at both MAIN sites and (recommended) both SG sites, which then collapses the skiplocation guard; the altmethod_min machinery becomes dead code. **Blast radius measured across 47 config pairs at rtol 0:** only `adversarial_pval_*`/`adversarial_line_*` move, and only on esf + below-window + altmethod runs; every default-regime, non-esf, and plain-adversarial config is byte-identical (all weights, both treatment modes, all estimands, uncond/customvce, subgroups, five pathological shapes). No battery or golden config reaches the affected regime (`esf_adv` is plain adversarial, reference-inside), so no re-bless — but two new pinned configs should be born with the fix (the silent-value cell and the crash-to-rc0 cell). One residual: regime-(a) panels with `min_t == 1` still produce a 1-row fit — closed by the guard in the item below, never by this fix; **the two fixes are coupled and belong in one wave** (see the coupling note below).
- **Adjacent bugs found during assessment (recommend folding into the same wave; also logged below the line):** (A) `skiplocation` is built with `glevelsof` over every 2-gap, so on gap-year data (dta_path_relyrgap class) it comes back multi-valued and the `qui replace t = 4 7 - 1` aborts the whole run r(198) — this hits PLAIN `adversarial` + esf too, no altmethod needed: a realistic crash. (B) On the same shapes, skiplocation can misidentify a DATA gap as the reference gap and plant the b=0 reference row at the wrong horizon — a silently wrong fitted line. (C) The SG blocks' non-empty-skiplocation branch cuts at `skiplocation-1`, which matches the MAIN cut only at esrelativeto(-1); at esrelativeto(-2) the SG fit drops relyr -1 where MAIN keeps it (measured sign-flip in `adversarial_line_att_s1`).

### NEW (found during skiplocation verification, 2026-07-22): espre=1 + esf + adversarial + altmethod -> r(2001) in the MAIN adversarial block
With `esprelength(1)` and a fixed reference below the window, the MAIN block's `espre + r(min) - 1` restriction selects a single pre-period row, so `qui reg b t` runs on 1 observation and aborts r(2001) ("insufficient observations"). Present identically in parent and the skiplocation fix (MAIN and SG), so pre-existing and NOT introduced by the fix; `espre==0` is unreachable (adversarial is cleared when espre==0, L785). Two skiplocation assessors/verifiers flagged it independently. Low priority, narrow trigger; a capture/guard on the MAIN adversarial `reg` (mirroring how the SG crash was handled) would close it someday.

**ASSESSMENT (2026-07-26, two-agent wave at 0c901c4; assess-only — nothing implemented): STILL AN ISSUE; the trigger set is WIDER than recorded; priority raised to MEDIUM.** Reproduced verbatim at HEAD (MAIN binary `qui reg b t` at L1628; continuous twin L2106; both SG fallbacks share the expression).
- **The crash condition is exactly `espre == 1` under esf + altmethod** (perfect separation over 60+ cells: fit-set size = espre in both reference regimes). Two entry points: (1) option-driven, `esrelativeto(-2) esprelength(1)` — note the reported event-study window is IDENTICAL to `esprelength(2) esrelativeto(-2)`, which runs fine; one spelling crashes, the other doesn't; (2) **data-driven at the DEFAULT esrelativeto(-1)**: any panel whose treated units have exactly one usable pre-period (`min_t == 1`) crashes at ANY esprelength() — no unusual option needed (with espostlength(1) the existing finaladvtest<=2 gate happens to rescue it). The crash is insensitive to everything else (verified across weights, estimators incl. poisson, estimands, jointtests, safety(off), if, plots, contreat variants, and seven data shapes); only impvalues and fastnose pre-empt it (r198). `espre == 0` unreachable: confirmed empirically. Post-crash session state is clean (B22 epilogue verified from a non-default frame; user estimates survive; healthy re-run and replay fine) — the damage is the bare, uninformative abort of an otherwise-successful estimation by a failed OPTIONAL diagnostic.
- **Guard threshold = 2 fit rows, NOT 3** (measured: the block consumes only e(b) and `predict, xb`; the chain test plugs xb in as fixed constants, so the line's sampling error never enters the p at any k; a 2-point fit is exact but functional, and the pooled esprelength(1) path already ships it today — a floor of 3 would silently kill that working path). Note a failed `reg` leaves e() STALE; any cap-based guard must not read e(b) on the failure branch.
- **Recommended guard (two layers):** PRIMARY — one line beside the existing espre==0/espl==0 gates (~L521): under esf + altmethod + `espre <= 1`, clear `adversarial` and print a C12 Caution naming the cause and the remedies (raise esprelength/espostlength or move esrelativeto closer). It fires before the adversarial frame or __marginsmodel exist (no unwind needed), covers all four twins with one edit, and is provably bit-neutral on every currently-succeeding run. BACKSTOP (optional, order-independent) — cap each of the four `reg b t` calls; on failure unwind (frame change/drop, cap est drop __marginsmodel — REQUIRED or the next horizitem's est sto fails — est restore ___wooldidfullmodel, clear adversarial+postinject) and issue the same Caution.
- **COUPLING with the B16 off-by-one above — the fixes interact and belong in ONE wave.** Under the corrected `t <= \`min_t'` cut, the esrelativeto(-2)+esprelength(1) crash becomes a working 2-row fit, but the min_t==1 route persists — the guard is needed regardless. Conversely, once the cut is fixed the `espre <= 1` closed form over-suppresses and must be re-derived to `esrelativeto == -1 & espre <= 1` (or rely on the cap backstop, which reads the actual fit set and survives either ordering).
- **Pre-registration for the eventual fix:** bit-neutral on all currently-succeeding runs; on currently-crashing configs e(adversarial_*) stays absent (the run dies today) and hidden e(wdid_adversarial) becomes empty; no baseline re-bless needed. The existing finaladvtest<=2 and espre==0 skips are SILENT — consider making all three noisy for consistency.
- **Adjacent finding (uninvestigated):** timetrends + esf + adversarial + altmethod returns `e(adversarial_pval_att) == .` in every cell that runs (rc 0 throughout) — consistent with the off-by-two zero-line fit found in the B16 assessment; fold into the same wave.







##### Please assess the following item for validity and the viability of implementing with harmless fixes, putting in extra work to check on margins documentation and examples to confirm that resuability would work and do not harm. Do not implement a fix yet. Please provide details on your assesment below. Focus a lot on feasibility and on what changes would need to be mdae to make things work.





------------------------------------------------------
------ DO NOT WORK ON THINGS BELOW THIS LINE ---------
--- ITEMS TO ADDRESS IN DUE TIME ---





### F6. Permutation / Conley–Taber inference for few *treated* clusters
WCB itself degrades when the number of *treated* clusters is very small (MacKinnon–Webb) — the classic 1-treated-state policy evaluation. A permutation option (reassign placebo treatment timing across control clusters, build the null distribution of the ATT estimate) covers this gap; Conley–Taber (2011) is the citable framework. Cost: full re-estimation per draw (reghdfe makes this tolerable; gate behind `reps()`), and design decisions about staggered-timing permutation schemes. 

Let's implement this, once again taking advantage of fast-mode style functionality (since only the coefficient). For permutation testing, let's offer a number of different approaches to doing randomization inference. 







### F3. ATT(g,t) building-block table
Under the Wooldridge approach the `i._i#i._t#1.targetflag` coefficients *are* the cohort-by-relative-period effects ATT(g,t) — with SEs already sitting in `e(b)`/`e(V)` at zero extra computation (linear case; xb scale for poisson). A `gteffects` option that extracts, labels (cohort × relyr), returns, and optionally displays this matrix would: enable arbitrary custom aggregations downstream, make output directly comparable to `csdid`/`xthdidregress` group-time tables, and give the `histogramcohorteffects` feature a rigorous sibling (the histogram currently uses `nose` margins). Mostly bookkeeping.



### **L1905/1909 df hardcode:** 
the binary att pre-trend joint tests hardcode df(df_r) while their continuous/SG siblings thread marginsdof, so under unconditionalse+jointtests the binary tests use the conditional VCE. Fix = df + conditional vce(unconditional), NOT a blind marginsdof swap (impvalues sets marginsdof to nose).


### **Duplicate SG colnames in e(b)/e(V)** 
(e.g. ES_ame_SG2_post1 x4 in a continuous+sg config, both builds): now that e(V) can exist on such runs (47defa0), _b[]/test can hit ambiguous-name errors where they previously failed for want of V. Own ticket.



### A list of minor items sits below. For these, please implement fixes for these issues. Since they are simpler, you can use just 1 subagent for tehse items, rather than 3. Feel free to add more scrutiny to these though if anything turns out to be non-innocuous or more complicated than things look. 
- **Poisson-path inference is z/chi2 by construction:** ppmlhdfe posts no e(df_r), wooldid passes df(.), and margins treats df(.) as no-df — so every margins p on the poisson path is a z/chi2, and e(pwsgs)'s hardcoded "t" column label is wrong there. Pre-existing; both joint-test engines inherit it identically. Can be solved by changing the column label in this case. 
- **Fast-mode display cosmetics:** removing the joint-p lines merges bracketing blank lines; the header still reads "Standard Errors: cluster(i)" above an all-"." table (one verifier considers a reword mandatory; the author's same-format directive won — flagging the disagreement).
- **fc0cadb commit-message omission (corrective):** the pre-leg prejoint blanking topology change was not called out; strictly corrective and unreachable in all battery configs.
- **Poisson + subgroups header lacks a space:** "Poisson RegressionResults for Subgroup ..." (the sg-loop poistext has no trailing space, unlike the main header). Fixing changes console bytes — needs a pre-registered golden-corpus re-bless; do not fold into another change silently.
- **sthlp / e() mismatches (doc):** the help documents e(joint_horizon_bysg), e(joint_horizon_i_byssg), e(joint_cont_horizon_bysg), e(joint_cont_horizon_i_bysg) as matrices — the code emits `ereturn local joint_horiz_s<sg>`-style macros instead.
- **Adversarial silently cancelled when finaladvtest trips:** the user asked for adversarial, gets no p-value and no message. A one-line Caution at the 4 cancellation sites would fix.
- **F7 partial-ereturn corner:** if wooldid_ptp failed mid-loop, earlier types' roth_* returns would persist with no marker (each individually cap'd; an e(roth_ok) flag would close it). Duplicate power levels (pretrendpower(80 80)) yield a duplicated row.
- **Internal ttre tempvar is float:** odd calendar codes >2^24 (YYYYMMDD style) round at the gen site (r(461); parent-identical). A gen double would change existing behavior for such data — wants its own careful wave.
- **contreatwithin silently overwrites the syntax local `safety`** (S1 pattern): F13's checks are shielded via the safetyinput copy, but every OTHER safety-gated check is still silently disabled under contreatwithin. Note - this is intended. Contreatwithin will always fail the safety checks originally turned off when safety(off) is specified, so by default it turns off the safety for you. No changes need be made, just flaggint this. 
- **Replay UX notes (author calls):** (1) after a run that fails at parse time, bare `wooldid` silently replays the PREVIOUS successful run (stock Stata semantics; documented in the sthlp) while a late failure gives the arity error — two different post-failure behaviors;  --- we want a failure at the parsing stage to trigger an error, not a replay (4) e(wdid_replayv) is presence-gated, not version-compared — any future change to the replay contract must bump and compare - please note this or address however would be appropriate.
- **saveplots prefix + subgroup-value filename concatenation is ambiguous in principle** (`sg`+`20` vs `sg2`+`0`) — never observed colliding. Maybe fix by always injecting a '_' character, since SG values will always be integers.


### Please run subagents to assess each of these items and provide greater detail on the issue, when it fires, and what impact it has on the program. Please also asses fixes.
- **Commensurability tier vs timetrends (one-t-unit band):** the tier's lower bound is min(t)+anticipation, but under timetrends the worker discards onsets at or below min(t)+1+anticipation (mintp1); a fractional onset in that one-unit band hard-errors r(459) where the worker would have dropped it. Conservative direction (error rather than silent continue), very narrow corner (s5_v2 confirmation); widen the bound by +1 under timetrends if it ever bites.
- **mint fragility:** the pre-drop mint is reused after row-level drops between its computation and the reference-forcing clause, a latent collinearity trigger (the real cause of the F13 wave's P6 r(322)); F13 now masks its only reachable trigger rather than fixing it. P6-at-espre(0) predicted silent, untested.



---
### S1. Treat syntax-parsed option locals as read-only; derive working copies
B7 (cfxplotnose/cfxlattice overwritten mid-loop, silently truncating later estimands' results) and B4 (`ll` reused as an IM CI bound, corrupting e(ll)) were both shared-mutable-local bugs — inevitable in a 3,500-line single program where every local is visible everywhere. A simple discipline eliminates the class: never assign to a local that `syntax` created; when a section needs a variant (an option rewritten for margins, a default filled in), derive a section-prefixed copy (the B7 fix's `cfxseopt`/`cfxlatticeplot` pattern). A one-pass audit of every `local <syntaxoption> ...` assignment left in the file would find any remaining instances. Longer-term, factoring major sections (cfx plots, IM inference, the adversarial machinery, results display) into subprograms would give each its own local namespace and make this structural rather than disciplinary — it would also enable F11's replay refactor.




### S6. A capture-wrapped worker with a cleanup epilogue (elevates B18/B22)
Multiple confirmed leak-on-error behaviors: `set emptycells drop` never restored, temp frames left as the current frame on mid-run failure, `__marginsmodel`/`___wooldidfullmodel` stranded in the estimates dir, `_ppmlhdfe_d` left in the data on aborts between estimation and cleanup. The standard structural fix is the wrapper pattern: `wooldid` becomes a thin shell that snapshots session state (emptycells, current frame), calls `capture noisily wooldid_worker ...`, then unconditionally restores state and drops temp estimates/frames before re-raising the error. One contained refactor retires the whole B18/B22 family and makes every future mid-run error safe by default.





---

## Appendix: verification environment & evidence

- Stata 19.0 MP (StataNow) via `/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e`; packages present: reghdfe, ppmlhdfe, gtools, ftools.
- Synthetic staggered panel: 45–200 units, 10–12 periods, two treated cohorts + never-treated group, aweights, subgroup variable, continuous dose.
- Feature-section checks: `boottest` and `honestdid` are installed here (`pretrends` port is not); `margins` returns `r(Jacobian)` with coefficient-named columns (basis for the F1 linear-combination route); `margins, post` returns the full V needed for F2.
- Confirmed by direct test: `syntax` returns `aweight` (not `aw`); `levelsof(x)` with parens is legal (not a bug); plain `poisson` posts no `e(df_r)` but `margins, df(.)` and `ereturn post, dof(.)` are both tolerated (so the poisson+customvce path works — z-statistics result, which matches the `z` column label); `matrix rownames A = eq:` is valid; margins `r(by#)` and `confirm matrix r(b)` work as the code assumes; long macros in `subinstr()` string literals are fine on modern Stata; `total, over()` column order matches the code's indexing; factor-variable values above 32,740 are accepted.
- Repro commands for the headline bugs (data as in `scratchpad/smoketest*.do`):
  - B1/B2: `wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats` → compare `e(summarystats)` rows.
  - B3: `wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends` → all SE/t/p/CI missing.
  - B4: `wooldid y i t ttre, att cluster(i) impvalues …` → `e(ll)` ≈ CI bound.
  - B5: `wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(untreatedZero) summarystats safety(off)` → r(198).
  - B6: `wooldid y i t ttre, att_itime cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests` → `e(joint_horiz_i_s1)` missing (present with `att_it`).
  - B7: `wooldid y i t ttre, att att_it cluster(i) contreat(dose) contreatcontrols(basic) lattice(all 3) makecfxplots safety(off) …` → `e(cfxatt_it_cfxpost)` has 1 column vs 6.
  - B12: any `poisson` run → `_ppmlhdfe_d` remains in dataset.
