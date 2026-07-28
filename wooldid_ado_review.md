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
 


### Help-file audit: wooldid.sthlp vs the current program (2026-07-28, three-auditor wave at HEAD f3b4602; full tables in scratchpad waves/sd_a1, sd_a2, sd_a3)

**Overall verdict: the help file is in strong shape.** All 65 user options are documented and every stated incompatibility is enforced by a real gate (18 probed empirically with the documented rc); the e() Stored Results section has ZERO wrong names, wrong categories, or wrong presence conditions, with the new surfaces (b_es/V_es, gteffects' 15 columns, the entire permute family incl. t-scale semantics, roth_*) matching 1:1; and every numeric claim in the recently written sections checks out to the digit (permute's 1,000/10,000/add-one/percentile arithmetic measured exactly; the gteffects weighted-mean identity to 4.4e-16). The discrepancies: one outright-wrong sentence, a handful of refuted or stale prose claims, and missing documentation for behavior added this campaign.

**WRONG (1 sentence):**
- ~L819-820, replay section: "if a wooldid run fails, e() still holds the previous successful estimation ... so a replay at that point redisplays the earlier run." Stale since d349bb1: a failed call now sets a session marker and bare `wooldid` REFUSES with r(198) ("the stored results predate that failed call") until a successful run clears it. The premise (e() unchanged) is still true; the stated consequence is the opposite of what happens.

**REFUTED / INACCURATE claims (probed):**
- ~L349-351: saveplots "pass a directory through this string" - the option is syntax type `name`; `saveplots(subdir/stem)` exits r(198) "invalid name". Only plain stems work.
- ~L901/L906: examples 8 and 9 carry a stray trailing ")" - example 8 fails verbatim with r(198) "option ) not allowed". (Also: no example #7; "quartiles" describes a 5-quantile lattice.)
- ~L70: cluster documented "[numerical]" but the syntax accepts any varlist; a string cluster variable dies with a bare r(109) type mismatch rather than an informative refusal (except under permute(probabilities), which refuses cleanly).
- ~L92: customweightmultiplier "[>=0 numerical]" - the >=0 domain is never enforced; negative multipliers pass silently.
- L68: summary-table typo `semi:elasticitiy`. Stale header stamp ("version 1.3 2024-07-31"); "reghdfe version 5.7.3" stated but the code checks presence only.

**Undocumented behavior (silent or new this campaign):**
- esrelativeto > -1 under esf silently clamped to -1; contreatelasticitytype outside dydx/eyex/eydx/dyex silently coerced to dydx.
- adversarial requires BOTH esprelength>=1 and espostlength>=1: the parse-level clear is SILENT (measured: espre(2) + espost 0 -> rc 0, no message, no pval); mid-run/per-SG cancellations DO print Cautions (the tranche-12 guards) but the help never says adversarial can be cancelled, and the promised 2-point-fit power note was never added.
- permute(): the never-treated-must-be-plain-missing refusal (extended missings r459); the permreps>=20 parse floor; the completion-floor refusal being skipped under permexact.
- espostlength(-999) (the post-period trip that wooldid's own Caution text recommends) is absent from the help, which says int>=0.
- Under default safety, cfx plots for non-att estimands are not produced (gate on att-or-safety(off)); the help documents their naming without that requirement.
- impvalues x adversarial and impvalues x customvce refusals appear in neither entry; makecfxplots' hard lattice()/cfxlattice() requirement is only implicit.
- e() gaps: e(roth_ok) (added deliberately as the user-facing completion marker) missing from Stored Results; numpre per-subgroup variants; the makecfxplotsbysg per-SG cfx matrices; the joint-test entries read as unconditional but are gated on the jointtests option.

**Index/cosmetic:**
- The Options-by-Purpose table (~L66-93) omits gteffects, the entire permute family (7 options), and altmethod_adversarial; permreps/permexact/permstrata have no standalone entries; several newer options lack abbreviation highlighting (sub, perm, permpr, permts, permtc).
- Line-64's never-treated parenthetical omits the anticipation/oostreatedcontrols adjustments and the timetrends second-period band (simplification, right in spirit); the subgroup name-length note omits the 7-digit makecfxplotsbysg tightening; the wdidmainresults/pwsgs column layouts are listed nowhere (so the poisson z-relabel is invisible - nothing contradicts, no claim exists).

**Confirmed accurate (load-bearing samples):** verbose retention; the fastnose section end to end; all impvalues/pretrendpower/unconditionalse requirement lists; the permute and gteffects numeric surfaces; the saveplots file inventory (incl. the underscore pattern); the e(V) placeholder paragraph; the event_plot stubs and b_es repost recipe (SG rows do not collide).


### F6 permutation / randomization inference - COMPLETE (2026-07-27/28; see wooldid_solved_issues.md and review_process/f6_permutation_plan.md)
Landed in eight commits (125bbd0, 86d37cb, e6b608a, 9f31b47, b544c55, e350220, 533052c, 099d93b): permute(shuffle) with permreps()/permexact and permstrata(); permute(probabilities) with the permprobs() grammar; permtstat and the inverted-t permtci; every check, caution, count, and output from the brief. Three-verifier wave clean: the inference reconstructs bitwise by hand (exact, strata, probabilities, and t-stat modes), the feature is bit-additive with the observed results untouched, and every failure path refuses with the observed results intact. perm_basic pins the e() surface in the battery; replay of a permutation run is byte-identical.


------------------------------------------------------
------ DO NOT WORK ON THINGS BELOW THIS LINE ---------
--- ITEMS TO ADDRESS IN DUE TIME ---










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
