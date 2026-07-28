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
 


### Help-file update work order - COMPLETE (2026-07-28; commits d98dc73, e5445c4, 1fb8f01)
Every item on the graduated list is addressed. Two were code changes per the author's rulings: cluster() now enforces numeric variables in the syntax itself (string clusters refused cleanly at parse; d98dc73), and units split across subgroups over time are now refused under default safety with share weighting + a caution under safety(off) (e5445c4; pinned by the new test/test_sgsplit.do, battery 29/29 bit-identical). The remainder is the help-file overhaul (1fb8f01): the stale replay sentence rewritten; the saveplots directory claim corrected; examples 8/9 fixed (9's pre-treatment plots now actually materialize - verified) plus new examples 7/11/12 covering anticipation+pretrendpower, fastnose+gteffects, and permute; the Description's inference enumeration, the index table (9 options added), and every entry-level gap from the audit list (contreatwithin wording per the author, poisson/customvce incompatibilities, ppmlhdfe retry, _ppmlhdfe_d reserved name, cfxplottypes precondition, espostlength override, permute floors/cautions, gteffects moved beside verbose, ttre double-storage advice, control-group requirement, Stored Results lead-in and b_es/fastnose condition, histogramestimates row-name note, Rambachan & Roth citation, anchors, version stamps, ALL-CAPS softening, typos). SMCL render verified clean; all examples run against test data.


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
