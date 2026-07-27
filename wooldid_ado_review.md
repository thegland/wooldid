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



### Adjacent findings from the 2026-07-26 B16-main assessment wave (logged, not acted on)
- **timetrends + adversarial returns a missing p in every cell** (rc 0) — RECLASSIFIED (2026-07-26, b16m_a3 run2 + confirmed unchanged by the landed cut/planting fix c0e17f8): the earlier off-by-two hypothesis was WRONG. Root cause: under `timetrends`, the final `qui test `testconstructor'` DROPS EVERY CONSTRAINT (r(drop)=1, r(df)=0), so F and p are missing. This happens in the POOLED regime too and with PLAIN adversarial (no altmethod, no esf, no reference planting involved) — `adversarial` is a silent no-op under `timetrends` in every regime measured (n09/n10/n11/n14/g03 probes). The cut fix did repair the fitted LINE under timetrends (n11 line (0,0) -> (-0.151, 0.824)); only the p is dead. Needs its own item: likely the projected-line constraints are collinear with the cohort trend terms the test conditions on — diagnose whether the test should target trend-adjusted contrasts or simply refuse with a Caution under timetrends.
- **2-point altmethod fits have no power against pre-period non-linearity** (zero residuals by construction) — intentional to keep allowing them; worth one sentence in the help when the fix wave touches the sthlp.

### Findings from the 2026-07-26 minor-items wave (logged, not acted on)
- **SG joint-test returns are macros while their pooled siblings are scalars** (`ereturn local joint_horiz_s<sg>` at ~L3491-3494 vs `ereturn scalar` for the pooled joint tests) — looks accidental; converting to scalars is trivial but changes the e() surface (battery dumps `e(macros)`), so it wants its own small pre-registered wave. The sthlp was patched to document current behavior (macros) in the meantime.
- **e(converged) is documented but never returned:** the worker harvests `converged` into a local (~L1134) but no `ereturn` ships it, on any path. One-word fix (add to the pass-through ereturn loop at ~L3458-3463) if wanted; the sthlp entry was removed for now — revert that doc hunk if the code fix is preferred.
- **The `o.` omission marker defeats the SG truerow match** (proximate trigger of the fixed label-shift bug): the strip handles `bn.<sgvar>`/`.<sgvar>` but not the omitted-marker forms, so an all-omitted per-sg continuous pre-margins round is treated as "no pre row". Fixing it would ADD a (currently missing) all-omitted pretreatment:ame row to the SG tables — a dimension/value change, own ticket.
- **poisson + subgroups + event-study windows is pathologically slow** (>12 min on the battery panel; pre-existing; fastnose is the workaround). Worth profiling someday.



### A list of minor items sits below. For these, please implement fixes for these issues. Since they are simpler, you can use just 1 subagent for tehse items, rather than 3. Feel free to add more scrutiny to these though if anything turns out to be non-innocuous or more complicated than things look. 
- **Internal ttre tempvar is float:** odd calendar codes >2^24 (YYYYMMDD style) round at the gen site (r(461); parent-identical). A gen double would change existing behavior for such data — wants its own careful wave. *(ASSESSED 2026-07-26 (m8_ttrefloat at 3f323a3): DEFERRED per the provably-inert mandate — and the item's premise is WRONG in a way that RAISES its severity to the major class (this is the ttre twin of the `_i` pooling bug fixed in 8c5bfd7). Measured: the r(461) fires only on sparse t grids; on dense YYYYMMDD-style grids the current build runs rc 0 and SILENTLY MERGES cohorts whose onsets float-round together (att -11% biased with a fabricated pre-trend spike; -64% under esf), and half-integer codes below 2^24 shift too; anticipation() can flip a wrong panel to accidentally-correct. A 4-site gen double prototype (L200 primary, L248 tp1, L670/L799 modified_ttre) is ready in waves/m8_ttrefloat: 17/17 float-exact configs bit-identical two-build, and large-code runs become byte-identical to their small-code recodes — but shipping MOVES rc-0 results (wrong -> correct), so it needs author graduation, plus a decision on the timetrends companion defect (raw t ~2e7 makes the i.i##c.t trend columns numerically inert on BOTH builds — not curable by storage type; center t or refuse loudly), plus optionally a validate-side loudness floor (float(ttre)!=ttre => Caution/error) as the do-no-harm minimum.)*


### STILL OPEN, awaiting graduation (from the mint assessment, rated above the item it was found under)
- **The esfixedbaseperiod-branch estimability drop (L429-431 class) can delete the very reference period a cohort was retained for, on CLEAN data with DEFAULT options** (realistic late-entering control group): safety on -> misdiagnosing r(322) ("...Specify contreatwithin..." on data with no continuous treatment); safety(off) -> silent att -1.01 where truth is +3.9. The pooled branch is immune (its controlcheck keys on _t==., which provably cannot remove a reference row). Suggested fix: re-verify refperiodobserved (or fold reference rows into controlcheck) after the L431 drop and remove affected cohorts with a message. Prototypes/evidence in waves/mint_a1; the landed mint fix verifier confirmed this branch is byte-identical pre/post 0c03010, so the finding stands exactly as assessed.







------------------------------------------------------
------ DO NOT WORK ON THINGS BELOW THIS LINE ---------
--- ITEMS TO ADDRESS IN DUE TIME ---






### F6. Permutation / Conley–Taber inference for few *treated* clusters
WCB itself degrades when the number of *treated* clusters is very small (MacKinnon–Webb) — the classic 1-treated-state policy evaluation. A permutation option (reassign placebo treatment timing across control clusters, build the null distribution of the ATT estimate) covers this gap; Conley–Taber (2011) is the citable framework. Cost: full re-estimation per draw (reghdfe makes this tolerable; gate behind `reps()`), and design decisions about staggered-timing permutation schemes. 

Let's implement this, once again taking advantage of fast-mode style functionality (since only the coefficient). For permutation testing, let's offer a number of different approaches to doing randomization inference. 



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
