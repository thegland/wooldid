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









### F9. Fast exploratory mode — RESOLVED (commit fc0cadb; 2026-07-26)
Implemented as FASTnose (abbreviation `fast`), per the author's spec: point estimates only — no SEs, CIs, p-values, or tests of any kind; the tables keep their exact format with "." in the SE-family cells; the output carries "SEs and related materials are not estimated in fast mode." once, replay-identically (hidden e(wdid_fastnose) macro). Five r(198) gates (unconditionalse, impvalues, adversarial, jointtests, pretrendpower); cfx matrices keep cfxplotnose's 1-column shape; the b-only-post caution is suppressed in fast mode. Off-state proven bit-neutral (21-config battery at rtol 0; a verifier additionally proved commit-isolated bit-identity on 23 configs). Speedups: 2-4x typical; a poisson event-study config that ran >22 minutes without finishing completes in ~5 s; with subgroups >81x. Suite: test/test_fastnose.do. Three verifiers: no defects. Details in wooldid_solved_issues.md.

### E1. Reuse margins results for the SG joint/pairwise tests — RESOLVED (commits cd8ce95, 2a1c9ca, 7ea1533, 184a6ff, ed780ac; 2026-07-26)
One posted pwcompare pass per (estimand x pre/post) x (binary/contreat) now serves both the pairwise table and the all-subgroups-equal joint test (chain `test` built by the new wooldid_sgchain subprogram from posted-e(b) names), replacing the separate contrast(overjoint) passes; restore-disciplined, pre leg capture-wrapped, pwsgs stacking decoupled from the joint-p local. A verifier found the one genuine equivalence break — on a rank-deficient posted VCE, `test` silently reported a residual-rank p where contrast reported missing — closed in ed780ac by a rank guard that defers to the original contrast engine whenever the chain test is not provably full rank, making the reported p identical to the pre-E1 build by construction; with the guard the 21-config battery is fully bit-identical to stock at rtol 0, so the author's abandon condition (same p-values up to ULP) is met exactly. Includes the separable df fix (7ea1533): the binary pre-trend joint tests now thread vce(unconditional) under unconditionalse (test/test_dffix.do). SG section ~1.5x faster; whole command ~20% at 48k obs. The battery gained three stock-blessed coverage configs and the stale baseline was re-blessed (cd8ce95). Three verifiers + coverage union: the single defect above, fixed and confirmed. Details in wooldid_solved_issues.md.







##### Please assess the following item for validity and the viability of implementing with harmless fixes, putting in extra work to check on margins documentation and examples to confirm that resuability would work and do not harm. Do not implement a fix yet. Please provide details on your assesment below. Focus a lot on feasibility and on what changes would need to be mdae to make things work.





------------------------------------------------------
------ DO NOT WORK ON THINGS BELOW THIS LINE ---------
--- ITEMS TO ADDRESS IN DUE TIME ---




### New below-the-line observations from the 2026-07-26 campaign (tranche 10: E1 + F9)
Logged per the PASS-with-observations rule; none acted on silently.
- **Poisson-path inference is z/chi2 by construction:** ppmlhdfe posts no e(df_r), wooldid passes df(.), and margins treats df(.) as no-df — so every margins p on the poisson path is a z/chi2, and e(pwsgs)'s hardcoded "t" column label is wrong there. Pre-existing; both joint-test engines inherit it identically.
- **poisson + jointtests + event-study horizons is CPU-bound** (>12 min at 1,400 obs; >25 min with subgroups) — an extension of the known ppmlhdfe+margins-over(_t) hang class. fastnose cures the fast-mode side; the full-mode hang stands.
- **The impvalues+subgroups r(198) gate is load-bearing for E1** (under impvalues marginsdof is nose; a chain test after a nose post is r321) — a warning comment now sits at the gate; do not relax the pairing without also emptying sgpost.
- **E1 widened the posted-margins abort window:** an abort between a posted pwcompare and its restore leaves e(cmd)="pwcompare" (pre-E1: the estimation model). No natural trigger found; replay correctly refuses afterwards. Capture-hardening is a cheap option if ever wanted.
- **wooldid_sgchain hand-call wart:** `do`-loaded copies expose it globally; after an unrelated estimation it silently writes a nonsense chain into two caller locals (namespaced in normal installs). Milder sibling of the guarded wooldid_esplots wart.
- **wooldid_sgchain `half` invariant:** the split assumes the targetflag posts exactly two realized levels (guaranteed today by the no-control-group r459 upstream); the o.-skip applies to the chain but not the pwsgs slice (an omitted over-level state could not be constructed).
- **update outranks every gate** (pre-existing): `wooldid, update <anything>` reinstalls before validation — including with fastnose. One-line move if gates should ever come first.
- **aw_uncond battery flap risk stands:** unconditionalse cells are not bit-reproducible run-to-run (same-build noise to ~1e-7 on contreat paths); uncond configs cannot be blessed at rtol 0 — the df-fix suite is tolerance-based for exactly this reason.
- **Fast-mode display cosmetics:** removing the joint-p lines merges bracketing blank lines; the header still reads "Standard Errors: cluster(i)" above an all-"." table (one verifier considers a reword mandatory; the author's same-format directive won — flagging the disagreement).
- **fc0cadb commit-message omission (corrective):** the pre-leg prejoint blanking topology change was not called out; strictly corrective and unreachable in all battery configs.
- **altmethod_adversarial is absent from the fastnose gate list** — inert by construction (all uses nest under adversarial, which is gated); confirmed rc 0.

### New below-the-line observations from the 2026-07-26 campaign (tranche 9: F11 + S5)
Logged per the PASS-with-observations rule; none acted on silently.
- **sthlp / e() mismatches (doc):** the help documents e(joint_horizon_bysg), e(joint_horizon_i_byssg), e(joint_cont_horizon_bysg), e(joint_cont_horizon_i_bysg) as matrices — the code emits `ereturn local joint_horiz_s<sg>`-style macros instead; the help also calls e(adversarial_pval_att) a local (it is a scalar).
- **Per-subgroup joint p-values are `ereturn local` while the whole-sample ones are `ereturn scalar`** — pre-existing inconsistency; replay handles both, but a future normalization would be tidier (and a battery re-bless).
- **Poisson + subgroups header lacks a space:** "Poisson RegressionResults for Subgroup ..." (the sg-loop poistext has no trailing space, unlike the main header). Fixing changes console bytes — needs a pre-registered golden-corpus re-bless; do not fold into another change silently.
- **"p-value for joint test that all SGs have the same ame: ." prints a bare "." (and ereturns a missing scalar) when the test is undefined** — pre-existing display choice, preserved by replay.
- **Replay UX notes (author calls):** (1) after a run that fails at parse time, bare `wooldid` silently replays the PREVIOUS successful run (stock Stata semantics; documented in the sthlp) while a late failure gives the arity error — two different post-failure behaviors; (2) `error 301` appends Stata's canned "last estimates not found" after wooldid's informative refusal; (3) a replay of a b-only (no e(V)) run shows blank SEs without the estimation-time caution — a replay-only one-liner would help but breaks strict est-vs-replay byte parity, so it needs pre-registration; (4) e(wdid_replayv) is presence-gated, not version-compared — any future change to the replay contract must bump and compare.
- **contreatwithin silently enables cfx plots for non-att estimands** (it sets `safety off` internally, and the cfx-production gate keys on that) — L432-era blast radius; the new S5 screens were re-gated onto safetyinput to escape it, but the cfx-enablement quirk itself remains.
- **esrelativeto(>=0) WITH esfixedbaseperiod is silently clamped to -1** (worker, pre-existing) — same silently-ignored family as the S5 screen for the no-esf case; promoting the clamp changes e(b_es_*) dimensions, so it needs its own careful wave.
- **Weight expressions are rejected:** [aw=2*w] fails r(111) "variable 2*w not found" (the worker treats exp as a variable name) — pre-existing on both builds; doc or feature candidate.
- **ccc_absorb() never reached a successful run on the battery data** (r(461)/r(198) across spec variants, both builds) — needs a working example or its own look.
- **altmethod_adversarial alone is inert** (it modifies the adversarial method and does nothing without `adversarial`) — doc candidate.
- **Plain poisson + clustered customvce prints "pseudo-R2 = ."** (that path posts no e(r2_p)) — pre-existing cosmetic.
- **A float USER subgroup variable silently pools 8-digit codes** (storage precision, invisible to wooldid) — a Caution when subgroups() is float-typed with large values would be cheap.
- **saveplots prefix + subgroup-value filename concatenation is ambiguous in principle** (`sg`+`20` vs `sg2`+`0`) — never observed colliding.
- **Message-vocabulary stragglers after C12:** the multicollinearity notice uses a fourth head-word ("Note:") and lowercase "warning"; one Caution continuation says "this warning message". Cosmetic.
- **Helper-program namespacing fact:** ado-loaded helpers are namespaced (program dir shows wooldid.wooldid_worker etc.), so wooldid_display/wooldid_esplots/wooldid_validate are NOT globally callable in normal installs — the hand-call guard matters only for `do`-loaded copies (i.e., development).
- **Unexercised corner:** a failure occurring while wooldid is inside its own histogram/IM scratch frames (B22-era shell frame-restore sub-case) resisted five deliberate attempts to trigger it.
- **Dev-surface corner:** hand-calling wooldid_esplots after a real run that produced no ES rows (espre=espl=0) still dies r(2000) inside its tempname frame (the e911347 guard checks presence of results, not ES-row count) — dev-only surface given helper namespacing; a `e(wdid_espre)`-based third guard would close it if ever wanted.
- **Commensurability tier vs timetrends (one-t-unit band):** the tier's lower bound is min(t)+anticipation, but under timetrends the worker discards onsets at or below min(t)+1+anticipation (mintp1); a fractional onset in that one-unit band hard-errors r(459) where the worker would have dropped it. Conservative direction (error rather than silent continue), very narrow corner (s5_v2 confirmation); widen the bound by +1 under timetrends if it ever bites.
- **One unreproduced batch-session death** during an assessor's ppmlhdfe run (no output, exit 0, ~150 program redefinitions in that session) — watch only; never reproduced on either build in fresh sessions.

### New below-the-line observations from the 2026-07-25/26 campaign (tranche 8)
Logged per the PASS-with-observations rule; none acted on silently.
- **Worker program-size ceiling (STANDING CONSTRAINT):** wooldid_worker is ~1.5 kB from Stata's per-program compiled-size limit (measured; inline additions beyond that make the ado fail to load r(1000)). All future features must ship as subprograms (F7/F4 helpers set the pattern; the agent-brief rename recipe was updated accordingly).
- **ppmlhdfe + margins-over(_t) performance:** poisson + esfixedbaseperiod on even the 960-obs battery data runs >10-17 minutes (hit independently by 4+ agents; the esf branch keeps a dummy per relative period, so i._i#i._t explodes under ppmlhdfe). Pre-existing scaling property, deserves its own look.
- **unconditionalse is not bit-reproducible run-to-run** (1-2 ulp drift on identical repeated commands, parent-reproduced). Practice note: rtol-0 comparisons on unconditionalse configs can flicker; the committed battery has not yet been bitten.
- **Internal ttre tempvar is float:** odd calendar codes >2^24 (YYYYMMDD style) round at the gen site (r(461); parent-identical). A gen double would change existing behavior for such data — wants its own careful wave.
- **L1905/1909 df hardcode:** the binary att pre-trend joint tests hardcode df(df_r) while their continuous/SG siblings thread marginsdof, so under unconditionalse+jointtests the binary tests use the conditional VCE. Fix = df + conditional vce(unconditional), NOT a blind marginsdof swap (impvalues sets marginsdof to nose).
- **Adversarial silently cancelled when finaladvtest trips:** the user asked for adversarial, gets no p-value and no message. A one-line Caution at the 4 cancellation sites would fix.
- **Verbose-replay feature path:** retaining _ppmlhdfe_d is necessary but not sufficient — the stored model's regressors are tempvars. A working prototype promoting them to stable names under verbose exists in the session scratchpad (waves/b22_a3); open questions: continuous-mode extras, name collisions, dataset mutation under verbose.
- **Duplicate SG colnames in e(b)/e(V)** (e.g. ES_ame_SG2_post1 x4 in a continuous+sg config, both builds): now that e(V) can exist on such runs (47defa0), _b[]/test can hit ambiguous-name errors where they previously failed for want of V. Own ticket.
- **Plain-continuous *_ame missing-SE shape:** some continuous configs carry missing SEs on all *_ame rows (whole e(V) silently dropped pre-47defa0, now rescued+loud). Whether those SEs should be missing at all is a separate pre-existing question.
- **{good, zero, missing} V-diagonal corner:** would post under 47defa0 where the parent was b-only (unconstructible in testing; zero-diagonals already post today when nothing else is missing).
- **mint fragility:** the pre-drop mint is reused after row-level drops between its computation and the reference-forcing clause, a latent collinearity trigger (the real cause of the F13 wave's P6 r(322)); F13 now masks its only reachable trigger rather than fixing it. P6-at-espre(0) predicted silent, untested.
- **contreatwithin silently overwrites the syntax local `safety`** (S1 pattern): F13's checks are shielded via the safetyinput copy, but every OTHER safety-gated check is still silently disabled under contreatwithin.
- **evertreated tempvar polarity flip** (built as ttre!=. in one block, ttre==. in another; only used as i.evertreated so the fit is identical) — benign, confusing.
- **F2 b_es labels ignore oldsyntax** (post0, never contemp) while escomp and the main matrices honor it — cosmetic inconsistency.
- **pretreatment_att has ~zero power against a linear trend** when the pre window is symmetric about the reference (its weighted mean relyr deviation is ~0) — display-note candidate.
- **F7 partial-ereturn corner:** if wooldid_ptp failed mid-loop, earlier types' roth_* returns would persist with no marker (each individually cap'd; an e(roth_ok) flag would close it). Duplicate power levels (pretrendpower(80 80)) yield a duplicated row.
- **Post-error leftover e() shows tempvar-name shifts between builds** on rc-322 exits (baseline-irrelevant since all battery configs are rc 0; targeted post-error dumps need name normalization).
- **Helper programs are globally visible** (wooldid_worker, wooldid_anticgrid/anticesf/ptpcheck/ptp); a direct worker call leaves the WOOLDID_CLEANUP global set; a user global of that name is cleared at wooldid entry. Namespaced, negligible.
- **anticgrid off-grid caution reads surviving rows only** — a period present in raw data only for units an earlier screen dropped is treated as unobserved (suppression-only; low realism). The esf timetrends second reference screen's anticipation drops remain unattributed (disclosed in e73287d).




### B16 (main-block half) — the SG-crash half was RESOLVED in `4483758` (see wooldid_solved_issues.md)
The two MAIN (non-SG) `altmethod_adversarial`+`esfixedbaseperiod` blocks compute the pre-trend-fit restriction as `espre + r(min) - 1`. When the fixed reference sits below the reported window (`esrelativeto < -1` with `esprelength <= |esrelativeto|-1`, so the ref row becomes the frame minimum), this excludes the relyr -1 coefficient from the fit and tests it as if post-treatment — a possible off-by-one. Deferred per the author's directive to fix only the rare SG-crash scenario in the skiplocation item; the MAIN blocks are the common path and were left untouched. Not bit-neutral if changed (would shift `adversarial_pval_*` on affected runs) -> needs its own careful wave + a re-blessed baseline config.

### NEW (found during skiplocation verification, 2026-07-22): espre=1 + esf + adversarial + altmethod -> r(2001) in the MAIN adversarial block
With `esprelength(1)` and a fixed reference below the window, the MAIN block's `espre + r(min) - 1` restriction selects a single pre-period row, so `qui reg b t` runs on 1 observation and aborts r(2001) ("insufficient observations"). Present identically in parent and the skiplocation fix (MAIN and SG), so pre-existing and NOT introduced by the fix; `espre==0` is unreachable (adversarial is cleared when espre==0, L785). Two skiplocation assessors/verifiers flagged it independently. Low priority, narrow trigger; a capture/guard on the MAIN adversarial `reg` (mirroring how the SG crash was handled) would close it someday.


 


## ➕ Feature suggestions (possible additions)

Update IMPvalues to take advantage of fast mode. 



### F3. ATT(g,t) building-block table
Under the Wooldridge approach the `i._i#i._t#1.targetflag` coefficients *are* the cohort-by-relative-period effects ATT(g,t) — with SEs already sitting in `e(b)`/`e(V)` at zero extra computation (linear case; xb scale for poisson). A `gteffects` option that extracts, labels (cohort × relyr), returns, and optionally displays this matrix would: enable arbitrary custom aggregations downstream, make output directly comparable to `csdid`/`xthdidregress` group-time tables, and give the `histogramcohorteffects` feature a rigorous sibling (the histogram currently uses `nose` margins). Mostly bookkeeping.


### F6. Permutation / Conley–Taber inference for few *treated* clusters
WCB itself degrades when the number of *treated* clusters is very small (MacKinnon–Webb) — the classic 1-treated-state policy evaluation. A permutation option (reassign placebo treatment timing across control clusters, build the null distribution of the ATT estimate) covers this gap; Conley–Taber (2011) is the citable framework. Cost: full re-estimation per draw (reghdfe makes this tolerable; gate behind `reps()`), and design decisions about staggered-timing permutation schemes. Highest-effort item here, but it is the main inference regime wooldid currently has no answer for (IM p-values need many clusters too).




---

## 🏗️ Structural suggestions from the 2026-07 editing campaign (added 2026-07-20)

Lessons distilled from fixing B1–B23: most bugs were not isolated typos but instances of a small number of recurring architectural patterns. Each suggestion below names the pattern, the bugs it produced, and a sketch of the remedy. Ordered roughly by reliability-payoff per effort.

### S1. Treat syntax-parsed option locals as read-only; derive working copies
B7 (cfxplotnose/cfxlattice overwritten mid-loop, silently truncating later estimands' results) and B4 (`ll` reused as an IM CI bound, corrupting e(ll)) were both shared-mutable-local bugs — inevitable in a 3,500-line single program where every local is visible everywhere. A simple discipline eliminates the class: never assign to a local that `syntax` created; when a section needs a variant (an option rewritten for margins, a default filled in), derive a section-prefixed copy (the B7 fix's `cfxseopt`/`cfxlatticeplot` pattern). A one-pass audit of every `local <syntaxoption> ...` assignment left in the file would find any remaining instances. Longer-term, factoring major sections (cfx plots, IM inference, the adversarial machinery, results display) into subprograms would give each its own local namespace and make this structural rather than disciplinary — it would also enable F11's replay refactor.

### S2. Name the sample you mean — a single vocabulary for the four sample concepts
B9, B10, B20, and B23 were all "right check, wrong sample" bugs. The program actually has four distinct sample notions: the user's if-sample (`(if')==1`), the screened pre-regression sample (`ifvar==1`/`ifstatement`), the realized estimation sample (`esamp`=e(sample)), and per-margins subpops. Suggest a short comment block near the top defining these four and which kinds of computation must use which (input validation → if-sample; estimability/construction → ifstatement; weights/summaries/anything feeding results → esamp; display slices → subpop), plus an audit pass tagging every data-touching line with the sample it uses. Most remaining below-the-line items (cfxlattice validation, the pre-period sg gap) reduce to applying this table.

### S3. Stop injecting raw data values into Stata names
The single largest remaining crash surface: subgroup values and relative years are pasted into local/matrix/ereturn/graph names (`wdidsgresults_<sg>`, `_joint_Choriz_i_s<sg>`, `cfx<item>_<nametext>_sg<sg>`, `<est>_sg<sg>`). This produced B20's crash, constrains values to non-negative integers, and still overflows Stata's 32-char name limit at ~9-digit subgroup labels (realistic administrative IDs) — a limit the author flags as genuinely pernicious. Structural remedy: map each subgroup level (and, where relevant, relative year) once to a sequential internal index (`sg1, sg2, ...`), use only the index in every constructed name, and return a crosswalk (e.g. `e(sgmap)` matrix or labeled rownames). This eliminates the whole class — name-length, negative-value hyphens, and non-integer naming crashes — in one refactor, and would let subgroups() accept arbitrary numeric labels. Touches many sites, so it wants its own careful campaign, but it is the highest-leverage robustness change available.


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
