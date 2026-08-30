# Working process for wooldid changes

The coordinator (Claude) owns this file; agents follow `agent_brief.md`. Updated
2026-08-18: the fixed per-issue lifecycle (3 assessors + 3 verifiers + full battery for
every change, regardless of size) is retired. Scale scrutiny to the risk and blast radius
of the change. Below: the review ladder, then the invariants that hold at every level —
most of them are lessons this campaign paid for empirically; do not drop them casually.

## 1. Scale review to risk
- **Every issue gets the coordinator's own written analysis first**, locked in writing
  BEFORE reading any reviewer's report (verdict, fix shape, risks, expected test
  footprint). This is cheap and has repeatedly caught anchoring; keep doing it.
- **Low risk** (doc fixes, message wording, dead-code removal, mechanical one-liners with
  obvious blast radius): own analysis + targeted checks, then commit. No subagent review
  required.
- **Medium risk**: add at least one independent Opus-class adversarial reviewer, tasked
  with probing the assessment (if rejecting an issue) or the proposed change (if fixing).
  Update the plan based on what it finds.
- **High risk**: after updating the plan from the first review, spin up a fresh
  adversarial reviewer with no knowledge of the prior review's content; iterate until
  coordinator and reviewer agree the course of action is correct.
- Full multi-assessor + multi-verifier waves remain available for the hardest structural
  changes — use them when independent lenses genuinely add information, not by default.
- Empirical testing (in Stata, on purpose-built data) is encouraged at every level and
  expected for medium+ issues where static analysis is unsuitable. Reviewers and the
  coordinator alike should test assumptions about Stata behavior (help files, online
  docs, or quick probe scripts) rather than argue from memory.

## 2. Invariants — coordination and code safety
- **Edits must not clash.** Pin the sha you analyzed; re-locate code before editing (line
  numbers drift); one issue per commit; serialize commits (never two agents editing the
  repo at once — implementation happens in the coordinator's hands, agents' repos are
  read-only).
- Never overwrite a syntax-parsed option local (derive prefixed copies — S1). Note
  validate NORMALIZES many option locals; harvest post-normalization values, and beware
  mid-line `if COND local x ...` mutations that line-anchored scans miss.
- State explicitly which of the four sample concepts each new data-touching line uses
  (user-if / screened ifvar / realized e(sample) / margins subpop — S2).
- New early-exit/validation checks get the standard probe: enumerate the full
  weights × mode × estimator lattice for inputs that currently SUCCEED (the B15 pweight
  lesson — a guard that rejects working calls is a regression, not a safety net).
- New features go in subprograms: `wooldid_worker` runs close to Stata's per-program
  size ceiling (~133.5 kB of leading-whitespace-stripped bytes).

## 3. Invariants — evidence and adjudication
- **Pre-register the expected footprint of any comparison before running it** ("expect
  identical except <exact cells>"). Any deviation from the pre-registered signature is an
  automatic stop-and-investigate.
- **Evidence beats vote**: when reviewers disagree on design, adopt the design backed by
  an empirical refutation of the alternatives (B14/B17/B20/B23 precedent), not majority.
- Disagreement with no safe workable path → leave the item unfixed with written
  commentary; that is a success, not a failure.
- Commit-message claims are testable assertions: every empirical statement traces to a
  test that actually ran; soften anything unverified. (Verifiers audit messages and have
  caught overstatements.)

## 4. Testing
- **Per-fix testing = relevant labeled suite components + purpose-built probes**, not the
  full suite. The full testing-and-assessment suite (see `test/README.md`) is run once
  after a large battery of changes, campaign-style — it is too hefty for per-fix use.
- Battery baselines are re-blessed only with an intentional-signature commit that
  documents the expected diff.
- Old-vs-new comparisons use rtol 1e-10 as the project indifference threshold; the
  clustered plain-reg/poisson ULP class (`unconditionalse`, `customvce(cluster)`) is
  never bit-reproducible — see the brief's tolerance note.
- **Pathological data answers a different question than well-behaved data** (author
  doctrine 2026-08-29; canonical statement in `test/suite/README.md` "Role of
  pathological data", classes in `test/suite/data_classes.csv`). A pathological-data
  test asks "correct results or observable failure?"; feature-quality evidence comes
  from well-behaved data, and features should tend to be tested on both. When a defect
  surfaces on pathological data, probe the same feature on well-behaved input before
  diagnosing — combine the two signal streams.

## 5. Follow-ups and close-out
- Narrow change requested BY a reviewer (one-clause guard tweak, doc line) → implement,
  test, commit, then incremental confirmation by the same reviewer. Anything with new
  scope → fresh review at the appropriate risk level.
- Reviewer PASS-with-observations → log the observations in the review file; do not act
  on them silently.
- Move fixed issues + a **RESOLVED** paragraph (commit sha, what actually shipped and why
  it differs from the suggestion if it does, verification summary) to
  `wooldid_solved_issues.md`. Rejected/illusory issues get a written rationale for the
  author. Update the coordinator's memory checkpoint so protocol state survives
  interruptions.

## Standing environment facts
- Stata binary: `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do file.do`
  (the PATH symlink is broken).
- **Compute cap: at most 8 concurrent Stata instances in total** (each agent limited
  to 1). There is NO cap on the number of concurrent subagents — only on how many may do
  serious Stata computation at once, via the 8-instance rule.
- gtools is a hard requirement since `40414c7`; reghdfe/ftools/gtools are checked at
  entry (r198 if missing).
- Author's edge-case calibration: realistic messy data and option combinations matter;
  deliberately absurd inputs don't — except name-length (32-char) hazards, which are real.
