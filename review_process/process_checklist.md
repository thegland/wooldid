# Per-issue lifecycle checklist (wooldid editing campaign)

The full protocol is defined in the game-plan section of `wooldid_ado_review.md`; this file
codifies the operational details adopted 2026-07-21 after tranches 1–2 (B1–B23). The
coordinator (Claude) owns this checklist; agents follow `agent_brief.md`.

## 1. Assessment wave (3 independent opus-class agents)
- [ ] Each assessor gets: the issue text verbatim, the pinned HEAD sha, a distinct
      emphasis/lens, the brief path, and its own working directory.
- [ ] Prompts require the coverage grid (`assessment_grid.md`) as the report's closing
      section, with untested-relevant cells (`U`) listed honestly.
- [ ] At most 6 Stata instances across ALL concurrent agents (each agent limited to 1);
      at most ~18 agents in flight.

## 2. Coordinator's independent analysis
- [ ] Completed and **locked in writing before reading any assessor report** (state the
      verdict, fix shape, risks, and predicted battery signature).

## 3. Adjudication
- [ ] Unanimous GENUINE + convergent design → proceed.
- [ ] Design disagreement → **evidence beats vote**: adopt the design backed by an
      empirical refutation of the alternatives (B14/B17/B20/B23 precedent), not majority.
- [ ] Disagreement with no safe workable path → abandon, leave the item in the review file
      with commentary explaining why (per the game plan, this is a success, not a failure).
- [ ] Issue judged illusory by the wave and own analysis → leave in the review file with
      commentary; do not fix.

## 4. Implementation
- [ ] Never overwrite a syntax-parsed option local (derive prefixed copies — S1).
- [ ] State explicitly which of the four sample concepts each new data-touching line uses
      (user-if / ifvar / esamp / subpop — S2).
- [ ] **Pre-register the battery signature in writing before running it** ("expect 18/18
      IDENTICAL except <exact cells>"). Any deviation from the pre-registered signature is
      an automatic stop-and-investigate.

## 5. Testing (coordinator's own, before any commit)
- [ ] Full battery (`test/`) vs the current baseline → must match the pre-registered
      signature exactly.
- [ ] Targeted suite for the fix itself: the bug corrected, the do-no-harm cells preserved
      (bit-compare against the parent build and/or physically-subset references), relevant
      pathological shapes (`test/data/dta_path_*`).

## 6. Commit
- [ ] Commit message claims are testable assertions: every empirical statement traces to a
      test that actually ran; soften anything unverified. (Verifiers audit messages and
      have caught overstatements — B14, B17, B20.)
- [ ] One issue per commit; battery baselines updated only with an intentional-signature
      commit.

## 7. Verification wave (3 fresh opus-class agents)
- [ ] Standard lenses: adversarial refutation / independent behavior matrix on fresh seeds /
      structural audit + commit-message fidelity.
- [ ] One verifier is explicitly assigned **the union of the assessment wave's `U` cells**
      (this is the single highest-value slot: both campaign escapes lived there).
- [ ] Verifiers pin to the commit sha; warned that the working tree drifts.
- [ ] Any new early-exit/validation check gets the standard probe: enumerate the full
      weights × mode × estimator lattice for inputs that currently SUCCEED (the B15
      pweight lesson).

## 8. Follow-ups
- [ ] Narrow change requested BY the verification wave itself (one-clause guard tweak, doc
      line) → implement, test, commit, then **incremental confirmation** by the same
      verifiers via a follow-up message (B14/B15 precedent).
- [ ] Anything with new scope → full fresh wave.
- [ ] Verifier PASS-with-observations → log observations below-the-line in the review
      file; do not act on them silently.

## 9. Close-out
- [ ] Move the issue's text + a **RESOLVED** paragraph (commit sha, what actually shipped
      and why it differs from the suggestion if it does, verification summary) to
      `wooldid_solved_issues.md`; delete it from the review file's work section.
- [ ] Abandoned/illusory items stay in the review file with commentary (see step 3).
- [ ] Update the coordinator's memory checkpoint (protocol state survives interruptions).

## Standing environment facts
- Stata binary: `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do file.do`
  (the PATH symlink is broken). Compute cap: 6 concurrent Stata instances total.
- gtools is a hard requirement since `40414c7`; the non-gtools code path no longer
  exists. The old D11 "gtools-neutralized" leg is now a dependency-check negative test.
- Author's edge-case calibration: realistic messy data and option combinations matter;
  deliberately absurd inputs don't — except name-length (32-char) hazards, which are real.
