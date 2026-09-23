# Change audit, 2026-09-22 - the record for the author's rulings

**Status line, dated: opened 2026-09-22 late evening; the audit itself is running.** This file
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
| R-A | Cards #29 (impvalues with a finer or cross-cutting cluster) and #214 (bootstrap(bayes) with cell-weighted estimands): you unticked "Silenced by safety(off)" in round 3 | build both (the untick is the ruling; test and help changes ride along: test_fastcollapse F15, test_bootstrap_gates G6, test_bootstrap_std:77, test_bootstrap_joint:152/229, test_bootstrap_bayes:332, help ~241/250/268/357, guide ~177; #214 also kills caution #215) / keep the escapes (safety(off) stays the documented way past both) | HELD; both cards carry a feedback badge |
| R-B | Subgroup containment guards S1 and S2 (the worker's split-unit and within-cell errors): the ticket says harden all three; your page edits hardened S3 only and kept "specify safety(off)" in S1/S2's text | harden S1/S2 too (their escape sentences go) / leave them liftable (ticket closes as S3-only) | HELD (open_ticket_sgshareweight) |
| R-C | Save failure, Q2: each saver re-runs the directory-part r(198) and file-exists r(602) checks after the loop; they are still fatal, so a file that appears mid-run still throws the finished run away | convert the r(602) leg to the same caution / delete the duplicates as redundant / leave | HELD (open_ticket_savefail_caution) |
| R-D | Save failure, Q3: the caution prints above the deferred display and scrolls off | move it below the tables / leave | HELD |
| R-E | Card #409: a caution when a bootstrap run completes more than 10,000 replications and e(bootdist) is silently not posted, the twin of #408 on the permutation side | approve the drafted wording (then it is built) / edit it / delete the card and leave the bootstrap side silent | HELD; wording is mine |

## 2b. Rulings given 2026-09-23 (multiple choice, in session)

| item | ruling | state |
|---|---|---|
| R-A #29 impvalues finer/crosscut cluster; #214 bootstrap(bayes) with cell-weighted estimands | **Build both**: safety(off) no longer lifts either; caution #215 goes | to build |
| R-B subgroup guards S1, S2 | **Leave them liftable**; ticket closes as S3-only | done (ticket text) |
| R-C save failure Q2, the post-loop duplicate file-exists check | **Leave it fatal** | no change |
| R-D save failure Q3, caution position | **Move below the tables** | to build |
| R-E card #409, bootstrap draw-cap caution | **Approve Claude's wording** | to build |
| esband-normal-law-vs-t-pointwise | **Multivariate t with df_r, floored at the pointwise t quantile**; re-bless band configs | to build |
| ptp-nis-tcrit-normal-law (blue-sky PTP-2) | **Roth: normal throughout**; help adds "under the normal reference" | to build |
| permtstat-cfxplotnose-blank | **Override cfxplotnose**, caution approved verbatim: "Caution: cfxplotnose was dropped; permtstat needs the plot cells' standard errors." | to build |
| sgdeg-inert-under-suppresssg | **Make the guard work** under suppresssgprimaryeffects | to build |
| amedetect-onerow-cell | **Exclude one-row cells** from the containment predicate | to build |
| ccscreen-contreatwithin-amevoid | **Run the void screen under contreatwithin** | to build |
| ccref-nevergroup-false-refusal | **Keep the refusal** | no change |
| perm-window-labels | **Recode late onsets in the permutation engine, and caution about it** (wording to the author first) | to build |
| val-permanalyticse-fastnose | **Refuse the pair** (widened message to the author first) | to build |
| fcol-lattice-cfx-key-liftable | **Lift the refusal**: the bin joins the collapse key | to build |
| perm-window-labels caution text | approved verbatim: "Caution: `n' cluster(s) with onsets after the last period are permuted as never treated." | to build |
| val-permanalyticse-fastnose text | approved: "Error: fastnose is not compatible with permtstat, permtci or permanalyticse." | to build |
| perm-save-comma | **Add bootsave's refusal**: "Error: permsave() should only contain a file name." | to build |
| permctl-partial-label-no-twin | **Add the screen, unliftable**, text approved: "Error: the subgroup variable `subgroups' is missing for some in-sample row(s) of units of `i' that carry a subgroup value on their other rows. This plausibly is an error, and is not compatible with non-att estimands." | to build |
| perm-esttab-analytic-under-ase | **Always permutation**: e(p)/e(ci_l)/e(ci_u) posted from the kernel on every permute() run | to build |
| perm-sgeq-uncentred | **List flagged contrast rows in e(perm_offcentre)** (existing footer names them); the equality test's p = 1 goes to research | to build + research |
| saveplots-no-precheck | **Keep Stata's own r(602) text** for the early refusal | to build |
| fcol-customvce-holdout-liftable | **Measure first** (one settling run), bring the result back | research |

## 3. The audit

(filled in when the workflow returns)
