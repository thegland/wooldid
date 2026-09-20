# hv-failed-ereturn — NOTES (state of play)

VERDICT: the claim is **TRUE**, mechanism as stated, on FOUR independent
failure paths and both estimation engines.  One sub-claim (`predict`
silently predicts) is **REFUTED** — predict errors r(111).

## Files
- `probe01.do/.log` — baseline Stata eclass semantics (built-in vs ado wrapper
  vs `_estimates hold`).  PART2 of probe01 is DEAD (my DGP left y missing for
  never-treated units -> wooldid r(459) correctly).  Superseded by probe02.
- `diag01–diag03.do/.log` — diagnosing that r(459); it was my data, not wooldid.
  diag03 V1..V7 all rc=0 once y is non-missing for never-treated.
- `probe02.do/.log` — the FAILURE-CLASS MAP (A..J).
- `probe03.do/.log` — consequence surface + peer attempts.
- `probe04.do/.log` — residue (frame/vars/estimates dir/emptycells), the
  discriminating replay pair, reproducibility of the next good run, peers.
- `repro_failed_ereturn.do/.log` — **the self-contained repro**.  Grep `^R|`.

## Measured facts
Baseline (probe01):
- built-in `regress price nosuchvar` (r111): e() UNTOUCHED, still `regress price mpg`.
- bare eclass ado that runs `regress` then `error 322`: e() = the INTERNAL model.
- same ado with `_estimates hold `h', restore nullok`: e() correctly restored.
- `grep -rl "_estimates hold" /Applications/Stata/ado/base` -> **119 files**.
  Protecting e() is the standard Stata idiom, and wooldid uses it nowhere
  (`grep -n "ereturn clear\|_estimates\|estimates hold" wooldid.ado` -> 0 hits).

Failure-class map (probe02, 60 units x 10 yrs, 12 never-treated):
| id | spec | rc | e(cmd) after | colsof(e(b)) |
|----|------|----|--------------|--------------|
| A | `attt` (option typo, parse) | 198 | **wooldid** (intact) | 1 |
| B | `contreat(dosep) contreatcontrols(untreatedZero)`, i=unit | 322 | **reghdfe** | 291 |
| C | `ccc_absorb(x1)`, i=unit (raised inside `margins`) | 461 | **reghdfe** | 2 |
| D | `fe(sy)` sy = state#year (saturating) | 322 | **reghdfe** | 146 |
| E | `controls(ccol)` collinear w/ one cell | 0 | wooldid | 1 |
| F | `poisson ccc_absorb(x1)` | 322 | **ppmlhdfe** | 291 |
| G | `unconditionalse ccc_absorb(x1)` | 198 (parse) | wooldid (intact) | 1 |
| H,I,J | subgroups / interactivecontrols / cluster(year) | 0 | wooldid | - |

Consequence surface (probe03 / repro):
- `e(cmdline)` = `reghdfe y   i.__000005#i.__000006#i.__000007#c.dosep ...`
- `e(sample)` present and populated; under C it counts **456**, not the user's 600.
- scripted harvest after `cap wooldid ...`: `e(b)[1,1]` = 0, `e(N)` = 600 (or 456
  under C) — **silent**, and both look like plausible wooldid quantities.
- `lincom _b[1.__000005#1.__000006#1.__000007#c.dosep]` -> rc 0, b = .74796608.
  A live numeric postestimation surface on the internal model.
- `predict` -> r(111) "variable __000005 not found".  **REFUTES** the recon
  sub-claim that predict silently predicts.
- `estimates table, keep(main_att)` -> r(198) "coefficient main_att does not
  occur in any of the models" (loud); with a stored good model alongside, the
  failed column is simply **blank** (silent hole).
- `est store` BEFORE the failing call + `est restore` fully recovers
  (`recovered_att` == `truth` == 2.3371429).  That is the only workaround.
- No other residue: frame restored to `default`, `c(emptycells)` back to `keep`,
  no stray variables, no stray `___wooldidfullmodel` in `estimates dir`,
  no `_ppmlhdfe_d` left behind.  The epilogue handles everything **except** e().
- A subsequent good run is bit-identical with and without an intervening failure
  (`+1.2b277f6afeef2X+001` both ways) — the leftover e() does not poison wooldid itself.

Replay pair (probe04 / repro §4):
- after PARSE failure: `wooldid` -> "Error: the last wooldid call failed with
  r(198); the stored results predate that failed call. Refit wooldid before
  replaying."  (the designed NFRPL-01/02 message)
- after POST-REGRESSION failure: `wooldid` -> "Error: command line input should
  feature all of y i t ttre; too few variables specified." (r198)
  The guard at ado:26 is gated on `"`e(cmd)'" == "wooldid"`, which the failure
  itself has just falsified.  The guard is unreachable for the class it exists for.

Peers: jwdid and did_imputation failures I could raise were parse-stage, so they
do not discriminate.  The decisive baselines are the built-in/ado/hold triple and
the 119-file grep.
