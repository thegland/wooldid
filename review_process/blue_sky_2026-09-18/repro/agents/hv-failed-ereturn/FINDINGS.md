# hv-failed-ereturn — FINDINGS

Repro for both findings: `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-failed-ereturn/repro_failed_ereturn.do`
Run: `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_failed_ereturn.do` then grep the log for `^R|`.

---

## HV-EPOST-01 — a wooldid failure after the internal regression leaves reghdfe's (or ppmlhdfe's) saturated model in e(), destroying the user's previous results

**Verdict on the recon claim: CONFIRMED**, mechanism exactly as stated, on four
independent failure paths and both estimation engines.  One stated consequence
(`predict` silently predicting) is **refuted** — see below.

**Claim.** wooldid never clears or holds e().  `grep -n "ereturn clear\|_estimates\|estimates hold"`
over the pinned `wooldid.ado` returns **zero** hits.  The internal regression at
`ado:1476` (`reghdfe`) and `ado:1462,1468` (`ppmlhdfe`) posts to e() immediately;
(the third engine, plain `reg` at `ado:1485` under `unconditionalse`/`customvce`,
is structurally identical but I never reached a post-regression failure on it --
see UNTESTED.) wooldid's own `ereturn post` is at `ado:4051`.  Any
error between those two points leaves the internal model as the session's
current estimation results.  The shell's rc-gated epilogue (`ado:79-97`) restores
the frame, `emptycells`, `level`, and drops the scratch estimates — everything
except e().

**Observed** (60 units x 10 years, 3 treated cohorts, 12 never-treated):

| trigger | rc | e(cmd) after | colsof(e(b)) | e(N) |
|---|---|---|---|---|
| good run | 0 | `wooldid` | 1 | 600 |
| `contreat(dosep) contreatcontrols(untreatedZero)` with i = unit | 322 | **reghdfe** | **291** | 600 |
| `fe(sy)`, sy = state#year (saturating region-by-year FE) | 322 | **reghdfe** | 146 | 600 |
| `ccc_absorb(x1)` with i = unit (raised inside `margins`) | 461 | **reghdfe** | 2 | **456** |
| `poisson ccc_absorb(x1)` | 322 | **ppmlhdfe** | 291 | 600 |
| option typo `attt` (parse stage) | 198 | `wooldid` (intact) | 1 | 600 |

`e(cmdline)` becomes `reghdfe y   i.__000005#i.__000006#i.__000007#c.dosep ...`;
`e(b)` colnames are internal tempvar cell dummies; `e(sample)` is present and
populated (and under the r(461) trigger counts 456, not the user's 600).

**Consequence surface (measured, not asserted):**
- *loud*: `predict` -> r(111) "variable __000005 not found" (the regressors are
  dropped tempvars).  `estimates table, keep(main_att)` -> r(198) "coefficient
  main_att does not occur in any of the models".  A bare `esttab`/`coefplot`
  would render 291 rows.
- *silent*: a `capture`-guarded scripted harvest — the ordinary way an applied
  economist loops over outcomes or subsamples — reads `e(b)[1,1]` = **0** and
  `e(N)` = **600** (or 456), both of which look like plausible wooldid
  quantities.  Nothing distinguishes the failed row from a real one except the
  `_rc` the script may or may not have kept.
- *silent*: `lincom _b[1.__000005#1.__000006#1.__000007#c.dosep]` returns rc 0
  and b = .74796608 with a t and a CI, computed on the internal saturated model.
- *silent hole*: with a good model stored alongside, `estimates table good failed,
  keep(main_att)` prints a **blank column** for the failed run, with no warning.
- *recovery*: `est store` BEFORE the failing call, then `est restore`, fully
  recovers (recovered att == 2.3371429 == truth).  Without a pre-emptive store
  the prior results are gone.

**Expected.**  Protecting prior e() across a failure is the standard Stata idiom:
`grep -rl "_estimates hold" /Applications/Stata/ado/base` matches **119** files.
Measured directly in `repro_failed_ereturn.do` §5:
- built-in `regress price nosuchvar` (r111) leaves e() = `regress price mpg`;
- a bare eclass ado that runs `regress` then `error 322` leaves e() = the
  internal model (this is wooldid's situation);
- the same ado with `_estimates hold `h', restore nullok` leaves e() = `regress price mpg`.
So the one-line fix is demonstrated to work in this exact setting.

**Not residue anywhere else.**  After a failure the frame is `default`,
`c(emptycells)` is back to `keep`, there are no stray variables, no stray
`___wooldidfullmodel`, no `_ppmlhdfe_d`, and the next good run is bit-identical
(`+1.2b277f6afeef2X+001`) with or without an intervening failure.  e() is the
single unhandled piece of state.

**Realistic:** yes.  All three binary/continuous triggers are ordinary
mis-specifications (a dose variable with unit-level `i`; saturating region-by-year
fixed effects; `ccc_absorb()` with unit-level `i`), and the pattern "run a
baseline, then try a variant that errors" is the normal interactive workflow.

**Severity: medium.**  It does not corrupt a successful run's numbers.  It costs a
user their previous results with no warning, and it feeds wrong numbers to
`capture`-guarded scripts and to `lincom`-style postestimation.

**Suggested fix (NOT verified by me):** `_estimates hold `h', restore nullok`
before `capture noisily wooldid_worker` at `ado:71`, then branch on `wdsh_rc`:
restore the held results when rc != 0, and **DISCARD the hold when rc == 0**
(`_estimates unhold `h', not`, or drop the `restore` suboption and unhold
explicitly only on the failure branch).  A bare `restore` hold auto-restores at
shell exit and would erase the worker's e() on every SUCCESSFUL run -- the
success branch must cancel it.

---

## HV-EPOST-02 — the replay-staleness guard is unreachable for exactly the failure class it was written for

**Found by me while verifying HV-EPOST-01.  This CONTRADICTS the design stated in
`ado:26-46` (NFRPL-01/02).**

The guard at `ado:26` is `if `"`0'"' == "" & "`e(cmd)'" == "wooldid"`.  Its
twenty-line comment says it exists so that "while [WOOLDID_LASTRC] is set, the
stored e() results predate the user's latest (failed) attempt, so redisplaying
them silently would be misleading."  But a post-regression failure sets
`e(cmd)` = `reghdfe`, which falsifies the guard's own gate.

Measured (repro §3a vs §4, same rc, different diagnosis):
- after a **parse-stage** failure, bare `wooldid` gives the designed message:
  `Error: the last wooldid call failed with r(198); the stored results predate
  that failed call. Refit wooldid before replaying.`
- after a **post-regression** failure, bare `wooldid` falls through to the worker
  and gives `Error: command line input should feature all of y i t ttre; too few
  variables specified.` — r(198), about a command line the user did not type.

So the guard fires only for failures that did **not** damage e(), and is silent
for every failure that did.  A user who types `wooldid` to redisplay is told the
wrong thing precisely when the truth matters most.

**Fix shape:** gate the staleness branch on `"${WOOLDID_LASTRC}" != ""` alone (it
is already a shell-level global set on every failed invocation), rather than on
`e(cmd)` surviving the failure.  Fixing HV-EPOST-01 would also fix this, since
`_estimates hold ... , restore` would put `e(cmd)` == `wooldid` back.

**Severity: medium.**  Uninformative error on realistic input; no wrong number.

---

## REFUTED sub-claim
Recon asserted "`predict` silently predicts from the internal regression."  It does
not: `predict double xbhat` after the failure errors **r(111), "variable __000005
not found"**, because the internal regressors are tempvars that no longer exist.
The claim's other stated consequences all reproduce.

## Not a discrepancy
Recon reported `e(b)` as 1x1083.  I measure 291 / 146 / 2 depending on the
trigger; the dimension is just the internal cell-dummy count for the panel used.

---

## Addendum: the author already applies this protection one level up

`grep -c "ereturn clear" wooldid.ado` -> **0**.
`grep -c "_estimates hold\|_est hold" wooldid.ado` -> **0**.

But the shell DOES protect e() around both replication engines:
- `ado:168` `qui est sto `wdpobs'` / `ado:172` `qui est restore `wdpobs'` — around
  `wooldid_permrun`, run unconditionally so a rep-loop failure cannot leave a
  replication's e() behind (this is the NFRPL-01 exception, and the code backs it).
- `ado:367` / `ado:396` — the same pair around the bootstrap rep loop.
- `ado:1494` `qui est sto `wdfm'` inside the worker, restored at 1797, 3092, 3132,
  3220, 3432, 3478 …

So the store/restore idiom is used routinely throughout wooldid.  The single place
it is missing is the one that matters to the user: the `capture noisily
wooldid_worker `0'` call at `ado:71`.  That asymmetry is the cleanest statement of
the defect — the permutation engine's temporary e() is protected from the user,
and the user's e() is not protected from the worker.


---

## UNTESTED
- The plain-`regress` engine (`ado:1485`, `unconditionalse`/`customvce`): structurally
  identical, but my attempt (`unconditionalse ccc_absorb(x1)`) was refused at parse
  (r198), so I never reached a post-regression failure there.
- The permute/bootstrap rep-loop protection (`ado:168/172`, `367/396`) is a CODE READ,
  not a measurement: I did not force a rep-loop failure and confirm e() is restored.
- A peer command failing AFTER its internal regression: both jwdid attempts and the
  csdid attempt were parse-stage, so the peer comparison does not discriminate.  The
  built-in baseline (`regress price nosuchvar`) is also parse-stage; the argument that
  protecting e() is standard practice rests on the 119-file `_estimates hold` grep and
  on the measured `zz_withHold` demonstration, not on that run.
- `coefplot`/`outreg2` after the failure: inferred loud from the 291-row e(b), not run.
