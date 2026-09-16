# Guard-cost investigation probes (2026-09-16)

Supporting material for `review_process/guard_cost_2026-09-16.md`. Each `.do` file expects an
instrumented copy of `wooldid.ado` beside it under the name it `do`-es; the instrumentation is
generated, not stored, because it is mechanical and would otherwise rot against the real file.

| file | what it does | how to regenerate its ado |
|---|---|---|
| `guardmap.py` | walks `wooldidinstall/wooldid.ado` tracking brace depth and program boundaries, tags every expensive gtools op with its INNERMOST enclosing condition | run from the repo root; reads the real file |
| `count.do` | counts how many times the two `safety` guard blocks execute under analytic / bootstrap / permute | insert `global F13N`/`SI8N` increments at the top of the blocks at `ado:762` and `ado:1341` |
| `fire.do` | counts how many times those blocks actually FIRE, tagged by `wdrepmode` | insert the same counters plus a `global F13FIRE`/`SI8FIRE` stamp inside each failure branch |
| `calib.do` | the control for `fire.do`: breaks one unit's onset inside a cohort so F13 genuinely must fire, and checks `safety(off)` silences it | same instrumented build as `fire.do` |
| `si8.do` | prints `si8cells/ilev/tlev/demand` against `c(maxvar)` on the fixture that raised Stata's own r(907) | insert one `dis` after `local si8demand` at `ado:1355`; no logic change |

The generator used for all of these is a short Python insert-by-line-number script; the anchors
are asserted before insertion so it fails loudly rather than instrumenting the wrong line.

`calib.do` is the important one. `fire.do` on its own reports "the check never fired in a
replication", which is indistinguishable from a dead probe - the failure shape this slate hit
four separate times. `calib.do` establishes that the same instrument fires when it should.
