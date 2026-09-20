# Coordinator repro do-files

These are the five runs the coordinator did personally, rather than delegating. Each is
self-contained apart from one line: `adopath ++ "<path>"` points at the session scratchpad, which
is temporary. Repoint it at a checkout of the pinned build `c1d8979`:

    git show c1d8979:wooldidinstall/wooldid.ado > <dir>/wooldid.ado

Then run e.g. `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do rem.do`.

| file | what it establishes |
|---|---|
| `nz.do` | README:141 / sthlp:164's "normalized zero in the earliest reported period" is false under the pooled default |
| `rem.do` | coarsening `i` to reach the documented remedy geometry moves the point estimate on a ragged panel; `fe(unit)` repairs it exactly |
| `femc.do` | the 200-draw MC showing `fe(unit)` preserves the coverage gain (0.970) and cuts sampling SD 24% |
| `uz.do` | `contreatcontrols(untreatedZero)` + `poisson` on the plain default engine returns ame -0.0284 where truth is +0.25 |
| `rk.do` | bootstrap `e(V)` rank at reps 8 / 40 — the measurement that reframed the rank claim |

Agent repro files are NOT copied here; they live in the session scratchpad and are named in each
finding's `repro` field in `findings_raw.json`. They do not survive the session.
