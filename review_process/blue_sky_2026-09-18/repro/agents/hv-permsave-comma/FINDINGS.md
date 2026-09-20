# hv-permsave-comma — FINDINGS

Build asserted in every log: `which wooldid` →
`/private/.../scratchpad/bluesky/build/wooldid.ado`, `*! version 1.4 2026-08-19`.
All data generated in-file; no repo data, no repo generator, no test suite.

---

## F1 — permsave() has no comma check; `permsave(f, replace)` writes "f, replace.dta" at rc 0
**VERDICT ON THE SEEDED CLAIM: TRUE, mechanism exactly as stated.**  (seeded)

* repro: `p1_comma.do` → `p1_comma.log`
* severity **medium**, realistic **true**, confidence **certain**, surface `permsave()` / `bootsave()`

Measured:

| test | command tail | rc | on disk |
|---|---|---|---|
| T1  | `permute(shuffle) permreps(20) permsave(draws, replace)` | **0** | `draws, replace.dta` (56 213 B, 21 obs, real draws) |
| T1b | the same line again | **602** | — |
| T2  | `bootstrap(standard) reps(20) bootsave(bd2, replace)`     | **198** | nothing |
| T6  | `use draws`                                               | **601** | — |
| E   (p4) | `permsave(d1.dta, replace)`                          | **0** | `d1.dta, replace.dta`; no `d1.dta` |

T1 prints, in plain text, at rc 0:
```
permsave(): the permutation draws were written to draws, replace.dta (20 replication rows
plus the _rep = 0 observed row; the onset_* columns hold the onset assigned to each cluster).
```
T2's twin prints the helpful refusal RF-6 added on 2026-09-13 (ado:6281-6289):
```
Error: bootsave() takes a plain file name and no comma; replace is a separate option of
wooldid, so write bootsave(bd2  replace) replace.
```
T1b is the UX punchline — the user who *wrote* `replace` is told to add it:
```
Error: permsave(draws, replace) - draws, replace.dta already exists in the working
directory; add the replace option to overwrite it.
```

Mechanism (verified in source): `wooldid_permctl` ado:5868-5880 tests only the `.dta`
suffix, `"/"`, `"\"` and the empty string. `wooldid_bootctl` ado:6279-6298 runs the same
three tests **plus** `strpos("`wdbfn'", ",")`. The comment at ado:5866 asserts
"permsave(): the name rules are bootsave()'s" — that assertion is now false. The writer
programs (`wooldid_permsave` ado:8845-8854, `wooldid_bootsave` ado:6532-6541) both lack
the comma test; only the `ctl` pre-flight has it, and `wooldid_bootctl` is called
unconditionally (ado:4941), so there is no reachable gap on the boot side.

Everything else about the two savers is symmetric and correct: `.dta` stripping, the
directory-part refusal, the pre-flight existence check, and the `replace` semantics
(T3a/b/c = 0/602/0 for permsave; T4a/b = 0/602 for bootsave).

Consequence surface: `saving(f, replace)` is the universal Stata idiom, and `permsave()`
is the one option in wooldid that silently accepts it. The cost is not a wrong number —
the file's contents are correct — but a full permutation run (the expensive kind: every
replication is a complete re-estimation) lands under a name containing a comma and a
space, `use draws` fails r(601), and re-running the same script fails r(602) telling the
user to add an option they already typed. `saveplots()` is declared `SAVEplots(name)`, so
the same slip is a parse error there; `permsave()` is declared `PERMSAve(string)` and is
the only exposed instance.

---

## F2 — under `boottstat` the studentized joint/equality statistics are not posted, though the help says they are  (self-found)

* repro: `f2_jointwobs.do` → `f2_jointwobs.log` (also p3_replay.do §B)
* severity **low**, realistic **true**, confidence **certain** (measurement), surface
  `bootstrap()` + `boottstat` + `jointtests` + `subgroups()`

Measured, same data and same options apart from `boottstat`:

| e() matrix | studentized (`boottstat`) | fixedmetric (no `boottstat`) | `permute()` twin |
|---|---|---|---|
| `boot_jointW`    | present | present | `perm_jointW` present |
| `boot_jointWobs` | **absent** | present | `perm_jointWobs` present |
| `boot_sgW`       | **absent** | present | `perm_sgW` present |
| `boot_sgWobs`    | **absent** | present | `perm_sgWobs` present |
| `boot_advWobs`   | **never posted on any path** | n/a | `perm_advWobs` present |

The studentized run nonetheless **reports** the tests it is not storing the statistics
for: `p-value for joint test that all SGs have the same att: 0.2500`.

Source: ado:6871-6872 and ado:7152-7153 post all four on the fixedmetric branch;
ado:6946 and ado:7007 post `boot_jointW` and `boot_advW` alone on the studentized branch.
`grep 'ereturn matrix boot_'` shows no `boot_advWobs` anywhere. The permutation engine
posts all six companions (ado:8443-8444, 8570-8571, 8638-8639).

What the help says:
* `sthlp:1340` (reference list) carves out **only** the `*obs` pair — "the observed ones
  in e(boot_jointWobs) / e(boot_sgWobs) (the latter two without boottstat)" — while still
  listing `e(boot_sgW)` as where the replication statistics are. `e(boot_sgW)` is absent
  under `boottstat`. **This is the crisp contradiction.**
* `sthlp:1340` also claims a `boottstat` run stores enough "so every reported studentized
  number is reconstructable from e() by name". The studentized joint and adversarial
  p-values are add-one ranks of W_r against W_obs, and W_obs is nowhere in e() by name
  (it can be recomputed from the preserved analytic covariance, but not read).
* `sthlp:1391` promises, with no carve-out, "under bootstrap() it is an add-one rank p of
  the replication statistics, of the family named by e(boot_jt_family) (studentized or
  fixedmetric), the statistics themselves in e(boot_jointW)/e(boot_jointWobs)".
* `sthlp:532` is ambiguous on the same point.

**The four absences are TWO different gaps, which matters for the fix:**
1. `boot_jointWobs` / `boot_advWobs` — the observed statistics DO exist at run time as the
   tempnames `wdbjobs` / `wdbaobs` (ado:365, handed to `wooldid_bootsave` at ado:437, and
   measured in row 0 of `jointW_att`). Not ereturning them is a two-line omission.
2. `boot_sgW` / `boot_sgWobs` — `wooldid_bootsgw` (Mata, ado:15380-15450) takes only
   `outloc` and `nloc` as outputs and returns the p-value and the completed count. It
   builds `Wobs` and the `Wr` colvector internally and **discards them**: there is no
   matrix argument to write them to. Under `boottstat` the studentized equality statistics
   are computed and thrown away, so the help's `e(boot_sgW)` claim cannot be met without
   changing that function's signature. This is a design gap, not a posting oversight.

**The reading of sthlp:1340 does not matter.** "(the latter two without boottstat)" most
naturally means the `*obs` pair, in which case `e(boot_sgW)` is promised and absent; under
the alternative reading (the two `sg*` matrices) `e(boot_jointWobs)` is promised and
absent. Either way one promised matrix is missing. sthlp:1391 carries no carve-out at all.

Where the number does survive: `bootsave()` writes the observed W into row `_rep == 0` of
`jointW_*` (measured: `jointW_att` row 0 = .2704913 in `e2.dta`, p4 §D), read from the
`jobs()` tempname rather than from `e()`. So the observed statistic exists at run time and
is simply not ereturned on that branch.

---

## F3 — `dots(-5)` is r(198) under `bootstrap()` and silently accepted under `permute()`  (self-found; the ado admits it in a comment)

* repro: `f3_dots.do` → `f3_dots.log`
* severity **low**, realistic **false**, confidence **certain**, surface `dots()`

Measured: `bootstrap(standard) reps(20) dots(-5)` → **r(198)**
"Error: dots() must be a nonnegative integer (0 suppresses the progress display)."
`permute(shuffle) permreps(20) dots(-5)` → **rc 0**, run completes, no progress display.

Mechanism: the range check lives in `wooldid_bootctl` (ado:6272), below the
stray-companion screen that `exit 0`s whenever `bootstrap()` is absent (ado:6225). Under
`permute()` the check is never reached, every `if `dots' > 0` test in the engine fails and
the value behaves as `dots(0)`. The ado's own comment at ado:6268-6271 states this.
Filed only because it is a measured twin asymmetry; nobody types `dots(-5)`.

---

## NEGATIVE RESULTS (attacked hard, came back clean)

1. **e() read from inside a new frame.** `wooldid_bootsave` passes literal `e(boot_jointW)`
   / `e(boot_advW)` strings into Mata that are resolved **inside** `frame `bsf' { ... }`,
   and `wooldid_permsave` does the same with `e(perm_ase_obs)`. Hypothesis: frame-local
   e() would silently blank those columns. **REFUTED** (p2_twins.do §E1-E3):
   `jointW_att` and `advW` columns are present and row 0 carries the observed statistic;
   `se_main_att` row 0 = .26416984, non-missing. e() survives the frame switch.
2. **The documented rngstate replay guarantee** (sthlp:542: "a run without rseed() can
   still be replayed exactly by set rngstate to that value before repeating the call").
   **HOLDS on both engines** (p3_replay.do): `set rngstate e(perm_rngstate)` then the same
   permute call → `cf _all` rc 0 against the first run's `permsave` file; same for
   `e(boot_rngstate)` + `bootsave`. `rseed(4242)` twice → `cf _all` rc 0. If either engine
   consumed the RNG between the recorded state and the first draw, cf would have flagged
   every column.
3. **`permsave` + `permonsetcohorts`** (§7, one day old at the time of probing):
   60 `onset_*` columns for 60 clusters, no D-10 caution, rc 0 (p2_twins.do §E6).
4. **Save-file layout parity** (p2_twins.do §E4): permsave and bootsave produce the same
   `_rep` / coefficient / `se_*` block under `permtstat` / `boottstat`, the same
   `main_SG0_att` / `main_SG1_att` naming under `subgroups()`, and each adds its own
   design block (`onset_<clus>` vs `w_<clus>`). No layout divergence found.
5. **Forced `w_` / `onset_` name collision** (p4_dropcol.do §C). `wooldid_bootsavemk`
   skips a colliding `w_`/`jointW_`/`advW` column with a bare `continue` (ado:14912,
   14921, 14929) and prints nothing, while `wooldid_permsave` counts the analogous
   `onset_` drops and prints the D-10 cautions (ado:8877-8880) — i.e. the D-10 ruling of
   2026-09-11 looks like it shipped on one twin only. I could **not** make the collision
   fire: cluster ids `.1` and `-1` render as `onset__1` and
   `onset__1000000000000000<...>` (the double `.1` is not `.1` under `%20.0g`), so 60
   clusters still gave 60 columns on both savers. Left as an unverified read-level
   observation, not a finding.
