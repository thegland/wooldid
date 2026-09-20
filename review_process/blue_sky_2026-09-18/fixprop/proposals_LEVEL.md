# Fix proposals — mechanism `LEVEL` (`level()` / `c(level)`)

Findings covered: **AR-4** and **SKEP-LEVELSTATE**. They are one defect with one fix site
cluster, so this is one proposal (L1) with three edits, plus one optional companion (L2)
and one "do not fix" (L3).

Everything below was **tested**: the three edits were applied to a private copy of the
pinned build at
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-LEVEL/patched/wooldid.ado`
and run against the pinned build on the same do-files. Logs are in
`/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/prop-LEVEL/`.

---

## 0. Mechanism, restated — with one correction to the findings' scope

`set level` parses the **literal text** it is handed and refuses more than two decimals.
wooldid hands it text produced by Stata's 16-significant-digit rendering of a double, at
two sites:

| site | code | what it renders |
|---|---|---|
| `ado:5507` (validate) + `ado:576` (worker) | `if \`level' == 0 local level = c(level)` … `set level \`level'` | the requested level: from `syntax , Level(real 0)`, or filled from `c(level)` |
| `ado:575` (worker) + `ado:70` (shell) | `global WOOLDID_LEVEL0 = c(level)` … `set level $WOOLDID_LEVEL0` | the caller's level, for the restore |

`level(99.9)` therefore reaches `set level` as `99.90000000000001` and is refused r(198).

**Correction to both findings (measured, `pB_scope2.log`).** The prober's rule
`{L : %.16g(L) != %.15g(L)}` — 1152 of 9000, floor 64.01 — is right for the **`level()`
option** limb only. The **session-level** limb is strictly wider: **1768 of 9000, floor
32.22**, and it is a superset. The reason is a fact neither finding has: `set level` does
not always store the double nearest the level you named.

```
set level 32.22 :  c(level) = +1.01c28f5c28f5dX+005      <- one ulp ABOVE
        literal :     32.22 = +1.01c28f5c28f5cX+005
```

So `= c(level)` renders digits the value never had, and one more decimal-heavy level in
nine falls over than the `level()` path suggests. Of the one-decimal levels, 40 are hit
(64.1, 64.4, …, 99.9). `set level 32.22` on its own is of course fine — the round trip is
not. Confirmed end to end: `set level 32.22` + a plain `wooldid y i t ttre, att robust`
(no `level()` anywhere) dies r(198) on the pinned build (`pC_before.log`, case G).

**`set level` is NOT program-scoped** (`pA_scope_scan.log`: a program sets 90, the caller
still sees 90 after two returns). So the `WOOLDID_LEVEL0` machinery is load-bearing and
"delete it, Stata restores the level itself" is *not* available as a fix. I checked this
first precisely because the ado's own `emptycells` comment records a similar assumption as
"verified illusory".

**One further consequence, in neither finding** (`pD_before.log`, section N). When the
restore at `ado:70` fails, `$WOOLDID_LEVEL0` is left set, and the **next** call's shell
consumes it. A later call that fails at the parse stage — before the worker ever reaches
`ado:575` — hits the stale global on its way out and returns **r(198) instead of its own
r(111)**, and skips its own epilogue too, so its `WOOLDID_LASTRC` is not recorded either.
One bad session level thus corrupts the return code and the replay marker of unrelated
later calls.

**The replay gap reproduces, but only from a clean marker.** In `pD_before.log`: good run
(clears `WOOLDID_LASTRC`) → `set level 99.9` → a call that fails r(198) → marker still
empty → bare `wooldid` **replays the pre-failure results, rc 0**. In `pC_before.log` the
same bare `wooldid` refused, because an earlier failure had left the marker at 198 and the
level failure could neither set nor clear it. Both behaviours are the same bug: the
epilogue does not run. Worth knowing when writing the regression test — the sequence
matters, and a test that starts from a dirty marker will pass for the wrong reason.

---

## L1. Render every level as the two-decimal text `set level` parses, and stop letting the
## restore abort the command  — **RECOMMENDED**

### Edit 1 of 3 (required) — `wooldid_validate`, after `ado:5507`

```stata
  * bootstrap() / permute(). Resolve it once here; every consumer sees a number
  if `level' == 0 local level = c(level)
```
becomes
```stata
  * bootstrap() / permute(). Resolve it once here; every consumer sees a number
  if `level' == 0 local level = c(level)
  * ...and canonicalize it to the two-decimal text `set level' parses. Both sources
  * render a double at 16 significant digits: `syntax , Level(real 0)' turns level(99.9)
  * into 99.90000000000001, and `= c(level)' is worse, because `set level' stores some
  * levels one ulp off the nearest double (32.22 -> +1.01c28f5c28f5dX+005), so the text
  * carries digits the value never had. `set level' then refuses it - "level() can have
  * at most two digits after the decimal point" - for 1152 of the 9000 legal two-decimal
  * levels asked for with level(), and 1768 of them taken from the session setting.
  * %9.0g renders every legal level (10.00-99.99, at most four significant digits)
  * exactly, and round-trips through `set level' bit-identically on all 9000 (measured);
  * where the 16-digit text was accepted it is byte-identical to this one, so no level
  * that works today changes value
  local level = strofreal(`level', "%9.0g")
```

`%9.0g` is not a new convention here: `ado:579` already builds `wdlvltxt` with
`strofreal(\`wdlvl',"%9.0g")`, which is why the column names `lb_92p5ci` have always been
right while the `set level` next to them was wrong. A legal level has at most four
significant digits, so `%9.0g` is exact, and `strofreal` returns it unpadded (measured:
`strofreal(99.9,"%9.0g")` = `99.9`, length 4).

Because this sits before the S5 `c_local` blanket (`level` is already in the list), every
downstream consumer — `wooldid_permctl`, `wooldid_bootctl`, `permwire`, `bootwire`,
`permjnt`, `bootjnt`, the cell rewirings, `ptp` — receives the canonical text.

### Edit 2 of 3 (required) — `wooldid_worker`, `ado:574-577`

```stata
  if `level' != 0 & "`wdrepmode'" == "" {
    global WOOLDID_LEVEL0 = c(level)
    set level `level'
  }
```
becomes
```stata
  * the guard compares CANONICAL TEXT, not numbers: an omitted level() was resolved to
  * c(level) in validate, so the set was already a no-op there, and skipping it keeps
  * wooldid from touching the session level (or the restore global) on any run that did
  * not ask for a different one. `= c(level)' renders the double at 16 significant
  * digits, which `set level' refuses, so the stored value is rendered the same way the
  * level text is built two lines below. (the comment above is stale in one detail:
  * replications do NOT inherit this setting - the shell restores the caller's level
  * before the permute()/bootstrap() rep loop runs, and every rep runs at the caller's
  * level. harmless, since a rep is harvested for b and se only)
  if `level' != 0 & "`wdrepmode'" == "" & "`level'" != strofreal(c(level), "%9.0g") {
    global WOOLDID_LEVEL0 = strofreal(c(level), "%9.0g")
    set level `level'
  }
```

Three changes, two of them in one line. The `strofreal` on the **global** is the one that matters (it is
what `ado:70` will feed back to `set level`). The added clause is belt-and-braces with a
real payoff: since validate fills an omitted `level()` from `c(level)`, the `set level` was
always a no-op on those runs, and wooldid now touches neither `c(level)` nor a global on
any run that did not ask for a different level. That closes the whole no-`level()` limb a
second time, independently of the rendering.

The third change is a comment correction while these lines are open. `ado:572-573` says
replications "inherit the outer call's setting and never touch it". They do not: the shell
restores `c(level)` at `ado:69`, **before** the `permute()`/`bootstrap()` rep loop at
`ado:112` onward, so every replication runs at the caller's level, not the requested one.
It is harmless — a replication is harvested for `b` and `se`, and the interval machinery is
handed `level(\`level')` explicitly — but the comment is wrong and the author tracks
doc-vs-code drift.

### Edit 3 of 3 (strongly recommended, structural) — the shell, `ado:68-73`

```stata
  * F054: the worker sets level() for the duration of the call; restore the caller's
  if "$WOOLDID_LEVEL0" != "" {
    set level $WOOLDID_LEVEL0
    macro drop WOOLDID_LEVEL0
  }
  nobreak {
```
becomes
```stata
  nobreak {
    * F054: the worker sets level() for the duration of the call; restore the caller's.
    * captured, and moved inside the nobreak, because an un-captured statement here sits
    * between the worker and the rc-gated epilogue: a restore that errors aborted the
    * command with THAT rc even when the worker had succeeded and posted e(), skipped
    * WOOLDID_LASTRC (so a bare wooldid then replayed pre-failure results) and skipped
    * the epilogue's frame and scratch-estimate cleanup entirely
    if "$WOOLDID_LEVEL0" != "" {
      capture set level $WOOLDID_LEVEL0
      if _rc dis as text "Caution: wooldid could not restore the confidence level in force before this call (level $WOOLDID_LEVEL0); c(level) is now `=c(level)'."
      macro drop WOOLDID_LEVEL0
    }
```

This is the edit that deserves the author's attention even if he dislikes the other two.
Edits 1-2 fix the rendering that fires the restore failure today; edit 3 fixes the reason a
restore failure was able to cost four separate guarantees. `local wdsh_rc = _rc` is read at
`ado:67`, before this block, so the added `capture` cannot disturb it. Moving the restore
inside `nobreak` also means a Ctrl-C arriving between the worker and the epilogue can no
longer strand the session level. Brace balance is unchanged (the block's `}` already
closes at the old `nobreak`'s `}`). The notice is `as text` / "Caution:" to match the
file's convention (`as error` / "Error:" is reserved for fatal messages here), and it did
not fire in any test after the patch — with the rendering fixed there is no level Stata
accepted that it cannot give back.

On the leak wording: `$WOOLDID_CLEANUP` left set was **measured** (`CLEANUP=[fullmodel]`
after the half-working case). A surviving scratch `___wooldidfullmodel` was **not**
reproduced — on the mid-worker failure the flag was still empty, and on the half-working
case the worker had succeeded, so there was nothing to clean. The skipped cleanup is real;
treat the scratch-estimate leak as a consequence that follows from the code path, not one
these logs show.

### Does it move a reported number?

**No — and this is proven over the whole legal domain, not sampled.**

The invariant (`pH_textid.log`): for **every** level the current code accepts — 7232 of the
9000 session levels and 7848 of the 9000 `level()` values — `strofreal(.,"%9.0g")` returns
the **byte-identical text** to today's 16-digit rendering (`text_differs = 0` in both
limbs). Identical text into `set level` gives an identical `c(level)`, hence identical
`e(level)`, quantiles, interval endpoints and column names. The levels where the text
differs are exactly the levels that fail today, and a run that fails has no numbers.

Empirically, on the same data and seeds, pinned vs patched:
- `e(wdidmainresults)` at `%21x` for `level(95)`, `level(90)`, `level(92.5)`: identical.
- an option sweep (`pF_swb.txt` vs `pF_swa.txt`), dumping `e(b)`, `e(wdidmainresults)`,
  `e(gteffects)` and `e(pwsgs)` at `%21x`: **byte-identical** across all four estimands,
  ES + `jointtests` + `gteffects`, `subgroups()`, `controls()` + `level(90)`, `poisson` +
  `level(92.5)`, `bootstrap(standard) reps(20) level(90)`, and a session-90 run with no
  `level()`. Two of that file's ten combos were **mis-specified by me** and refused
  identically on both builds (`unconditionalse` without `robust`; `permute()` without
  `cluster()`), so they contributed no numbers — they were re-run correctly in
  `pJ_jb.txt` vs `pJ_ja.txt`: `permute(shuffle) permreps(20) cluster(i)` at the session
  level and at `level(90)`, and `unconditionalse robust` at the session level and at
  `level(92.5)`, all four **byte-identical**. That is the empirical leg for the `permute()`
  consumers of the canonicalized text (`permctl`'s range check, `permwire`, `permjnt`,
  `permcellsx`). The single difference across both sweep files is the intended one:
  `permute(...) level(99.9)` is r(198) on the pinned build and rc 0 with `lb_99p9ci` on the
  patched one.

So there is **no re-bless**. What changes is *which runs run at all*: 1768 session levels
and 1152 `level()` values move from r(198) to a result. Those results are new numbers in
the sense that they did not exist before, and they are correct — the patched
`level(99.9)` interval has implied quantile 3.32016169, which is `invt(331, .9995)` to
every digit printed, and its columns are `lb_99p9ci`/`ub_99p9ci` as `sthlp:540` promises.

One refusal is relaxed, and it is the only behaviour change on input that fails today:
`level(95.000000001)` now runs at 95 instead of erroring, because `%9.0g` keeps eight
significant digits. `level(95.123)`, `level(5)`, `level(99.995)`, `level(-1)`, `level(100)`
are still refused on both builds (`pG_*.log`). If that bothers the author, L2 restores the
refusal with a better message.

### Cost

Zero on data. Per call: one `strofreal` and one string comparison at parse time, one
`strofreal` and one `capture` in the shell. No loop, no pass over the data, nothing that
scales with N, T or cell count. Nothing here is a guard that gets more expensive on large
datasets.

### Blast radius

- `set level` appears at exactly two sites in the ado (`ado:70`, `ado:576`); both are in
  this patch. `$WOOLDID_LEVEL0` is read at exactly one site (`ado:70`) and written at one
  (`ado:575`).
- `level` as a macro flows through the S5 `c_local` blanket to ~10 subprograms. All of them
  use it numerically or through `%9.0g` for names; none is text-compared to a 16-digit
  form. The `if \`level' != 95` tests at `ado:7376`/`ado:8709` still behave identically
  (`95` renders as `95` either way).
- Neither S5 `syntax` line is touched, so `test/check_validate_sync.py` is unaffected.
- Replications are untouched: `wdrepmode` still short-circuits the `set level`.
- Not covered, deliberately: `level(0)` continues to mean "unspecified" on both builds,
  since `Level(real 0)` uses 0 as its sentinel. Pre-existing, unchanged, arguably a help
  footnote at most.

### Alternatives rejected

1. **Round the level instead of re-rendering it** (`round(\`level',.01)`). Same effect on
   the good cases but silently accepts `level(95.123)` as 95.12 — it destroys the ability
   to refuse genuinely over-precise input, which `%9.0g` preserves.
2. **Pass `level()` to each `margins` call instead of `set level`.** Reverses author ruling
   F054, touches 16+ harvest sites, and does not even fix it: `margins , level()` parses a
   real the same way.
3. **`capture set level` at `ado:576` and carry on at the session level.** Silently gives
   the user a different confidence level than they asked for; the whole point of F054 is
   that `level()` reaches the analytic tables.
4. **A `wdlvlgiven` flag c_local'd from validate** (the cleaner-reading version of edit 2's
   gate). Rejected because it means editing the S5 blanket list, which the file's own
   comments flag as a sync hazard, to buy nothing the text comparison does not already buy.
5. **Delete the `WOOLDID_LEVEL0` machinery and let Stata restore the level.** Measured
   false: `set level` is not program-scoped.

### The test that would prove it correct

Three, none needing an oracle:

1. **Exhaustive text identity** (the strong one, ~20 s, no wooldid needed — `pH_textid.do`).
   For all 9000 legal two-decimal levels, in both limbs: where the current rendering is
   accepted by `set level`, assert the new rendering is the identical string. This is what
   licenses "moves no number" without a re-bless. Result today: 7232 and 7848 accepted,
   0 text differences.
2. **Round-trip closure** (`pA_scope_scan.do`): for all 9000, `set level T`; render;
   `set level 95`; re-`set level` the rendering; assert rc 0 and `c(level)` bit-identical
   to the first. Result: 0 failures, 0 mismatches (the current rendering: 1768 failures).
3. **State invariants around one call** (`pC_ba.do`, cheap enough for the cert suite): with
   the session at 99.9 and `level(95)`, assert `rc == 0`, `e(level) == 95`,
   `c(level) == 99.9` afterwards, `"$WOOLDID_LEVEL0" == ""`, `"$WOOLDID_CLEANUP" == ""`;
   and with the session at 99.9 and no `level()`, assert `rc == 0` and
   `e(level) == c(level)`. Plus the replay leg from a **clean** marker: good run → failing
   run → bare `wooldid` must refuse. The closed-form leg: on a `level(99.9)` run,
   `(ub - b)/se == invt(e(df_r), .9995)`.

### Measured before → after

| case | pinned | patched |
|---|---|---|
| `level(99.9)` | r(198) | rc 0, `e(level)=99.9`, `lb_99p9ci` |
| `level(99.99)` (the help's own bound) | r(198) | rc 0, `e(level)=99.99` |
| `set level 99.9`, no `level()` | r(198) | rc 0, `e(level)=99.9`, `c(level)` still 99.9 |
| `set level 32.22`, no `level()` | r(198) | rc 0, `e(level)=32.22` |
| `set level 99.9` + `level(95)` | **rc 198 with `e(N)=600` posted**; `c(level)` silently left at 95; `$WOOLDID_LEVEL0` leaked | rc 0; `c(level)` back to 99.9; no leak |
| good run → level failure → bare `wooldid` | replays pre-failure results, rc 0 | n/a (the call succeeds); genuine-failure control still refuses |
| later parse failure after a leak | returns r(198), true rc 111 lost | r(111), `LASTRC=[111]` |
| `bootstrap(standard) reps(20) level(99.9)` | r(198) | rc 0, `lb_99p9ci` |
| `permute(shuffle) permreps(20) cluster(i) level(99.9)` | r(198) | rc 0, `lb_99p9ci` |
| option sweep at `%21x` (12 runs incl. permute, bootstrap, poisson, subgroups, ES) | — | byte-identical |

---

## L2. Give wooldid its own words for an illegal `level()` — **OPTIONAL, author's call**

The plain analytic path has **no** `level()` validation at all: the only range checks live
in `wooldid_permctl` (`ado:5862`) and `wooldid_bootctl` (`ado:6262`), so on an ordinary run
`level(5)` and `level(99.995)` are refused by Stata's `set level`, in Stata's words. After
L1 that is harmless and accurate, so this is a polish item, not a fix.

Insert after the canonicalization in edit 1:

```stata
  if `level' < 10 | `level' > 99.99 | abs(`level'*100 - round(`level'*100)) > 1e-9 {
    dis as error "Error: level() must be between 10 and 99.99, with at most two digits after the decimal point."
    error 198
  }
```

The tolerance is measured, not guessed (`pI_tol.log`): over all 9000 legal levels the worst
`|L*100 - round(L*100)|` is **1.819e-12** (through `c(level)`; 9.095e-13 through the
literal, at 64.01), so 1e-9 clears every legal value by three orders of magnitude while
still refusing `level(95.000000001)`, whose deviation is 1e-7. A 1e-6 tolerance would
*not* refuse it — worth knowing if this is written by hand later.

Cost: three floating-point operations at parse time. Moves no number: every value it
refuses is refused today, and it re-refuses the one input L1 alone would start accepting
(`level(95.000000001)` as 95). Gains: the message is wooldid's own and matches `sthlp:540`
word for word, and the `permctl`/`bootctl` copies become redundant. **The edit itself was
not run** — only its tolerance was measured — unlike everything in L1.

Not recommended on its own: without L1 it makes nothing better, because every legal level
passes this check and `set level` still refuses 1768 of them afterwards.

---

## L3. What should NOT be changed

- **No help edit is needed.** `sthlp:540` ("between 10 and 99.99", "a non-integer level
  writes its decimal point as p") describes the behaviour L1 delivers. Documenting the
  refused levels instead of fixing them would be the wrong trade: the set is
  `{L : the 16-digit rendering of L has three or more decimals}`, which is not a sentence
  anyone can act on.
- **Do not add a `capture` around `ado:576`.** A `level()` that genuinely cannot be honoured
  should still fail loudly; L1 leaves that path intact (`level(95.123)` still r(198)).
- **`level(0)` meaning "unspecified"** is a consequence of the `real 0` sentinel and is not
  worth a code change.
