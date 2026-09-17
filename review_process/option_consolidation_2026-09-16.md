# Ticket: consolidate fastcollapse / fullprecisionslow into nofastcollapse

**Status 2026-09-16: COMPLETE and shipped (`38f0362`), after the T1/T2/T3 merge cleared the block
on the two syntax lines.** `fastcollapse` is removed, the silent ineligibility skip is promoted to
a caution, and `fullprecisionslow` is a clean break with no alias, all as ruled. The option surface
check now lists `nofastcollapse` where `fastcollapse` and `fullprecisionslow` used to appear, and
the full suite of 2026-09-16 passed every component with it in place.

This header previously read "OPEN, ruled, NOT started" and was stale; the work had shipped the
same day. Corrected 2026-09-16.

## What the audit found before the ruling

The proposal as raised was "fastcollapse is always on by default, so fold both options into one
`nofastcollapse`". Half of that is right. Read against `ado:5096-5117`:

```stata
if "`fastcollapse'" != "" {
    if "`fcbad'" != "" { ...enumerated diagnostic... ; error 198 }
}
local fcon 0
local fcskip `fcbad'
if "`fullprecisionslow'" == "" & "`fcbad'" == "" local fcon 1
```

`fastcollapse` does **not** appear in the computation of `fcon`, the realized collapse decision.
On an eligible request it is a **no-op**. Its only effect is to turn the silent skip of an
INELIGIBLE request into an error. It is an assertion, not an on-switch, exactly as
`sthlp:307` says: "this option is never needed to obtain it: specifying it asks wooldid to
{bf:insist} on the collapse".

So the surface carries three states on two axes:

| state | how | axis |
|---|---|---|
| collapse on, silent skip when ineligible | default | - |
| collapse on, ERROR when ineligible | `fastcollapse` | strictness |
| collapse off | `fullprecisionslow` | on/off |

`fullprecisionslow` is the off switch and genuinely is the proposed `nofastcollapse`. Naming the
pair `fastcollapse`/`nofastcollapse` would have been worse than the status quo: they read as
opposites and are not.

Suite evidence that the assertion is nearly unexercised: of six configs naming `fastcollapse`,
five (`fcol_base`, `fcol_fn`, `fcol_sg`, `fcol_pois`) sit on eligible requests where the flag
does nothing and would behave identically without it. Only `fcol_gate` (rc_new 198) reaches the
insist path.

## The ruling (author, 2026-09-16)

1. **One option.** `nofastcollapse` is the only option. `fastcollapse` is removed outright.
2. **The skip stops being silent.** When a request is ineligible and the user did NOT ask for
   `nofastcollapse`, print a caution naming the reason. The assertion's whole value was "tell me
   when the collapse did not engage", and a caution delivers that without an option. The existing
   error text is an unusually good diagnostic - it enumerates every incompatible option with
   measured magnitudes - and moves into the caution rather than being discarded.
3. **Clean break, no alias.** `fullprecisionslow` is not kept as a hidden spelling. An archived
   do-file naming it will error with r(198), and that is the intended, honest signal.

## Implementation notes, gathered while auditing

- **Both syntax lines**, `ado:546` (worker) and `ado:4547` (validate), must change together and
  stay byte-identical; `test/check_validate_sync.py` enforces this and runs inside cert.
- **Verify the declaration form before writing it.** Stata's `syntax` has native `no`-prefix
  handling (`noLOG` yields a local `log` holding "nolog"), which would collide confusingly with
  the existing local named `fastcollapse`. A plain option declared `NOFASTcollapse`, yielding a
  local `nofastcollapse`, is likely cleaner. Test both forms rather than assuming which Stata
  accepts and what it names the macro.
- `fcon`'s guard becomes `if "`nofastcollapse'" == "" & "`fcbad'" == ""`.
- Remove the contradiction screen at `ado:5096` (it has nothing left to contradict).
- **Replication forwarding**: `ado:140` and `ado:333` forward `fastcollapse` into permute and
  bootstrap replication command lines. Those go; confirm the off-switch is forwarded in its place
  (`fullprecisionslow` is documented as forwarded - verify the new name inherits that).
- The c_local blanket list at `ado:5480` names both options and must be updated.
- Several entries say an explicit `fastcollapse` alongside them is an error (`unconditionalse`
  at `sthlp:268`, `customvce()` at `:273`). Those screens and their prose go away with the option.

## Suite work this forces

- `fcol_gate` currently pins rc 198 on the insist path. That path is gone; re-point it at the new
  caution, so the disclosure has a host rather than losing the only config that covered this.
- `fcol_base`, `fcol_fn`, `fcol_sg`, `fcol_pois` drop a no-op flag. Keep the configs - they cover
  the collapse on the robust, fastnose, subgroup and poisson paths - but they no longer differ
  from the default path, so check they are not now duplicates of existing rows.
- `fcol_slow` takes the new spelling.
- COVERAGE: `coverage_check.py` asserts the option surface, so `nofastcollapse` needs at least one
  config and the two retired names must leave the assertion cleanly.
- The oldnew deck WILL move on these six configs if either option name reaches `e()`; check
  whether the `wdv_ret` blanket posts them before predicting the diff.

## Help work this forces

One entry replaces two. Cross-references to update, all measured 2026-09-16: `sthlp:79` (the
synopt option list), `:268`, `:273`, `:299-310` (the `fastcollapse` entry itself), `:322`, `:344`,
`:353`, `:389`, `:411`, `:414-418` (the `fullprecisionslow` entry itself), `:416`, `:448`,
`:1270`, `:1456`. Register measured against `2a2da17` as usual.
