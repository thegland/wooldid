# The .mlib conversion: findings and the decision it needs (2026-08-31)

> ## SUPERSEDING FINDING (same day, after the research below)
>
> **The bug this chunk exists to fix does not affect normally-installed users.** It is an
> artifact of loading wooldid with `do "<path>/wooldid.ado"` — which is what every test
> suite, harness and agent brief in this project does, and what all the SI-17/SI-18
> evidence was gathered under.
>
> Measured, identical commands, two loading methods:
>
> | loading method | gteffects+cleanmatrices, one run | bootstrap after cleanmatrices |
> |---|---|---|
> | `do "…/wooldid.ado"` (our harness) | r3499 | r3499 |
> | adopath (a real install) | **rc 0** | **rc 0** |
>
> Mechanism, isolated in a ten-line reproduction: when a mata function is missing, Stata
> falls back to **re-reading the defining ado from the adopath**, which re-runs its mata
> block. A `do`-loaded program has no such association and cannot recover. The minimal
> case — an ado that clears mata mid-program and then calls its own mata function —
> returns rc 0 when adopath-resolved and r3499 when `do`-loaded.
>
> **Consequences.**
> 1. The `.mlib` would be solving a problem users do not have. Its remaining value is for
>    developers running working copies by path — real, but a much smaller prize than the
>    chunk assumed.
> 2. **SI-17's refusal is a user-facing regression.** `permute() + cleanmatrices`
>    returns rc 0 on the pre-SI-17 build when adopath-loaded (measured). The refusal I
>    shipped removes a combination that worked. The same applies to the
>    `bootstrap() + cleanmatrices` refusal.
> 3. My error, and worth naming: the SI-17/SI-18 evidence was gathered under the dev
>    harness's loading method and generalised to users without checking the other one.
>
> **Revised recommendation: do not ship a library.** Instead reconsider the two refusals
> (they block working user configurations), and if the developer-facing breakage is worth
> addressing at all, address it in the test harness — which could load wooldid via
> `adopath ++` instead of `do`, matching how users actually run it.
>
> Everything below remains accurate as mechanics research and is retained for the record.


Chunk: *"Create a .mlib file containing all mata functions used by wooldid, so that
`clear mata` does not break wooldid's ability to use them … I want a story for [dev vs
shipped] that can't leave the dev build and the shipped build silently different."*

**Status: researched and prototyped, not shipped.** The fix works — I proved it end to
end — but the research turned up three complications that make "ship a prebuilt library"
a decision rather than a formality, and one of them lands directly on the requirement
quoted above. This document is the evidence and the fork.

---

## 1. The fix works

Built `lwooldid.mlib` from wooldid.ado's own mata tail (16 functions, 341 kB) and put it
on the adopath. **Every previously-failing case now completes**, with `wooldid.ado`
itself unchanged:

| case | before | after |
|---|---|---|
| `gteffects` + `cleanmatrices`, one run | r3499 | **0** |
| `adversarial` + `cleanmatrices`, one run | r3499 | **0** |
| `bootstrap` after a `cleanmatrices` run | r3499 | **0** |
| `permute` after a `cleanmatrices` run | r3499 | **0** |
| plain run, no mata options (control) | 0 | 0 |

Mechanics confirmed: library functions survive `clear mata` indefinitely and reload
transparently; a net-installed library in `PLUS/l/` is found by a fresh session with no
index call and no ado-load-time code; the packaging change is literally one line
(`f lwooldid.mlib`); the `l` prefix is mandatory and enforced at build time; members
cannot be replaced individually, so each release must recreate the library.

**The shape has to be "inline *and* library", not "library only."** Deleting the mata
from the ado would break the assessment suite's historical build staging (`git show
<sha>:wooldidinstall/wooldid.ado` yields one self-contained file today), both documented
two-build comparison methods, ~25 test suites that load the package by explicit path, and
the colleagues' external validation harness, which installs a single supplied file.
Inline definitions win while they are in memory, so the library is pure insurance that
only engages after `clear mata`.

---

## 2. The three complications

### 2a. A library built here is silently ignored by older Stata — and wooldid declares 16.0

Measured: an `.mlib` records, per function, the Stata version that compiled it. A Stata
older than that **ignores the member and errors r(3499) "not found"** — the same return
code as a genuinely missing function — and the explanatory note is **completely
suppressed under `quietly`**, which is how an ado calls its own mata. There is no
`version()` option on `mata mlib create`/`add`/`lmbuild` (r198, and absent from the help
syntax), and `version 16` in any position — source file, command prefix, whole build file
— changes only the compile-time field, never the built-by stamp. Twelve builds under
`version 9` through `19` all produced identical stamps.

Under the inline+library shape this is **not harmful**: Stata 16 users keep working
normally off the inline definitions and simply do not get the insurance. But it does mean
a library built on this machine delivers the fix to Stata 19 users only, and the
package's own `.pkg` says `Requires: STATA 16`.

The ecosystem agrees this is the sharp edge. Every installed package that ships a
prebuilt library built it at roughly its own declared minimum (gtools 13.1, boottest
13.1, honestdid 14.2, pretrends 15.1). **reghdfe and ftools ship no prebuilt library at
all** — they carry mata source and recompile whenever `c(stata_version)` changes.

### 2b. A stale library shadows silently — the exact failure the chunk asked to prevent

Measured, and this is the important one. Inline definitions win while they are in memory,
so a wrong library is invisible during normal use. But **the moment `clear mata` runs,
resolution falls to whatever `lwooldid.mlib` the adopath holds — with rc 0 and no
message.** A developer running a working copy by path, with an older released wooldid
installed in `PLUS`, would after `cleanmatrices` silently execute the *old* release's
mata against the *new* ado. That is precisely the "dev build and shipped build silently
different" outcome the chunk text set out to avoid, and shipping a prebuilt binary is
what creates the opportunity for it.

It is mitigable — a version-stamp function inside the library that the ado checks after a
clear, plus a build-time hash guard in `check_validate_sync.py` — but the mitigation is
part of the cost, not an afterthought.

### 2c. A binary in the repo

The library is a 341 kB opaque artifact. Committing it makes review diffs uninformative
and re-stages a stale copy for every historical sha. Not disqualifying, but it argues for
generating rather than committing wherever possible.

---

## 3. The three viable designs

**A — Ship a prebuilt library (the chunk as written).** One `.pkg` line; a release build
step; generated from the ado's own mata block so it cannot drift by accident; plus a
version-stamp guard for 2b and a sync check in the cert component. *Costs:* a binary in
git; the fix reaches only users on the build Stata or newer, so the build should happen on
the oldest supported Stata (not available on this machine); the 2b hazard has to be
actively guarded rather than assumed away.

**B — Build the library at runtime from the ado's own definitions (the ftools pattern,
adapted).** When `cleanmatrices` is about to clear mata, first write a library from the
functions already in memory, then clear. Afterwards, calls resolve from a library that was
built by *this* Stata, from *this* ado, moments earlier. *This dissolves all three
complications at once:* no version gate (built by the running Stata), no staleness
possible (built from the running source), no binary in git, no packaging change. *Costs:*
needs a writable directory with a fallback chain, adds a few seconds to the first
`cleanmatrices` run, and puts new machinery inside a shared option.

**C — Re-source a mata-only file after clearing (SI-18 option 2).** Simplest and fully
portable, but it moves the mata source into a second shipped file, which the single-file
staging workflow would then have to carry.

**Recommendation: B**, with A as a later addition if prebuilt distribution is wanted for
other reasons. B is the only one that satisfies the chunk's own stated requirement — a
story that *cannot* leave dev and shipped silently different — by construction rather than
by discipline. A satisfies it only if the guards are built and maintained.

---

## 4. What is already built and ready

- `build_mlib.py` — generates the library from wooldid.ado's own mata tail, asserts its
  extraction boundary (fails loudly if a program ever appears below the mata section), and
  verifies all 16 functions landed. Runs clean; produces the 341 kB library above.
- `test_mlib.do` — discovery, contents, survival across `clear mata`, the four cure cases,
  a plain-run control, and a **bisectability check** that hides the library and confirms
  the cure cases fail again, so they cannot pass vacuously.
Both are in the session scratchpad, unshipped, pending the choice above.

## 5. If a prebuilt library is chosen anyway

Build on the oldest Stata that can be reached, not this one. Ship it as
`wooldidinstall/lwooldid.mlib` with `f lwooldid.mlib` added to the `.pkg` and the
distribution date bumped. Add the version-stamp check for 2b and a mata-region hash to
`check_validate_sync.py` so source/library drift fails the cert component. Keep the
`permute`/`bootstrap` + `cleanmatrices` refusals: on any Stata older than the build
machine the library is inert, so the refusals remain the portable guarantee — the library
is the belt to their braces, not a replacement for them.

Note that `gteffects` and `adversarial` have **no** refusal today and fail in the same
run, so under design A they would remain broken for old-Stata users. Design B fixes them
for everyone.
