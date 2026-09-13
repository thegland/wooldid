# Reproductions of record — 2026-09 full assessment run, Stage B proof sections

This directory holds the fourteen reproduction do-files (and the helpers they need) that the six
Stage B proof sections — `legacy`, `si20`, `contsg`, `adversarial`, `inference`, `intended` — cite
as their evidence. They were drafted in a session scratchpad with absolute paths baked in; every
one of those paths is now a **positional argument** (do-files) or an **environment variable**
(shell and python helpers), so the whole set runs from a clean checkout on any machine.

Nothing here was re-run during the rewrite. The transformation was exact string replacement,
reviewed statically; behaviour is unchanged. **Verify a leg's first run against the section's
quoted decisive lines before trusting it.**

---

## 1. Why the cited dumps are a *regenerable* artifact, not a stored one

Every section cites `test/suite/out/run_old/<config>.csv` and `run_new/<config>.csv` — the %21x
per-config dumps of the dual-build run — as the evidence of record, and several tie their own
reproduction dumps to them byte-for-byte at `--rtol 0`.

`test/suite/out/` is **gitignored** (`test/.gitignore` line 5: `suite/out/`). Nothing in the
repository preserves any run's dumps, by design: they are outputs, and `run_full.sh`'s `do_oldnew`
begins by `rm -rf`-ing `run_old`/`run_new` on every full run. The 2026-09-11 dumps the sections
were written against were overwritten in place by the 2026-09-13 assessment run and no copy of
them exists.

That is not a loss of evidence, because the dumps are **deterministic functions of three tracked
inputs**: the two builds of `wooldidinstall/wooldid.ado`, the suite definition
(`configs.csv`, `configs_ext.csv`, `preps.json`, `mkrun.py`, the fixtures under `test/data` and
`test/suite/data`), and the fixed seeds `mkrun.py` writes into every shard. All three are in git.
Section 2 gives the command; section 3 states precisely what has and has not moved since 09-11.

## 2. How to regenerate the dumps the sections cite

```sh
cd <repo>/test/suite
./run_full.sh oldnew          # ~25 min at STATA_JOBS=6; writes out/run_old and out/run_new
```

`do_oldnew` stages both sides from git itself, never from the working tree:

| leg | staged from | dump provenance line |
|---|---|---|
| old | `git show main:wooldidinstall/wooldid.ado` | `prov\|old\|<short sha of main>` |
| new | `git show HEAD:wooldidinstall/wooldid.ado` | `prov\|new\|<short sha of HEAD>` |

It writes `out/run_old/<config>.csv` and `out/run_new/<config>.csv` (one file per config, every
`e()` scalar / macro / matrix cell in `%21x`, plus the `rc|`, `varlist|`, `files|` and `matnames|`
lines), the shard do-files under `out/shards_<side>/`, and the judged comparison at
`out/oldnew_compare.txt`. `STATA_BIN` defaults to
`/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp` — the direct binary, because the `PATH`
symlink `stata-mp` is broken on this machine; export `STATA_BIN` elsewhere.

`./run_full.sh compare` re-judges dumps that already exist without re-running Stata, which is the
cheap way to re-tie a claim.

## 3. Which builds the two sides are, and what has moved since 2026-09-11

**old** is published wooldid v1.3 = `git main` = **`301b706`** (3,450 lines). It has not moved and
is the old side of every run cited anywhere in the report.

**new** is git `HEAD` at the moment of the run, so the label depends on when you run it:

| run | new-side sha | ado content |
|---|---|---|
| 2026-09-11 (the run the sections were written against) | `14c30d3` | — |
| the sections' pin | `08ec0b4` | **byte-identical to `14c30d3` except** 9 comment lines and a shortened `dis as text` permutation-caution string (`git diff 14c30d3 08ec0b4 -- wooldidinstall/wooldid.ado`: 10 insertions, 1 deletion, no estimation code) |
| 2026-09-13 assessment run | `6a05127` | **byte-identical to `08ec0b4`** (`git diff 08ec0b4 6a05127 -- wooldidinstall/wooldid.ado` is empty) |
| current `HEAD` `21db885` | `21db885` | **byte-identical to `08ec0b4`** |

So no stored `e()` object moves between the 09-11 build and the build the 09-13 run used: **every
value cell the sections quote is reproduced by the 09-13 run.** That is the claim on which the
regenerability of section 1 rests, and it is checkable in one command.

**What did move, and is a caveat the reader is owed.** Commit `4280f3b` (2026-09-13 11:23, "Update
the assessment suite for the 2026-09 inference rework") changed the *harness*, hours after the
sections were drafted. Measured against `26ce421` (the suite state the sections read):

- `configs.csv` 154 → 191 rows: **37 added, 0 removed, and not one option string on a common row
  changed.** The only edits to pre-existing rows are the `note` field (6 rows: `adv_r2001`,
  `boot_tdf`, `cov_ccc_absorb`, `fcol_gate`, `fcol_sg`, `trig_jtclus`) and the `slow` tier column
  (3 rows: `boot_bca`, `boot_obs`, `boot_tstat`).
- `configs_ext.csv` 358 → 362 rows: **4 added, 0 removed**; on common rows only `rc_new`
  (`ext_contract10_nonpositive_time`), `msg_old` and `slow` (`ext_contract10_rcs_factor_fe`).
- `preps.json`: additions only, no prep removed or rewritten.
- `mkrun.py`: `dumpres`'s `svpre` becomes a *list* of file prefixes (so `save()`, `permsave()` and
  `bootsave()` output is pinned, not just `saveplots()`), and an `e()` matrix the dump cannot read
  now writes `matnames|<name>|UNREADABLE|rc=<rc>` instead of vanishing silently. Shard tiering and
  timeouts changed.

Consequence, stated exactly: **every scalar, macro and matrix VALUE cell for the 507 configs the
09-11 dumps cover is expected identical at `--rtol 0`.** A regenerated dump may nevertheless differ
from a 09-11 dump in its `files|` line (on configs using `save`/`permsave`/`bootsave`) or gain a
`matnames|…|UNREADABLE` line. Those are additive harness records, not estimates. The honest form of
the sections' "bit-identical to `run_old` at rtol 0" claims is therefore: *re-tie them with
`test/compare.py --rtol 0` on the regenerated dumps and record the result.* **That re-tie has not
been performed** — the 09-13 run was still writing `run_old` when this directory was created — and
it is listed as owed in `LEDGER_DELTA.md`.

## 4. The verbose re-judge (`compare_verbose_0911.txt`) is likewise a recipe, not a file

Four sections cite `compare_verbose_0911.txt` — the per-config verbose re-judge of the 09-11 dumps
— and two of them spell it `scratchpad/compare_verbose_0911.txt`, which reads as repo-relative. It
is not: there is no `scratchpad/` in the repository and the file lives only in a session tmp tree.
It is 1.1 MB of derived text and is not copied here.

It is a pure function of the dumps, the ledger and the comparer, and `compare_oldnew.py` takes
both of the latter as arguments, so any historical judgement can be rebuilt exactly:

```sh
cd <repo>/test/suite
python3 compare_oldnew.py out/run_old out/run_new --rtol 1e-10 \
    --ledger ledger.csv --configs configs.csv --configs configs_ext.csv --verbose
```

To rebuild the *26ce421* judgement the sections quote (399 EXPLAINED), point `--ledger` and
`--configs` at that commit's files instead, extracted to a scratch path:

```sh
for f in ledger.csv configs.csv configs_ext.csv; do
  git -C <repo> show 26ce421:test/suite/$f > /tmp/26ce421_$f
done
python3 compare_oldnew.py out/run_old out/run_new --rtol 1e-10 \
    --ledger /tmp/26ce421_ledger.csv --configs /tmp/26ce421_configs.csv \
    --configs /tmp/26ce421_configs_ext.csv --verbose
```

Two judgements of the same 09-11 dumps are quoted across the report and they differ **only** in
which ledger and comparer produced them, not in any measurement. Which one the report judges
against, and the full row-level delta between them, is settled in the appendix **`LEDGER_DELTA.md`**
carried alongside the six sections.

## 5. Running the reproductions

**One-time setup per section.** Each do-file expects a *scratch workdir* containing the extracted
builds as subdirectories. Extract them from git — never copy a working tree:

```sh
W=/tmp/stageB/legacy            # any path; one per section
mkdir -p "$W"/old "$W"/new
git -C <repo> show 301b706:wooldidinstall/wooldid.ado > "$W"/old/wooldid.ado    # old = v1.3
git -C <repo> show 08ec0b4:wooldidinstall/wooldid.ado > "$W"/new/wooldid.ado    # new
```

`adversarial` needs a third build, the BS-12 isolation tree `7e30e85` (= `25727a6^`):

```sh
git -C <repo> show 7e30e85:wooldidinstall/wooldid.ado > "$W"/mid/wooldid.ado
```

`si20` runs **instrumented** copies. Build them with the helper shipped here, which inserts the
read-only `REPRO_SI20` probe (52 inserted lines, 0 deletions per build) and is proven inert by
`check_inert.py`:

```sh
#                           <src>              <dst>                <expect_f2> <expect_f1>
python3 si20/instrument.py "$W"/old/wooldid.ado "$W"/old_i/wooldid.ado 1773 2148
python3 si20/instrument.py "$W"/new/wooldid.ado "$W"/new_i/wooldid.ado 2170 2662
```

The two integers are asserted probe-site line numbers, **in the order the script reads them**
(`src, dst, expect_f2, expect_f1 = sys.argv[1:5]`): **`expect_f2` is the binary-`att` site**
(old 1773, new 2170) and **`expect_f1` the continuous-`ame` site** (old 2148, new 2662). The
script `assert`s each and stops if the build does not match, so a wrong pair fails loudly rather
than instrumenting the wrong line.

`contsg` additionally needs its shared helper inside the workdir, because both of its do-files
`include` it from there (`contsg/repro_run.sh` does this for you):

```sh
cp contsg/repro_helpers.doh "$W"/
```

**Paths must not contain spaces** — the do-files pass them through Stata locals.

### Argument contract, file by file

| file | arguments, in order |
|---|---|
| `legacy/repro_A1A2.do` | `<old\|new> <repo-root> <workdir>` |
| `legacy/repro_A12TTB.do` | `<old\|new> <repo-root> <workdir>` |
| `si20/repro_si20.do` | `<old\|new> <repo-root> <workdir>` (sources `<workdir>/<build>_i/wooldid.ado`) |
| `inference/repro_inference.do` | `<old\|new> <repo-root> <workdir>` (run from `<workdir>/<build>`) |
| `intended/repro_intended.do` | `<old\|new> <repo-root> <workdir>` (run from `<workdir>/<build>`) |
| `intended/repro_svpsg.do` | `<old\|new> <repo-root> <workdir>` (run from `<workdir>/<build>`) |
| `contsg/repro_old.do` | `<repo-root> <workdir>` (build fixed by filename) |
| `contsg/repro_new.do` | `<repo-root> <workdir>` (build fixed by filename) |
| `adversarial/old/repro_adv_old.do` | `<repo-root> <workdir>` |
| `adversarial/mid/repro_adv_mid.do` | `<repo-root> <workdir>` |
| `adversarial/new/repro_adv_new.do` | `<repo-root> <workdir>` |
| `adversarial/repro_template.do` | `<repo-root> <workdir>`, **after** substituting `__BUILD__` |
| `legacy/harness_template.do` | not run directly — the pattern the others follow; `$REPO` in its comments |
| `contsg/harness_template.do` | not run directly — same template |

Every do-file now refuses with `r(198)` and a one-line usage message if an argument is missing,
and each carries its exact invocation in a header comment. Each prints `REPRO_<id>_DONE` sentinels;
a leg that does not print all of them did not finish.

**Shell drivers** take their roots from the environment and pass them on:

```sh
REPO=<repo> WD=$W zsh legacy/run_repro.sh              # optionally: ... run_repro.sh A1A2
REPO=<repo> W=$W  zsh contsg/repro_run.sh
```

**Python helpers** read `WOOLDID_REPO` (checkout root) and, where they need it, `WOOLDID_SCRATCH`
(the section workdir); each exits with a named message if the variable is unset.
`legacy/repro_cross.py`, `contsg/repro_summarize.py`, `contsg/gen_repro.py`, `contsg/dumpparse.py`,
`si20/check_inert.py`, `si20/summarize.py`, `si20/named_cells.py`. `si20/instrument.py` and
`adversarial/hexdec.py` are verbatim copies — they never held a path.

Both drivers honour the project's hard cap of 8 concurrent `stata-mp` instances. **Do not start a
leg while a full assessment run is live.**

### What is deliberately NOT here

The bulk artifacts — per-config `%21x` dumps, Stata logs, the copied `wooldid.ado` builds, the
1.1 MB verbose re-judge, and the summarised `.txt` outputs — are 15 MB and every one of them is
regenerable from the commands above. Copying them would have stored a snapshot whose provenance
decays; the recipe does not.

---

*This directory is untracked. It was written while a full assessment run held `test/` and
`wooldidinstall/`; no Stata was run and nothing under either tree was modified. See
`EDITS.md` alongside the sections for the per-edit log.*
