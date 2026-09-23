# TICKET — make a failed bootsave()/permsave() write a caution, not a fatal error

**Opened 2026-09-21 on the author's instruction,** left as a reviewer note on the message review page:
*"Please explore the feasibility of converting this from an error to a caution, so that if the rest of
the run succeeded, users can still look at the results."*

**Status: BUILT.** The conversion itself shipped in `ed6a0da` (2026-09-22) for BOTH savers, on the
author's page ruling: the write failure is a caution, the success line sits in an else arm,
permsave's content cautions are gated onto success, and the return code is posted as
`e(boot_save_rc)` / `e(perm_save_rc)` so a quiet script can still detect the loss (s5's residual
cost, accepted as posted-not-mitigated). The test fallout in s8 followed on 2026-09-22: the bare
`use` reads after a saver run in `test_bootstrap_gates`, `test_permute3` and the Monte Carlo
harness now confirm the file first, and a forced-failure gate exists for each saver (read-only
target, `replace`, rc 0, the caution text, the posted return code, and a bare `wooldid` replay).
The help and guide entries for both options carry the failure clause.

**Still for the author (s7, undecided, carried to the ruling document):** Q2, the two in-program
pre-flight duplicates (the directory-part r(198) and file-exists r(602) checks each saver re-runs
after the loop) are still fatal, so a file that appears mid-run still throws the finished run
away; Q3, the caution prints above the deferred display. Line numbers below are as of `a9f9997`.

## The finding that raises the stakes

The author's motive understates what the refusal costs today. Under `bootstrap()` and `permute()` the
display is **deferred**: the worker skips the tables and the shell prints them after the post-loop
rewiring. The save call sits above that.

    :410  global WOOLDID_LASTRC postloop        <- replay marker set
    :462  if "`bootsave'" != "" wooldid_bootsave , ...     <- BARE, no cap
    :466  wooldid_display , part(1)             <- the tables
    :537  global WOOLDID_LASTRC                 <- marker cleared

So a failed write does three things, not one:

1. **No results are printed at all.** Not once. The abort happens three lines above the first table.
2. **Replay is then refused.** The marker is still `postloop`, so bare `wooldid` answers *"the last
   wooldid call failed; refit wooldid before replaying."*
3. **Three further post-loop steps are skipped**: `e(boot_rngstate)`, the cell-surface rewiring, and
   the deferred plots.

The coefficient inference in `e()` is complete and correct at that moment, and it survives the error.
It is simply unreachable by any route the user has. Converting the branch to a caution dissolves all
three by itself, because control then falls through to the display and reaches the marker clear. **No
separate marker change is needed.**

## Feasibility: the case against does not survive

An adversarial pass found no inconsistent state left by a failed save.

- **No stranded frame.** `frame X { }` restores the caller's frame even when the block exits by error
  under `capture`, so the `cap frame drop` on the next line succeeds. Probed in Stata 19.
- **User data and working directory untouched.** The build and the save happen entirely inside a
  tempname frame, and neither program contains a `cd`.
- **Nothing downstream reads the file.** Lines 463 to 537 contain no reference to the option, the
  `.dta`, or the program's locals.
- **No half-written file.** On a full volume the save left nothing at all, in both the new-name and
  overwrite cases. Measured on one filesystem, so write it as measured rather than as a guarantee.

## House precedent is squarely for it

Twelve lines below the block in question, the analytic cfx-data write already does exactly this, and
states the principle the author is invoking:

> *"cap, but not silent: a failure to write the data file must not lose an estimation that has already
> succeeded and printed ... and equally must not be swallowed."*

The precedent also carries its own discriminator: a failure that loses only a **data file** gets a
caution, while one that also loses a **drawn figure** does not. A save failure is the first kind.

Second precedent, inside the very program to be changed: `wooldid_permsave` already announces its
assignment-block failures as cautions rather than errors, which was author ruling D-10.

## The residual cost, to be accepted rather than mitigated

**A stale file can be read as this run's draws.** Reproduced: with a read-only target, `save, replace`
fails r(608) and leaves the previous file **byte-intact**. Under a caution, a script that runs
`bootsave(d) replace` quietly and then `use d.dta` silently reads the *previous* run's draws. Today it
dies instead.

No wording fixes this, because `dis as text` is suppressed by `quietly` while a nonzero return code is
not. Two partial mitigations, both worth doing:

- Word the caution to say the file may be one that **predates this run**, name the file, and print the
  return code verbatim.
- Post a machine-readable marker so a scripted user can detect it. `wooldid_bootsave` is not eclass,
  so the code must return by `c_local` and be posted at the call site, the way `boot_reid` already is.

## Implementation traps

1. **Do not touch the `_rc` capture ordering.** `local wdbsrc = _rc` must stay on the line immediately
   after the `cap noi frame` block and `cap frame drop` must stay below it. The in-code comment says
   swapping them turns a failed save into rc 0 and prints the success message for a file that was
   never created.
2. **The success line must move into an `else` arm.** `wdbsn` is assigned inside the frame block, so on
   the failure path it is empty and the line prints *"-1 replication rows"* for a nonexistent file.
   Verified.
3. **Re-gate permsave's two content cautions on the save succeeding.** They are set before the save, so
   left alone they would describe a file that was never written.
4. **Do not reuse the cfx precedent's stated cause.** Its "a file of that name exists and replace was
   not specified" case is pre-checked here. The real causes are filesystem (read-only, full, locked) or
   a wooldid layout defect.

## Scope decisions for the author

**Q1. Both programs, or one?** `wooldid_bootsave` and `wooldid_permsave` are byte-identical in shape
and the source comment says so. Converting one leaves neighbouring options with opposite failure
semantics and no test would notice. **Recommend both.**

**Q2. What about the two in-program pre-flight duplicates?** Each program re-runs the directory-part
check (r198) and the file-exists check (r602) *after* the loop. The real pre-flight lives in
`wooldid_bootctl` / `wooldid_permctl` and runs before any estimation; those should stay errors. The
duplicates are near-unreachable, but the 602 leg survives a file appearing mid-run, and as written it
throws away a finished run exactly like the write failure does. Convert, delete as redundant, or leave.

**Q3. Where should the caution print?** At the call site it prints **above** the entire deferred
display, so it scrolls off the top. The block's other cautions print after the last table. Moving the
call below the display would read better, but first check that nothing between the two positions
reposts the accumulators the call reads.

## Test fallout

**Zero existing assertions break.** Nothing in the suite provokes this path, nothing reads the message,
and every asserted `bootsave()`/`permsave()` refusal targets the pre-loop checks, which this does not
touch. Two tests exist specifically to keep those refusals ahead of the loop.

The cost is prospective and needs work in the same commit:

- **Eleven downstream reads depend on the return code**, nine of them a bare `use` inside `preserve`
  with no `cap confirm file`. Today a write failure gives a named FAIL; afterwards the do-file dies at
  r(601) mid-run and every later gate in that file is lost. Convert these to the `cap confirm file`
  plus FAIL pattern that three gates already use.
- **Two harnesses use the return code as a stop gate.** The Monte Carlo one is the worse case: it would
  proceed to analyse a file that was never written.
- **There is no forced-failure gate anywhere.** Add one for both programs: run once, make the file
  read-only, run again with `replace` and a different seed, then assert rc 0, the caution text, the
  posted return-code marker, an empty replay marker, and that bare `wooldid` replays.

## Documentation

The help file documents both options as simply writing the file and says nothing about failure. Three
passages need a clause. Keep the register test in mind: new help prose is measured against `2a2da17`.

## Verification plan

1. Convert both sites, keeping the `_rc` ordering and moving the success line into an `else`.
2. Post the return code at both call sites; re-gate permsave's content cautions.
3. Rewrite the eleven downstream reads to the guarded pattern; add the forced-failure gate.
4. `check_validate_sync.py`, then the bootstrap and permute suites.
5. Full battery, expecting **IDENTICAL**: no passing run changes behaviour.
