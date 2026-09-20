# SKEPTIC notes — hv-permsave-comma-refute

Role: refute (or confirm) three findings from agent `hv-permsave-comma`.
Rule followed: every verdict rests on MY OWN do-file, MY OWN data, run in MY OWN
subdirectory. I read their do-files only to know what to attack; I never ran them.

Build asserted in every run: `which wooldid` -> .../scratchpad/bluesky/build/wooldid.ado

## Files
- `r1_permsave_comma.do` / `.log`, workdir `w1/` — finding 01
- `r2_dots.do` / `.log`, workdir `w2/` — finding 03 (+ two new dots() sentinel probes)
- `r3_bootW.do` / `.log`, workdir `w3/` — finding 02
- `r4_permsave_contents.do` / `.log`, workdir `w4/` — NEW: is the permsave() FILE correct?
- `r5_strata_both.do` / `.log`, workdir `w5/` — NEW: permsave under permstrata();
  permute()+bootstrap() given together

## STATUS

### 01 permsave(f, replace) — **CONFIRMED** on my own repro (r1)
- `permsave(pd, replace)` -> rc 0, prints the success line, writes `pd, replace.dta`
  (80,453 bytes, 21 obs, `_rep main_att onset_1 ... onset_90`) into the cwd.
- `bootsave(bd, replace)` on the SAME data -> r(198)
  "Error: bootsave() takes a plain file name and no comma; replace is a separate option
   of wooldid, so write bootsave(bd  replace) replace." No file written.
- Rerunning the permsave line verbatim -> r(602) "permsave(pd, replace) - pd, replace.dta
  already exists in the working directory; add the replace option to overwrite it."
  (the user typed `replace`).
- `permsave(pd) replace` (the correct form) -> rc 0, writes `pd.dta`.
Mechanism as claimed: `wooldid_permctl` ado:5868-5880 has the .dta/"/"/"\"/empty tests
but no `strpos(...,",")`; `wooldid_bootctl` ado:6279-6298 has all four, the comma one
dated RF-6 2026-09-13. The ado's own comment at 5866 ("the name rules are bootsave()'s")
is falsified by this.
Severity: keep MEDIUM. No number is wrong and no result is lost, but the author already
judged this exact slip worth a dedicated refusal on the twin, and the 602 message on the
retry is actively misleading.

### 03 dots(-5) — **CONFIRMED** on my own repro (r2)
- `bootstrap(standard) reps(20) dots(-5)` -> r(198) "dots() must be a nonnegative integer".
- `permute(shuffle) permreps(20) dots(-5)` -> rc 0, and the run prints NO progress line,
  i.e. it behaves as dots(0). Compare: dots(0) prints nothing, dots(5) prints
  "Permutation replications (20): . (~14 s projected for 20 replications)... done".
- Control: the rseed() range check IS hoisted — `permute(shuffle) rseed(-3)` -> r(198).
Severity: keep LOW.

### 02 boottstat W matrices — see r3 (running)

## NEW (mine, not seeded)
- `dots(-1)` collides with the internal "not given" sentinel:
  * `bootstrap(standard) reps(20) dots(-1)` -> rc 0 AND the progress display PRINTS,
    i.e. the explicit request is silently reinterpreted as the default dots(1).
  * `dots(-1)` with NO engine -> rc 0 (the stray-companion screen cannot see it), while
    `dots(0)` and `dots(7)` with no engine -> r(198) "dots() require(s) the bootstrap()
    option". So exactly one value of dots() escapes the stray screen.
  Cosmetic; reported LOW.

## FINAL STATE (all probes done)
- 01 CONFIRMED (r1), 03 CONFIRMED (r2), 02 PARTIAL (r7: absences real, mechanism corrected).
- r3_bootW.do (heavy arm of 02) hit its 15-min deadline and was killed; r7_W_lean.do is the
  lean replacement and settles the same question.
- NEW N1 (high): bootsave() under boottstat writes the observed COEFFICIENT, not the observed
  Wald statistic, into _rep=0 of jointW_* and advW. r8 found it, r9 is the clean repro.
- NEW N2 (low): dots(-1) sentinel collision.
- Negative results: permsave file contents correct (r4); permstrata onset block correct (r5);
  name-rule surface otherwise symmetric (r6 S4); permute()+bootstrap() refused r(198) (r7 S3).
- FINDINGS.md holds the full write-up.
