# NOTES — agent: anticipation() / replay / level()  (method: exact equivalence testing)

Build asserted via `adopath ++ .../bluesky/build`; every probe prints `which wooldid`.

## Code-reading map (pinned build, sha c1d8979)

### anticipation()
- ado:751-752 — the single shift:
  `qui gen double ttre = ttre_asgiven <shift> if ttre_asgiven <shift> <= maxyearobserved + oostreatedcontrols & ifstatement`
  Both the VALUE and the out-of-window screen move together, so it matches
  `gen ttre2 = ttre - k` exactly at that line.
- ado:5347-5348 (validate, non-commensurate-grid tier) pre-adds `anticipation` to both
  window bounds, so that screen is also shift-equivalent.
- The only anticipation-conditional code paths are CAUTIONS:
  `wooldid_anticesf` (ado:861/5537-5559), `wooldid_anticgrid` (ado:945/5519-5535),
  reltime ado:12036-12047, and the three "Note:" lines at ado:919/960/12173.
  None of them changes the sample; the drops themselves come from the shared
  `nopre` rule computed on the SHIFTED `ttre`.
- Therefore a uniform shift cannot create a cohort COLLISION (it is injective), and the
  seeded "collision" arm is vacuous. Dropped.
- Rep option strings carry `anticipation(k)` separately (ado:163, ado:364) — that is the
  one non-obvious place a divergence could hide.

### level()
- worker ado:574-579 — `global WOOLDID_LEVEL0 = c(level)` then `set level level`, only when
  `level != 0 & wdrepmode == ""`. `wdlvl = c(level)` read AFTER the set.
- shell ado:68-71 — restores and drops the global, **immediately after the worker returns
  and BEFORE the permute/bootstrap blocks**. Those blocks pass `level(`level')` explicitly,
  so that ordering is (so far) benign.
- Range check `"Error: level() must be between 10 and 99.99."` exists ONLY at ado:5863
  (wooldid_permctl) and ado:6263 (wooldid_bootctl). **The analytic path has no range check**
  — it falls through to Stata's own `set level`.  <-- probe E in p1.
- Label locals: ado:579 (`wdlvltxt`), 7377 (`wdbplvl`), 8036/8710 (`wdplvltxt`),
  9722 (`wdblvltxt`) — all the same expression; the boot/perm ones are guarded by
  `if level != 95`, which is a no-op rename when the source column is already correct.
- Roth block: live prints `level'%` (ado:10619); replay prints
  `cond(e(level)<., strofreal(e(level)), "95")%` (ado:11486). Two different sources.

### replay
- ado:24-53. Gated on a strictly empty command line + `e(cmd)=="wooldid"`;
  refuses while `$WOOLDID_LASTRC` is set. Prints part(1), part(roth), part(2),
  part(perm), part(boot).
- Live analytic order: worker part(1) @4491, ptp @4498, part(2) @4500 — same order.
- Live bootstrap order: shell part(1) @441, ptp @449, part(2) @486, part(boot) @507 — same.
- `permute() + pretrendpower()` is REFUSED at wooldid_permctl, so the advisor's
  permute-order hypothesis is unreachable. `hbaronly` (worker @4499, bootstrap leg)
  `exit`s before any display (ado:10421), so it prints nothing live.

## Probes
- `p1_level.do`  — level() identity on every analytic surface + restoration + range. RUNNING.
- `p2_antic.do`  — 16-config exact-equivalence sweep, %21x on every visible e(). WRITTEN.
- `p3_replay.do` — 12-config live-vs-replay log diff + staleness. WRITTEN.
- helpers: `mkdata.do` (generator), `dumpe.do` (e() -> hex text), `diffe.py`, `diffrep.py`.

## Results so far
- p1 A-block (level(80), robust, att/att_it/att_i/att_itime + subgroups + gteffects):
  `e(wdidmainresults)`, `e(wdidsgresults_1)`, `e(pwsgs)` all carry `lb_80ci/ub_80ci`
  and satisfy `lb == b - invttail(611,.10)*se` to <3e-15 relative. `e(gteffects)` colnames
  include `lb_80ci ub_80ci`. `e(level)==80`. Session `c(level)` restored to 95 after.
  => level() is honoured on the analytic tables. No finding yet.

## 2026-09-17, progress log

### p1 (level, wave 1) COMPLETE — clean
- level(80): main / subgroup / pairwise tables all `lb_80ci`, identity `lb = b - invttail(df,.10)*se`
  to <3e-15 rel.  `e(gteffects)` carries `lb_80ci`, and its cells DO move with level
  (`GTEFFECTS_IDENTICAL_ACROSS_LEVELS=0`).
- level(92.5) -> `lb_92p5ci`, e(level)=92.5, identity holds.
- level(92.567) -> refused r(198) "level() can have at most two digits after the decimal point"
  (Stata's own `set level` message); session level restored; $WOOLDID_LEVEL0 dropped.
- level(5)/level(100)/level(-1)/level(0.5) -> refused r(198)
  "level() must be between 10 and 99.99 inclusive"; session level restored; global dropped.
- restoration after a FAILED run (r459 after the set fired): c(level) back to 90, global empty.
- back-to-back level(80), level(99), then no level() with session at 90: e(level) 80/99/90,
  c(level) stays 90 throughout, labels lb_90ci on the third. No leak.
=> NO FINDING on the level() core contract.

### tooling
- `dumpe.do` had to use `local X : e(scalars)`; `` `e(scalars)' `` returns nothing.
- positive control confirmed: diffe.py flags 38 differences between two genuinely
  different configs, and CLEAN on the anticipation pair c01.

### p2 (anticipation exact equivalence) RUNNING
- c01 (att robust, k=2): CLEAN across 14 scalars / 17 macros / 3 matrices.
  Only ignored entries are e(anticipation) (documented) and e(cmdline).

### Structural finding candidate found by reading (to be confirmed by p5)
- 21 `Caution:` sites live in wooldid_worker, 15 in wooldid_validate, 5 in the `wooldid`
  shell itself, and only 8 in wooldid_display. sthlp:1458 documents that replay "does not
  repeat estimation-time Caution messages". BUT ado:299 and ado:511 print, from the SHELL,
  a caution that describes the RESULTS rather than the estimation run:
  "Caution: histogram outputs are drawn from the ANALYTIC covariance ...; the reported
   tables ... are bootstrap-based." That one is lost on replay.

### p6 (anticipation follow-ups) COMPLETE
- SEMANTIC check: with anticipation(2) and recorded onset 2008, e(gteffects) puts relyr=0 at
  caltime=2006 and reports onset=2006. 841 rows, zero rows where relyr==0 and caltime!=onset.
  Exactly what sthlp:821 promises.
- k=0 is an EXACT no-op (out_k0_A vs out_k0_B CLEAN, 21 scalars / 18 macros / 6 matrices);
  e(anticipation) is not posted at k=0 (ado:4202 gates on >0), matching sthlp:826.
- anticipation(-1) refused r(198) with a named, informative message.
- anticipation(2.5) refused r(198) "option anticipation() incorrectly specified" (integer option).
- subgroups(ttre) under anticipation(2) keeps the RECORDED onsets as labels: e(wdid_sglist)
  = "2005 2008 2010", not "2003 2006 2008". That is sthlp:826's own carve-out, and it is the
  one literal exception to the same paragraph's "exactly equivalent to replacing ttre with
  ttre - #". Documented; NOT filed.

### p3 (replay) partial
- r01 r02 r04 r06 r07: every printed table, diagnostic and footer is BYTE-IDENTICAL.
  The only divergence is the block of leading estimation-time cautions. -> AR-1.
- r08 (permute): same, PLUS an extra blank line the replay prints before the
  "Randomization Inference" block that the live run does not (ado:48
  `if "`e(perm_method)'" != "" dis as result " "`). Cosmetic. -> AR-2.
- r03/r05/r09 were my own option errors (contreat needs contreatcontrols(),
  pretrendpower() needs esfixedbaseperiod); redone in p7.

### p2 (anticipation) partial: c01-c09, c15 all CLEAN. c06/c14 were my option errors.

### p5 (replay, adversarial) COMPLETE
- h01 (bootstrap + histogramcohorteffects) and h02 (permute + hist): the shell-only histogram
  caution IS lost on replay. CONFIRMED -> AR-1 item 6.
- h03 (analytic + hist): REPLAY-EXACT (35 lines). No caution fires there, so nothing to lose.
- b01: replay with NO DATASET IN MEMORY is byte-identical to the live display (228 vs 231
  lines; the 3 extra lines are do-file echoes). Matches sthlp:1458. CLEAN.
- c01: estimates store -> another estimation -> estimates restore -> replay: byte-identical.
- c03: with TWO stored wooldid models, restoring the first and replaying shows the first
  model's own output byte-identically. CLEAN.
- e01: two consecutive replays are byte-identical (idempotent). CLEAN.
- C02 TRAP CONFIRMED -> AR-3:
    e(cmd)==regress; `wooldid'  ->  "Error: command line input should feature all of
      y i t ttre; too few variables specified." r(198)
    $WOOLDID_LASTRC becomes 198
    estimates restore WD2 (e(cmd)==wooldid); `wooldid'  ->  refused r(198) with
      "the last wooldid call failed with r(198); the stored results predate that failed
       call. Refit wooldid before replaying."
  The "failed call" was the user's own replay attempt; nothing is stale.
- F block: fe(nosuchvar) fails r(111); e() still holds the PREVIOUS good run (e(N)=1200,
  e(cmd)=wooldid) and replay is refused r(198) -- documented. `macro drop WOOLDID_LASTRC`
  then replays the pre-failure results with no marker at all (recon already knew this; the
  global is documented as reserved, so not filed).

### p4 (level wave 2) partial
- L1: level(80) is EXACTLY `set level 80` -- 40 scalars, 19 macros, 22 matrices all
  bit-identical (att/att_it/att_i/att_itime + subgroups + gteffects + jointtests).

### p2 (anticipation) c17 (permute(shuffle) permreps(25) rseed(4242)) CLEAN:
  37 scalars, 23 macros, 23 matrices bit-identical. The replication engine carries
  anticipation(k) in its rep option string (ado:163) and the result is the same as
  permuting the already-shifted onsets.

### p2 + p7 ANTICIPATION SWEEP COMPLETE -- 26 of 26 estimable configs BIT-IDENTICAL
CLEAN (%21x on every visible e() scalar/macro/matrix, ignoring only e(anticipation) and
e(cmdline)): c01 att robust; c02 four estimands + cluster(i) + jointtests; c03 esfixedbaseperiod
+ esrelativeto(-2); c04 oostreatedcontrols(3); c05 subgroups; c07 poisson; c08 timetrends;
c09 gteffects + summarystats; c10 adversarial; c11 controls + fe(i.sg); c12 [aw=w] k=1;
c13 k=5; c15 if sg==1; c16 att_it impvalues cluster(i); c17 permute(shuffle) permreps(25);
c18/c19 bootstrap(standard) reps(25) (att and att_i legs); c20 BIENNIAL with ODD k=1 (off-grid
effective onsets); c21 biennial even k; c24/c25 unbalanced administrative panel (late entry,
early exit, 4% interior gaps) with and without esfixedbaseperiod; c26 k=6 (pushes the earliest
cohort's effective onset to the first period); c27 subgroups(cohort); f06 contreat +
contreatcontrols(basic); f14/f14b pretrendpower + esfixedbaseperiod.
EQUAL FAILURE in both arms: c22/c23 (no never-treated units at all) r(459) in BOTH arms;
f06b r(907) maxvar in BOTH arms.  Not divergences.
k=0 exact no-op; anticipation(-1) and anticipation(2.5) refused informatively.

### p4 LEVEL WAVE 2
- level(L) == `set level L' exactly, on L1 (analytic 4 estimands + subgroups + gteffects),
  L2 (bootstrap), L3 (permute + permtstat + permtci), L4 (pretrendpower), L5 (poisson).
- LEVEL INVARIANCE: with a fixed rseed, level(80) vs level(95) is bit-identical on everything
  except interval columns.  I1 (analytic) fully invariant.  I2/I3 differ ONLY in
  e(boot_ptails) rows cp_lo/cp_hi and e(perm_se_shift) -- both documented as level-derived
  (sthlp:1345 "the interval-implied standard error").  So level() perturbs no point estimate,
  standard error, statistic or p-value, and does not disturb the RNG stream.
- bootstrap at level(92.5): e(wdidmainresults), e(pwsgs) and e(wdidsgresults_1) all carry
  lb_92p5ci -- the ado:7377 `wdbplvl' path is correct.
- impvalues at level(80) (p7): lb_80ci, and lb == b - invttail(74,.10)*se to 0.0e+00.

### p8-p15: the level() boundary (AR-4) -- the big one
- Label round-trip: nbad=0 over all 9000 two-decimal levels (worker's `wdlvltxt' vs the
  gteffects Mata `gtlt' round-trip). Clean.
- gteffects CI identity at L = 80, 92.5, 10, 33.33, 66.67: worst rel dev 1.2e-13. Clean.
- level(99.99) REFUSED r(198). Mapping (p11): refused 99.9 99.93 99.96 99.99 98.9 79.9;
  accepted 99.91 99.92 99.94 99.95 99.97 99.98 99.5 92.5 10.01 88.8 44.4 12.3.
  The split is EXACTLY {L : %.16g(L) != %.15g(L)} -- 1152 of 9000 (12.8%), all >= 64.01.
- Stata's own `set level' accepts every one of them (p10/repro_AR4 block A).
- Fires with NO level() option via validate ado:5507 (p12).
- A SUCCESSFUL run is reported as a failure, the caller's level is left changed, the
  epilogue is skipped and the replay marker is not armed (p15, p14).
- p14 is the documented-guarantee contradiction: bare `wooldid' after the failed call
  replays the pre-failure results rc 0 instead of refusing.

### FINAL STATE: FINDINGS.md holds AR-1 (medium), AR-2 (low), AR-3 (medium), AR-4 (high).
Repro files: repro_AR1.do repro_AR2.do repro_AR3.do repro_AR4.do -- all verified to run.
