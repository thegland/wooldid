# architecture recon — NOTES (live)

Pinned ado: .../bluesky/build/wooldid.ado (17021 lines, sha c1d8979)

## Established facts (verified by reading + one Stata run)
- Program layout: `wooldid` shell = L19-528; `wooldid_worker` = L531-4554; `wooldid_validate` = L4568-5512.
  43 subprograms; mata block from ~L13500 to EOF.
- Worker stage banners: L583 init/parse, L648 ifvar, L728 event study+dummies, L1003 lattice,
  L1127 controls, L1333 ESTIMATE, L1572 weights, L1744 gteffects, L1779 sumstats, L1856 histograms+FCOL,
  L1951 main effects/horizons, L2289 cfx, L2787 by-subgroup, L3504 results matrix, L3808 ES composition,
  L3991 present/return, L4506 ES plots, L4529 cleanup.
- Replication loops live in the SHELL, after the worker returns (perm block L~100-303, boot block L~305-513).
  They re-run `qui wooldid_validate `0'` to recover normalized option locals, then call
  wooldid_permrun / wooldid_bootrun, which call `wooldid` (the SHELL) recursively per replication.
- SIZE CEILING MEASURED (szprobe.do/szprobe2.do in this dir): padding wooldid_worker with
  +18 kB of code lines LOADS FINE; +19 kB fails r(1000) at file load. Worker code size
  (leading-whitespace-stripped, comment/blank lines removed) = 114,308 bytes.
  => ceiling ~= 132.7-133.8 kB, HEADROOM ~= 18.4-19.5 kB (~16%), NOT "1.5 kB".
- validate syntax line (L4572) is BYTE-IDENTICAL to worker's (L555). 1511 chars. No drift.
- validate c_local blanket (L5508) covers EVERY option in the syntax line. No stranded option.
- Four sample concepts: raw `if` used at only 5 sites (L615 forward to shell, L657 gen ifvar,
  L676 ifrule, L679/680 dropcount, L1503 comment). `ifvar`/`ifstatement` (L709 `ifvar'==1`)
  is what every margins subpop uses. `esamp` = e(sample) (L1511) drives the WEIGHTS and
  the fast-collapse keep, and is what is posted as e(sample) (L4051).
- w_att = w*cwm is NOT zeroed off e(sample); w_att_it = w*cwm/w_per_it is NOT zeroed off
  e(sample) either (w_per_it counts only esamp rows but is positive for any cell with
  some esamp row). w_att_i and w_att_itime ARE zeroed (their denominators are 0 off-sample).
  This is SAFE only because margins' own default restricts to e(sample) on the full-data path
  and because wooldid_fcol `keep if esamp==1` on the collapsed path. Any margins call that
  gains `noesample` on the FULL data would silently include reghdfe-dropped singletons.
- Fast collapse (default): worker `preserve` L1876, wooldid_fcol keeps esamp rows and
  collapses to one row per (_i,_t,targetflag,ifvar,t,ttre,subgroups[,cluster if impvalues][,i if contreat]).
  8 `restore, preserve` sites in the subgroup block (fcdeg fallback), final `restore` L3496.
  Ineligible list (fcbad, L5069-5123): contreatpoly>1, lattice(), makecfxplots, elasticity
  eyex/eydx, contreat-in-controls, interactivecontrols(), unconditionalse, semielasticity,
  poisexpresults, customvce().
- RNG: permrun L7471 `set seed`, bootrun L9015 `set seed` — NEITHER restores the caller's
  rngstate afterwards. wooldid_esband (L10892/10963) DOES save+restore, but only when
  mode != bootstrap.
- Globals: WOOLDID_CLEANUP, WOOLDID_LASTRC, WOOLDID_LEVEL0, WOOLDID_REPMODE, WOOLDID_BOOTNOSUB.
  REPMODE is nobreak-protected at all 6 setter sites.
- permute() and bootstrap() are mutually exclusive (L6398).

## Open leads to write up
See recon/architecture.md.

## STATUS: DELIVERABLE WRITTEN
/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/architecture.md

## Probes run in this dir (all with logs)
- szprobe.do / szprobe2.do + sz/w_*.ado  -> SIZE CEILING: +18kB code OK, +19kB r(1000).
  Worker code size 114,308 B; ceiling ~132.7-133.8 kB; headroom ~18.5 kB (~400 lines).
- sampleconcept.do -> CONFIRMED: validate's asserts use the RAW user if, not the markout.
  legs B/D/E all r(411) on rows the estimator would have dropped; leg F (add `if id>0`) rc=0.
- state.do -> sort order/vars/obs/level/emptycells/type/estdir all restored; RNG untouched
  on a plain run; permute() advances the stream (expected); WOOLDID_LASTRC survives a
  failed run and refuses replay; e(cmd) still "wooldid" after a failed run (stale e()).
- atomic.do -> post-regression failure leaves nothing behind (estdir empty); verbose LEAVES
  ___wooldidfullmodel in the estimates directory; default vs nofastcollapse bit-identical
  with reghdfe singleton drops.

## Top leads handed downstream (see recon file §9)
1. engines' `smp`/`insmp` = raw if + nonmissing cluster ONLY -> empty clusters counted in G
2. validate asserts on the wrong sample (CONFIRMED)
3. permute() never rewires e(V)/e(V_es_*)
4. stale e() addressable after a failed run
5. ppmlhdfe not checked -> bare r(199)
6. esband RNG/seed-from-name

## POST-ADVISOR CORRECTIONS (final)
- atomic2.do added: REAL poisson run measured. P1 plain poisson -> _ppmlhdfe_d ABSENT,
  estdir EMPTY (clean). P2 poisson+verbose -> _ppmlhdfe_d PRESENT in the user's data
  (nvars 7->8) AND ___wooldidfullmodel left in estdir. MEASURED, not code-read.
- P3: poisson y=-1 on a ZERO-WEIGHT row under [pw=wt] -> r(411) EVEN WITH safety(off)
  (the y>=0 assert at L4984 is not safety-gated). P4: binary id==0 on zero-weight row
  + safety(off) -> rc=0 (the i>0 assert at L4958 IS safety-gated).
- szprobe3.do: ceiling is BYTE-denominated, not line-denominated. 99-char pad lines put
  the boundary at the same 18kB/19kB place (186 lines ok, 196 lines r(1000)) as 28-char
  lines did (674 ok / 711 r(1000)).
- wooldid_bootctl L6460-6466 builds its OWN copy of the loose sample (`smp2`) with the
  same rule, so bootctl and bootrun agree with each other and both disagree with
  e(sample). The "<2 clusters in the estimation sample" message is literally false.
- 9.6 reframed: the esband seed-from-name is D-9 BY DESIGN; do not judge at rtol 1e-10
  (MC error at nd=10000 is ~1e-2). The reportable part is the UNCAPPED call site L4275.
