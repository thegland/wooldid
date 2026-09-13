# mc/ — statistical calibration Monte Carlo (deep-validation component)

Validates the STATISTICAL claims of the new-feature surface — not numerical
reproduction (that is `oldnew`/`battery`/`newfeat`'s job), but whether the p-values,
coverage and power claims are CORRECT under known truth. Five campaigns live here,
each with its own preregistration, driver(s), runner and analyzer.

## Inventory (2026-09-13)

**1. NEWFEAT calibration** — prereg `../PREREG_NEWFEAT_2026-08-24.md`; report
`MC_REPORT_2026-08-24.md`; runner `run_mc.sh`; analyzer `analyze_mc.py`; outputs `out/`
(gitignored).
- **RI size** (`mc_ri.do`, configs R1-R7): under the sharp null with the true
  assignment mechanism equal to the `permute()` scheme, `e(perm_p_att)` must be
  exactly sized on its achievable grid — including under AR(1) serial correlation and
  heteroskedasticity where cluster asymptotics fail (RI's selling point), across
  shuffle / strata / probabilities / permtstat+permtci / permexact / weights.
- **Adversarial line test size** (`mc_adv.do`, A1-A6; `mc_adv2.do`, A4 = the
  continuous-treatment arm under a sound VCE): flat and planted-line nulls, both
  methods, moderate and small G, AR(1) errors.
- **pretrendpower power** (`mc_ptp.do`, arms pilot/W/N/Z; `pilot_median.py` reads the
  planted trend off the pilot): planting a pre-trend equal to the reported minimum
  detectable trend must reject at ~the requested power; the reported NIS size must
  match the empirical any-pre-significant rate under the null.
- **Null-fidelity pretest** (`mc_pretest.do`): gates the planted-line MC — the REPORTED
  ES coefficients must lie on the planted line (expects `NULLFID_ALL_PASS`).

Design, bands, and seed protocol are PREREGISTERED (committed before the first run).
The analyzer hardcodes those integer count bands; deviations require a documented
amendment. Every rep records its seed — any row reproduces standalone.

**2. ADVBOUND — where the adversarial test is correctly sized** — prereg
`../PREREG_ADVBOUND_2026-08-25.md`; report `REPORT_ADVBOUND_2026-08-26.md`; figure
`adversarial_calibration.html`.
- `mc_advb.do` (configs AB*) with the grid, rep counts and LOCKED predictions in
  `advb_grid.json`; runner `run_advb.sh`; analyzer `analyze_advb.py`.

**3. AME-SE diagnosis grid** — prereg `../PREREG_AMESE_2026-08-24.md`.
- `mc_amese.do` (MC-truth SEs of the displayed ame quantities, per case) and its
  att-side twin `mc_amese_att.do`; shared case definitions `amese_cases.doh`;
  deterministic structure pass `amese_structure.do`; runner `run_amese.sh` (runs
  `mc_amese.do` only); analyzer `analyze_amese.py`.

**4. F058 permutation extension** — author ruling 2026-09-06; conclusions in
`out2/permutation/F058_perm_conclusions.md` (a few reference CSVs are tracked beside it).
- `mc_ri2.do` (four phases in priority order, D1 = sharp-null size with the step-down
  record; see its header); runner `run_mcri2.sh` (zsh; honours the Stata instance cap;
  usage in its header); analyzer `analyze_mcri2.py <outdir>` (size / FWER / coverage
  tables with Monte-Carlo standard errors).

**5. F058 bootstrap extension** — author ruling 2026-09-06; tables in
`out2/bootstrap/tables.md`.
- `mc_bootse2.do` + `mc_bootse2_dgp.doh` (the two designs): percentile / BC / BCa /
  percentile-t interval coverage, p-value sizes, the joint-test families and the att_it
  re-identification leg; runner `run_bootse2.sh <shard-name> CFG:FIRST:LAST ...`;
  analyzer `analyze_bootse2.py`; `sched_bootse2.sh` is the one-day wall-clock
  time-slicer of that run (no reuse value beyond it).
- `mc_bootse.do` is the ORIGINAL, self-checking Wald-family consistency harness
  (prints `MCB_CHECK` lines; NO analyzer consumes its `mc_bootse_*.csv`). It is NOT the
  bootstrap calibration guard — `mc_bootse2.do` + `analyze_bootse2.py` is.

**6. Stage C scaffolds (written 2026-09-13, NOT yet run in production)** — the three
measurements the 2026-09-13 assessment audit found owed. They are scaffolds: each parses,
runs and writes its CSV, but no campaign has been run and there is no analyzer yet, so
nothing here is evidence of anything until one is. Shared runner `run_stagec.sh`.
- `mc_adv_bs12.do` — the paired SIZE ladder for the BS-12 one-parameter adversarial null.
  Every adversarial size number on record (`MC_REPORT_2026-08-24.md`,
  `REPORT_ADVBOUND_2026-08-26.md`, the help) was measured on the free-intercept statistic
  BS-12 replaced on 2026-09-10, and the ruling's own prerequisite simulation was never run.
  Arms: analytic p, permutation p, bootstrap. `./run_stagec.sh bs12 <outdir> <nreps> <B>
  <slope> <geoms|-> <shards...>`; slope 0 is the null arm.
- `mc_bands.do` — coverage of the simultaneous sup-t band (`simultaneousbands`, PS-45) on
  the analytic-simulated and bootstrap-studentized engines, including a cluster strictly
  BELOW `i`, which no MC has ever used and which is the geometry the help recommends.
  `test_esband.do` pins the band's structure only, never its coverage.
  `./run_stagec.sh bands <outdir> <nreps> <B> <tau> <arms|-> <shards...>`.
- `mc_ri_ext.do` — the randomization-inference extensions with no MC coverage at HEAD:
  continuous treatment, `permside(right)`, `level(90)`, `permstrata()` and a `permtstat`
  null arm. `./run_stagec.sh riext <outdir> <nreps> <B> <tau> <phases|-> <shards...>`.

All three share one positional contract (`<ado> <shard> <nreps> <seedbase> <outdir> <B>
<arg7> <arg8>`), which is why they share a runner; `run_stagec.sh` with no arguments prints
that contract. `B` is the bootstrap/permutation replication count, and wooldid's floor for
`permreps()` is 20, so `B` below 20 makes the permutation arms refuse rather than run.

**Archived diagnostics** (no runner, no analyzer; kept for the record).
- `mc_advdiag.do` (D1-D4) and `mc_advdiag2.do` (V1, V2, D5): the 2026-08-24 A1
  over-rejection diagnostics. They diagnose the PRE-BS-12 adversarial statistic (BS-12,
  2026-09-10, made the adversarial null a one-parameter line through the reference
  period); they still run on HEAD (the safety screen they once needed `safety(off)` for
  is now the PS-64 caution), but their conclusions describe a statistic that no longer
  ships.
- `run_verify2.sh` runs `smoke2/verify2.do`. Both `smoke2/smoke2.do` and
  `smoke2/verify2.do` are tracked as of 2026-09-13 (the 2026-09-13 inventory check found
  the runner tracked but its target untracked, so a fresh checkout had nothing to run);
  their `.log` output stays ignored. The runner keeps a guard that fails loudly if the
  target ever goes missing again.

**Last production runs** (`git log -- test/suite/mc` is the authority) — all BEFORE the
2026-09-08..09-11 inference reworks (BS-12, F033, PS-45, PS-62, batch-3 stored-result
names): campaign 1 reported 2026-08-24 (`MC_REPORT_2026-08-24.md`); campaign 2 reported
2026-08-26, testing stopped by the author 2026-08-27 with the AB29 partial recorded;
campaign 3's grid results (26 cases x 256 reps) committed 2026-08-24; campaigns 4 and 5
recorded 2026-09-06 (`out2/`). Re-run a campaign before citing its numbers for HEAD.

## Usage

```bash
cd test/suite/mc
./run_mc.sh all          # pretest -> ri -> adv -> pilot -> arms -> analyze (~1.5 h)
./run_mc.sh pretest      # null-fidelity check only (fast; gates the planted-line MC)
./run_mc.sh ri analyze   # one family + judgment
STATA_JOBS=4 ./run_mc.sh all
./run_advb.sh [AB1 AB2 ...]        # campaign 2, then analyze_advb.py
./run_amese.sh                     # campaign 3 (26 cases x 256 reps), then analyze_amese.py
./run_mcri2.sh <outdir> <bD1> <bD2> <bD3> <bD4> <maxreps> <verify> <permreps> <shards...>   # campaign 4
./run_bootse2.sh <shard-name> CFG:FIRST:LAST ...                                              # campaign 5
./run_stagec.sh bands|riext|bs12 <outdir> <nreps> <B> <arg7> <arg8|-> <shards...>          # scaffolds (6)
```

NOT part of `run_full.sh all` (too costly for the standing campaign run). Run it when
the inference machinery changes (permute/bootstrap/advstat/ptp code paths) or at major
releases. Campaign 1 outputs land in `mc/out/` (gitignored): per-rep CSVs, shard logs,
`mc_report.md`; the F058 campaigns write `mc/out2/permutation/` and `mc/out2/bootstrap/`
(their conclusion documents and a few reference CSVs are tracked; logs are not).
The judged summary is copied into the campaign report when a run is a deliverable.

Verdicts: PASS (inside preregistered band) / MARGINAL (documented corridor) /
FAIL (real mis-sizing -> fix-or-document per the prereg's outcome classification).
Secondary flags (chi-square grid uniformity, null-band and permtci coverage,
joint_horiz sizing) are findings to investigate, never silent.
