# mc/ — statistical calibration Monte Carlo (deep-validation component)

Validates the STATISTICAL claims of the new-feature surface — not numerical
reproduction (that is `oldnew`/`battery`/`newfeat`'s job), but whether the p-values
and power claims are CORRECT under known truth:

- **RI size** (`mc_ri.do`, configs R1-R7): under the sharp null with the true
  assignment mechanism equal to the `permute()` scheme, `e(perm_p_att)` must be
  exactly sized on its achievable grid — including under AR(1) serial correlation and
  heteroskedasticity where cluster asymptotics fail (RI's selling point), across
  shuffle / strata / probabilities / permtstat+permtci / permexact / weights.
- **Adversarial line test size** (`mc_adv.do`, A1-A6): flat and planted-line nulls,
  both methods, moderate and small G, AR(1) errors.
- **pretrendpower power** (`mc_ptp.do`): planting a pre-trend equal to the reported
  minimum detectable trend must reject at ~the requested power; the reported NIS size
  must match the empirical any-pre-significant rate under the null.

Design, bands, and seed protocol are PREREGISTERED in
`../PREREG_NEWFEAT_2026-08-24.md` (committed before the first run). The analyzer
(`analyze_mc.py`) hardcodes those integer count bands; deviations require a documented
amendment. Every rep records its seed — any row reproduces standalone.

## Usage

```bash
cd test/suite/mc
./run_mc.sh all          # pretest -> ri -> adv -> pilot -> arms -> analyze (~1.5 h)
./run_mc.sh pretest      # null-fidelity check only (fast; gates the planted-line MC)
./run_mc.sh ri analyze   # one family + judgment
STATA_JOBS=4 ./run_mc.sh all
```

NOT part of `run_full.sh all` (too costly for the standing campaign run). Run it when
the inference machinery changes (permute/advstat/ptp code paths) or at major releases.
Outputs land in `mc/out/` (gitignored): per-rep CSVs, shard logs, `mc_report.md`.
The judged summary is copied into the campaign report when a run is a deliverable.

Verdicts: PASS (inside preregistered band) / MARGINAL (documented corridor) /
FAIL (real mis-sizing -> fix-or-document per the prereg's outcome classification).
Secondary flags (chi-square grid uniformity, null-band and permtci coverage,
joint_horiz sizing) are findings to investigate, never silent.
