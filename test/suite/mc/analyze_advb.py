#!/usr/bin/env python3
"""analyze_advb.py — judge the adversarial-boundary (AB) MC campaign against the
PREREGISTERED bands and the LOCKED theory predictions in
test/suite/PREREG_ADVBOUND_2026-08-25.md.

The grid (configs, rep counts, locked predictions) lives in advb_grid.json so this
analyzer is fixed before any rep fires. Band arithmetic, the 10-bin uniformity chi2
and the marginal corridors are the house M2 conventions, reused verbatim from
analyze_mc.py so results are comparable to the A1-A6/A4 record.

Writes out/advb_report.md. Exit 1 iff any config fails its own MEASUREMENT integrity
checks (rc!=0, missing p-values, wrong n) — a size FAIL is a finding, not a harness
error, so oversizing does NOT exit nonzero.

usage: analyze_advb.py [outdir] [gridjson]
"""
import csv
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "out")
GRID = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "advb_grid.json")
EPS = 1e-9
CHI2_9_CRIT = 21.666  # df=9, p=.01 (house convention)


def load(prefix):
    rows = []
    for f in sorted(glob.glob(f"{OUT}/{prefix}_s*.csv")):
        with open(f) as fh:
            rows.extend(csv.DictReader(fh))
    return rows


def num(row, k):
    v = row.get(k, ".")
    return None if v in (".", "", None) else float(v)


def chi2_grid(ps):
    """10 equal-probability bins on (0,1]; returns (stat, flag). House convention."""
    n = len(ps)
    if n == 0:
        return 0.0, False
    bins = [0] * 10
    for p in ps:
        b = min(9, max(0, math.ceil(p * 10 - EPS) - 1))
        bins[b] += 1
    e = n / 10.0
    stat = sum((c - e) ** 2 / e for c in bins)
    return stat, stat > CHI2_9_CRIT


def binom_region(n, p):
    """Central ~95% region, normal approx rounded outward (prereg convention)."""
    mu, sd = n * p, math.sqrt(n * p * (1 - p))
    return max(0, math.floor(mu - 1.96 * sd)), min(n, math.ceil(mu + 1.96 * sd))


def rate_verdict(r05, r10, m, bad_rc, missing_p, n):
    """PASS / MARGINAL / FAIL on the house corridors."""
    b05, b10 = binom_region(m, 0.05), binom_region(m, 0.10)
    clean = bad_rc == 0 and missing_p == 0 and n == m
    ok = clean and b05[0] <= r05 <= b05[1] and b10[0] <= r10 <= b10[1]
    rt05, rt10 = r05 / max(1, n), r10 / max(1, n)
    marg = clean and (0.025 <= rt05 <= 0.085) and (0.070 <= rt10 <= 0.135)
    return ("PASS" if ok else "MARGINAL" if marg else "FAIL"), b05, b10


def score_prediction(rate_pct, pred):
    """Did the locked prediction contain the measured rate?"""
    if pred is None or pred.get("band_low_pct") is None:
        return "n/a"
    lo, hi = pred["band_low_pct"], pred["band_high_pct"]
    return "HIT" if lo - EPS <= rate_pct <= hi + EPS else "MISS"


def main():
    with open(GRID) as fh:
        grid = json.load(fh)

    L = ["# Adversarial-boundary MC report (AB campaign)\n",
         "Preregistration: `test/suite/PREREG_ADVBOUND_2026-08-25.md` (grid, locked "
         "predictions and the keep-vs-abandon decision rule were committed before any rep fired).\n",
         "Bands are the house M2 convention: central ~95% binomial region, normal approx "
         "rounded outward; uniformity is a 10-bin chi2 (df=9, crit 21.666 at p=.01); the "
         "MARGINAL corridor is rate@.05 in [.025,.085] and rate@.10 in [.070,.135].\n"]

    integrity_fail = []
    results = {}

    for cfg in grid["configs"]:
        name = cfg["config"]
        m = cfg["M"]
        rows = load(name)
        n = len(rows)
        bad_rc = sum(1 for r in rows if r.get("rc") not in ("0", 0))
        ok_rows = [r for r in rows if r.get("rc") in ("0", 0)]

        L.append(f"\n## {name} — {cfg.get('label', '')}\n")
        L.append(f"`{cfg.get('spec', '')}`\n")
        L.append("| quantity | verdict | pred | detail |")
        L.append("|---|---|---|---|")

        # primary: adversarial p (att side for binary configs, att twin for contreat)
        for qty, col in cfg.get("quantities", [["adversarial", "padv"]]):
            ps = [p for p in (num(r, col) for r in ok_rows) if p is not None]
            missing_p = len(ok_rows) - len(ps)
            r05 = sum(1 for p in ps if p <= 0.05 + EPS)
            r10 = sum(1 for p in ps if p <= 0.10 + EPS)
            stat, flag = chi2_grid(ps)
            v, b05, b10 = rate_verdict(r05, r10, m, bad_rc, missing_p, n)
            rate_pct = 100.0 * r05 / max(1, len(ps))
            pred = (cfg.get("predictions") or {}).get(qty)
            hit = score_prediction(rate_pct, pred)
            predtxt = "n/a" if pred is None else f"{pred['point_pct']:.1f}% [{pred['band_low_pct']:.1f},{pred['band_high_pct']:.1f}] {hit}"
            det = (f"n={n} rc!=0:{bad_rc} p-miss:{missing_p} rej@.05={r05}/{m} "
                   f"({rate_pct:.2f}%) band{list(b05)} rej@.10={r10} band{list(b10)} "
                   f"chi2={stat:.1f}{'*' if flag else ''}")
            L.append(f"| {qty} | {v} | {predtxt} | {det} |")
            results[f"{name}:{qty}"] = {"verdict": v, "rate_pct": rate_pct,
                                        "chi2_flag": flag, "pred": hit, "n": n}
            if bad_rc or missing_p or n != m:
                integrity_fail.append(
                    f"{name}:{qty} integrity — n={n} (want {m}), rc!=0:{bad_rc}, p-missing:{missing_p}")

        # secondaries: joint pre-trend p's, informational
        for qty, col in cfg.get("secondaries", []):
            ps = [p for p in (num(r, col) for r in ok_rows) if p is not None]
            if not ps:
                L.append(f"| {qty} (sec) | INFO | n/a | no non-missing values (expected if a screen suppressed it) |")
                continue
            j05 = sum(1 for p in ps if p <= 0.05 + EPS)
            stat, flag = chi2_grid(ps)
            L.append(f"| {qty} (sec) | INFO | n/a | rej@.05={j05}/{len(ps)} "
                     f"({100.0*j05/len(ps):.2f}%) chi2={stat:.1f}{'*' if flag else ''} |")
            results[f"{name}:{qty}"] = {"verdict": "INFO", "rate_pct": 100.0 * j05 / len(ps),
                                        "chi2_flag": flag, "pred": "n/a", "n": len(ps)}

    # ---- the preregistered decision rule, applied mechanically ----
    L.append("\n## Decision rule (preregistered, applied mechanically)\n")
    rule = grid.get("decision_rule", {})
    certified = [k for k, v in results.items()
                 if v["verdict"] == "PASS" and not v["chi2_flag"] and k.endswith(":adversarial")]
    cluster_keys = set(rule.get("cluster_type_keys", []))
    cert_cluster = [k for k in certified if k in cluster_keys]
    L.append(f"- Certified adversarial configs (PASS band + no chi2 flag): "
             f"{', '.join(certified) if certified else 'NONE'}")
    L.append(f"- Of those, cluster-type VCE configs: "
             f"{', '.join(cert_cluster) if cert_cluster else 'NONE'}")
    L.append(f"\n**{rule.get('text', 'see the preregistration')}**\n")

    # ---- prediction scorecard ----
    L.append("\n## Locked-prediction scorecard\n")
    hits = sum(1 for v in results.values() if v["pred"] == "HIT")
    miss = sum(1 for v in results.values() if v["pred"] == "MISS")
    L.append(f"Consensus band contained the measured rate in {hits} of {hits+miss} scored quantities.\n")

    L.append("\n## Measurement integrity\n")
    if integrity_fail:
        L.extend(f"- FAIL {x}" for x in integrity_fail)
    else:
        L.append("- all configs complete, no nonzero rcs, no missing p-values")

    report = os.path.join(OUT, "advb_report.md")
    with open(report, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L[-14:]))
    print(f"\nwrote {report}")
    return 1 if integrity_fail else 0


if __name__ == "__main__":
    sys.exit(main())
