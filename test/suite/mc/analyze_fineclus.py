#!/usr/bin/env python3
"""analyze_fineclus.py - judge prereg F3 (P4/P5) from the mc_fineclus.do CSVs.

Judgment convention is INHERITED, not new (PREREG_FEATURE_FINECLUS_2026-09-16.md F0,
which reuses PREREG_STAGEC_2026-09-13.md C0): central ~95% binomial region, normal
approximation rounded outward, with the Wilson interval printed beside it; the MINN=30
not-judged gate; the MARGINAL corridor [.90, .98] for a .95 coverage claim.

Truths are the zero-noise readout of dgpE, verified before any run:
    att_i = 0.800000000000      att = 0.953846153846  (= 248/260)
att coverage is REPORTED, NEVER JUDGED - it is the control that the DGP's cohort
heterogeneity is live (a DGP where att and att_i coincide cannot detect a weighting
error, and these two differ by 0.1538).
"""
import csv
import glob
import math
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
TRUE_ATTI = 0.800000000000
TRUE_ATT = 0.953846153846
NOMINAL = 0.95
MINN = 30
MARG = (0.90, 0.98)
EPS = 1e-12


def num(row, k):
    v = row.get(k, ".")
    if v in (".", "", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def binom_region(n, p):
    """Central ~95% region, normal approx rounded outward (house prereg convention)."""
    mu, sd = n * p, math.sqrt(n * p * (1 - p))
    return max(0, math.floor(mu - 1.96 * sd)), min(n, math.ceil(mu + 1.96 * sd))


def wilson(h, n, z=1.959963985):
    if n <= 0:
        return float("nan"), float("nan")
    p = h / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - hw), min(1.0, c + hw)


def load(cfg):
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, f"fineclus_{cfg}_*.csv"))):
        with open(f) as fh:
            rows.extend(csv.DictReader(fh))
    # dedupe on m (a resumed shard can overlap)
    seen, out = set(), []
    for r in rows:
        if r["m"] in seen:
            continue
        seen.add(r["m"])
        out.append(r)
    return out


def cover(rows, lo, hi, truth):
    h = n = 0
    for r in rows:
        a, b = num(r, lo), num(r, hi)
        if a is None or b is None:
            continue
        n += 1
        if a - EPS <= truth <= b + EPS:
            h += 1
    return h, n


L = []
FAILS = []
L.append("# F3 - att_i coverage under finer-than-i clustering (prereg P4/P5)\n")
L.append(f"truths (zero-noise dgpE readout): att_i = {TRUE_ATTI}, att = {TRUE_ATT}\n")
L.append(f"nominal {NOMINAL}; central ~95% binomial region, normal approx rounded outward; "
         f"MINN={MINN}; MARGINAL corridor {MARG}\n")

for pred, cfg, what in (("P4", "FINE", "cluster(unit), STRICTLY FINER than i=state"),
                        ("P5", "NEST", "cluster(state) == i, the regression leg")):
    rows = load(cfg)
    ok = [r for r in rows if r["rc"] == "0"]
    bad = len(rows) - len(ok)
    L.append(f"\n**{pred} - {cfg}**: {what}. rows={len(rows)}, rc!=0: {bad}. ")
    if rows:
        nd = [num(r, "ndone") for r in ok if num(r, "ndone") is not None]
        nf = [num(r, "nfail") for r in ok if num(r, "nfail") is not None]
        ncl = sorted({int(num(r, "nclust")) for r in ok if num(r, "nclust") is not None})
        L.append(f"replications completed (min/median/max): "
                 f"{min(nd) if nd else '-'}/{sorted(nd)[len(nd)//2] if nd else '-'}/{max(nd) if nd else '-'}; "
                 f"failed total {int(sum(nf)) if nf else 0}; e(N_clust) values seen: {ncl}\n")
    L.append("| pred | cfg | family | estimand | covered/n | rate | band | Wilson | verdict |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for fam, lo, hi in (("percentile (displayed)", "pl", "pu"), ("normal/Wald", "nl", "nu")):
        for est, truth in (("att_i", TRUE_ATTI), ("att", TRUE_ATT)):
            # The PREREGISTERED readout is exactly one cell per config: att_i on the
            # DISPLAYED (percentile) family. att is the control that the heterogeneity
            # is live; the Wald family is a secondary report. Neither is judged.
            judged = (est == "att_i") and fam.startswith("percentile")
            h, n = cover(ok, f"{lo}_{'atti' if est == 'att_i' else 'att'}",
                         f"{hi}_{'atti' if est == 'att_i' else 'att'}", truth)
            band = binom_region(n, NOMINAL)
            wl, wu = wilson(h, n)
            rate = (h / n) if n else float("nan")
            if not judged:
                v = ("REPORTED (control: heterogeneity live)" if est == "att"
                     else "REPORTED (secondary family)")
            elif n < MINN:
                v = f"NOT JUDGED (n={n} < MINN={MINN})"
            elif band[0] <= h <= band[1]:
                v = "PASS"
            elif MARG[0] <= rate <= MARG[1]:
                v = "MARGINAL"
            else:
                v = "FAIL"
                FAILS.append(f"{pred}/{cfg}/{fam}/{est}: {h}/{n} = {rate:.4f}, band {band}")
            L.append(f"| {pred} | {cfg} | {fam} | {est} | {h}/{n} | "
                     f"{rate:.4f} | {band} | [{wl:.4f}, {wu:.4f}] | {v} |")

L.append("")
if FAILS:
    L.append("## FAILURES")
    for f in FAILS:
        L.append(f"- {f}")
else:
    L.append("No FAIL verdicts.")

txt = "\n".join(L)
print(txt)
with open(os.path.join(OUT, "fineclus_report.md"), "w") as fh:
    fh.write(txt + "\n")
sys.exit(1 if FAILS else 0)
