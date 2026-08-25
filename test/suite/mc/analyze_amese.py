#!/usr/bin/env python3
"""analyze_amese.py — judge the AME-SE diagnosis grid against the preregistered
ratio bands (PREREG_AMESE_2026-08-24.md). Per case and per quantity (main ame,
pretreatment ame, 8 ES horizons): ratio = mean(reported SE)/SD_true(estimate).
Bands: [0.85,1.15] OK; outside [0.70,1.40] WRONG; between FLAG. Also: SE-collapse
signature (distinct SE values per rep), missing rates, estimate degeneracy, eigen
range. Writes out/amese_report.md."""
import csv
import glob
import math
import statistics
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
CASES = """d1_basic_ci d1_basic_cu d1_basic_rb d1_bycoh_ci d1_bycoh_cu d1_bycoh_rb
d1_bytime_ci d1_bytime_cu d1_bytime_rb
d2_basic_ci d2_basic_cu d2_basic_rb d2_bycoh_ci d2_bycoh_cu d2_bycoh_rb
d2_bytime_ci d2_bytime_cu d2_bytime_rb
d4_within_ci d4_within_cu d4_within_rb d5_within_ci d5_within_cu d5_within_rb
d1c_basic_ci d1c_basic_rb""".split()
QTY = ["main", "pre"] + [f"h{j}" for j in range(1, 9)]  # b1..b10 map: 1=main 2=pre 3..10=ES

def num(v):
    return None if v in (None, ".", "") else float(v)

L = ["# AME-SE diagnosis grid — MC results vs preregistered bands\n"]
L.append("| case | n | rc!=0 | est-SD range (10 qty) | SE-missing | distinct-SE/rep | ratio range | per-qty ratios (main,pre,h1..h8) | verdict |")
L.append("|---|---|---|---|---|---|---|---|---|")

summary = {}
for case in CASES:
    rows = []
    for f in sorted(glob.glob(f"{OUT}/{case}_s*.csv")):
        rows.extend(csv.DictReader(open(f)))
    n = len(rows)
    ok = [r for r in rows if r.get("rc") == "0"]
    badrc = n - len(ok)
    ratios, sds, semiss = [], [], 0
    # distinct SE values within a rep (collapse signature): count over first ok rep
    distinct = None
    for r in ok[:1]:
        vals = [num(r.get(f"s{j}")) for j in range(1, 11)]
        vv = [v for v in vals if v is not None]
        distinct = len(set(round(v, 12) for v in vv)) if vv else 0
    pq = []
    for j in range(1, 11):
        bs = [num(r.get(f"b{j}")) for r in ok]
        ss = [num(r.get(f"s{j}")) for r in ok]
        bs = [b for b in bs if b is not None]
        ssv = [s for s in ss if s is not None]
        semiss += sum(1 for s in ss if s is None)
        if len(bs) >= 30:
            sd = statistics.stdev(bs)
            sds.append(sd)
            if ssv and sd > 0:
                ratio = (sum(ssv) / len(ssv)) / sd
                ratios.append(ratio)
                pq.append(f"{ratio:.2f}")
            else:
                pq.append("na")
        else:
            pq.append("-")
    if ratios:
        rmin, rmax = min(ratios), max(ratios)
        if all(.85 <= x <= 1.15 for x in ratios):
            verdict = "OK"
        elif any(x < .70 or x > 1.40 for x in ratios):
            verdict = "WRONG"
        else:
            verdict = "FLAG"
        rr = f"[{rmin:.2f},{rmax:.2f}]"
    else:
        verdict = "NO-SE"  # SEs never posted (or estimates degenerate)
        rr = "-"
    sdr = f"[{min(sds):.3f},{max(sds):.3f}]" if sds else "-"
    mtot = len(ok) * 10
    L.append(f"| {case} | {n} | {badrc} | {sdr} | {semiss}/{mtot} | {distinct} | {rr} | {' '.join(pq)} | {verdict} |")
    summary[case] = verdict

L.append("\n## Verdict groups\n")
for v in ["OK", "FLAG", "WRONG", "NO-SE"]:
    cs = [c for c, x in summary.items() if x == v]
    if cs:
        L.append(f"- **{v}**: {', '.join(cs)}")

with open(f"{OUT}/amese_report.md", "w") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
