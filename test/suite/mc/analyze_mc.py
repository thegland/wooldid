#!/usr/bin/env python3
"""analyze_mc.py — judge the MC calibration CSVs against the PREREGISTERED bands
(test/suite/PREREG_NEWFEAT_2026-08-24.md; the integer count bands there are
authoritative and hardcoded here). Writes out/mc_report.md. Exit 1 iff any FAIL
verdict (MARGINAL and secondary flags do not fail the component)."""
import csv
import glob
import math
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
EPS = 1e-9
CHI2_9_CRIT = 21.666  # df=9, p=.01

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
    """10 equal-probability bins on (0,1]; returns (stat, flag)."""
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

L = []          # report lines
FAILS = []

def verdict(name, ok, marginal=None, detail=""):
    if ok:
        v = "PASS"
    elif marginal:
        v = "MARGINAL"
    else:
        v = "FAIL"
        FAILS.append(f"{name}: {detail}")
    L.append(f"| {name} | {v} | {detail} |")
    return v

# ---------------- M1: RI size ----------------
RI_BAND_05 = (11, 29)     # prereg counts, M=400
RI_BAND_10 = (28, 52)
BAND_COV = (348, 372)     # comparator in [lb,ub], Binomial(400,.90)
TCI_COV = (371, 389)      # 0 in [tlb,tub], Binomial(400,.95)

L.append("# MC calibration report\n")
L.append("## M1 — RI size calibration (M=400 per config; preregistered count bands: "
         f"alpha=.05 -> {RI_BAND_05}, alpha=.10 -> {RI_BAND_10})\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")

for cfg in ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]:
    rows = load(cfg)
    n = len(rows)
    ok_rows = [r for r in rows if r["rc"] == "0"]
    bad_rc = n - len(ok_rows)
    ps = [num(r, "p") for r in ok_rows]
    ps = [p for p in ps if p is not None]
    r05 = sum(1 for p in ps if p <= .05 + EPS)
    r10 = sum(1 for p in ps if p <= .10 + EPS)
    nfail_pos = sum(1 for r in ok_rows if (num(r, "nfail") or 0) > 0)
    redraws = sum(1 for r in rows if r.get("attempts", "1") not in ("1", ".", ""))
    stat, flag = chi2_grid(ps)
    det = (f"n={n} rc!=0:{bad_rc} rej@.05={r05} rej@.10={r10} "
           f"chi2={stat:.1f}{'*' if flag else ''} permfailed>0:{nfail_pos} redraws:{redraws}")
    ok = (RI_BAND_05[0] <= r05 <= RI_BAND_05[1]) and (RI_BAND_10[0] <= r10 <= RI_BAND_10[1]) and n == 400 and bad_rc == 0
    if cfg == "R5":  # redraw conditioning is expected; bad_rc==0 still required post-redraw
        ok = (RI_BAND_05[0] <= r05 <= RI_BAND_05[1]) and (RI_BAND_10[0] <= r10 <= RI_BAND_10[1]) and n == 400 and bad_rc == 0
    verdict(f"{cfg} size", ok, detail=det)
    if cfg != "R5" and nfail_pos:
        L.append(f"| {cfg} perm_failed | INVESTIGATE | {nfail_pos} reps with perm_failed>0 |")
    # secondary: null band coverage
    cov = sum(1 for r in ok_rows
              if None not in (num(r, "comp"), num(r, "lb"), num(r, "ub"))
              and num(r, "lb") - EPS <= num(r, "comp") <= num(r, "ub") + EPS)
    cflag = not (BAND_COV[0] <= cov <= BAND_COV[1])
    L.append(f"| {cfg} null-band (sec) | {'FLAG' if cflag else 'ok'} | comparator in [lb,ub]: {cov}/400 (band {BAND_COV}) |")
    if cfg == "R3":
        tcov = sum(1 for r in ok_rows
                   if None not in (num(r, "tlb"), num(r, "tub"))
                   and num(r, "tlb") - EPS <= 0 <= num(r, "tub") + EPS)
        tflag = not (TCI_COV[0] <= tcov <= TCI_COV[1])
        L.append(f"| R3 permtci cover (sec) | {'FLAG' if tflag else 'ok'} | 0 in CI: {tcov}/400 (band {TCI_COV}) |")

# ---------------- M2: adversarial size ----------------
L.append("\n## M2 — adversarial size (M=1000: bands .05 -> [36,64], .10 -> [81,119]; "
         "marginal rate corridors [.025,.085] / [.070,.135]; A5 report-only)\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")

ADV = {"A1": 1000, "A2": 1000, "A3a": 1000, "A3b": 1000, "A6": 1000, "A5": 504}
for cfg, m in ADV.items():
    rows = load(cfg)
    n = len(rows)
    ok_rows = [r for r in rows if r["rc"] == "0"]
    bad_rc = n - len(ok_rows)
    ps = [num(r, "padv") for r in ok_rows]
    missing_p = sum(1 for p in ps if p is None)
    ps = [p for p in ps if p is not None]
    r05 = sum(1 for p in ps if p <= .05 + EPS)
    r10 = sum(1 for p in ps if p <= .10 + EPS)
    stat, flag = chi2_grid(ps)
    det = (f"n={n} rc!=0:{bad_rc} p-missing:{missing_p} rej@.05={r05} rej@.10={r10} "
           f"rate=({r05/max(1,len(ps)):.4f},{r10/max(1,len(ps)):.4f}) chi2={stat:.1f}{'*' if flag else ''}")
    if cfg == "A5":
        L.append(f"| A5 size (report-only) | INFO | {det} |")
        continue
    b05, b10 = (36, 64), (81, 119)
    ok = (b05[0] <= r05 <= b05[1]) and (b10[0] <= r10 <= b10[1]) and bad_rc == 0 and missing_p == 0
    rt05, rt10 = r05 / max(1, len(ps)), r10 / max(1, len(ps))
    marg = (.025 <= rt05 <= .085) and (.070 <= rt10 <= .135) and bad_rc == 0 and missing_p == 0
    verdict(f"{cfg} size", ok, marginal=marg, detail=det)
    pj = [num(r, "pjoint") for r in ok_rows]
    pj = [p for p in pj if p is not None]
    if pj:
        j05 = sum(1 for p in pj if p <= .05 + EPS)
        L.append(f"| {cfg} joint_horiz (sec) | INFO | rej@.05={j05}/{len(pj)} ({j05/len(pj):.4f}) |")

# ---------------- M3b: pretrendpower power ----------------
L.append("\n## M3b — pretrendpower power arms (bands: W/N rejection in [.72,.87] of 400; "
         "Z NIS-rate within binomial 95% of mean reported nissize)\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")

for arm, crit in [("W", "prep"), ("N", "maxt")]:
    rows = load(f"ptp_{arm}")
    ok_rows = [r for r in rows if r["rc"] == "0"]
    n = len(ok_rows)
    if arm == "W":
        rej = sum(1 for r in ok_rows if (num(r, "prep") or 1) <= .05 + EPS)
    else:
        rej = sum(1 for r in ok_rows if (num(r, "maxtpre") or 0) > 1.959963985)
    rate = rej / max(1, n)
    gam = num(rows[0], "gamma") if rows else None
    ok = (.72 <= rate <= .87) and n >= 390
    verdict(f"ptp arm {arm}", ok, detail=f"gamma={gam} n={n} rej={rej} rate={rate:.4f}")

# Z arm + pilot: NIS empirical size vs reported nissize
zrows = [r for r in load("ptp_Z") + load("ptp_pilot") if r["rc"] == "0"]
nz = len(zrows)
if nz:
    nis_emp = sum(1 for r in zrows if (num(r, "maxtpre") or 0) > 1.959963985)
    nis_rep = [num(r, "nissize") for r in zrows]
    nis_rep = [v for v in nis_rep if v is not None]
    mean_rep = sum(nis_rep) / max(1, len(nis_rep))
    lo, hi = binom_region(nz, mean_rep)
    ok = lo <= nis_emp <= hi
    verdict("ptp Z nissize", ok,
            detail=f"n={nz} empirical={nis_emp} ({nis_emp/nz:.4f}) reported mean={mean_rep:.4f} band=[{lo},{hi}]")
    wrej = sum(1 for r in zrows if (num(r, "prep") or 1) <= .05 + EPS)
    L.append(f"| ptp Z wald size (sec) | INFO | prep rej@.05 = {wrej}/{nz} ({wrej/nz:.4f}) |")

L.append("\n## Verdict summary\n")
L.append("FAILS: " + (", ".join(FAILS) if FAILS else "none"))

os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/mc_report.md", "w") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
sys.exit(1 if FAILS else 0)
