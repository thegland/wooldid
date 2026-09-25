#!/usr/bin/env python3
"""analyze_bootse2.py - tables for the F058 bootstrap calibration extension.

Reads out2/bootstrap/mcb2_<cfg>.csv (one row per arm run) and reports, per
(config, arm): run health, bias, the mean-SE/SD ratio, interval coverage of the
TRUE att with Monte-Carlo standard errors, the p-value sizes, and the joint-test
sizes. Every cell carries n and an MCSE; a cell is FLAGGED when the observed rate
is more than 2 MCSE from nominal.

The 90%-level column: wooldid's percentile and bias-corrected p-values are
LEVEL-INDEPENDENT and are the exact duals of their intervals (wooldid_bootci), so
under a tau=0 DGP the rejection rate at alpha IS the non-coverage rate of the truth
at level 100(1-alpha). That is how the 90 column is obtained for those families
without a second set of replications. The Wald family is recomputed analytically
from b and se. Under a heterogeneous alternative (design B, tau_s varying) no such
dual exists and the 90 cell is left blank by construction.
"""
import csv, glob, math, os, sys
from scipy import stats

OUT = sys.argv[1] if len(sys.argv) > 1 else "out2/bootstrap"
Z95 = stats.norm.ppf(0.975)
Z90 = stats.norm.ppf(0.95)

def num(x):
    if x is None or x in (".", "", "NA"):
        return None
    try:
        v = float(x)
    except ValueError:
        return None
    return None if v != v else v

rows = []
NBAD = 0
for f in sorted(glob.glob(os.path.join(OUT, "mcb2_*.csv"))):
    with open(f) as fh:
        rdr = csv.reader(fh)
        head = next(rdr)
        for raw in rdr:
            # a shard killed mid-write can leave a short final line; drop it rather
            # than let csv.DictReader pad it with None and score it as a result
            if len(raw) != len(head):
                NBAD += 1
                continue
            rows.append(dict(zip(head, raw)))

# de-duplicate on (cfg, draw, arm): a resumed shard may overlap a killed one
seen, keep = set(), []
for r in rows:
    k = (r["cfg"], r["draw"], r["arm"])
    if k in seen:
        continue
    seen.add(k)
    keep.append(r)
rows = keep

CFGDESC = {
    "AN": "design A (24 equal clusters, T=6), tau=0",
    "AA": "design A (24 equal clusters, T=6), tau=2",
    "AT": "design A, tau=0, boottstat",
    "BA": "design B (12 states of 1..12 units, T=8), tau_s=0.1*n_s",
    "BN": "design B, tau=0",
    "BI": "design B, tau_s=0.1*n_s, att_it re-identification leg",
    "BJ": "design B, tau_s=0.1*n_s, event-study window + jointtests",
    "BT": "design B, tau=0, boottstat",
}
ARMDESC = {
    "plain": "analytic clustered",
    "plainit": "analytic clustered (att_it)",
    "std": "bootstrap(standard)",
    "stdbca": "bootstrap(standard) bootdisplay(bca)",
    "bayes": "bootstrap(bayes)",
    "attit": "bootstrap(standard), att_it",
    "bts": "bootstrap(standard) boottstat",
}

MINN = int(os.environ.get("MCB2_MINN", "15"))   # below this a cell is reported, never flagged

def mcse(p, n):
    return math.sqrt(p * (1 - p) / n) if n > 0 else float("nan")

def fmtrate(h, n, nominal):
    """Report the observed rate with its own MCSE; FLAG on the MCSE under correct
    calibration, sqrt(nominal(1-nominal)/n), which is the right yardstick for
    'is this cell off by more than 2 MC standard errors' and does not collapse to
    zero when the observed rate is 0 or 1."""
    if n == 0:
        return "-", False
    p = h / n
    s = mcse(p, n)
    s0 = math.sqrt(nominal * (1 - nominal) / n)
    flag = (n >= MINN) and (abs(p - nominal) > 2 * s0)
    return f"{100*p:.1f} ({100*s:.1f}) n={n}", flag

def implied_df(c, tail=0.975):
    """The t degrees of freedom whose `tail` quantile is c (inf -> None)."""
    if c is None or c <= 0:
        return None
    if c <= stats.norm.ppf(tail) + 1e-6:
        return None
    lo, hi = 1.0, 5000.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if stats.t.ppf(tail, mid) > c:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

L = []
FLAGS = []
if NBAD:
    L.append(f"_{NBAD} malformed CSV line(s) dropped (a shard killed mid-write)._\n")

L.append("## Run inventory\n")
L.append("| cfg | arm | n rows | rc!=0 | mean secs | reps completed | reps failed | boot_N_V |")
L.append("|---|---|---|---|---|---|---|---|")
groups = {}
for r in rows:
    groups.setdefault((r["cfg"], r["arm"]), []).append(r)
order = [k for k in sorted(groups, key=lambda k: (list(CFGDESC).index(k[0]) if k[0] in CFGDESC else 99, k[1]))]
for k in order:
    g = groups[k]
    bad = sum(1 for r in g if r["rc"] != "0")
    secs = [num(r["secs"]) for r in g if num(r["secs"]) is not None]
    nc = [num(r["ncomp"]) for r in g if num(r["ncomp"]) is not None]
    nf = [num(r["nfail"]) for r in g if num(r["nfail"]) is not None]
    nv = [num(r["nv"]) for r in g if num(r["nv"]) is not None]
    def mm(v):
        return f"{sum(v)/len(v):.1f}" if v else "-"
    L.append(f"| {k[0]} | {k[1]} | {len(g)} | {bad} | {mm(secs)} | {mm(nc)} | {mm(nf)} | {mm(nv)} |")

L.append("\n## Point estimate and the standard-error ratio\n")
L.append("(ratio = mean reported SE / SD of the estimate across datasets; 1.00 = calibrated. "
         "The SD denominator is that arm's OWN estimates, because the event-study window "
         "changes the estimator.)\n")
L.append("| cfg | arm | n | truth | mean b | bias | bias t | SD(b) | mean SE | ratio (SE) |")
L.append("|---|---|---|---|---|---|---|---|---|---|")
for k in order:
    g = [r for r in groups[k] if r["rc"] == "0"]
    bs = [num(r["b"]) for r in g if num(r["b"]) is not None]
    ses = [num(r["se"]) for r in g if num(r["se"]) is not None]
    if len(bs) < 2:
        continue
    truth = num(g[0]["truth"])
    n = len(bs)
    mb = sum(bs) / n
    sd = math.sqrt(sum((x - mb) ** 2 for x in bs) / (n - 1))
    bias = mb - truth
    bt = bias / (sd / math.sqrt(n)) if sd > 0 else float("nan")
    if ses and len(ses) > 1:
        ms = sum(ses) / len(ses)
        sds = math.sqrt(sum((x - ms) ** 2 for x in ses) / (len(ses) - 1))
        ratio = ms / sd
        rse = ratio * math.sqrt((sds / math.sqrt(len(ses)) / ms) ** 2 + 1.0 / (2 * (n - 1)))
        rtxt = f"{ratio:.3f} ({rse:.3f})"
    else:
        ms, rtxt = float("nan"), "-"
    note = " *coefficient-draw SD" if k[1] == "bts" else ""
    L.append(f"| {k[0]} | {k[1]} | {n} | {truth:.5f} | {mb:.4f} | {bias:+.4f} | {bt:+.2f} | "
             f"{sd:.4f} | {ms:.4f}{note} | {rtxt} |")

FAMS = [("pct", "percentile"), ("bc", "bias-corrected"), ("bca", "BCa"),
        ("wald", "Wald (normal/t on the bootstrap SD)"), ("main", "displayed family")]

L.append("\n## Interval coverage of the true att, nominal 95 (percent, MCSE in brackets)\n")
L.append("| cfg | arm | percentile | bias-corrected | BCa | Wald | displayed |")
L.append("|---|---|---|---|---|---|---|")
for k in order:
    g = [r for r in groups[k] if r["rc"] == "0"]
    truth = num(g[0]["truth"])
    cells = []
    for key, _ in FAMS:
        h = n = 0
        for r in g:
            lo, hi = num(r[f"{key}_lb"]), num(r[f"{key}_ub"])
            if lo is None or hi is None:
                continue
            n += 1
            h += (lo <= truth <= hi)
        txt, fl = fmtrate(h, n, 0.95)
        if fl:
            FLAGS.append(f"{k[0]}/{k[1]} {key} coverage@95 = {100*h/n:.1f} (n={n}, MCSE {100*mcse(h/n,n):.1f})")
        cells.append(("**" + txt + "**") if fl else txt)
    L.append(f"| {k[0]} | {k[1]} | " + " | ".join(cells) + " |")

L.append("\n## Interval coverage of the true att, nominal 90 (percent, MCSE in brackets)\n")
L.append("(percentile / BC / BCa taken from the level-independent dual p-values, so this row "
         "exists only where the DGP is the null the p tests; Wald recomputed as b +/- z*se. "
         "On bts arms these three are the undisplayed coefficient-draw families - F3, 2026-09-13.)\n")
L.append("| cfg | arm | percentile | bias-corrected | BCa | Wald |")
L.append("|---|---|---|---|---|---|")
for k in order:
    g = [r for r in groups[k] if r["rc"] == "0"]
    isnull = g[0]["dgp"] == "null"
    cells = []
    for pk in ("p_pct", "p_bc", "p_bca"):
        if not isnull:
            cells.append("n/a")
            continue
        ps = [num(r[pk]) for r in g if num(r[pk]) is not None]
        if not ps:
            cells.append("-")
            continue
        h = sum(1 for p in ps if p > 0.10 + 1e-12)
        txt, fl = fmtrate(h, len(ps), 0.90)
        if fl:
            FLAGS.append(f"{k[0]}/{k[1]} {pk} coverage@90 = {100*h/len(ps):.1f} (n={len(ps)}, MCSE {100*mcse(h/len(ps),len(ps)):.1f})")
        cells.append(("**" + txt + "**") if fl else txt)
    h = n = 0
    truth = num(g[0]["truth"])
    for r in g:
        b, se = num(r["b"]), num(r["se"])
        if b is None or se is None or se <= 0 or k[1] in ("bts",):
            continue
        lo95 = num(r["wald_lb"]) if num(r["wald_lb"]) is not None else num(r["main_lb"])
        hi95 = num(r["wald_ub"]) if num(r["wald_ub"]) is not None else num(r["main_ub"])
        c = None
        if lo95 is not None and hi95 is not None:
            c95 = (hi95 - lo95) / (2 * se)
            dfi = implied_df(c95)
            c = stats.t.ppf(0.95, dfi) if dfi else Z90
        if c is None:
            df = num(r["df"])
            c = stats.t.ppf(0.95, df) if (df is not None and df > 0) else Z90
        n += 1
        h += (b - c * se <= truth <= b + c * se)
    txt, fl = fmtrate(h, n, 0.90)
    if fl and n:
        FLAGS.append(f"{k[0]}/{k[1]} Wald coverage@90 = {100*h/n:.1f} (n={n}, MCSE {100*mcse(h/n,n):.1f})")
    cells.append(("**" + txt + "**") if fl else txt)
    L.append(f"| {k[0]} | {k[1]} | " + " | ".join(cells) + " |")

L.append("\n## Mean interval width by family (same datasets, so widths are directly comparable)\n")
L.append("| cfg | arm | n | percentile | bias-corrected | BCa | Wald | displayed |")
L.append("|---|---|---|---|---|---|---|---|")
for k in order:
    g = [r for r in groups[k] if r["rc"] == "0"]
    cells, nn = [], 0
    for key, _ in FAMS:
        w = [num(r[f"{key}_ub"]) - num(r[f"{key}_lb"]) for r in g
             if num(r[f"{key}_lb"]) is not None and num(r[f"{key}_ub"]) is not None]
        nn = max(nn, len(w))
        cells.append(f"{sum(w)/len(w):.4f}" if w else "-")
    L.append(f"| {k[0]} | {k[1]} | {nn} | " + " | ".join(cells) + " |")

L.append("\n## Paired standard-error ratios (per dataset, so free of sampling noise in b)\n")
L.append("(each bootstrap arm's reported SE divided by the SAME dataset's analytic clustered SE, "
         "from the analytic arm run in the same config with the same event-study window.)\n")
L.append("| cfg | arm | n | mean se_boot/se_analytic (SE) | min | max |")
L.append("|---|---|---|---|---|---|")
for k in order:
    if k[1] in ("plain", "plainit"):
        continue
    g = [r for r in groups[k] if r["rc"] == "0"]
    rr = []
    for r in g:
        se, an = num(r["se"]), num(r["analse"])
        if se is None or an is None or an <= 0:
            continue
        rr.append(se / an)
    if len(rr) < 2:
        continue
    m = sum(rr) / len(rr)
    sd = math.sqrt(sum((x - m) ** 2 for x in rr) / (len(rr) - 1))
    note = " *coefficient-draw SD from the hidden Wald table, not the displayed studentized se" if k[1] == "bts" else ""
    L.append(f"| {k[0]} | {k[1]} | {len(rr)} | {m:.3f} ({sd/math.sqrt(len(rr)):.3f}){note} | {min(rr):.3f} | {max(rr):.3f} |")

L.append("\n## p-value size under the null DGP (rejection percent, MCSE in brackets)\n")
L.append("(F3, 2026-09-13: on bts arms p_percentile / p_bc / p_bca are the UNDISPLAYED coefficient-draw "
         "families, posted regardless of boottstat; the DISPLAYED studentized (percentile-t, add-one) p "
         "is the p_bts column, read from e(boot_tstud) row p_pctl. Rows written before 2026-09-13 have "
         "no p_bts column and show '-'.)\n")
L.append("| cfg | arm | p_percentile @5 | @10 | p_bc @5 | @10 | p_bca @5 | @10 | p_bts (studentized, displayed on bts) @5 | @10 |")
L.append("|---|---|---|---|---|---|---|---|---|---|")
for k in order:
    g = [r for r in groups[k] if r["rc"] == "0" and r["dgp"] == "null"]
    if not g:
        continue
    cells = []
    any_p = False
    for pk in ("p_pct", "p_bc", "p_bca", "p_bts"):
        ps = [num(r.get(pk)) for r in g if num(r.get(pk)) is not None]
        lab = pk + (" (coefficient-draw, undisplayed)" if (k[1] == "bts" and pk != "p_bts") else "")
        for a in (0.05, 0.10):
            if not ps:
                cells.append("-")
                continue
            any_p = True
            h = sum(1 for p in ps if p <= a + 1e-12)
            txt, fl = fmtrate(h, len(ps), a)
            if fl:
                FLAGS.append(f"{k[0]}/{k[1]} {lab} size@{a:.2f} = {100*h/len(ps):.1f} (n={len(ps)}, MCSE {100*mcse(h/len(ps),len(ps)):.1f})")
            cells.append(("**" + txt + "**") if fl else txt)
    if any_p:
        L.append(f"| {k[0]} | {k[1]} | " + " | ".join(cells) + " |")

L.append("\n## Joint pre-trend test (parallel pre-trends hold in both DGPs, so this is size)\n")
L.append("| cfg | arm | family | n | rej @5 | rej @10 |")
L.append("|---|---|---|---|---|---|")
for k in order:
    g = [r for r in groups[k] if r["rc"] == "0"]
    ps = [(num(r["jt_p"]), r.get("jt_fam", "")) for r in g if num(r["jt_p"]) is not None]
    if not ps:
        continue
    fam = sorted({f for _, f in ps if f}) or ["-"]
    vals = [p for p, _ in ps]
    c = []
    for a in (0.05, 0.10):
        h = sum(1 for p in vals if p <= a + 1e-12)
        txt, fl = fmtrate(h, len(vals), a)
        if fl:
            FLAGS.append(f"{k[0]}/{k[1]} joint_horiz ({'/'.join(fam)}) size@{a:.2f} = {100*h/len(vals):.1f} (n={len(vals)}, MCSE {100*mcse(h/len(vals),len(vals)):.1f})")
        c.append(("**" + txt + "**") if fl else txt)
    L.append(f"| {k[0]} | {k[1]} | {'/'.join(fam)} | {len(vals)} | " + " | ".join(c) + " |")

L.append("\n## Cells more than 2 MC standard errors from nominal\n")
if FLAGS:
    for f in FLAGS:
        L.append(f"- {f}")
else:
    L.append("- none")

txt = "\n".join(L) + "\n"
with open(os.path.join(OUT, "tables.md"), "w") as fh:
    fh.write(txt)
print(txt)
