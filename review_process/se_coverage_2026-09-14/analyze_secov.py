#!/usr/bin/env python3
"""Judge se_cov_mc.do: are wooldid's analytic SEs calibrated, with csdid and did_imputation
as controls on identical point estimates?

Part 8 C.6 measured wooldid covering at 0.838 / 0.800 on two arms where both external
packages covered correctly in the same draws - on the SAME point estimates, to 1e-14 and
1e-7, from the same data with the same cluster(i). That isolates the variance path exactly.
It rested on 30 draws at one cluster count. This reads the 500-draw run.

TWO measures, and they answer different questions:
  meanSE / SD  - is the reported standard error the right size? Needs no truth at all.
  coverage     - does the interval contain the DGP's true horizon effect?
The controls are what make either interpretable: if wooldid AND both external packages all
came out low, the DGP or the harness would be at fault, not wooldid's variance.

C.6's own flag rule, stated there before any of this ran: an arm fails if its meanSE/SD is
below 0.90 or its 95% coverage below 0.93 WHILE both controls in the same draws sit inside
those bands.

usage: analyze_secov.py <dir-of-se_cov_s*_G*.txt>
"""
import sys, os, glob, re, math, statistics, collections

MINN = 30
LOW_SE, LOW_COV = 0.90, 0.93

ARMS = {
 "esfeb":   "wooldid esfixedbaseperiod, i = unit id",
 "esfceb":  "wooldid esfixedbaseperiod, i = cohort id",
 "pooleb":  "wooldid pooled esprelength(0), i = unit id",
 "poolceb": "wooldid pooled esprelength(0), i = cohort id",
 "cs":      "csdid long2 estat event   [CONTROL]",
 "di":      "did_imputation horizons   [CONTROL]",
}
CONTROLS = {"cs", "di"}

def truth(kind, lab):
    """The DGP's true effect for a coefficient name. se_mc.do's DGP, unchanged."""
    if lab in ("main_att", "ATT", "Post_avg"):
        return 0.5 if kind == "hom" else 119/130
    m = re.match(r"(?:ES_att_)?(pre|post|Tm|Tp|tau)(\d+)$", lab)
    if not m:
        return None
    if m.group(1) in ("pre", "Tm"):
        return 0.0
    h = int(m.group(2))
    if kind == "hom":
        return 0.5
    return (0.7 + 0.05*h) if h <= 4 else (1.0 + 0.05*h)

d = sys.argv[1] if len(sys.argv) > 1 else "."
files = sorted(glob.glob(os.path.join(d, "se_cov_s*_G*.txt")))
if not files:
    sys.exit(f"no se_cov_s*_G*.txt under {d}")

# tag -> coefname -> draw -> (b, se)
D = collections.defaultdict(lambda: collections.defaultdict(dict))
errs = collections.Counter()
Gs = set()
for path in files:
    g = re.search(r"_G(\d+)\.txt$", path)
    if g: Gs.add(int(g.group(1)))
    for line in open(path):
        f = line.rstrip("\n").split("|")
        if f[0] == "err":
            errs[f[1]] += 1
        elif f[0] in ("progress",) or f[0].startswith("SECOV_DONE"):
            continue
        elif len(f) == 5:
            try:
                D[f[0]][f[2]][int(f[1])] = (float(f[3]), float(f[4]))
            except ValueError:
                pass
print(f"files: {len(files)}   G values: {sorted(Gs)}")
print(f"estimation errors by tag: {dict(errs) if errs else 'none'}")

for kind in ("hom", "stag"):
    print(f"\n{'='*78}\n{kind}  ({'common onset' if kind=='hom' else 'staggered, two onsets'})\n{'='*78}")
    summary = {}
    for pre, label in ARMS.items():
        tag = f"{pre}_{kind}"
        if tag not in D:
            continue
        print(f"\n-- {label}   [{tag}]")
        print(f"   {'coef':16s} {'n':>4s} {'SD(b)':>9s} {'meanSE':>9s} {'SE/SD':>7s} {'cov95':>7s}")
        ratios, covs = [], []
        for lab in sorted(D[tag], key=lambda c: (0 if re.search(r'(pre|Tm)\d', c) else 1, c)):
            rec = D[tag][lab]
            ds = sorted(rec)
            bs = [rec[x][0] for x in ds]
            ss = [rec[x][1] for x in ds if rec[x][1] == rec[x][1] and rec[x][1] > 0]
            if len(bs) < 3:
                continue
            sd = statistics.stdev(bs)
            ms = statistics.mean(ss) if ss else float("nan")
            r = ms/sd if sd > 0 else float("nan")
            th = truth(kind, lab)
            cov = float("nan")
            if th is not None:
                hit = [abs(rec[x][0]-th) <= 1.959964*rec[x][1]
                       for x in ds if rec[x][1] == rec[x][1] and rec[x][1] > 0]
                if hit: cov = sum(hit)/len(hit)
            if r == r: ratios.append(r)
            if cov == cov: covs.append(cov)
            print(f"   {lab:16s} {len(bs):4d} {sd:9.4f} {ms:9.4f} {r:7.3f} {cov:7.3f}")
        if ratios:
            n = max(len(D[tag][l]) for l in D[tag])
            summary[pre] = (statistics.median(ratios),
                            statistics.median(covs) if covs else float("nan"), n)
            gate = "" if n >= MINN else f"  (INSUFFICIENT_N {n} < {MINN}, not scored)"
            print(f"   -> median SE/SD {summary[pre][0]:.3f}   median cov95 "
                  f"{summary[pre][1]:.3f}   n={n}{gate}")

    # C.6's flag rule, applied
    ctrl_ok = all(pre in summary and summary[pre][0] >= LOW_SE and summary[pre][1] >= LOW_COV
                  for pre in CONTROLS if pre in summary)
    print(f"\n   CONTROLS ({', '.join(sorted(c for c in CONTROLS if c in summary))}) inside "
          f"[SE/SD >= {LOW_SE}, cov >= {LOW_COV}]: {ctrl_ok}")
    if not ctrl_ok:
        print("   -> C.6's flag rule does NOT apply on this DGP: a wooldid arm below band here "
              "is not evidence about wooldid, because a control is below band too.")
    for pre, (r, c, n) in summary.items():
        if pre in CONTROLS or n < MINN:
            continue
        low = r < LOW_SE or c < LOW_COV
        if low and ctrl_ok:
            print(f"   FLAG  {pre:8s} SE/SD {r:.3f}  cov {c:.3f}  <- below band while both "
                  f"controls are inside it")
        elif low:
            print(f"   note  {pre:8s} SE/SD {r:.3f}  cov {c:.3f}  (below band, but see above)")
        else:
            print(f"   ok    {pre:8s} SE/SD {r:.3f}  cov {c:.3f}")
