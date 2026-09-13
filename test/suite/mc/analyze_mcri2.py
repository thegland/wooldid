#!/usr/bin/env python3
"""analyze_mcri2.py - F058 PERMUTATION part: read the mc_ri2 shard CSVs and print the
size / FWER / coverage tables with Monte-Carlo standard errors and n per cell.

Usage: analyze_mcri2.py <outdir> [permreps]
  permreps = the run's permreps() (default 199, the production value); pass the smoke's value to
  judge a smaller run - the integrity and lattice checks are exact in it (F12, 2026-09-13).
Every table cell is a binomial share; the MC se printed is sqrt(phat(1-phat)/n) and the
DEVIATION test is |phat - nominal| > 2*sqrt(nominal(1-nominal)/n) (the se under the null,
which is the right yardstick for a calibration verdict)."""
import csv, glob, math, os, sys
from collections import defaultdict

OUT = sys.argv[1] if len(sys.argv) > 1 else "out2/permutation"
NDRAW = int(sys.argv[2]) if len(sys.argv) > 2 else 199   # permreps() of the run (production 199)
HALF = (NDRAW + 1) / 2.0                                  # equal-tailed p lattice step is 2/(NDRAW+1)
# theory of the p-value on the (NDRAW+1)-rank lattice, by enumeration of the observed rank R:
# equal-tailed p = min(1, 2 min(R, NDRAW+2-R)/(NDRAW+1)); one-sided share = R/(NDRAW+1)
_R = range(1, NDRAW + 2)
_peq = [min(1.0, 2.0 * min(R, NDRAW + 2 - R) / (NDRAW + 1)) for R in _R]
_plo = [R / (NDRAW + 1) for R in _R]
MEAN_EQ = sum(_peq) / len(_peq)
SD_EQ = math.sqrt(sum(p * p for p in _peq) / len(_peq) - MEAN_EQ ** 2)
MEAN_LO = sum(_plo) / len(_plo)
EPS = 1e-9
FLAGS = []

def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, "mc_ri2_s*.csv"))):
        with open(f) as fh:
            rows.extend(csv.DictReader(fh))
    return rows

def num(r, k):
    v = r.get(k, ".")
    if v in (".", "", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None

def share(hits, n):
    if n == 0:
        return None, None
    p = hits / n
    return p, math.sqrt(p * (1 - p) / n) if n else None

def line(label, hits, n, nominal=None, extra=""):
    p, se = share(hits, n)
    if p is None:
        print(f"| {label} | - | - | 0 | {extra} |")
        return
    dev = ""
    if nominal is not None:
        se0 = math.sqrt(nominal * (1 - nominal) / n)
        z = (p - nominal) / se0 if se0 > 0 else 0.0
        dev = f"{z:+.2f}"
        if abs(z) > 2:
            FLAGS.append(f"{label}: {p:.4f} vs nominal {nominal:.3f} ({z:+.2f} MC se, n={n})")
    print(f"| {label} | {p:.4f} | {se:.4f} | {n} | {dev} | {extra} |")

def hdr(cols):
    print("| " + " | ".join(cols) + " |")
    print("|" + "|".join("---" for _ in cols) + "|")

rows = load()
if not rows:
    sys.exit("no CSV rows found in " + OUT)

by = defaultdict(list)
for r in rows:
    by[(r["phase"], r["arm"])].append(r)

print("# F058 - permutation calibration (mc_ri2)\n")
print(f"rows read: {len(rows)}  from {len(glob.glob(os.path.join(OUT,'mc_ri2_s*.csv')))} shard files; judged at permreps = {NDRAW}\n")
for k in sorted(by):
    ok = sum(1 for r in by[k] if r["rc"] == "0")
    print(f"- phase {k[0]} arm {k[1]}: {len(by[k])} rows, {ok} with rc=0")
print()

D1 = [r for r in by[("D1", "main")] if r["rc"] == "0"]

# ------------------------------------------------------------------ integrity
print("## 0. Integrity of the run\n")
bad_names = sum(1 for r in D1 if num(r, "namesok") != 1)
incomplete = sum(1 for r in D1 if (num(r, "ncomp") or 0) != NDRAW)
nfailpos = sum(1 for r in D1 if (num(r, "nfail") or 0) > 0)
ndr = [num(r, "ndraw") for r in D1 if num(r, "ndraw") is not None]
print(f"- D1 reps with rc=0: {len(D1)}")
print(f"- column-name lookups that fell back to position: {bad_names}")
print(f"- reps with e(perm_completed) != {NDRAW}: {incomplete}; with e(perm_failed)>0: {nfailpos}")
if ndr:
    print(f"- usable draws per rep for main_att: min {min(ndr):.0f} max {max(ndr):.0f}")
    if max(ndr) != NDRAW:
        print(f"- WARNING: the CSV's ndraw column peaks at {max(ndr):.0f} but this run is being judged "
              f"at permreps = {NDRAW}; pass the run's permreps as the second argument")
        FLAGS.append(f"ndraw column max {max(ndr):.0f} != judged permreps {NDRAW}")
# self-consistency: the driver's recomputed equal-tailed p must equal e(perm_p)
mism = [(r["shard"], r["rep"], num(r, "p_main"), num(r, "cpeq")) for r in D1
        if num(r, "p_main") is not None and num(r, "cpeq") is not None
        and abs(num(r, "p_main") - num(r, "cpeq")) > 1e-12]
print(f"- recomputed equal-tailed p vs e(perm_p) main_att: {len(mism)} mismatches of {len(D1)}")
if mism[:3]:
    print(f"  first mismatches: {mism[:3]}")
# permside(left) identity: e(perm_p_lo) must be the permside(left) p
mism2 = [1 for r in D1 if num(r, "plo_main") is not None and num(r, "cpleft") is not None
         and abs(num(r, "plo_main") - num(r, "cpleft")) > 1e-12]
print(f"- recomputed left-tail p vs e(perm_p_lo) main_att: {len(mism2)} mismatches of {len(D1)}")

# lattice check: with all NDRAW draws complete for a column, the equal-tailed p must be an
# exact multiple of 2/(NDRAW+1) (.01 at 199; or exactly 1, the cap) and the one-sided shares multiples of 1/(NDRAW+1). A column
# whose accumulator lost draws (an unidentified horizon in some reassignment) lands off that
# lattice, which is exactly the selection that would break the exactness argument in s4.
print("\n### lattice check (a p off the 1/200 lattice means that column lost draws)\n")
hdr(["element", "reps on lattice", "n", "off-lattice examples"])
for nm, col in [("main_att", "p_main"), ("pretreatment_att", "p_pre"), ("ES_att_pre2", "p_pre2"),
                ("ES_att_pre1", "p_pre1"), ("ES_att_post0", "p_post0"), ("ES_att_post1", "p_post1")]:
    vals = [num(r, col) for r in D1 if num(r, col) is not None]
    # p = 1 is the min(1, .) cap of the equal-tailed rule and is attainable at an even NDRAW
    off = [v for v in vals if v < 1 - EPS and abs(v*HALF - round(v*HALF)) > 1e-6]
    print(f"| {nm} | {len(vals)-len(off)} | {len(vals)} | {sorted(set(off))[:3]} |")
    if off:
        FLAGS.append(f"{nm}: {len(off)} of {len(vals)} reps have a p off the {NDRAW}-draw lattice "
                     f"(that column lost replications; the exactness argument is conditional)")
# e(joint_horiz) is the one-sided add-one p = Q/(NDRAW+1), so it must be a multiple of 1/(NDRAW+1).
# Off-lattice would mean the metric saw fewer than NDRAW complete draws - or that the ANALYTIC
# joint p survived under the permutation header.
_jv = [num(r, "jointp") for r in D1 if num(r, "jointp") is not None]
_jo = [v for v in _jv if abs(v*(NDRAW+1) - round(v*(NDRAW+1))) > 1e-6]
print(f"| e(joint_horiz) | {len(_jv)-len(_jo)} | {len(_jv)} | {sorted(set(_jo))[:3]} |")
if _jo:
    FLAGS.append(f"e(joint_horiz): {len(_jo)} of {len(_jv)} reps off the 1/{NDRAW+1} lattice "
                 f"(fewer than {NDRAW} complete draws in the metric, or an analytic value survived)")

# ------------------------------------------------------------------ package verification
print("\n## 1. Package verification of the alternative p-value conventions\n")
print("Every VERIFY-th D1 rep was re-run through the package with permabsp / permrawp /")
print("permside(left) on the SAME rseed (identical draws). The package p is compared with")
print("the value the driver computed from e(permdist_b) + e(perm_comparator_b).\n")
mainbykey = {(r["shard"], r["rep"]): r for r in D1}
hdr(["arm", "n verified", "exact matches", "max abs diff"])
for arm, comp in (("absp", "cpabs"), ("rawp", "cpraw"), ("left", "cpleft")):
    vr = [r for r in by[("D1", arm)] if r["rc"] == "0"]
    nmatch, mx, nn = 0, 0.0, 0
    for r in vr:
        m = mainbykey.get((r["shard"], r["rep"]))
        if not m:
            continue
        a, b = num(r, "p_main"), num(m, comp)
        if a is None or b is None:
            continue
        nn += 1
        d = abs(a - b)
        mx = max(mx, d)
        if d <= 1e-12:
            nmatch += 1
    print(f"| {arm} | {nn} | {nmatch} | {mx:.3e} |")
    if nn and nmatch != nn:
        FLAGS.append(f"convention {arm}: package/driver disagreement on {nn-nmatch} of {nn} verified reps")

# ------------------------------------------------------------------ D1 sizes
print("\n## 2. D1 - size under the sharp null (tau = 0)\n")
print(f"Nominal size is attained EXACTLY on the {NDRAW}-replication grid: with the observed")
print(f"statistic ranked R uniform on 1..{NDRAW+1}, the equal-tailed p is 2*min(R,{NDRAW+2}-R)/{NDRAW+1}, so")
print(f"P(p<=.05)={sum(1 for p in _peq if p <= .05 + EPS)}/{NDRAW+1}={sum(1 for p in _peq if p <= .05 + EPS)/(NDRAW+1):.4f} "
      f"and P(p<=.10)={sum(1 for p in _peq if p <= .10 + EPS)}/{NDRAW+1}={sum(1 for p in _peq if p <= .10 + EPS)/(NDRAW+1):.4f}; the same holds for permabsp,")
print("permrawp and permside(). Any departure below is a real calibration finding.\n")
ELEMS = [("main_att", "p_main"), ("pretreatment_att", "p_pre"), ("ES_att_pre2", "p_pre2"),
         ("ES_att_pre1", "p_pre1"), ("ES_att_post0", "p_post0"), ("ES_att_post1", "p_post1")]
for alpha in (0.05, 0.10):
    print(f"\n### equal-tailed p, alpha = {alpha:.2f}\n")
    hdr(["element", "size", "MC se", "n", "dev (MC se)", ""])
    for nm, col in ELEMS:
        ps = [num(r, col) for r in D1]
        ps = [p for p in ps if p is not None]
        line(f"{nm} @{alpha:.2f}", sum(1 for p in ps if p <= alpha + EPS), len(ps), alpha)

print("\n### one-sided shares e(perm_p_lo) / e(perm_p_hi), alpha = .05\n")
hdr(["element / side", "size", "MC se", "n", "dev (MC se)", ""])
for nm, base in [("main_att", "main"), ("pretreatment_att", "pre"), ("ES_att_pre2", "pre2"),
                 ("ES_att_pre1", "pre1"), ("ES_att_post0", "post0"), ("ES_att_post1", "post1")]:
    for side in ("lo", "hi"):
        ps = [num(r, f"p{side}_{base}") for r in D1]
        ps = [p for p in ps if p is not None]
        line(f"{nm} p_{side} @.05", sum(1 for p in ps if p <= 0.05 + EPS), len(ps), 0.05)

print("\n### alternative conventions, main_att (computed on every D1 rep)\n")
hdr(["convention", "size", "MC se", "n", "dev (MC se)", ""])
for alpha in (0.05, 0.10):
    for nm, col in (("equal-tailed (default)", "p_main"), ("permabsp", "cpabs"),
                    ("permrawp", "cpraw"), ("permside(left)", "cpleft")):
        ps = [num(r, col) for r in D1]
        ps = [p for p in ps if p is not None]
        line(f"{nm} @{alpha:.2f}", sum(1 for p in ps if p <= alpha + EPS), len(ps), alpha)

# mean of the p-value: on the NDRAW-draw lattice the equal-tailed p has mean (NDRAW+3)/(2(NDRAW+1))
# (.505 at 199) and SD SD_EQ (.2887 at 199), the one-sided share mean (NDRAW+2)/(2(NDRAW+1)) - a
# global calibration check that uses the whole distribution rather than only the 5% tail.
print(f"\n### mean p-value (global check; theory at permreps {NDRAW}: equal-tailed {MEAN_EQ:.4f}, one-sided {MEAN_LO:.4f}, SD {SD_EQ:.4f})\n")
hdr(["element", "mean equal-tailed p", "MC se", "dev (MC se)", "mean p_lo", "n"])
for nm, col, lo in [("main_att", "p_main", "plo_main"), ("pretreatment_att", "p_pre", "plo_pre"),
                    ("ES_att_pre2", "p_pre2", "plo_pre2"), ("ES_att_pre1", "p_pre1", "plo_pre1"),
                    ("ES_att_post0", "p_post0", "plo_post0"), ("ES_att_post1", "p_post1", "plo_post1")]:
    v = [num(r, col) for r in D1 if num(r, col) is not None]
    w = [num(r, lo) for r in D1 if num(r, lo) is not None]
    if not v:
        continue
    m = sum(v)/len(v)
    se0 = SD_EQ/math.sqrt(len(v))
    z = (m - MEAN_EQ)/se0
    mw = sum(w)/len(w) if w else float("nan")
    print(f"| {nm} | {m:.4f} | {se0:.4f} | {z:+.2f} | {mw:.4f} | {len(v)} |")
    if abs(z) > 2:
        FLAGS.append(f"{nm}: mean equal-tailed p {m:.4f} vs theory {MEAN_EQ:.4f} ({z:+.2f} MC se, n={len(v)})")

print("\n### fixed-metric joint pre-trend test e(joint_horiz) (jointtests)\n")
hdr(["test", "size", "MC se", "n", "dev (MC se)", ""])
for alpha in (0.05, 0.10):
    ps = [num(r, "jointp") for r in D1]
    ps = [p for p in ps if p is not None]
    line(f"joint_horiz @{alpha:.2f}", sum(1 for p in ps if p <= alpha + EPS), len(ps), alpha)

# ------------------------------------------------------------------ step-down
print("\n## 3. D1 - multiplicity over the post-treatment horizons (e(perm_stepdown))\n")
print("Two post horizons (post0, post1). FWER = share of reps rejecting ANY horizon.\n")
have = [r for r in D1 if num(r, "psd0") is not None and num(r, "psd1") is not None]
if not have:
    print("FAIL: no step-down rows in D1 - was permstepdown passed to mc_ri2.do?")
    sys.exit(1)
for alpha in (0.05, 0.10):
    hdr(["quantity", "share", "MC se", "n", "dev (MC se)", ""])
    line(f"FWER step-down min(p_stepdown)<={alpha:.2f}",
         sum(1 for r in have if min(num(r, "psd0"), num(r, "psd1")) <= alpha + EPS), len(have), alpha)
    line(f"FWER single-step min(p_singlestep)<={alpha:.2f}",
         sum(1 for r in have if min(num(r, "pss0"), num(r, "pss1")) <= alpha + EPS), len(have), alpha)
    line(f"UNCORRECTED min(p_raw)<={alpha:.2f}",
         sum(1 for r in have if min(num(r, "praw0"), num(r, "praw1")) <= alpha + EPS), len(have), None,
         "no nominal: this is the multiplicity inflation being corrected")
    line(f"per-horizon p_raw post0<={alpha:.2f}",
         sum(1 for r in have if num(r, "praw0") <= alpha + EPS), len(have), alpha)
    line(f"per-horizon p_raw post1<={alpha:.2f}",
         sum(1 for r in have if num(r, "praw1") <= alpha + EPS), len(have), alpha)
    line(f"per-horizon p_stepdown post0<={alpha:.2f}",
         sum(1 for r in have if num(r, "psd0") <= alpha + EPS), len(have), None, "per-horizon adjusted")
    line(f"per-horizon p_stepdown post1<={alpha:.2f}",
         sum(1 for r in have if num(r, "psd1") <= alpha + EPS), len(have), None, "per-horizon adjusted")
    print()
# cross-check: studentizing every draw of a horizon by ONE constant (the SD of that
# horizon's null draws) preserves ranks, so the step-down's p_raw must equal the package's
# permabsp p for the same horizon. The permabsp verification arms give that for free.
print("Cross-check: p_raw (step-down, |t| share) vs the package permabsp p for the same")
print("horizon on the verification reps - studentizing by one constant preserves ranks,")
print("so these must agree exactly.\n")
hdr(["horizon", "n verified", "exact matches", "max abs diff"])
for base, col in (("ES_att_post0", ("praw0", "p_post0")), ("ES_att_post1", ("praw1", "p_post1"))):
    nn, nm2, mx = 0, 0, 0.0
    for r in by[("D1", "absp")]:
        if r["rc"] != "0":
            continue
        m = mainbykey.get((r["shard"], r["rep"]))
        if not m:
            continue
        a, b = num(m, col[0]), num(r, col[1])
        if a is None or b is None:
            continue
        nn += 1
        d = abs(a - b)
        mx = max(mx, d)
        if d <= 1e-12:
            nm2 += 1
    print(f"| {base} | {nn} | {nm2} | {mx:.3e} |")
    if nn and nm2 != nn:
        FLAGS.append(f"step-down p_raw vs permabsp on {base}: {nn-nm2} of {nn} disagree")

# uniformity of the p-value on the attainable grid (secondary diagnostic)
print("\n### uniformity of the equal-tailed p on 10 equal-probability bins (chi-square, df=9)\n")
hdr(["element", "chi2", "n", "p>.01 crit 21.67?"])
for nm, col in ELEMS:
    ps = [num(r, col) for r in D1 if num(r, col) is not None]
    if len(ps) < 10:
        continue
    bins = [0]*10
    for p_ in ps:
        b = min(9, max(0, math.ceil(p_*10 - EPS) - 1))
        bins[b] += 1
    e = len(ps)/10.0
    st = sum((c-e)**2/e for c in bins)
    print(f"| {nm} | {st:.2f} | {len(ps)} | {'ok' if st <= 21.666 else 'EXCEEDS'} |")
    if st > 21.666:
        FLAGS.append(f"{nm}: p-value grid chi-square {st:.1f} > 21.67 (df=9, .01) on n={len(ps)}")

# ------------------------------------------------------------------ shift interval
def covtable(phase, arm, tau, label):
    R = [r for r in by[(phase, arm)] if r["rc"] == "0"]
    if not R:
        print(f"\n### {label}: no rows\n")
        return {}
    print(f"\n### {label}  (n reps = {len(R)}, true effect on every post horizon = {tau})\n")
    hdr(["column", "coverage of tau", "MC se", "n non-blank", "dev (MC se)", "mean width; blank share"])
    widths = {}
    for nm, base, truth in (("main_att", "main", tau), ("ES_att_post0", "post0", tau),
                            ("ES_att_post1", "post1", tau)):
        live = [r for r in R if num(r, f"p_{base}") is not None]
        nb = [r for r in live if num(r, f"ll_{base}") is not None and num(r, f"ul_{base}") is not None]
        blank = len(live) - len(nb)
        hits = sum(1 for r in nb if num(r, f"ll_{base}") - EPS <= truth <= num(r, f"ul_{base}") + EPS)
        w = [num(r, f"ul_{base}") - num(r, f"ll_{base}") for r in nb]
        mw = sum(w) / len(w) if w else float("nan")
        widths[base] = {(r["shard"], r["rep"]): (num(r, f"ul_{base}") - num(r, f"ll_{base}")) for r in nb}
        bs = blank / len(live) if live else 0.0
        line(f"{label} {nm}", hits, len(nb), 0.95, f"width {mw:.3f}; blank {blank}/{len(live)} = {bs:.3f}")
        if blank:
            FLAGS.append(f"{label} {nm}: shift interval BLANK on {blank}/{len(live)} reps "
                         f"(null band excludes 0; e(perm_kernel) offcentre)")
    return widths

print("\n## 4. Shift-interval coverage e(perm_ci) at level(95)\n")
w0 = covtable("D1", "main", 0.0, "tau=0 (D1)")
w1 = covtable("D2", "main", 1.0, "tau=1 (D2, default)")
w3 = covtable("D3", "main", 3.0, "tau=3 (D3, default)")
wt = covtable("D4", "tstat", 1.0, "tau=1 (D4, permtstat)")

print("\n### paired width growth (same DGP seed and rseed across tau; ratio per rep)\n")
hdr(["comparison", "column", "n paired", "mean ratio", "median ratio"])
def pairw(wa, wb, la, lb):
    for base in ("main", "post0", "post1"):
        A, B = wa.get(base, {}), wb.get(base, {})
        keys = sorted(set(A) & set(B))
        if not keys:
            continue
        rr = sorted(A[k] / B[k] for k in keys if B[k] and B[k] > 0)
        if not rr:
            continue
        med = rr[len(rr) // 2] if len(rr) % 2 else 0.5 * (rr[len(rr)//2 - 1] + rr[len(rr)//2])
        print(f"| {la} / {lb} | {base} | {len(rr)} | {sum(rr)/len(rr):.3f} | {med:.3f} |")
pairw(w1, w0, "tau=1", "tau=0")
pairw(w3, w0, "tau=3", "tau=0")
pairw(wt, w0, "tau=1 permtstat", "tau=0")
pairw(wt, w1, "tau=1 permtstat", "tau=1 default")

# scale of the null distribution. On the default (coefficient) arms e(perm_sd) is in outcome
# units; under permtstat it is the SD of the null t-ratios, which should sit near 1 if the
# studentization is doing its job. The shift interval's width is (q_hi - q_lo) times the
# observed analytic se, so a null-t SD far from 1 is exactly what makes the permtstat
# interval wider or narrower than the coefficient-scale one.
print("\n### scale of the permutation null distribution e(perm_sd)\n")
hdr(["phase / arm", "mean sd main_att", "mean sd post0", "mean sd post1", "n"])
for phase, arm in (("D1", "main"), ("D2", "main"), ("D3", "main"), ("D4", "tstat")):
    R = [r for r in by[(phase, arm)] if r["rc"] == "0"]
    if not R:
        continue
    def msd(c):
        v = [num(r, c) for r in R if num(r, c) is not None]
        return sum(v)/len(v) if v else float("nan")
    print(f"| {phase} / {arm} | {msd('sd_main'):.3f} | {msd('sd_post0'):.3f} | {msd('sd_post1'):.3f} | {len(R)} |")
print()
print("D4's row is on the t-ratio scale (permtstat); the others are in outcome units.")

# permtstat health: were analytic SEs available on the post rows?
D4 = [r for r in by[("D4", "tstat")] if r["rc"] == "0"]
if D4:
    print("\n### permtstat health under cluster(i)\n")
    for base in ("main", "post0", "post1"):
        miss = sum(1 for r in D4 if num(r, f"ase_{base}") is None or (num(r, f"ase_{base}") or 0) <= 0)
        print(f"- {base}: analytic SE missing or non-positive on {miss} of {len(D4)} reps")

# pre-treatment p drift under a planted effect (documented widening)
print("\n## 5. Pre-treatment p-values under a planted post-treatment effect\n")
hdr(["phase (tau)", "mean p ES_att_pre2", "mean p ES_att_pre1", "mean p pretreatment_att", "n"])
for phase, arm, tau in (("D1", "main", 0), ("D2", "main", 1), ("D3", "main", 3)):
    R = [r for r in by[(phase, arm)] if r["rc"] == "0"]
    if not R:
        continue
    def mn(c):
        v = [num(r, c) for r in R if num(r, c) is not None]
        return sum(v) / len(v) if v else float("nan")
    print(f"| {phase} (tau={tau}) | {mn('p_pre2'):.3f} | {mn('p_pre1'):.3f} | {mn('p_pre'):.3f} | {len(R)} |")

print("\n## 6. Surfaces off nominal by more than 2 MC standard errors\n")
if FLAGS:
    for f in FLAGS:
        print(f"- {f}")
else:
    print("- none")
