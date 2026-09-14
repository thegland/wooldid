#!/usr/bin/env python3
"""analyze_advbs12.py - judge the mc_adv_bs12.do campaign: the paired SIZE ladder for the
BS-12 one-parameter adversarial null.

WHY THIS ONE MATTERS. BS-12 (2026-09-10) made the adversarial null the one-parameter line
THROUGH the reference period. Every adversarial size number on record - MC_REPORT_2026-08-24,
REPORT_ADVBOUND_2026-08-26, and the help - was measured on the free-intercept statistic
BS-12 replaced, and the ruling's own prerequisite simulation was never run. Until this
campaign is judged, the shipped size claims describe a statistic that no longer ships.

THE CLAIM. On the null arm (slope == 0) the adversarial p must be correctly sized: the
rejection count at .05 and .10 inside the central ~95% region of Binomial(n, attained
alpha), per config and per arm. The pm (permutation) arm is exact under the sharp null and
is the paired reference the other two are read against on the SAME datasets.

SLOPE != 0 IS NOT POWER. mc_adv_bs12.do plants SLOPE*(relyr+1) on treated units: zero at
the reference period, linear in relyr, i.e. a line THROUGH the reference period, which is
INSIDE the BS-12 null for every SLOPE. Those rows are therefore a SECOND NULL - the null
the old free-intercept statistic could not distinguish from the new one - not an
alternative. They are reported in their own section and are NOT scored, because the
campaign brief calls them the power arm while the driver's own geometry makes them size,
and the two readings imply opposite verdicts. That discrepancy is an author ruling, flagged
in the report rather than silently decided here.

Bands and thresholds are CONSTANTS at the top with their provenance. Conventions (load /
num / binom_region / chi2 grid / MARGINAL corridors / PASS-MARGINAL-FAIL / n beside every
rate) are the house ones from analyze_mc.py and analyze_advb.py.

usage: analyze_advbs12.py [outdir] [B] [--out <report.md>] [--M <reps-per-cell>]
  outdir  directory holding mc_adv_bs12_s*.csv (default: ./out next to this file)
  B       permreps()/reps() of the run; default: inferred from the ncomp column
  --out   write the report here instead of <outdir>/bs12_report.md
  --M     the planned reps per (config, arm); without it the observed count is used
Writes <outdir>/bs12_report.md, prints the report, and ends with the sentinel
ADVBS12_ANALYZE_OK / ADVBS12_ANALYZE_FAIL=<n> at column 1. Exit 1 iff any check FAILs.
"""
import csv
import glob
import math
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 1e-9

# ---------------------------------------------------------------- constants and provenance
# mc_adv_bs12.do's ACCEPTANCE block (2026-09-13) is the source; that block itself says "to be
# preregistered before the production run" and no PREREG_*.md carries it, so these are NOT
# PREREGISTERED in the sense campaigns 1-2 use the word - they are the driver author's
# written-in-advance acceptance, and are treated as binding here.
ALPHAS = (0.05, 0.10)   # driver ACCEPTANCE: "rejection COUNT at .05 and .10 inside
                        # Binomial(M, alpha)'s region" - with the ATTAINED lattice alpha
                        # substituted on the resampling arms (analyze_mc.py DEVIATION D1).
RC_MAX = 0.0            # driver ACCEPTANCE: "integrity: rc = 0 on every row".
ADVABSENT_MAX = 0       # driver ACCEPTANCE: "advabsent = 0 everywhere". A missing
                        # adversarial p at rc 0 means the statistic never formed on that rep.
ADVW_BIG_MAX_PM = 0     # driver ACCEPTANCE: "advw_big = 0 on PM". advw_big counts
                        # replication W statistics above 1e10 - KF-2's escaped 3.58e25 draw.
# House MARGINAL corridors for SIZE, verbatim from analyze_mc.py's M2 block.
MARG_05 = (0.025, 0.085)
MARG_10 = (0.070, 0.135)
# The driver's own negative control: an on AB12 is cluster(i) containment, ruled
# UNCALIBRATED on 2026-08-26 (REPORT_ADVBOUND_2026-08-26.md) and expected to fail. It is
# reported like analyze_mc.py's A5 - INFO, unscored - so a documented defect does not hold
# the sentinel down forever. If it unexpectedly passes, that is itself a finding and is said.
EXPECTED_UNCALIBRATED = {("AB12", "an"): "cluster(i) containment, ruled uncalibrated "
                                         "2026-08-26; driver header: 'expected to FAIL'"}
# Below MINN reps a cell is REPORTED, never judged - a 1-rep smoke must say "insufficient
# n", not "PASS". NOT PREREGISTERED; house analogue analyze_bootse2.py's MCB2_MINN = 15;
# 30 is used because binom_region(n, .05) cannot discriminate below roughly that.
MINN = int(os.environ.get("MCBS12_MINN", "30"))
CHI2_9_CRIT = 21.666    # df=9, p=.01 (house uniformity convention; secondary flag only)

# The p-value lattice, per arm. an is the analytic p (continuous, so the attained alpha IS
# alpha); pm is the permutation add-one p and bt the bootstrapped studentized p, both
# right-tailed on the positive W statistic, so both take values R/(B+1), R = 1..B+1.
# NOT PREREGISTERED and, for bt, NOT OBSERVED: the smoke has no bt row, so the add-one form
# is an assumption. It matters - at B = 199 an add-one p attains .05 exactly while a k/B p
# attains 9/199 = .0452, and the count band moves with it. Section 0's "lattice observed"
# FINDING measures the convention off the recorded p's and says so out loud, so a wrong
# assumption announces itself instead of silently mis-banding a verdict.
ARM_CONV = {"an": "continuous", "pm": "right", "bt": "right"}
ARM_WHAT = {"an": "analytic p, e(adversarial_pval_att)",
            "pm": "permutation p (sharp null, cohort level) - the exact reference",
            "bt": "boottstat p (bootstrapped W, containment cell)"}
REQUIRED = ["config", "arm", "shard", "rep", "seed", "slope", "rc", "secs", "p", "pjoint",
            "n", "ncomp", "nfail", "advabsent", "advw_n", "advw_big", "advw_max", "tstud_n"]
# Live on every row; all-missing here means the family is VACUOUS, not zero.
LIVE_ALWAYS = ["p", "pjoint", "n", "advabsent"]
# Live only on the resampling arms (an posts no replications and no W accumulator).
LIVE_RESAMPLE = ["ncomp", "nfail", "advw_n", "advw_big", "advw_max"]
LIVE_BT = ["tstud_n"]   # e(boot_tstud) row n for main_att - bootstrap-studentized only


def _parse_args(argv):
    """positional [outdir] [B] plus --out <path> and --M <n> (analyze_advb.py's parser)."""
    pos, out, M, i = [], None, None, 0
    while i < len(argv):
        a = argv[i]
        if a in ("--out", "--M"):
            v = argv[i + 1] if i + 1 < len(argv) else None
            out, M = (v, M) if a == "--out" else (out, v)
            i += 2
            continue
        if a.startswith("--out="):
            out = a[len("--out="):]
            i += 1
            continue
        if a.startswith("--M="):
            M = a[len("--M="):]
            i += 1
            continue
        pos.append(a)
        i += 1
    return pos, out, (int(M) if M else None)


_POS, REPORT_OUT, M_TARGET = _parse_args(sys.argv[1:])
OUT = _POS[0] if len(_POS) > 0 else os.path.join(HERE, "out")
B_ARG = int(_POS[1]) if len(_POS) > 1 else None

L = []
FAILS = []
SHORTN = []
NOTJUDGED = []


def num(row, k):
    """Stata writes missing as '.': missing is None, never 0."""
    v = row.get(k, ".")
    if v in (".", "", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def load(stem):
    """Concatenate the shard CSVs; drop a short final line from a shard killed mid-write
    (analyze_bootse2.py) and de-duplicate on (config, arm, slope, shard, rep). SLOPE is part
    of the key because the driver APPENDS to <outdir>/mc_adv_bs12_s<shard>.csv and its header
    says to "run BOTH" slopes: a second run at another slope into the same outdir writes rows
    with the same (config, arm, shard, rep) - the seed is slope-independent - and must NOT be
    discarded as duplicates, which would empty the sloped-null section silently."""
    files = sorted(glob.glob(os.path.join(OUT, f"{stem}_s*.csv")))
    rows, nbad, head = [], 0, None
    for f in files:
        with open(f) as fh:
            rdr = csv.reader(fh)
            try:
                h = next(rdr)
            except StopIteration:
                continue
            head = h if head is None else head
            for raw in rdr:
                if len(raw) != len(h):
                    nbad += 1
                    continue
                rows.append(dict(zip(h, raw)))
    seen, keep, ndup = set(), [], 0
    for r in rows:
        k = (r.get("config"), r.get("arm"), r.get("slope"), r.get("shard"), r.get("rep"))
        if k in seen:
            ndup += 1
            continue
        seen.add(k)
        keep.append(r)
    return keep, head, files, nbad, ndup


def binom_region(n, p):
    """Central ~95% region, normal approx rounded outward (house prereg convention)."""
    mu, sd = n * p, math.sqrt(n * p * (1 - p))
    return max(0, math.floor(mu - 1.96 * sd)), min(n, math.ceil(mu + 1.96 * sd))


def wilson(h, n, z=1.959963985):
    """95% Wilson score interval for a binomial share (stdlib, no scipy)."""
    if n <= 0:
        return float("nan"), float("nan")
    p = h / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - hw), min(1.0, c + hw)


def rate(h, n):
    """A rate is NEVER printed without its n and its interval."""
    if n <= 0:
        return "-/0 (no reps)"
    lo, hi = wilson(h, n)
    return f"{h}/{n} = {h/n:.4f} [{lo:.4f},{hi:.4f}]"


def cnt(h, n):
    """A defect COUNT: h of n. No interval - it is not an estimate of anything."""
    return f"{h}/{n}"


def chi2_grid(ps):
    """10 equal-probability bins on (0,1]; (stat, flag). House convention, secondary."""
    n = len(ps)
    if n == 0:
        return 0.0, False
    bins = [0] * 10
    for p in ps:
        bins[min(9, max(0, math.ceil(p * 10 - EPS) - 1))] += 1
    e = n / 10.0
    stat = sum((c - e) ** 2 / e for c in bins)
    return stat, stat > CHI2_9_CRIT


def attained(alpha, B, conv):
    """P(p <= alpha) the test can actually attain. Continuous p: alpha. Right-tailed
    add-one p on B replications: floor(alpha(B+1))/(B+1) - at B = 199 that is .05 and .10
    exactly; at small B it can be 0, i.e. the claim is UNATTAINABLE and is not judged
    (analyze_mc.py's R6 rule)."""
    if conv == "continuous" or B is None:
        return alpha
    return math.floor(alpha * (B + 1) + EPS) / (B + 1.0)


def observed_lattice(vals, B):
    """Which right-tailed lattice do these p's actually lie on (analyze_riext.py's
    on_lattice, run in reverse)? Returns the labels consistent with EVERY value: 'add-one
    k/(B+1)', 'plain k/B', or neither. Reported, never scored - it is evidence about the
    engine's convention, which decides the attained alpha this analyzer bands against."""
    if not vals or B is None or B <= 0:
        return []
    out = []
    for label, den in (("add-one k/(B+1)", B + 1.0), ("plain k/B", float(B))):
        if all(abs(v * den - round(v * den)) < 1e-6 for v in vals):
            out.append(label)
    return out


def verdict(name, ok, marginal=None, detail="", n=None):
    """PASS / MARGINAL / FAIL, plus INSUFFICIENT_N below MINN. Only FAIL moves the
    sentinel."""
    if n is not None and n < MINN:
        SHORTN.append(f"{name}: n={n} < MINN={MINN}")
        L.append(f"| {name} | INSUFFICIENT_N | {detail} |")
        return "INSUFFICIENT_N"
    if ok:
        v = "PASS"
    elif marginal:
        v = "MARGINAL"
    else:
        v = "FAIL"
        FAILS.append(f"{name}: {detail}")
    L.append(f"| {name} | {v} | {detail} |")
    return v


def fail(name, detail):
    FAILS.append(f"{name}: {detail}")
    L.append(f"| {name} | FAIL | {detail} |")


def finish():
    txt = "\n".join(L) + "\n"
    report = REPORT_OUT if REPORT_OUT else os.path.join(OUT, "bs12_report.md")
    try:
        if os.path.dirname(report):
            os.makedirs(os.path.dirname(report), exist_ok=True)
        with open(report, "w") as fh:
            fh.write(txt)
        wrote = f"\nwrote {report}"
    except OSError as e:
        wrote = f"\ncould not write the report ({e})"
    print(txt + wrote)
    if NOTJUDGED:
        print(f"NOT_JUDGED={len(NOTJUDGED)} (alpha unattainable on this run's p lattice)")
    if SHORTN:
        print(f"INSUFFICIENT_N={len(SHORTN)} (checks not judged: too few reps; MINN={MINN})")
    if FAILS:
        print(f"ADVBS12_ANALYZE_FAIL={len(FAILS)}")
        sys.exit(1)
    print("ADVBS12_ANALYZE_OK")
    sys.exit(0)


# ------------------------------------------------------------------ load
rows, head, files, nbad, ndup = load("mc_adv_bs12")
L.append("# BS-12 adversarial size ladder (mc_adv_bs12)\n")
L.append(f"outdir `{OUT}`; {len(files)} shard file(s); {len(rows)} rows after de-duplication "
         f"on (config, arm, slope, shard, rep); {nbad} malformed line(s) dropped; {ndup} "
         "duplicate "
         "row(s) dropped.\n")
L.append("Every adversarial size number on record was measured on the FREE-INTERCEPT "
         "statistic that BS-12 replaced on 2026-09-10 (MC_REPORT_2026-08-24, "
         "REPORT_ADVBOUND_2026-08-26, the help). This campaign is the re-measurement; until "
         "it is judged, those numbers describe a statistic that no longer ships.\n")
L.append("Headline: the rejection rate of `p` at nominal .05 and .10 on the NULL arm "
         "(slope == 0), per config and per arm, against the central ~95% binomial region. "
         "The `pm` arm is exact under the sharp null and is the paired reference.\n")
L.append(f"Cells with fewer than MINN={MINN} reps are reported, never judged.\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")

if not rows:
    fail("rows present", f"no rows in {OUT}/mc_adv_bs12_s*.csv - the campaign produced "
                         f"nothing (searched {len(files)} file(s))")
    finish()

missing_cols = [c for c in REQUIRED if head is None or c not in head]
if missing_cols:
    fail("required columns", "the CSV header is missing " + ", ".join(missing_cols) +
         " - this analyzer cannot judge a file it cannot read; the driver's header changed")
    finish()

bad_rc_rows = [r for r in rows if r.get("rc") not in ("0", 0)]
verdict("rc == 0 on every row", len(bad_rc_rows) / len(rows) <= RC_MAX + EPS,
        detail=f"rc!=0 on {cnt(len(bad_rc_rows), len(rows))}; threshold {RC_MAX:.2f} "
               "(driver ACCEPTANCE: rc = 0 on every row)")
ok_rows = [r for r in rows if r.get("rc") in ("0", 0)]
if not ok_rows:
    fail("usable rows", "every row has rc != 0: no adversarial statistic ever formed")
    finish()

NULL = [r for r in ok_rows if (num(r, "slope") or 0) == 0]
SLOPED = [r for r in ok_rows if (num(r, "slope") or 0) != 0]
cells = defaultdict(list)
for r in NULL:
    cells[(r.get("config", "?"), r.get("arm", "?"))].append(r)
CELLS = sorted(cells)

ncomps = [int(v) for v in (num(r, "ncomp") for r in ok_rows) if v is not None]
B_MODE = Counter(ncomps).most_common(1)[0][0] if ncomps else None
B = B_ARG if B_ARG is not None else B_MODE
if B_ARG is not None and B_MODE is not None and B_ARG != B_MODE:
    L.append(f"| B argument vs data | WARN | B given as {B_ARG} but the modal ncomp is "
             f"{B_MODE}; judging the resampling lattices at {B_ARG} |")
if B is None:
    L.append("| replication count | NOTE | no ncomp values and no B argument: the pm/bt "
             "lattices cannot be established, so their attained alpha falls back to the "
             "nominal alpha and the completeness check is skipped |")
else:
    L.append(f"| replication count | NOTE | B = {B} "
             f"({'argument' if B_ARG is not None else 'modal ncomp'}); attained alpha on the "
             "right-tailed add-one lattice " +
             ", ".join(f"{a:.2f}->{attained(a, B, 'right'):.4f}" for a in ALPHAS) + " |")

if not NULL:
    fail("null arm present", f"no rows with slope == 0 in {len(ok_rows)} rc-0 rows - this "
                             "campaign has no SIZE arm at all, only the sloped null "
                             "(slope values: "
                             f"{sorted({r.get('slope') for r in SLOPED})}). The size ladder "
                             "is the point of the driver; re-run with slope 0")

# all-missing columns, per cell (the vacuous-family trap)
for cfg, arm in CELLS:
    g = cells[(cfg, arm)]
    expect = list(LIVE_ALWAYS)
    if arm in ("pm", "bt"):
        expect += LIVE_RESAMPLE
    if arm == "bt":
        expect += LIVE_BT
    dead = [c for c in expect if all(num(r, c) is None for r in g)]
    if dead:
        fail(f"{cfg}/{arm} columns all-missing",
             "EVERY value is missing in " + ", ".join(dead) + f" on {len(g)} rc-0 null rows "
             "- the driver wrote nothing for this family, so any rate computed from it "
             "would be VACUOUS, not zero. Check the e() surface the arm reads")
    if arm == "an":
        unexpected = [c for c in LIVE_RESAMPLE + LIVE_BT
                      if any(num(r, c) is not None for r in g)]
        if unexpected:
            L.append(f"| {cfg}/{arm} unexpected values | NOTE | " + ", ".join(unexpected) +
                     " carries values on the analytic arm, which runs no replications |")

# ------------------------------------------------------------------ inventory
L.append("\n## Run inventory\n")
L.append("| config | arm | what | slope | rows (rc=0) | rc!=0 | N (modal) | ncomp (modal) | "
         "nfail>0 | advabsent | mean secs |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
inv = defaultdict(list)
for r in ok_rows:
    inv[(r.get("config", "?"), r.get("arm", "?"), r.get("slope", "?"))].append(r)
for (cfg, arm, sl) in sorted(inv):
    g = inv[(cfg, arm, sl)]
    badrc = sum(1 for r in rows if r.get("config") == cfg and r.get("arm") == arm
                and r.get("slope") == sl and r.get("rc") not in ("0", 0))
    nn = [int(v) for v in (num(r, "n") for r in g) if v is not None]
    nc = [int(v) for v in (num(r, "ncomp") for r in g) if v is not None]
    secs = [num(r, "secs") for r in g if num(r, "secs") is not None]
    L.append(f"| {cfg} | {arm} | {ARM_WHAT.get(arm, '?')} | {sl} | {len(g)} | {badrc} | "
             f"{Counter(nn).most_common(1)[0][0] if nn else '-'} | "
             f"{Counter(nc).most_common(1)[0][0] if nc else '-'} | "
             f"{sum(1 for r in g if (num(r, 'nfail') or 0) > 0)} | "
             f"{sum(1 for r in g if (num(r, 'advabsent') or 0) > 0)} | "
             f"{(sum(secs)/len(secs) if secs else float('nan')):.2f} |")
missing_cells = [k for k in [("AB1", "an"), ("AB13", "an"), ("AB12", "an"), ("AB12", "bt"),
                             ("PM", "pm"), ("AB1C", "an")] if k not in cells]
if missing_cells:
    L.append("\nRungs of the ladder with NO null rows in this outdir: " +
             ", ".join(f"{c}/{a}" for c, a in missing_cells) +
             " (a partial run, or the geometry list did not include them - not judged).\n")

# ------------------------------------------------------------------ 0. integrity
L.append("\n## 0. Integrity (per-row facts; judged at any n - these are not statistical)\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")
for cfg, arm in CELLS:
    g = cells[(cfg, arm)]
    n = len(g)
    absent = sum(1 for r in g if (num(r, "advabsent") or 0) > 0)
    verdict(f"{cfg}/{arm} advabsent == 0", absent <= ADVABSENT_MAX,
            detail=f"advabsent>0 on {cnt(absent, n)} (the adversarial p is missing at rc 0, "
                   "or e(boot_adv_absent) is set: the statistic never formed on that rep, so "
                   "it contributes no test - a rate computed without it is conditional)")
    oob = [num(r, "p") for r in g
           if num(r, "p") is not None and not (-EPS <= num(r, "p") <= 1 + EPS)]
    verdict(f"{cfg}/{arm} p in [0,1]", not oob,
            detail=f"{len(oob)} p value(s) outside [0,1], e.g. {oob[:3]}" if oob
                   else "no p outside [0,1]")
    joob = [num(r, "pjoint") for r in g
            if num(r, "pjoint") is not None and not (-EPS <= num(r, "pjoint") <= 1 + EPS)]
    verdict(f"{cfg}/{arm} pjoint in [0,1]", not joob,
            detail=(f"{len(joob)} e(joint_horiz) value(s) outside [0,1], e.g. {joob[:3]} - a "
                    "'p' above 1 is a broken statistic, not a calibration finding. NOTE: "
                    "the driver header's report of 3.34 on AB12 at slope .15 was REFUTED "
                    "2026-09-13 - the value is 3.34387233470e-07 and the exponent was "
                    "dropped reading a %18.12g field by eye; no code path can emit a p "
                    "above 1. Keep this check anyway: it is cheap and it is exactly the "
                    "kind of claim that should be verified by a parser, not an eye.") if joob
                   else "no e(joint_horiz) outside [0,1]")
    if arm in ("pm", "bt"):
        if B is not None:
            inc = sum(1 for r in g if (num(r, "ncomp") or -1) != B)
            verdict(f"{cfg}/{arm} ncomp == B ({B})", inc == 0,
                    detail=f"ncomp != {B} on {cnt(inc, n)} (lost replications put that rep's "
                           "p on a coarser lattice than the one judged here)")
        nf = sum(1 for r in g if (num(r, "nfail") or 0) > 0)
        verdict(f"{cfg}/{arm} nfail == 0", nf == 0,
                detail=f"nfail>0 on {cnt(nf, n)} (a failed replication means the estimate "
                       "never formed on that draw)")
        ps_ = [num(r, "p") for r in g if num(r, "p") is not None]
        lat = observed_lattice(ps_, B)
        L.append(f"| {cfg}/{arm} lattice observed | FINDING | {len(ps_)} p value(s) are "
                 f"consistent with: {', '.join(lat) if lat else 'NEITHER k/(B+1) nor k/B'} "
                 f"at B={B}. This analyzer bands against the add-one form "
                 f"({ARM_CONV.get(arm)}); if the only consistent label is 'plain k/B', the "
                 "attained alpha is floor(alpha*B)/B and the count bands in section 1 are "
                 "off by one rank - change ARM_CONV and re-run. Reported, never scored |")
        big = sum(1 for r in g if (num(r, "advw_big") or 0) > 0)
        wmax = max([num(r, "advw_max") for r in g if num(r, "advw_max") is not None],
                   default=float("nan"))
        if cfg == "PM":
            verdict(f"{cfg}/{arm} advw_big == 0", big <= ADVW_BIG_MAX_PM,
                    detail=f"reps with a replication W above 1e10: {cnt(big, n)}; max W seen "
                           f"{wmax:.4g} (driver ACCEPTANCE: advw_big = 0 on PM)")
        else:
            L.append(f"| {cfg}/{arm} advw_big | FINDING | reps with a replication W above "
                     f"1e10: {cnt(big, n)}; max W seen {wmax:.4g}. On a containment rung this "
                     "is REPORTED, not scored (driver ACCEPTANCE) - it is KF-2's escape "
                     "route, where a single degenerate draw carries the whole distribution |")
    if M_TARGET is not None:
        verdict(f"{cfg}/{arm} row count == M ({M_TARGET})", n == M_TARGET,
                detail=f"{n} rc-0 null rows of {M_TARGET} planned")
    else:
        L.append(f"| {cfg}/{arm} row count | NOTE | {n} rc-0 null rows; no --M target given, "
                 "so completeness of the campaign is not checked |")

# ------------------------------------------------------------------ 1. size ladder
L.append("\n## 1. SIZE on the null arm (slope == 0) - the headline\n")
L.append("(`attained` is P(p <= alpha) the arm's lattice can reach: alpha itself for the "
         "analytic p, floor(alpha(B+1))/(B+1) for the add-one permutation and bootstrap p. "
         "chi2 is the 10-bin uniformity statistic, df=9, crit 21.666 at p=.01 - secondary.)\n")
L.append("| config | arm | conv | alpha | attained | rej (n, Wilson 95%) | count band | chi2 |")
L.append("|---|---|---|---|---|---|---|---|")
claims = []
for cfg, arm in CELLS:
    g = cells[(cfg, arm)]
    conv = ARM_CONV.get(arm, "continuous")
    ps = [num(r, "p") for r in g if num(r, "p") is not None]
    n = len(ps)
    stat, flag = chi2_grid(ps) if n >= 10 else (float("nan"), False)
    cc = []
    for a in ALPHAS:
        att = attained(a, B, conv)
        h = sum(1 for p in ps if p <= a + EPS)
        lo, hi = binom_region(n, att)
        cc.append((a, att, h, n, lo, hi))
        L.append(f"| {cfg} | {arm} | {conv} | {a:.2f} | {att:.4f} | {rate(h, n)} | "
                 f"[{lo},{hi}] | {stat:.1f}{'*' if flag else ''} |")
    claims.append((cfg, arm, conv, cc, stat, flag, n))

L.append("\n| check | verdict | detail |")
L.append("|---|---|---|")
for cfg, arm, conv, cc, stat, flag, n in claims:
    bad, unjud, det = [], [], []
    for a, att, h, nn, lo, hi in cc:
        det.append(f"alpha {a:.2f} (attained {att:.4f}): {rate(h, nn)} band [{lo},{hi}]")
        if att <= EPS:
            unjud.append(f"alpha {a:.2f} is UNATTAINABLE at B={B} on the {conv} lattice")
            continue
        if not (lo <= h <= hi):
            bad.append(f"alpha {a:.2f}")
    if unjud:
        NOTJUDGED.extend(f"{cfg}/{arm}: {u}" for u in unjud)
    detail = "; ".join(det + unjud) + (f"; chi2={stat:.1f}{'*' if flag else ''}"
                                       if stat == stat else "")
    if (cfg, arm) in EXPECTED_UNCALIBRATED:
        # analyze_mc.py's A5 treatment: report-only, so a documented defect cannot hold the
        # sentinel down. An unexpected PASS here is itself a finding and is said out loud.
        L.append(f"| {cfg}/{arm} size (report-only) | INFO | {detail}. NOT SCORED: "
                 f"{EXPECTED_UNCALIBRATED[(cfg, arm)]} |")
        if not bad and n >= MINN:
            L.append(f"| {cfg}/{arm} unexpectedly in band | FINDING | the negative control "
                     "sits INSIDE its binomial region at every alpha. Either BS-12 fixed the "
                     "containment cell the 2026-08-26 ruling called uncalibrated, or the "
                     "geometry is not exercising it - worth an author look either way |")
        continue
    rt = {a: (h / nn if nn else float("nan")) for a, att, h, nn, lo, hi in cc}
    marg = (MARG_05[0] <= rt.get(0.05, -1) <= MARG_05[1]
            and MARG_10[0] <= rt.get(0.10, -1) <= MARG_10[1])
    verdict(f"{cfg}/{arm} size", not bad, marginal=marg, detail=detail, n=n)
    if flag:
        L.append(f"| {cfg}/{arm} uniformity (sec) | FLAG | 10-bin chi2 = {stat:.1f} > "
                 f"{CHI2_9_CRIT} on n={n}: p is not uniform even where the tail rate is in "
                 "band - investigate, never silent |")

# ------------------------------------------------------------------ 2. paired ladder
L.append("\n## 2. Paired against the exact reference (same datasets, same reps)\n")
L.append("(PM/pm is the permutation p under the sharp null: exact by construction. Every "
         "other cell is run on the SAME dataset per rep, so the comparison is paired and "
         "free of dataset noise. Discordance is where the two arms disagree at alpha = .05.)\n")
ref = {(r.get("shard"), r.get("rep")): r for r in cells.get(("PM", "pm"), [])}
if not ref:
    L.append("No PM/pm null rows in this outdir: the paired reference did not run, so the "
             "ladder below is unpaired and only the absolute rates in section 1 apply.\n")
else:
    L.append("| config | arm | n paired | arm rej@.05 | pm rej@.05 | arm only | pm only |")
    L.append("|---|---|---|---|---|---|---|")
    for cfg, arm in CELLS:
        if (cfg, arm) == ("PM", "pm"):
            continue
        g = cells[(cfg, arm)]
        pair = [(r, ref[(r.get("shard"), r.get("rep"))]) for r in g
                if (r.get("shard"), r.get("rep")) in ref]
        pair = [(a, b) for a, b in pair if num(a, "p") is not None and num(b, "p") is not None]
        if not pair:
            L.append(f"| {cfg} | {arm} | 0 | - | - | - | - |")
            continue
        ha = sum(1 for a, b in pair if num(a, "p") <= 0.05 + EPS)
        hb = sum(1 for a, b in pair if num(b, "p") <= 0.05 + EPS)
        ao = sum(1 for a, b in pair if num(a, "p") <= 0.05 + EPS and num(b, "p") > 0.05 + EPS)
        bo = sum(1 for a, b in pair if num(b, "p") <= 0.05 + EPS and num(a, "p") > 0.05 + EPS)
        L.append(f"| {cfg} | {arm} | {len(pair)} | {rate(ha, len(pair))} | "
                 f"{rate(hb, len(pair))} | {ao} | {bo} |")
    L.append("\nThe pairing is the point of the ladder: an analytic arm that rejects where "
             "the exact permutation p does not is over-rejecting on THESE datasets, not on "
             "an incomparable draw.\n")

# ------------------------------------------------------------------ 3. the sloped null
L.append("\n## 3. slope != 0 - reported separately, NOT scored\n")
if not SLOPED:
    L.append("No slope != 0 rows in this outdir.\n")
else:
    L.append("**Read this before reading the table.** mc_adv_bs12.do plants "
             "SLOPE*(relyr+1) on treated units. That is zero at the reference period and "
             "linear in relyr, i.e. a line THROUGH the reference period - which is INSIDE "
             "the BS-12 null for every SLOPE. So these rows are a SECOND NULL (the null the "
             "old free-intercept statistic could not distinguish from the new one), not an "
             "alternative, and their rejection rates are SIZE numbers that belong in the "
             "same bands as section 1. The campaign brief calls this the power arm. The two "
             "readings imply opposite verdicts on the same numbers, so this section is "
             "INFO-only and the discrepancy is an author ruling, not a thing this analyzer "
             "decides.\n")
    L.append("| config | arm | slope | alpha | attained | rej (n, Wilson 95%) | "
             "count band if read as SIZE |")
    L.append("|---|---|---|---|---|---|---|")
    sl_cells = defaultdict(list)
    for r in SLOPED:
        sl_cells[(r.get("config", "?"), r.get("arm", "?"), r.get("slope", "?"))].append(r)
    for k in sorted(sl_cells):
        cfg, arm, sl = k
        g = sl_cells[k]
        conv = ARM_CONV.get(arm, "continuous")
        ps = [num(r, "p") for r in g if num(r, "p") is not None]
        for a in ALPHAS:
            att = attained(a, B, conv)
            h = sum(1 for p in ps if p <= a + EPS)
            lo, hi = binom_region(len(ps), att)
            L.append(f"| {cfg} | {arm} | {sl} | {a:.2f} | {att:.4f} | {rate(h, len(ps))} | "
                     f"[{lo},{hi}] |")

# ------------------------------------------------------------------ 4. secondary
L.append("\n## 4. Secondary (information only, never fails the component)\n")
L.append("| config | arm | e(joint_horiz) rej@.05 | mean p | advw_n (modal) | max advw_max | "
         "tstud_n (modal) |")
L.append("|---|---|---|---|---|---|---|")
for cfg, arm in CELLS:
    g = cells[(cfg, arm)]
    pj = [num(r, "pjoint") for r in g if num(r, "pjoint") is not None]
    j05 = sum(1 for p in pj if p <= 0.05 + EPS)
    ps = [num(r, "p") for r in g if num(r, "p") is not None]
    wn = [int(v) for v in (num(r, "advw_n") for r in g) if v is not None]
    wm = [num(r, "advw_max") for r in g if num(r, "advw_max") is not None]
    tn = [int(v) for v in (num(r, "tstud_n") for r in g) if v is not None]
    L.append(f"| {cfg} | {arm} | {rate(j05, len(pj))} | "
             f"{(sum(ps)/len(ps) if ps else float('nan')):.4f} | "
             f"{Counter(wn).most_common(1)[0][0] if wn else '-'} | "
             f"{max(wm) if wm else float('nan'):.4g} | "
             f"{Counter(tn).most_common(1)[0][0] if tn else '-'} |")
L.append("\n`advw_n` is the count of usable replication W statistics behind that rep's p and "
         "`advw_max` their maximum; `tstud_n` is e(boot_tstud)'s n for main_att. A mean p "
         "far from the lattice mean, or an advw_max in the 1e10s, says the tail is one draw.\n")

L.append("\n## Verdict summary\n")
L.append(f"- checks FAILED: {len(FAILS)}")
for f_ in FAILS:
    L.append(f"  - {f_}")
L.append(f"- claims NOT JUDGED (alpha unattainable at B={B}): {len(NOTJUDGED)}")
for x in NOTJUDGED:
    L.append(f"  - {x}")
L.append(f"- checks not judged for want of reps (MINN={MINN}): {len(SHORTN)}")
for s_ in SHORTN:
    L.append(f"  - {s_}")

finish()
