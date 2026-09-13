#!/usr/bin/env python3
"""analyze_bands.py - judge the mc_bands.do campaign: COVERAGE of the simultaneous
sup-t band (simultaneousbands, PS-45) on the analytic-simulated and the bootstrap
engines, including a cluster strictly BELOW i.

THE CLAIM. The band promises SIMULTANEOUS coverage: with probability level()/100 EVERY
reported horizon's truth lies inside its band. That event is the CSV's band_all. The
pointwise intervals of the same run (pw_*) are the CONTROL: each horizon separately
should cover at ~nominal, and the JOINT pointwise event (pw_all) must sit BELOW the
simultaneous one - that gap is the multiplicity the band exists to close. mc_bands.do
plants a constant tau, so with esfixedbaseperiod every reported horizon's truth is known
exactly (0 on pre horizons, tau on post), which is what makes a hit computable per rep.

Bands and thresholds are CONSTANTS at the top with their provenance; read them before
reading a verdict. Conventions (load / num / binom_region / chi2 grid / PASS-MARGINAL-FAIL
/ n beside every rate) are the house ones from analyze_mc.py, analyze_advb.py and
analyze_bootse2.py.

usage: analyze_bands.py [outdir] [B] [--out <report.md>] [--M <reps-per-arm>]
  outdir  directory holding mc_bands_s*.csv (default: ./out next to this file)
  B       bootstrap reps() of the run; default: inferred from the ncomp column
  --out   write the report here instead of <outdir>/bands_report.md
  --M     the planned reps per arm; without it the observed row count is used and the
          completeness check only warns
Writes <outdir>/bands_report.md, prints the report, and ends with the sentinel
BANDS_ANALYZE_OK / BANDS_ANALYZE_FAIL=<n> at column 1. Exit 1 iff any check FAILs.
"""
import csv
import glob
import math
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 1e-9

# ---------------------------------------------------------------- preregistered constants
# PROVENANCE. mc_bands.do's ACCEPTANCE block (2026-09-13) states these; that block itself
# says "to be preregistered before the production run" and NO PREREG_*.md carries them
# (PREREG_2026-09-13.md does not mention the band). So: NOT PREREGISTERED in the strict
# sense of campaigns 1-2; they are the driver author's written-in-advance acceptance and
# are treated as binding here. Any change is an amendment that belongs in this file.
NOMINAL = 0.95          # mc_bands.do passes no level(), so wooldid's default level(95)
                        # is the band's nominal simultaneous coverage. NOT PREREGISTERED
                        # (driver ACCEPTANCE: "the band's nominal claim").
NH_EXPECTED = 8         # esprelength(4) espostlength(4) esfixedbaseperiod -> relyr
                        # -5..-2 and 0..3. Driver ACCEPTANCE: "nh = 8 (the AB1 window)".
RC_MAX = 0.0            # share of rows with rc != 0 tolerated. Driver ACCEPTANCE:
                        # "integrity: rc = 0 ... on every row" -> zero tolerance.
BANDABSENT_MAX = 0      # count of rc-0 rows with no e(esband_att). Driver ACCEPTANCE:
                        # "band_absent = 0 on every row".
# Judgment band for a count: the house central ~95% binomial region, normal approx rounded
# outward (analyze_mc.py binom_region, "prereg convention"). The reported interval beside
# every rate is a Wilson 95% score interval (exact-ish, closed form, stdlib).
# MARGINAL corridor for a COVERAGE claim: nominal -/+ 5 points. NOT PREREGISTERED - the
# house corridors ([.025,.085] / [.070,.135]) are written for SIZE at .05/.10 and do not
# transfer to a .95 coverage claim. 5 points is the coarsest departure this campaign would
# still call "documented, not broken"; a reader who disagrees should move this constant.
COV_MARGINAL = (NOMINAL - 0.05, min(1.0, NOMINAL + 0.03))
# Below MINN reps a cell is REPORTED, never judged: a 1-rep smoke must say "insufficient n",
# not "PASS". NOT PREREGISTERED; the house analogue is analyze_bootse2.py's MCB2_MINN = 15.
# 30 is used here because binom_region(n, .95) is [n, n] or wider than the whole unit
# interval below ~25, i.e. it cannot discriminate. Override with MCBANDS_MINN.
MINN = int(os.environ.get("MCBANDS_MINN", "30"))
CHI2_9_CRIT = 21.666    # df=9, p=.01 (house uniformity convention; secondary only)

HORIZONS = ["m5", "m4", "m3", "m2", "0", "1", "2", "3"]
BASECOLS = ["arm", "shard", "rep", "seed", "tau", "rc", "secs", "engine",
            "crit", "nband", "nh", "band_absent", "band_all", "band_hits",
            "pw_all", "pw_hits", "pw_n", "ncomp", "nfail"]
REQUIRED = BASECOLS + [f"bd{h}" for h in HORIZONS] + [f"pw{h}" for h in HORIZONS]

# Which columns must carry values, by arm. A column that is entirely missing where it is
# expected live means the driver wrote nothing and the family is VACUOUS - that is the
# defect that emptied one existing driver's family-wise phase, so it is a FAILURE of the
# run, never a zero result. Columns expected DEAD are dead by construction:
#   nband/ncomp/nfail on `an`  - e(esband_n_att)/e(boot_*) exist only on bootstrap engines
EXPECT_LIVE_ALL = (["crit", "nh", "band_absent", "band_all", "band_hits",
                    "pw_all", "pw_hits", "pw_n"]
                   + [f"bd{h}" for h in HORIZONS] + [f"pw{h}" for h in HORIZONS])
EXPECT_LIVE_BOOT = ["nband", "ncomp", "nfail"]
EXPECT_DEAD = {"an": ["nband", "ncomp", "nfail"]}


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

L = []          # report lines
FAILS = []      # every FAIL verdict, in order
SHORTN = []     # checks not judged for want of reps


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
    """Concatenate the shard CSVs. A shard killed mid-write can leave a short final line;
    drop it rather than let DictReader pad it with None and score it (analyze_bootse2.py).
    De-duplicate on (arm, tau, shard, rep) so a resumed shard overlapping a killed one does
    not double-count. TAU is part of the key because the driver APPENDS to
    <outdir>/mc_bands_s<shard>.csv: a second run at another tau into the same outdir writes
    rows with the same (arm, shard, rep) and must NOT be discarded as duplicates."""
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
        k = (r.get("arm"), r.get("tau"), r.get("shard"), r.get("rep"))
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


def verdict(name, ok, marginal=None, detail="", n=None):
    """PASS / MARGINAL / FAIL, plus INSUFFICIENT_N when the cell is too small to judge.
    Only FAIL moves the sentinel."""
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
    report = REPORT_OUT if REPORT_OUT else os.path.join(OUT, "bands_report.md")
    try:
        if os.path.dirname(report):
            os.makedirs(os.path.dirname(report), exist_ok=True)
        with open(report, "w") as fh:
            fh.write(txt)
        wrote = f"\nwrote {report}"
    except OSError as e:
        wrote = f"\ncould not write the report ({e})"
    print(txt + wrote)
    if SHORTN:
        print(f"INSUFFICIENT_N={len(SHORTN)} (checks not judged: too few reps; MINN={MINN})")
    if FAILS:
        print(f"BANDS_ANALYZE_FAIL={len(FAILS)}")
        sys.exit(1)
    print("BANDS_ANALYZE_OK")
    sys.exit(0)


# ------------------------------------------------------------------ load and integrity
rows, head, files, nbad, ndup = load("mc_bands")
L.append("# Simultaneous sup-t band coverage (mc_bands, PS-45)\n")
L.append(f"outdir `{OUT}`; {len(files)} shard file(s); {len(rows)} rows after de-duplication "
         f"on (arm, tau, shard, rep); {nbad} malformed line(s) dropped; {ndup} duplicate "
         "row(s) dropped.\n")
L.append("The band's claim is SIMULTANEOUS: `band_all` = 1 iff EVERY reported horizon's truth "
         f"is inside the band on that rep. Nominal {NOMINAL:.2f} (mc_bands.do passes no level(), "
         "so wooldid's default level(95) is in force). Pointwise coverage is the CONTROL: each "
         "horizon near nominal, and the joint pointwise event `pw_all` BELOW `band_all`.\n")
L.append(f"Judgment band = central ~95% binomial region, normal approx rounded outward "
         f"(house convention, analyze_mc.py); the interval printed beside every rate is a Wilson "
         f"95% score interval. Cells with fewer than MINN={MINN} reps are reported, never judged.\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")

if not rows:
    fail("rows present", f"no rows in {OUT}/mc_bands_s*.csv - the campaign produced nothing "
                         f"(searched {len(files)} file(s))")
    finish()

missing_cols = [c for c in REQUIRED if head is None or c not in head]
if missing_cols:
    fail("required columns", "the CSV header is missing " + ", ".join(missing_cols) +
         " - this analyzer cannot judge a file it cannot read; the driver's header changed")
    finish()

bad_rc_rows = [r for r in rows if r.get("rc") not in ("0", 0)]
rc_rate = len(bad_rc_rows) / len(rows)
verdict("rc == 0 on every row", rc_rate <= RC_MAX + EPS,
        detail=f"rc!=0 on {cnt(len(bad_rc_rows), len(rows))}; threshold {RC_MAX:.2f} "
               f"(driver ACCEPTANCE: rc = 0 on every row)")
ok_rows = [r for r in rows if r.get("rc") in ("0", 0)]
if not ok_rows:
    fail("usable rows", "every row has rc != 0: no estimate ever formed")
    finish()

by_arm = defaultdict(list)
for r in ok_rows:
    by_arm[r.get("arm", "?")].append(r)
ARMS = sorted(by_arm)

# infer the bootstrap replication count from ncomp (the `an` arm posts none)
ncomps = [int(v) for v in (num(r, "ncomp") for r in ok_rows) if v is not None]
B_MODE = Counter(ncomps).most_common(1)[0][0] if ncomps else None
B = B_ARG if B_ARG is not None else B_MODE
if B_ARG is not None and B_MODE is not None and B_ARG != B_MODE:
    L.append(f"| B argument vs data | WARN | B given as {B_ARG} but the modal ncomp is {B_MODE}; "
             f"judging completeness at {B_ARG} |")

# all-missing columns, per arm (the vacuous-family trap)
for arm in ARMS:
    g = by_arm[arm]
    expect = list(EXPECT_LIVE_ALL) + ([] if arm == "an" else EXPECT_LIVE_BOOT)
    dead = [c for c in expect if all(num(r, c) is None for r in g)]
    if dead:
        fail(f"{arm} columns all-missing",
             "EVERY value is missing in " + ", ".join(dead) +
             f" on {len(g)} rc-0 rows - the driver wrote nothing for this family, so any "
             "rate computed from it would be vacuous, not zero. Check the arm's e() surface "
             "(and, for bd*/pw*, that the run's geometry really reports those horizons)")
    live_but_dead = [c for c in EXPECT_DEAD.get(arm, [])
                     if any(num(r, c) is not None for r in g)]
    if live_but_dead:
        L.append(f"| {arm} unexpected values | NOTE | " + ", ".join(live_but_dead) +
                 " carries values on an arm where it is dead by construction (analytic engine) |")

# ------------------------------------------------------------------ inventory
L.append("\n## Run inventory\n")
L.append("| arm | engine(s) | rows (rc=0) | rc!=0 | band_absent | mean secs | mean crit | "
         "ncomp (modal) | nfail>0 |")
L.append("|---|---|---|---|---|---|---|---|---|")
for arm in ARMS:
    g = by_arm[arm]
    eng = sorted({r.get("engine", "") for r in g if r.get("engine")})
    badrc = sum(1 for r in rows if r.get("arm") == arm and r.get("rc") not in ("0", 0))
    absent = sum(1 for r in g if (num(r, "band_absent") or 0) > 0)
    secs = [num(r, "secs") for r in g if num(r, "secs") is not None]
    crit = [num(r, "crit") for r in g if num(r, "crit") is not None]
    nc = [num(r, "ncomp") for r in g if num(r, "ncomp") is not None]
    nfpos = sum(1 for r in g if (num(r, "nfail") or 0) > 0)
    mc = Counter(int(x) for x in nc).most_common(1)[0][0] if nc else "-"
    L.append(f"| {arm} | {'/'.join(eng) if eng else '-'} | {len(g)} | {badrc} | {absent} | "
             f"{(sum(secs)/len(secs) if secs else float('nan')):.2f} | "
             f"{(sum(crit)/len(crit) if crit else float('nan')):.4f} | {mc} | {nfpos} |")

# ------------------------------------------------------------------ integrity, per arm
L.append("\n## 0. Integrity (per-row facts; judged at any n - these are not statistical)\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")
for arm in ARMS:
    g = by_arm[arm]
    absent = sum(1 for r in g if (num(r, "band_absent") or 0) > 0)
    verdict(f"{arm} band present on every rep", absent <= BANDABSENT_MAX,
            detail=f"band_absent>0 on {cnt(absent, len(g))} (a rep with no e(esband_att) means "
                   "the kernel refused the whole family on one degenerate horizon, so that rep "
                   "contributes NO band event; driver ACCEPTANCE: band_absent = 0)")
    nh_bad = [r for r in g if (num(r, "nh") or -1) != NH_EXPECTED]
    verdict(f"{arm} nh == {NH_EXPECTED}", not nh_bad,
            detail=f"nh != {NH_EXPECTED} on {cnt(len(nh_bad), len(g))} "
                   f"(observed {sorted({num(r,'nh') for r in g})}); the window is "
                   "esprelength(4) espostlength(4) esfixedbaseperiod")
    pwn_bad = [r for r in g if num(r, "pw_n") is not None and num(r, "nh") is not None
               and num(r, "pw_n") != num(r, "nh")]
    verdict(f"{arm} pw_n == nh", not pwn_bad,
            detail=f"pw_n != nh on {cnt(len(pwn_bad), len(g))} (the pointwise control would then "
                   "be measured on a different horizon set from the band)")
    # arithmetic self-consistency: the summary counts must equal the per-horizon flags
    bd_bad = pw_bad = 0
    for r in g:
        bds = [num(r, f"bd{h}") for h in HORIZONS]
        pws = [num(r, f"pw{h}") for h in HORIZONS]
        if num(r, "band_hits") is not None and all(x is not None for x in bds):
            bd_bad += (abs(sum(bds) - num(r, "band_hits")) > EPS)
        if num(r, "pw_hits") is not None and all(x is not None for x in pws):
            pw_bad += (abs(sum(pws) - num(r, "pw_hits")) > EPS)
    verdict(f"{arm} band_hits == sum(bd*)", bd_bad == 0,
            detail=f"mismatch on {cnt(bd_bad, len(g))}")
    verdict(f"{arm} pw_hits == sum(pw*)", pw_bad == 0,
            detail=f"mismatch on {cnt(pw_bad, len(g))}")
    # band_all must be the indicator of band_hits == nh
    allbad = sum(1 for r in g
                 if None not in (num(r, "band_all"), num(r, "band_hits"), num(r, "nh"))
                 and num(r, "band_all") != float(num(r, "band_hits") == num(r, "nh")))
    verdict(f"{arm} band_all == 1[band_hits == nh]", allbad == 0,
            detail=f"mismatch on {cnt(allbad, len(g))} (band_all IS the simultaneous event)")
    if arm != "an" and B is not None:
        inc = sum(1 for r in g if (num(r, "ncomp") or -1) != B)
        verdict(f"{arm} ncomp == B ({B})", inc == 0,
                detail=f"ncomp != {B} on {cnt(inc, len(g))} (a lost replication means the critical "
                       "value was simulated from fewer draws than requested)")
        nf = sum(1 for r in g if (num(r, "nfail") or 0) > 0)
        verdict(f"{arm} nfail == 0", nf == 0,
                detail=f"nfail>0 on {cnt(nf, len(g))}")
    elif arm != "an":
        L.append(f"| {arm} ncomp == B | NOTE | no ncomp values and no B argument: "
                 "completeness not checked |")
    if M_TARGET is not None:
        verdict(f"{arm} row count == M ({M_TARGET})", len(by_arm[arm]) == M_TARGET,
                detail=f"{len(by_arm[arm])} rc-0 rows of {M_TARGET} planned")
    else:
        L.append(f"| {arm} row count | NOTE | {len(by_arm[arm])} rc-0 rows; no --M target "
                 "given, so completeness of the campaign is not checked |")

# ------------------------------------------------------------------ 1. simultaneous coverage
L.append(f"\n## 1. SIMULTANEOUS coverage - the headline (nominal {NOMINAL:.2f})\n")
L.append("| arm | engine | band_all rate (n, Wilson 95%) | count band |")
L.append("|---|---|---|---|")
simul = {}
band_rows = []
for arm in ARMS:
    g = [r for r in by_arm[arm] if num(r, "band_all") is not None
         and (num(r, "band_absent") or 0) == 0]
    n = len(g)
    h = sum(1 for r in g if num(r, "band_all") == 1)
    lo, hi = binom_region(n, NOMINAL)
    eng = "/".join(sorted({r.get("engine", "") for r in g if r.get("engine")})) or "-"
    simul[arm] = (h, n)
    band_rows.append((arm, eng, h, n, lo, hi))
    L.append(f"| {arm} | {eng} | {rate(h, n)} | [{lo},{hi}] |")

L.append("\n| check | verdict | detail |")
L.append("|---|---|---|")
for arm, eng, h, n, lo, hi in band_rows:
    r_ = h / n if n else float("nan")
    verdict(f"{arm} ({eng}) simultaneous coverage", lo <= h <= hi,
            marginal=(n > 0 and COV_MARGINAL[0] <= r_ <= COV_MARGINAL[1]),
            detail=f"band_all {rate(h, n)} vs nominal {NOMINAL:.2f}, count band [{lo},{hi}], "
                   f"MARGINAL corridor [{COV_MARGINAL[0]:.2f},{COV_MARGINAL[1]:.2f}]", n=n)

# ------------------------------------------------------------------ 2. pointwise control
L.append(f"\n## 2. POINTWISE coverage per horizon - the control (nominal {NOMINAL:.2f} each)\n")
L.append("| arm | " + " | ".join(f"rel {h.replace('m','-')}" for h in HORIZONS) + " | n |")
L.append("|---" * (len(HORIZONS) + 2) + "|")
pw_offenders = {}
for arm in ARMS:
    g = by_arm[arm]
    cells, bad = [], []
    nn = 0
    for hz in HORIZONS:
        vals = [num(r, f"pw{hz}") for r in g if num(r, f"pw{hz}") is not None]
        n = len(vals)
        nn = max(nn, n)
        h = sum(1 for v in vals if v == 1)
        lo, hi = binom_region(n, NOMINAL)
        out = not (lo <= h <= hi)
        if n >= MINN and out:
            bad.append(f"rel {hz.replace('m','-')}: {rate(h, n)} outside [{lo},{hi}]")
        cells.append(f"{h}/{n}" + ("*" if (n >= MINN and out) else ""))
    pw_offenders[arm] = (bad, nn)
    L.append(f"| {arm} | " + " | ".join(cells) + f" | {nn} |")
L.append("\n(* = outside the central ~95% binomial region for that horizon.)\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")
for arm in ARMS:
    bad, nn = pw_offenders[arm]
    det = ("not judged: fewer than MINN reps per horizon" if nn < MINN
           else "all horizons inside their central ~95% binomial region" if not bad
           else "; ".join(bad))
    verdict(f"{arm} per-horizon pointwise coverage ({len(HORIZONS)} horizons)", not bad,
            detail=det, n=nn)

# ------------------------------------------------------------------ 3. the multiplicity gap
L.append("\n## 3. The multiplicity gap - joint POINTWISE coverage must sit BELOW the band's\n")
L.append("(With 8 horizons, intervals that only claim pointwise coverage cover ALL of them "
         "well under 95% unless the horizons are nearly perfectly correlated. If pw_all is "
         "not below band_all the band is buying nothing.)\n")
L.append("| arm | pw_all rate (n, Wilson) | band_all rate (n, Wilson) | gap |")
L.append("|---|---|---|---|")
gaps = {}
for arm in ARMS:
    g = [r for r in by_arm[arm] if num(r, "pw_all") is not None]
    n = len(g)
    h = sum(1 for r in g if num(r, "pw_all") == 1)
    bh, bn = simul[arm]
    gap = (bh / bn if bn else float("nan")) - (h / n if n else float("nan"))
    gaps[arm] = (h, n, bh, bn, gap)
    L.append(f"| {arm} | {rate(h, n)} | {rate(bh, bn)} | {gap:+.4f} |")
L.append("\n| check | verdict | detail |")
L.append("|---|---|---|")
for arm in ARMS:
    h, n, bh, bn, gap = gaps[arm]
    nmin = min(n, bn)
    verdict(f"{arm} pw_all share below band_all share", gap > 0,
            detail=f"pw_all {rate(h, n)}, band_all {rate(bh, bn)}, gap {gap:+.4f} "
                   "(directional check; no band - the sign is the claim)", n=nmin)

# ------------------------------------------------------------------ 4. secondary information
L.append("\n## 4. Secondary (information only, never fails the component)\n")
L.append("| arm | mean crit | sd crit | mean band_hits/nh | mean pw_hits/pw_n | mean secs |")
L.append("|---|---|---|---|---|---|")
for arm in ARMS:
    g = by_arm[arm]

    def mean(c):
        v = [num(r, c) for r in g if num(r, c) is not None]
        return sum(v) / len(v) if v else float("nan")

    cv = [num(r, "crit") for r in g if num(r, "crit") is not None]
    m = sum(cv) / len(cv) if cv else float("nan")
    sd = (math.sqrt(sum((x - m) ** 2 for x in cv) / (len(cv) - 1)) if len(cv) > 1
          else float("nan"))
    bh = [num(r, "band_hits") / num(r, "nh") for r in g
          if num(r, "band_hits") is not None and (num(r, "nh") or 0) > 0]
    ph = [num(r, "pw_hits") / num(r, "pw_n") for r in g
          if num(r, "pw_hits") is not None and (num(r, "pw_n") or 0) > 0]
    L.append(f"| {arm} | {m:.4f} | {sd:.4f} | "
             f"{(sum(bh)/len(bh) if bh else float('nan')):.4f} | "
             f"{(sum(ph)/len(ph) if ph else float('nan')):.4f} | {mean('secs'):.2f} |")
L.append("\nThe critical value is the band's whole content: `crit` well above the pointwise "
         "1.96 (or the engine's t) is what widens the band, and a `crit` that does NOT move "
         "with the engine would mean the simultaneity is nominal only.\n")

L.append("\n## Verdict summary\n")
L.append(f"- checks FAILED: {len(FAILS)}")
for f_ in FAILS:
    L.append(f"  - {f_}")
L.append(f"- checks not judged for want of reps (MINN={MINN}): {len(SHORTN)}")
for s_ in SHORTN:
    L.append(f"  - {s_}")

finish()
