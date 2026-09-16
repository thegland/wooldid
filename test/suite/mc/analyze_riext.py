#!/usr/bin/env python3
"""analyze_riext.py - judge the mc_ri_ext.do campaign: the randomization-inference
extensions with no MC coverage at HEAD - continuous treatment (E1), permside(right) (E2),
level(90) (E3), permstrata() (E4) and a permtstat null arm (E5).

THE CLAIM. On every tau = 0 phase the permutation p-values are EXACT on their lattice, so
each recorded p is uniform on that lattice and the rejection rate at a nominal alpha must
equal the ATTAINABLE alpha - not alpha itself. At permreps B the equal-tailed add-one p
takes values 2*min(R, B+2-R)/(B+1) and the one-sided add-one p takes R/(B+1), R = 1..B+1
uniform; this analyzer enumerates that lattice for the run's own B and judges the count
against Binomial(n, attained). That is the analyze_mc.py DEVIATION D1/D2 arithmetic (the
equal-tailed default attains .04 at B = 99 and cannot reach .05 at all at B = 20), applied
to whichever B the campaign actually ran. E3 is not a size arm: it is COVERAGE of the
planted tau by the shift interval e(perm_ci) at level(90).

Bands and thresholds are CONSTANTS at the top with their provenance. Conventions (load /
num / binom_region / chi2 grid / lattice check / PASS-MARGINAL-FAIL / n beside every rate)
are the house ones from analyze_mc.py, analyze_advb.py and analyze_mcri2.py.

usage: analyze_riext.py [outdir] [B] [--out <report.md>] [--M <reps-per-phase>]
  outdir  directory holding mc_ri_ext_s*.csv (default: ./out next to this file)
  B       permreps() of the run; default: inferred from the ncomp column (mcri2's
          convention is a required argument; here the CSV carries it, so it is inferred
          and a disagreement with the argument is reported)
  --out   write the report here instead of <outdir>/riext_report.md
  --M     the planned reps per phase; without it the observed row count is used
Writes <outdir>/riext_report.md, prints the report, and ends with the sentinel
RIEXT_ANALYZE_OK / RIEXT_ANALYZE_FAIL=<n> at column 1. Exit 1 iff any check FAILs.
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
# mc_ri_ext.do's ACCEPTANCE block (2026-09-13) is the source for the structure of these
# checks; that block itself says "to be preregistered before the production run" and no
# PREREG_*.md carries it, so the bands below are NOT PREREGISTERED except where marked.
ALPHAS = (0.05, 0.10)   # driver ACCEPTANCE: "rejection COUNT at .05 and .10 inside the
                        # central ~95% region of Binomial(M, alpha)" - with the attained
                        # lattice alpha substituted for alpha (analyze_mc.py DEVIATION D1).
LEVEL_COV = 0.90        # E3 runs level(90); driver ACCEPTANCE: "coverage COUNT of tau by
                        # [ll_main, ul_main] inside the Binomial(M, .90) region".
NPOS_MIN = 3000         # PREREGISTERED: PREREG_NEWFEAT_2026-08-24.md - "assignment spaces
                        # chosen >= ~3000 so tie prob <= ~3%". A small assignment space puts
                        # a FLOOR under the achievable p (min p = 2/npos-ish), so a rate at
                        # alpha can be unattainable for reasons that have nothing to do with
                        # the estimator.
RC_MAX = 0.0            # driver ACCEPTANCE: "integrity: rc = 0 ... on every row".
# House MARGINAL corridors for SIZE, verbatim from analyze_mc.py's M2 block.
MARG_05 = (0.025, 0.085)
MARG_10 = (0.070, 0.135)
# MARGINAL corridor for the E3 COVERAGE claim: nominal -/+ 5 / +4 points. NOT PREREGISTERED
# (the house corridors are written for size at .05/.10 and do not transfer to a .90 claim).
COV_MARGINAL = (LEVEL_COV - 0.05, min(1.0, LEVEL_COV + 0.04))
# Below MINN reps a cell is REPORTED, never judged - a 1-rep smoke must say "insufficient
# n", not "PASS". NOT PREREGISTERED; house analogue analyze_bootse2.py's MCB2_MINN = 15;
# 30 is used because binom_region(n, .05) is [0, k] below ~25 and cannot discriminate.
MINN = int(os.environ.get("MCRIX_MINN", "30"))
CHI2_9_CRIT = 21.666    # df=9, p=.01 (house uniformity convention; secondary flag only)

# Per phase: the arm label the driver writes, the null status, the expected level and
# p-value convention, and which p columns are LIVE. The ame family and the C-family joint p
# exist only where contreat() runs (E1); everywhere else they are '.' BY CONSTRUCTION and
# their being empty is not a defect. Everything else missing everywhere IS a defect.
ATT_PS = [("p_main", "main_att"), ("p_pre", "pretreatment_att"),
          ("p_post0", "ES_att_post0"), ("p_post1", "ES_att_post1")]
AME_PS = [("p_ame", "main_ame"), ("p_amepost0", "ES_ame_post0"),
          ("p_amepost1", "ES_ame_post1")]
PHASES = {
    "E1": {"arm": "contreat", "null": True, "level": 95, "conv": "equaltailed",
           "ps": ATT_PS + AME_PS, "joint": [("jointp", "e(joint_horiz)"),
                                            ("jointCp", "e(joint_Choriz)")],
           "extra_live": ["ll_ame", "ul_ame", "b_ame"],
           "what": "contreat(dose) + contreatcontrols(basic), sharp null"},
    "E2": {"arm": "right", "null": True, "level": 95, "conv": "right",
           "ps": ATT_PS, "joint": [("jointp", "e(joint_horiz)")], "extra_live": [],
           "what": "permside(right): one-sided add-one p"},
    "E3": {"arm": "level90", "null": False, "level": 90, "conv": "equaltailed",
           "ps": ATT_PS, "joint": [("jointp", "e(joint_horiz)")], "extra_live": [],
           "what": "level(90): COVERAGE of the planted tau by e(perm_ci)"},
    "E4": {"arm": "strata", "null": True, "level": 95, "conv": "equaltailed",
           "ps": ATT_PS, "joint": [("jointp", "e(joint_horiz)")], "extra_live": [],
           "what": "permstrata(str): exact within strata"},
    "E5": {"arm": "tstat0", "null": True, "level": 95, "conv": "equaltailed",
           "ps": ATT_PS, "joint": [("jointp", "e(joint_horiz)")], "extra_live": [],
           "what": "permtstat at tau = 0 (the size companion of mc_ri2 D4)"},
    # E6 added 2026-09-16 for ticket T1 (PREREG_FEATURE_FINECLUS_2026-09-16.md F2/P1). It is
    # a REGISTRY ENTRY ONLY: same null status, same convention, same p columns and therefore
    # the same lattice rule, the same binomial region and the same MARGINAL corridors as
    # E2/E4/E5. No band is derived for it - the prereg forbids that, and if the arm could not
    # be judged by the existing rule that would itself be the finding.
    "E6": {"arm": "fineclus", "null": True, "level": 95, "conv": "equaltailed",
           "ps": ATT_PS, "joint": [("jointp", "e(joint_horiz)")], "extra_live": [],
           "what": "cluster() strictly finer than i under permonsetcohorts (T1)"},
}
# Live on EVERY phase; all-missing here is a vacuous family, not a zero result.
LIVE_ALWAYS = ["ncomp", "nfail", "npos", "namesok", "p_main", "p_pre", "p_post0",
               "p_post1", "plo_main", "phi_main", "jointp", "ll_main", "ul_main",
               "sd_main", "b_main", "ndraw"]
REQUIRED = (["phase", "arm", "tau", "shard", "rep", "dgpseed", "rseed", "rc", "secs",
             "pconv", "level"] + LIVE_ALWAYS
            + ["p_ame", "p_amepost0", "p_amepost1", "jointCp", "ll_ame", "ul_ame", "b_ame"])


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
NOTJUDGED = []      # claims the lattice cannot express at this B


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
    (analyze_bootse2.py) and de-duplicate on (phase, arm, tau, shard, rep). TAU is part of
    the key because the driver APPENDS to <outdir>/mc_ri_ext_s<shard>.csv: a second run at
    another TAU into the same outdir writes E3 rows with the same (phase, arm, shard, rep)
    and must NOT be discarded as duplicates."""
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
        k = (r.get("phase"), r.get("arm"), r.get("tau"), r.get("shard"), r.get("rep"))
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


def lattice(B, conv):
    """The p-value lattice at permreps B, by enumeration of the observed rank R = 1..B+1
    (analyze_mcri2.py's construction): equal-tailed add-one p = min(1, 2 min(R, B+2-R)/(B+1)),
    one-sided add-one p = R/(B+1). Under the sharp null R is uniform on 1..B+1."""
    if conv == "right":
        return [R / (B + 1.0) for R in range(1, B + 2)]
    return [min(1.0, 2.0 * min(R, B + 2 - R) / (B + 1.0)) for R in range(1, B + 2)]


def attained(alpha, B, conv):
    """P(p <= alpha) on that lattice: the size the test can actually attain at this B.
    Judging against alpha itself when the lattice cannot reach it is the error
    analyze_mc.py's DEVIATION D1/D2 corrects (R6 is judged at .10 only for exactly this
    reason). Returns 0.0 when alpha is UNATTAINABLE."""
    lat = lattice(B, conv)
    return sum(1 for p in lat if p <= alpha + EPS) / len(lat)


def on_lattice(v, B, conv):
    """Is this p on the B-draw lattice? A p off it means that column lost replications
    (an unidentified horizon in some reassignment), so its exactness argument is
    conditional - analyze_mcri2.py's lattice check."""
    if v is None:
        return True
    if v >= 1 - EPS:
        return True                      # the min(1, .) cap
    step = (2.0 if conv != "right" else 1.0) / (B + 1.0)
    return abs(v / step - round(v / step)) < 1e-6


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
    report = REPORT_OUT if REPORT_OUT else os.path.join(OUT, "riext_report.md")
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
        print(f"RIEXT_ANALYZE_FAIL={len(FAILS)}")
        sys.exit(1)
    print("RIEXT_ANALYZE_OK")
    sys.exit(0)


# ------------------------------------------------------------------ load
rows, head, files, nbad, ndup = load("mc_ri_ext")
L.append("# Randomization-inference extensions (mc_ri_ext)\n")
L.append(f"outdir `{OUT}`; {len(files)} shard file(s); {len(rows)} rows after de-duplication "
         f"on (phase, arm, tau, shard, rep); {nbad} malformed line(s) dropped; {ndup} "
         "duplicate "
         "row(s) dropped.\n")
L.append("The headline is SIZE: on every tau = 0 phase each p column must be uniform on its "
         "lattice, so the rejection count at a nominal alpha sits in the central ~95% region "
         "of Binomial(n, ATTAINED alpha). The attained alpha is enumerated from the run's own "
         "permreps (analyze_mc.py's DEVIATION D1/D2 arithmetic): at B = 199 the equal-tailed "
         "default attains .05 and .10 exactly; at small B it may not reach .05 at all, and a "
         "claim it cannot express is NOT JUDGED rather than failed.\n")
L.append("`npos` is the size of the assignment space: a small space puts a floor under the "
         "achievable p. `pconv` names the p-value convention in force. `namesok` = 0 means a "
         "surface was read by POSITION after a name lookup failed - the numbers may be the "
         "wrong column's.\n")
L.append(f"Cells with fewer than MINN={MINN} reps are reported, never judged.\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")

if not rows:
    fail("rows present", f"no rows in {OUT}/mc_ri_ext_s*.csv - the campaign produced nothing "
                         f"(searched {len(files)} file(s))")
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
    fail("usable rows", "every row has rc != 0: no permutation ever ran")
    finish()

by_phase = defaultdict(list)
for r in ok_rows:
    by_phase[r.get("phase", "?")].append(r)
SEEN = [p for p in PHASES if p in by_phase]
UNKNOWN = [p for p in sorted(by_phase) if p not in PHASES]
if UNKNOWN:
    L.append(f"| unknown phases | NOTE | rows carry phase(s) {', '.join(UNKNOWN)} which this "
             "analyzer has no specification for; they are listed in the inventory and "
             "otherwise ignored |")

ncomps = [int(v) for v in (num(r, "ncomp") for r in ok_rows) if v is not None]
B_MODE = Counter(ncomps).most_common(1)[0][0] if ncomps else None
B = B_ARG if B_ARG is not None else B_MODE
if B is None:
    fail("permreps known", "no ncomp values and no B argument: the p lattice - and so the "
                           "attainable alpha - cannot be established, so no size claim can "
                           "be judged")
    finish()
if B_ARG is not None and B_MODE is not None and B_ARG != B_MODE:
    L.append(f"| B argument vs data | WARN | B given as {B_ARG} but the modal ncomp is "
             f"{B_MODE}; judging the lattice at {B_ARG} |")
L.append(f"| permreps in force | NOTE | B = {B} "
         f"({'argument' if B_ARG is not None else 'modal ncomp'}); attained alpha "
         + ", ".join(f"{a:.2f}->{attained(a, B, 'equaltailed'):.4f} (eq) / "
                     f"{attained(a, B, 'right'):.4f} (right)" for a in ALPHAS) + " |")

# all-missing columns, per phase (the vacuous-family trap)
for ph in SEEN:
    g = by_phase[ph]
    spec = PHASES[ph]
    expect = list(LIVE_ALWAYS) + [c for c, _ in spec["ps"] if c not in LIVE_ALWAYS] \
        + [c for c, _ in spec["joint"] if c not in LIVE_ALWAYS] + spec["extra_live"]
    dead = [c for c in dict.fromkeys(expect) if all(num(r, c) is None for r in g)]
    if dead:
        fail(f"{ph} columns all-missing",
             "EVERY value is missing in " + ", ".join(dead) + f" on {len(g)} rc-0 rows - the "
             "driver wrote nothing for this family, so any rate computed from it would be "
             "VACUOUS, not zero. Check the e() surface names the driver reads")
    if ph != "E1":
        unexpected = [c for c in ["p_ame", "p_amepost0", "p_amepost1", "jointCp",
                                  "ll_ame", "ul_ame", "b_ame"]
                      if any(num(r, c) is not None for r in g)]
        if unexpected:
            L.append(f"| {ph} unexpected values | NOTE | " + ", ".join(unexpected) +
                     " carries values outside E1, where contreat() does not run |")

# ------------------------------------------------------------------ inventory
L.append("\n## Run inventory\n")
L.append("| phase | arm | what | rows (rc=0) | rc!=0 | tau | pconv | level | npos | "
         "ncomp (modal) | nfail>0 | namesok=0 | mean secs |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for ph in sorted(by_phase):
    g = by_phase[ph]
    spec = PHASES.get(ph, {})
    badrc = sum(1 for r in rows if r.get("phase") == ph and r.get("rc") not in ("0", 0))
    taus = sorted({r.get("tau") for r in g})
    conv = sorted({r.get("pconv", "") for r in g})
    lev = sorted({r.get("level", "") for r in g})
    npos = sorted({int(v) for v in (num(r, "npos") for r in g) if v is not None})
    nc = [int(v) for v in (num(r, "ncomp") for r in g) if v is not None]
    secs = [num(r, "secs") for r in g if num(r, "secs") is not None]
    L.append(f"| {ph} | {'/'.join(sorted({r.get('arm','') for r in g}))} | "
             f"{spec.get('what', '(unspecified phase)')} | {len(g)} | {badrc} | "
             f"{'/'.join(str(t) for t in taus)} | {'/'.join(conv)} | {'/'.join(lev)} | "
             f"{'/'.join(str(x) for x in npos) if npos else '-'} | "
             f"{Counter(nc).most_common(1)[0][0] if nc else '-'} | "
             f"{sum(1 for r in g if (num(r, 'nfail') or 0) > 0)} | "
             f"{sum(1 for r in g if num(r, 'namesok') != 1)} | "
             f"{(sum(secs)/len(secs) if secs else float('nan')):.2f} |")

# ------------------------------------------------------------------ 0. integrity
L.append("\n## 0. Integrity (per-row facts; judged at any n - these are not statistical)\n")
L.append("| check | verdict | detail |")
L.append("|---|---|---|")
for ph in SEEN:
    g = by_phase[ph]
    spec = PHASES[ph]
    n = len(g)
    inc = sum(1 for r in g if (num(r, "ncomp") or -1) != B)
    verdict(f"{ph} ncomp == B ({B})", inc == 0,
            detail=f"ncomp != {B} on {cnt(inc, n)} (a lost replication puts that rep's p on a "
                   "coarser lattice than the one judged here)")
    nf = sum(1 for r in g if (num(r, "nfail") or 0) > 0)
    verdict(f"{ph} nfail == 0", nf == 0,
            detail=f"nfail>0 on {cnt(nf, n)} (a failed replication means the estimate never formed "
                   "on that draw; driver ACCEPTANCE: nfail = 0 on every row)")
    nd = sum(1 for r in g if (num(r, "ndraw") or -1) != B)
    verdict(f"{ph} ndraw == B ({B}) for main_att", nd == 0,
            detail=f"ndraw != {B} on {cnt(nd, n)} (main_att's accumulator lost draws)")
    nk = sum(1 for r in g if num(r, "namesok") != 1)
    verdict(f"{ph} namesok == 1", nk == 0,
            detail=f"namesok = 0 on {cnt(nk, n)} (a surface read by POSITION after a name lookup "
                   "failed, so the recorded p may be a different column's)")
    small = sorted({int(v) for v in (num(r, "npos") for r in g)
                    if v is not None and v < NPOS_MIN})
    verdict(f"{ph} npos >= {NPOS_MIN}", not small,
            detail=("assignment space " +
                    (f"{sorted({int(v) for v in (num(r,'npos') for r in g) if v is not None})}"
                     if any(num(r, "npos") is not None for r in g) else "unknown") +
                    f"; the PREREGISTERED floor is ~{NPOS_MIN} (PREREG_NEWFEAT_2026-08-24.md: "
                    "spaces >= ~3000 so tie probability <= ~3%). Below it the achievable p "
                    "has a floor that is not the estimator's fault"))
    want = spec["conv"]
    badc = sorted({r.get("pconv", "") for r in g if r.get("pconv", "") != want})
    verdict(f"{ph} pconv == '{want}'", not badc,
            detail=f"pconv values off '{want}': {badc if badc else 'none'} (the convention decides "
                   "the lattice and therefore the attainable alpha)")
    badl = sorted({r.get("level", "") for r in g
                   if str(r.get("level", "")).strip() != str(spec["level"])})
    verdict(f"{ph} level == {spec['level']}", not badl,
            detail=f"level values off {spec['level']}: {badl if badl else 'none'}")
    if spec["null"]:
        badt = sorted({r.get("tau") for r in g if (num(r, "tau") or 0) != 0})
        verdict(f"{ph} tau == 0 (it is a SIZE arm)", not badt,
                detail=f"tau values off 0: {badt if badt else 'none'} (a nonzero tau here would make "
                       "every 'size' number below a power number instead)")
    # a p outside [0,1] is a broken statistic, not a calibration finding
    oob = []
    for col, nm in spec["ps"] + spec["joint"]:
        bad = [num(r, col) for r in g
               if num(r, col) is not None and not (-EPS <= num(r, col) <= 1 + EPS)]
        if bad:
            oob.append(f"{col} ({nm}): {len(bad)} value(s), e.g. {bad[:3]}")
    verdict(f"{ph} every p in [0,1]", not oob,
            detail="; ".join(oob) if oob else "no p outside [0,1]")
    if M_TARGET is not None:
        verdict(f"{ph} row count == M ({M_TARGET})", n == M_TARGET,
                detail=f"{n} rc-0 rows of {M_TARGET} planned")
    else:
        L.append(f"| {ph} row count | NOTE | {n} rc-0 rows; no --M target given, so "
                 "completeness of the campaign is not checked |")

# ------------------------------------------------------------------ 1. size
L.append("\n## 1. SIZE under tau = 0 - the headline\n")
L.append("(att and, on E1, ame surfaces plus the joint p's. `attained` is P(p <= alpha) on "
         "this run's lattice; the count band is the central ~95% region of Binomial(n, "
         "attained). chi2 is the 10-bin uniformity statistic, df=9, crit 21.666 at p=.01 - "
         "a secondary flag, never a verdict.)\n")
L.append("| phase | column | surface | conv | alpha | attained | rej (n, Wilson 95%) | "
         "count band | chi2 |")
L.append("|---|---|---|---|---|---|---|---|---|")
size_claims = []
for ph in SEEN:
    spec = PHASES[ph]
    if not spec["null"]:
        continue
    g = by_phase[ph]
    for col, nm in spec["ps"] + spec["joint"]:
        # the joint p's are the one-sided add-one p (analyze_mcri2.py), whatever permside
        # the arm ran; the per-column p's follow the arm's pconv
        conv = "right" if col.startswith("joint") else spec["conv"]
        ps = [num(r, col) for r in g if num(r, col) is not None]
        n = len(ps)
        stat, flag = chi2_grid(ps) if n >= 10 else (float("nan"), False)
        cells = []
        for a in ALPHAS:
            att = attained(a, B, conv)
            h = sum(1 for p in ps if p <= a + EPS)
            lo, hi = binom_region(n, att)
            cells.append((a, att, h, n, lo, hi))
            L.append(f"| {ph} | {col} | {nm} | {conv} | {a:.2f} | {att:.4f} | {rate(h, n)} | "
                     f"[{lo},{hi}] | {stat:.1f}{'*' if flag else ''} |")
        size_claims.append((ph, col, nm, conv, cells, stat, flag, n))

L.append("\n| check | verdict | detail |")
L.append("|---|---|---|")
for ph, col, nm, conv, cells, stat, flag, n in size_claims:
    bad, unjud, det = [], [], []
    for a, att, h, nn, lo, hi in cells:
        det.append(f"alpha {a:.2f} (attained {att:.4f}): {rate(h, nn)} band [{lo},{hi}]")
        if att <= EPS:
            unjud.append(f"alpha {a:.2f} is UNATTAINABLE at B={B} on the {conv} lattice "
                         "(analyze_mc.py's R6 rule: judged at the other alpha only)")
            continue
        if not (lo <= h <= hi):
            bad.append(f"alpha {a:.2f}")
    if unjud:
        NOTJUDGED.extend(f"{ph}/{col}: {u}" for u in unjud)
    rt = {a: (h / nn if nn else float("nan")) for a, att, h, nn, lo, hi in cells}
    marg = (MARG_05[0] <= rt.get(0.05, -1) <= MARG_05[1]
            and MARG_10[0] <= rt.get(0.10, -1) <= MARG_10[1])
    verdict(f"{ph} {col} ({nm}) size", not bad, marginal=marg,
            detail="; ".join(det + unjud) + (f"; chi2={stat:.1f}{'*' if flag else ''}"
                                             if stat == stat else ""), n=n)
    if flag:
        L.append(f"| {ph} {col} uniformity (sec) | FLAG | 10-bin chi2 = {stat:.1f} > "
                 f"{CHI2_9_CRIT} on n={n}: the p is not uniform even if the tail rate is in "
                 "band - investigate, never silent |")

# ------------------------------------------------------------------ 2. lattice
L.append("\n## 2. Lattice check (a p off the lattice means that column lost draws)\n")
L.append("(analyze_mcri2.py's check: with all B draws complete an equal-tailed p is an exact "
         f"multiple of 2/{B+1} and a one-sided p of 1/{B+1}. Off-lattice on main_att is a "
         "FAILURE - the exactness argument does not hold for the headline column; off-lattice "
         "elsewhere is a FINDING, because a horizon can be unidentified in some reassignments.)\n")
L.append("| phase | column | on lattice | n | off-lattice examples |")
L.append("|---|---|---|---|---|")
lat_claims = []
for ph in SEEN:
    spec = PHASES[ph]
    g = by_phase[ph]
    for col, nm in spec["ps"] + spec["joint"]:
        conv = "right" if col.startswith("joint") else spec["conv"]
        vals = [num(r, col) for r in g if num(r, col) is not None]
        off = [v for v in vals if not on_lattice(v, B, conv)]
        L.append(f"| {ph} | {col} | {len(vals)-len(off)} | {len(vals)} | "
                 f"{sorted(set(off))[:3] if off else '-'} |")
        lat_claims.append((ph, col, len(off), len(vals)))
L.append("\n| check | verdict | detail |")
L.append("|---|---|---|")
for ph, col, noff, nv in lat_claims:
    if col == "p_main":
        verdict(f"{ph} p_main on the B={B} lattice", noff == 0,
                detail=f"off-lattice on {cnt(noff, nv)} reps (that column lost replications, so its p "
                       "is not exact and the size verdict above is conditional)")
    elif noff:
        L.append(f"| {ph} {col} lattice (sec) | FLAG | off-lattice on {cnt(noff, nv)} reps: "
                 "that column lost replications (an unidentified horizon in some "
                 "reassignment); its exactness argument is conditional |")

# ------------------------------------------------------------------ 3. E3 coverage
L.append(f"\n## 3. E3 - COVERAGE of the planted tau by e(perm_ci) at level(90)\n")
L.append("(The only level() != 95 arm in any MC. A blank interval is the off-centre case "
         "mc_ri2 documents: the null band excludes 0, so no shift interval exists - those "
         "reps are excluded from the coverage denominator and reported separately.)\n")
if "E3" in by_phase:
    g = by_phase["E3"]
    taus = {num(r, "tau") for r in g if num(r, "tau") is not None}
    tau = sorted(taus)[0] if taus else None
    L.append("| quantity | covered (n, Wilson 95%) | count band | mean width | blank |")
    L.append("|---|---|---|---|---|")
    live = [r for r in g if num(r, "p_main") is not None]
    nb = [r for r in live if num(r, "ll_main") is not None and num(r, "ul_main") is not None]
    blank = len(live) - len(nb)
    hits = sum(1 for r in nb if tau is not None
               and num(r, "ll_main") - EPS <= tau <= num(r, "ul_main") + EPS)
    w = [num(r, "ul_main") - num(r, "ll_main") for r in nb]
    lo, hi = binom_region(len(nb), LEVEL_COV)
    L.append(f"| main_att, tau = {tau} | {rate(hits, len(nb))} | [{lo},{hi}] | "
             f"{(sum(w)/len(w) if w else float('nan')):.4f} | {cnt(blank, len(live))} |")
    L.append("\n| check | verdict | detail |")
    L.append("|---|---|---|")
    r_ = hits / len(nb) if nb else float("nan")
    verdict("E3 shift-interval coverage of tau at level(90)", lo <= hits <= hi,
            marginal=(COV_MARGINAL[0] <= r_ <= COV_MARGINAL[1]),
            detail=f"tau = {tau}: {rate(hits, len(nb))} vs nominal {LEVEL_COV:.2f}, count band "
                   f"[{lo},{hi}], MARGINAL corridor "
                   f"[{COV_MARGINAL[0]:.2f},{COV_MARGINAL[1]:.2f}]", n=len(nb))
    if blank:
        L.append(f"| E3 blank intervals (sec) | FLAG | {cnt(blank, len(live))} reps have no "
                 "shift interval (null band excludes 0; e(perm_kernel) off-centre). They are "
                 "out of the coverage denominator, so the rate above is CONDITIONAL on an "
                 "interval existing |")
    L.append("\nE3's p columns are NOT scored: with a planted tau they are power, not size, "
             "and mc_ri2 s5 documents pre-treatment p drift under a planted effect. Their "
             "means are in the secondary table below.")
else:
    L.append("No E3 rows: the level(90) coverage arm did not run.")
    L.append("| check | verdict | detail |")
    L.append("|---|---|---|")
    L.append("| E3 present | NOTE | no rows for the level(90) arm in this outdir |")

# ------------------------------------------------------------------ 4. secondary
L.append("\n## 4. Secondary (information only, never fails the component)\n")
L.append("| phase | mean p_main | mean p_pre | mean p_post0 | mean p_post1 | mean sd_main | "
         "mean b_main | n |")
L.append("|---|---|---|---|---|---|---|---|")
for ph in sorted(by_phase):
    g = by_phase[ph]

    def mn(c):
        v = [num(r, c) for r in g if num(r, c) is not None]
        return sum(v) / len(v) if v else float("nan")

    L.append(f"| {ph} | {mn('p_main'):.4f} | {mn('p_pre'):.4f} | {mn('p_post0'):.4f} | "
             f"{mn('p_post1'):.4f} | {mn('sd_main'):.4f} | {mn('b_main'):.4f} | {len(g)} |")
L.append("\nUnder a correct null the mean equal-tailed p is not .5 but the lattice mean "
         f"({sum(lattice(B, 'equaltailed'))/(B+1):.4f} at B = {B}); the one-sided mean is "
         f"{sum(lattice(B, 'right'))/(B+1):.4f}. Read the means against those, not against .5.\n")

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
