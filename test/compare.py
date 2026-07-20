#!/usr/bin/env python3
"""Compare two wooldid battery output dirs. Usage: compare.py BASE_DIR NEW_DIR [--rtol 1e-9]

Battery dumps store every e() value in Stata %21x format (exact hexadecimal doubles,
e.g. +1.921fb54442d18X+001 meaning 0x1.921fb54442d18 * 2^0x001). NOTE: the exponent
after X is HEXADECIMAL. Values are compared numerically with a relative tolerance;
strings that fail to parse fall back to exact string comparison.
"""
import sys, os, glob

def hex2f(s):
    """Parse Stata %21x output to a float. Returns the raw string for missing values
    (., .a, ...) and anything unparseable, so those compare by string equality."""
    s = s.strip()
    if s.startswith(".") and len(s) <= 2:
        return s  # missing value sentinel (., .a, ...)
    try:
        t = s.replace(" ", "")
        neg = t.startswith("-")
        if t.startswith("+") or t.startswith("-"):
            t = t[1:]
        mant, exp = t.split("X")
        e = int(exp, 16)  # exponent is hex in %21x
        v = float.fromhex("0x" + mant + "p" + str(e))
        return -v if neg else v
    except Exception:
        return s

def _selftest():
    cases = [
        ("+1.0000000000000X+000", 1.0),
        ("+1.8000000000000X+001", 3.0),
        ("-1.0000000000000X-001", -0.5),
        ("+1.4000000000000X+002", 5.0),
        ("+1.0000000000000X-3ff", 2.0 ** -1023),
        ("+0.0000000000000X-3ff", 0.0),
    ]
    for s, want in cases:
        got = hex2f(s)
        assert isinstance(got, float) and abs(got - want) <= 1e-300 + 1e-15 * abs(want), (s, got, want)
    assert hex2f(".") == "." and hex2f(".a") == ".a"

def load(path):
    d = {"scalar": {}, "macro": {}, "mat": {}, "matnames": {}, "rc": None, "varlist": None}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if not line:
            continue
        p = line.split("|")
        if p[0] == "rc":
            d["rc"] = p[1] if len(p) > 1 else ""
        elif p[0] == "varlist":
            d["varlist"] = p[1] if len(p) > 1 else ""
        elif p[0] == "scalar":
            d["scalar"][p[1]] = hex2f(p[2])
        elif p[0] == "macro":
            d["macro"][p[1]] = "|".join(p[2:])
        elif p[0] == "matnames":
            d["matnames"][p[1]] = "|".join(p[2:])
        elif p[0] == "mat":
            d["mat"][(p[1], int(p[2]), int(p[3]))] = hex2f(p[4])
    return d

def close(a, b, rtol):
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    if a == b:
        return True
    return abs(a - b) <= rtol * max(abs(a), abs(b), 1e-300) or abs(a - b) < 1e-12

def main():
    _selftest()
    base, new = sys.argv[1], sys.argv[2]
    rtol = 1e-9
    if "--rtol" in sys.argv:
        rtol = float(sys.argv[sys.argv.index("--rtol") + 1])
    bcfg = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(base, "*.csv"))}
    ncfg = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(new, "*.csv"))}
    ndiffs = 0
    for miss in sorted(bcfg - ncfg):
        print(f"[{miss}] MISSING in new dir"); ndiffs += 1
    for add in sorted(ncfg - bcfg):
        print(f"[{add}] NEW config (no baseline)"); ndiffs += 1
    for cfg in sorted(bcfg & ncfg):
        B = load(os.path.join(base, cfg + ".csv"))
        N = load(os.path.join(new, cfg + ".csv"))
        diffs = []
        if B["rc"] != N["rc"]:
            diffs.append(f"rc: {B['rc']} -> {N['rc']}")
        if B["varlist"] != N["varlist"]:
            bs, ns = set((B["varlist"] or "").split()), set((N["varlist"] or "").split())
            diffs.append(f"varlist: removed={sorted(bs-ns)} added={sorted(ns-bs)}")
        for kind in ("scalar", "macro", "matnames"):
            keys = set(B[kind]) | set(N[kind])
            for k in sorted(keys):
                if k not in B[kind]:
                    diffs.append(f"{kind} {k}: ADDED = {N[kind][k]}")
                elif k not in N[kind]:
                    diffs.append(f"{kind} {k}: REMOVED (was {B[kind][k]})")
                elif kind == "scalar":
                    if not close(B[kind][k], N[kind][k], rtol):
                        diffs.append(f"scalar {k}: {B[kind][k]} -> {N[kind][k]}")
                elif B[kind][k] != N[kind][k]:
                    diffs.append(f"{kind} {k} changed: [{B[kind][k]}] -> [{N[kind][k]}]")
        cells = set(B["mat"]) | set(N["mat"])
        matdiff = {}
        for k in cells:
            if k not in B["mat"] or k not in N["mat"] or not close(B["mat"][k], N["mat"][k], rtol):
                matdiff.setdefault(k[0], []).append(k)
        for m, ks in sorted(matdiff.items()):
            ex = sorted(ks)[:4]
            exs = "; ".join(
                f"[{r},{c}] {B['mat'].get((m,r,c),'ABSENT')} -> {N['mat'].get((m,r,c),'ABSENT')}"
                for (_, r, c) in ex)
            diffs.append(f"matrix {m}: {len(ks)} cell diff(s), e.g. {exs}")
        if diffs:
            ndiffs += len(diffs)
            print(f"[{cfg}] {len(diffs)} difference group(s):")
            for dd in diffs:
                print(f"    {dd}")
        else:
            print(f"[{cfg}] IDENTICAL")
    print(f"TOTAL_DIFF_GROUPS={ndiffs}")
    sys.exit(0 if ndiffs == 0 else 1)

main()
