#!/usr/bin/env python3
"""Cross-build summary of the reproduction logs. Parses the REPRO_<id>_<build> KEY=<%21x>
lines, decodes them with the repo's own %21x parser (test/compare.py), prints the
old-wrong / new-right triangles, and ties the old/new legs byte-for-byte to the
2026-09-11 suite dumps (test/suite/out/run_old|run_new)."""
import sys, re, os
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
REPO = _REPO
def hex2f(s):
    """%21x -> float; same rule as test/compare.py (hex exponent; missing stays a string)."""
    s = s.strip()
    if s.startswith(".") and len(s) <= 2:
        return s
    try:
        t = s.replace(" ", "")
        neg = t.startswith("-")
        if t.startswith("+") or t.startswith("-"):
            t = t[1:]
        mant, exp = t.split("X")
        v = float.fromhex("0x" + mant + "p" + str(int(exp, 16)))
        return -v if neg else v
    except Exception:
        return s
def same(a, b, tol=1e-12):
    """equal to the project's indifference class (rtol 1e-10 is the ledger threshold;
    1e-12 here so an ulp-level summation-order difference is reported, not hidden)"""
    fa, fb = hex2f(a), hex2f(b)
    if isinstance(fa, float) and isinstance(fb, float):
        r = abs(fa - fb) / max(abs(fa), abs(fb), 1e-300)
        return f"{'YES' if r <= tol else 'NO'} (reldif {r:.1e})"
    return "YES" if a == b else "NO"
def dumptie(build, cfg, mine_dir):
    """every scalar/mat cell in my dump that also exists in the 0911 suite dump: byte-compare"""
    mine = {}
    for line in open(os.path.join(mine_dir, f"{cfg}.csv"), encoding="utf-8", errors="replace"):
        parts = line.rstrip("\n").split("|")
        if parts[0] in ("mat", "scalar"): mine[tuple(parts[:-1])] = parts[-1]
    theirs = {}
    for line in open(os.path.join(REPO, "test/suite/out", f"run_{build}", f"{cfg}.csv"), encoding="utf-8", errors="replace"):
        parts = line.rstrip("\n").split("|")
        if parts[0] in ("mat", "scalar"): theirs[tuple(parts[:-1])] = parts[-1]
    common = [k for k in mine if k in theirs]
    diff = [k for k in common if mine[k] != theirs[k]]
    print(f"  dump-tie [{build}] {cfg}: {len(common)} scalar/mat cells in common with the 0911 dump, {len(diff)} differ at rtol 0" + (f": {diff[:6]}" if diff else ""))
def reldif(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300) if isinstance(a, float) and isinstance(b, float) else "n/a"
def load(idn, build):
    d = {}
    try:
        for line in open(f"repro_{idn}_{build}.log", encoding="utf-8", errors="replace"):
            m = re.match(rf"REPRO_{idn}_{build} (\S+)=(.*)$", line.strip())
            if m: d[m.group(1)] = m.group(2)
    except FileNotFoundError:
        pass
    return d
def dump(build, cfg):
    d = {}
    p = os.path.join(REPO, "test/suite/out", f"run_{build}", f"{cfg}.csv")
    for line in open(p, encoding="utf-8", errors="replace"):
        parts = line.rstrip("\n").split("|")
        if parts[0] == "mat": d[("mat", parts[1], parts[2], parts[3])] = parts[4]
        elif parts[0] == "scalar": d[("scalar", parts[1])] = parts[2]
    return d
def tie(label, mine, theirs):
    print(f"  byte-tie {label}: repro {mine} vs 0911 dump {theirs} -> {'SAME' if mine == theirs else 'DIFFERENT'}")
for idn in sys.argv[1:]:
    o, n = load(idn, "old"), load(idn, "new")
    print(f"==== {idn}: cross-build summary ====")
    if not o or not n:
        print("  missing log(s): old" if not o else "", "new" if not n else ""); continue
    if idn == "A1A2":
        for b, d in (("old", o), ("new", n)):
            f = {k: hex2f(v) for k, v in d.items()}
            print(f"  [{b}] control_mean={f['SS_control_mean']:.16g} pretreat_mean={f['SS_pretreat_mean']:.16g} TRUTH_control_mean={f['TRUTH_control_mean']:.16g}")
            print(f"  [{b}]   control==pretreat: {same(d['SS_control_mean'], d['SS_pretreat_mean'])};  control==TRUTH: {same(d['SS_control_mean'], d['TRUTH_control_mean'])}")
            print(f"  [{b}] insample_iqr={f['SS_insample_iqr']:.16g} pre_iqr={f['SS_pretreat_iqr']:.16g} post_iqr={f['SS_posttreat_iqr']:.16g} control_iqr={f['SS_control_iqr']:.16g}")
            print(f"  [{b}]   TRUTH pre_iqr={f['TRUTH_pre_iqr']:.16g} post_iqr={f['TRUTH_post_iqr']:.16g} control_iqr={f['TRUTH_control_iqr']:.16g}")
            print(f"  [{b}]   pre_iqr==insample_iqr: {same(d['SS_pretreat_iqr'], d['SS_insample_iqr'])};  pre_iqr==TRUTH: {same(d['SS_pretreat_iqr'], d['TRUTH_pre_iqr'])};  post_iqr==TRUTH: {same(d['SS_posttreat_iqr'], d['TRUTH_post_iqr'])};  control_iqr==TRUTH: {same(d['SS_control_iqr'], d['TRUTH_control_iqr'])};  e(ypre_iqr)==SS[2,4]: {same(d['E_ypre_iqr'], d['SS_pretreat_iqr'])}")
            dd = dump(b, "basic")
            tie(f"[{b}] summarystats[4,1]", d['SS_control_mean'], dd[("mat","summarystats","4","1")])
            tie(f"[{b}] summarystats[2,4]", d['SS_pretreat_iqr'], dd[("mat","summarystats","2","4")])
            tie(f"[{b}] e(ypre_iqr)", d['E_ypre_iqr'], dd[("scalar","ypre_iqr")])
            dumptie(b, "basic", f"dumps_{b}")
    if idn == "A12TTB":
        fo = {k: hex2f(v) for k, v in o.items()}; fn = {k: hex2f(v) for k, v in n.items()}
        for b, f, d in (("old", fo, o), ("new", fn, n)):
            print(f"  [{b}] att TT_RAW={f['B_TT_RAW']:.16g}  NOTT_RAW={f['B_NOTT_RAW']:.16g}  TT_RECODED={f['B_TT_RECODED']:.16g}  props TT_RAW={d['PROPS_TT_RAW']}")
            print(f"  [{b}] rss TT_RAW={f['RSS_TT_RAW']:.16g}  NOTT_RAW={f['RSS_NOTT_RAW']:.16g}  TT_RECODED={f['RSS_TT_RECODED']:.16g}")
            print(f"  [{b}]   reldif(att TT_RAW, NOTT_RAW)={reldif(f['B_TT_RAW'], f['B_NOTT_RAW']):.3e}  reldif(rss TT_RAW, NOTT_RAW)={reldif(f['RSS_TT_RAW'], f['RSS_NOTT_RAW']):.3e}   [0 = trends INERT]")
            print(f"  [{b}]   reldif(att TT_RAW, TT_RECODED)={reldif(f['B_TT_RAW'], f['B_TT_RECODED']):.3e}  reldif(rss TT_RAW, TT_RECODED)={reldif(f['RSS_TT_RAW'], f['RSS_TT_RECODED']):.3e}   [~0 = shift-invariant]")
            dd = dump(b, "trig_tt_big")
            tie(f"[{b}] b[1,1] (leg A)", d['B_TT_RAW'], dd[("mat","b","1","1")])
            tie(f"[{b}] e(rss) (leg A)", d['RSS_TT_RAW'], dd[("scalar","rss")])
            dumptie(b, "trig_tt_big", f"dumps_{b}")
        print(f"  [cross] reldif(att old TT_RECODED benchmark, new TT_RAW)={reldif(fo['B_TT_RECODED'], fn['B_TT_RAW']):.3e}  rss: {reldif(fo['RSS_TT_RECODED'], fn['RSS_TT_RAW']):.3e}")
        print(f"  [cross] reldif(att old TT_RAW, new TT_RAW)={reldif(fo['B_TT_RAW'], fn['B_TT_RAW']):.3e}   (the A12-ttb-m ledger movement)")
        print(f"  [cross] reldif(att old TT_RECODED, new TT_RECODED)={reldif(fo['B_TT_RECODED'], fn['B_TT_RECODED']):.3e}   (both builds agree on small codes)")
