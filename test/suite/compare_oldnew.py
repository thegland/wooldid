#!/usr/bin/env python3
"""Old-vs-new wooldid comparison with intersection semantics and a divergence ledger.

Usage:
  compare_oldnew.py OLD_DIR NEW_DIR [--rtol 1e-10] [--ledger ledger.csv]
                    [--configs configs.csv] [--verbose]

Both dirs hold %21x dumps written by the suite runner (same format as test/runbattery.do:
rc|, prov|, varlist|, scalar|name|%21x, macro|name|val, matnames|m|RxC|rn|cn, mat|m|r|c|%21x).

Comparison semantics (differs from the bit-exact regression compare in test/compare.py):
  - rc: if the two builds' rcs differ, OR either rc is nonzero, the config is judged on rc
    ALONE (an erroring old build leaves residual internal margins e() state — content is
    meaningless). The rc pair must match the per-config expectations in configs.csv
    (columns rc_old / rc_new, blank = 0); an expected pair is EXPLAINED, anything else is
    UNEXPLAINED. Without --configs, only rc_old==rc_new==0 proceeds to content.
  - Content (both rcs 0): objects present in BOTH dumps compared at --rtol (default 1e-10,
    the project indifference threshold; NO absolute floor). Objects only in NEW are
    informational ADDITIONS (new features). Objects only in OLD are REMOVALS = failures
    unless ledgered (the old e() surface is a strict subset of the new one by audit).
  - Ledger rows explain differences: a row matches by config glob + kind + object-name
    glob, and — when its rtol field is nonempty — only BOUNDS the difference (|a-b| <=
    row_rtol * max(|a|,|b|)). A matching diff is EXPLAINED[id]. Rows with empty rtol
    explain any difference of that object (reserve for documented LARGE fixes).
  - STALE ledger rows: a row whose config glob matched at least one compared config but
    which never fired is reported STALE and fails the run (keeps the ledger honest).

Ledger CSV columns:
  id,config,kind,pattern,rtol,reason
    id      - short fix id (e.g. A1-B1)
    config  - config-name glob (fnmatch)
    kind    - scalar | macro | matnames | mat | varlist
    pattern - object-name glob (for mat rows: the matrix name)
    rtol    - optional bound for numeric diffs (empty = any difference matches)
    reason  - one line: what fix, which commit

Exit 0 iff no UNEXPLAINED differences and no STALE ledger rows.
"""
import sys, os, glob, csv, fnmatch, re

TEMPVAR = re.compile(r"__0000[0-9A-Za-z]+")

def canon_tempvars(s):
    """Replace Stata tempvar tokens (__0000XX) with order-of-appearance placeholders.

    cfx matrices embed tempvar names in their row/colnames; the two builds allocate
    different counters, so raw names always differ. Distinct tempvars within one string
    stay distinguishable (__TMP1, __TMP2, ...)."""
    seen = {}
    def sub(m):
        t = m.group(0)
        if t not in seen:
            seen[t] = f"__TMP{len(seen) + 1}"
        return seen[t]
    return TEMPVAR.sub(sub, s)

def hex2f(s):
    """Parse Stata %21x to float; missing-value sentinels and unparseables stay strings."""
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

def close(a, b, rtol):
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    if a == b:
        return True
    if rtol == 0:
        return False
    return abs(a - b) <= rtol * max(abs(a), abs(b))

def _selftest():
    assert hex2f("+1.8000000000000X+001") == 3.0
    assert hex2f("-1.0000000000000X-001") == -0.5
    assert hex2f(".") == "." and hex2f(".a") == ".a"
    assert close(1.0, 1.0 + 1e-11, 1e-10) and not close(1.0, 1.0 + 1e-9, 1e-10)
    assert not close(1e-14, 2e-14, 1e-10), "no absolute floor: small magnitudes stay strict"
    assert not close(0.0, 1e-300, 1e-10)

def load(path):
    d = {"scalar": {}, "macro": {}, "mat": {}, "matnames": {}, "rc": None,
         "varlist": None, "prov": ""}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if not line:
            continue
        p = line.split("|")
        if p[0] == "rc": d["rc"] = p[1] if len(p) > 1 else ""
        elif p[0] == "prov": d["prov"] = "|".join(p[1:])
        elif p[0] == "varlist": d["varlist"] = p[1] if len(p) > 1 else ""
        elif p[0] == "scalar": d["scalar"][p[1]] = hex2f(p[2])
        elif p[0] == "macro": d["macro"][p[1]] = "|".join(p[2:])
        elif p[0] == "matnames": d["matnames"][p[1]] = canon_tempvars("|".join(p[2:]))
        elif p[0] == "mat": d["mat"][(p[1], int(p[2]), int(p[3]))] = hex2f(p[4])
    return d

def load_ledger(path):
    rows = []
    with open(path, newline="") as f:
        for i, r in enumerate(csv.DictReader(f)):
            if not r.get("id") or r["id"].startswith("#"):
                continue
            r["_rtol"] = float(r["rtol"]) if r.get("rtol", "").strip() else None
            r["_fired"] = False
            r["_scoped"] = False
            rows.append(r)
    return rows

def load_rc_expect(path):
    exp = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("config") or r["config"].startswith("#"):
                continue
            exp[r["config"]] = (r.get("rc_old", "").strip() or "0",
                                r.get("rc_new", "").strip() or "0")
    return exp

def main():
    _selftest()
    args = sys.argv[1:]
    def optval(flag, default=None):
        return args[args.index(flag) + 1] if flag in args else default
    old_dir, new_dir = args[0], args[1]
    rtol = float(optval("--rtol", "1e-10"))
    ledger = load_ledger(optval("--ledger")) if "--ledger" in args else []
    rc_expect = load_rc_expect(optval("--configs")) if "--configs" in args else {}
    verbose = "--verbose" in args
    stale_check = "--no-stale-check" not in args  # disable on SUBSET runs only

    def ledger_find(cfg, kind, name, reldiff):
        """First ledger row matching this diff; bounded rows must bound the diff."""
        for r in ledger:
            if not fnmatch.fnmatch(cfg, r["config"]): continue
            if r["kind"] != kind: continue
            if not fnmatch.fnmatch(name, r["pattern"]): continue
            if r["_rtol"] is not None and (reldiff is None or reldiff > r["_rtol"]):
                continue  # bounded rows only explain numeric diffs within their bound
            r["_fired"] = True
            return r
        return None

    ocfg = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(old_dir, "*.csv"))}
    ncfg = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(new_dir, "*.csv"))}
    unexplained = explained = additions = 0
    verdicts = {}
    for cfg in sorted(ocfg | ncfg):
        # scope ledger staleness to configs that actually ran
        for r in ledger:
            if fnmatch.fnmatch(cfg, r["config"]):
                r["_scoped"] = True
        if cfg not in ocfg or cfg not in ncfg:
            print(f"[{cfg}] MISSING on {'new' if cfg in ocfg else 'old'} side")
            unexplained += 1
            verdicts[cfg] = "MISSING"
            continue
        O = load(os.path.join(old_dir, cfg + ".csv"))
        N = load(os.path.join(new_dir, cfg + ".csv"))
        lines, adds = [], []

        def classify(kind, name, msg, reldiff=None):
            nonlocal unexplained, explained
            m = ledger_find(cfg, kind, name, reldiff)
            if m:
                explained += 1
                lines.append(("E", f"EXPLAINED[{m['id']}] {msg}"))
            else:
                unexplained += 1
                lines.append(("U", f"UNEXPLAINED {msg}"))

        # --- rc gate ---
        exp_old, exp_new = rc_expect.get(cfg, ("0", "0"))
        rc_o, rc_n = O["rc"] or "", N["rc"] or ""
        rc_ok = (rc_o == exp_old and rc_n == exp_new)
        if rc_o != rc_n or rc_o != "0" or rc_n != "0":
            # rc-only judgement (erroring builds leave residual e() state)
            if rc_ok:
                if rc_o != rc_n:
                    explained += 1
                    lines.append(("E", f"EXPLAINED[rc-by-design] rc: {rc_o} -> {rc_n} "
                                       f"(expected {exp_old} -> {exp_new})"))
                else:
                    lines.append(("E", f"rc {rc_o} on both builds, as expected"))
            else:
                unexplained += 1
                lines.append(("U", f"UNEXPECTED rc pair: old {rc_o} (expected {exp_old}), "
                                   f"new {rc_n} (expected {exp_new})"))
            statmap = {True: "RC_EXPLAINED", False: "RC_UNEXPLAINED"}
            verdicts[cfg] = statmap[rc_ok]
            print(f"[{cfg}] {verdicts[cfg]}")
            for _, l in lines:
                print("    " + l)
            continue
        if not rc_ok:
            classify("rc", "rc", f"rc 0/0 but expectations say {exp_old}/{exp_new} "
                                 f"(stale configs.csv expectation?)")

        # --- content compare (both rc 0) ---
        for kind in ("scalar", "macro", "matnames"):
            for k in sorted(set(O[kind]) | set(N[kind])):
                if k not in O[kind]:
                    adds.append(f"{kind} {k}")
                elif k not in N[kind]:
                    classify(kind, k, f"{kind} {k}: REMOVED in new (was {O[kind][k]})")
                elif kind == "scalar":
                    a, b = O[kind][k], N[kind][k]
                    if not close(a, b, rtol):
                        rd = (abs(a - b) / max(abs(a), abs(b))
                              if not isinstance(a, str) and not isinstance(b, str)
                              and max(abs(a), abs(b)) > 0 else None)
                        classify(kind, k, f"scalar {k}: {a} -> {b}", rd)
                elif O[kind][k] != N[kind][k]:
                    classify(kind, k, f"{kind} {k}: [{O[kind][k]}] -> [{N[kind][k]}]")
        if (O["varlist"] or "") != (N["varlist"] or ""):
            bs, ns = set((O["varlist"] or "").split()), set((N["varlist"] or "").split())
            for v in sorted(ns - bs):
                adds.append(f"varlist {v}")
            for v in sorted(bs - ns):
                classify("varlist", v, f"varlist var REMOVED in new: {v}")
        omats = {k[0] for k in O["mat"]}
        nmats = {k[0] for k in N["mat"]}
        for m in sorted(omats - nmats):
            classify("mat", m, f"matrix {m}: REMOVED in new")
        for m in sorted(nmats - omats):
            adds.append(f"matrix {m}")
        bad = {}
        for k in set(O["mat"]) & set(N["mat"]):
            if not close(O["mat"][k], N["mat"][k], rtol):
                bad.setdefault(k[0], []).append(k)
        for k in set(O["mat"]) ^ set(N["mat"]):
            if k[0] in (omats & nmats):
                bad.setdefault(k[0], []).append(k)  # shape change within a shared matrix
        for m, ks in sorted(bad.items()):
            worst = 0.0
            has_str = False
            for kk in ks:
                a, b = O["mat"].get(kk), N["mat"].get(kk)
                if a is None or b is None or isinstance(a, str) or isinstance(b, str):
                    has_str = True
                elif max(abs(a), abs(b)) > 0:
                    worst = max(worst, abs(a - b) / max(abs(a), abs(b)))
            ex = "; ".join(f"[{r},{c}] {O['mat'].get((m,r,c),'ABSENT')} -> "
                           f"{N['mat'].get((m,r,c),'ABSENT')}" for (_, r, c) in sorted(ks)[:3])
            classify("mat", m,
                     f"matrix {m}: {len(ks)} differing cell(s), worst reldiff "
                     f"{'n/a(str)' if has_str else f'{worst:.2e}'}, e.g. {ex}",
                     None if has_str else worst)

        n_un = sum(1 for t, _ in lines if t == "U")
        if not lines:
            verdicts[cfg] = "MATCH"
        elif n_un == 0:
            verdicts[cfg] = "EXPLAINED"
        else:
            verdicts[cfg] = "DIFFERS"
        additions += len(adds)
        print(f"[{cfg}] {verdicts[cfg]}" + (f" (+{len(adds)} new-only objects)" if adds else ""))
        for t, l in lines:
            if t == "U" or verbose:
                print("    " + l)
        if verbose and adds:
            print("    ADDITIONS: " + ", ".join(adds))

    stale = [r for r in ledger if r["_scoped"] and not r["_fired"]] if stale_check else []
    for r in stale:
        print(f"STALE ledger row {r['id']} ({r['config']} {r['kind']} {r['pattern']}): "
              f"never fired on any matching config")
    counts = {}
    for v in verdicts.values():
        counts[v] = counts.get(v, 0) + 1
    print("VERDICTS " + " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"SUMMARY unexplained={unexplained} explained={explained} "
          f"additions={additions} stale_ledger={len(stale)}")
    sys.exit(0 if unexplained == 0 and not stale else 1)

main()
