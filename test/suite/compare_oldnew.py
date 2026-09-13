#!/usr/bin/env python3
"""Old-vs-new wooldid comparison with intersection semantics and a divergence ledger.

Usage:
  compare_oldnew.py OLD_DIR NEW_DIR [--rtol 1e-10] [--ledger ledger.csv]
                    [--configs configs.csv]... [--only comp1,comp2]
                    [--no-stale-check] [--verbose]

Both dirs hold %21x dumps written by the suite runner (mkrun.py shards): rc|, prov|,
varlist|, scalar|name|%21x, macro|name|val, matnames|m|RxC|rn|cn, mat|m|r|c|%21x,
and a trailing end|<config> sentinel.

Comparison semantics (differs from the bit-exact regression compare in test/compare.py):
  - UNIVERSE: the set of configs to compare comes from the --configs CSVs (restricted by
    --only components if given). A universe config missing from either dir is a failure;
    a dump present for a config NOT in the universe is a failure (stale run dir).
  - INTEGRITY: dumps without the end| sentinel are MALFORMED failures (truncated by a
    timeout kill). Each side's dumps must carry ONE provenance (label|sha) — mixed
    vintages are a failure. The prov pair is printed in the summary.
  - rc: if the two builds' rcs differ, OR either rc is nonzero, the config is judged on
    the rc pair ALONE against the per-config rc_old/rc_new expectations (blank = 0); an
    expected pair is EXPLAINED, anything else UNEXPLAINED. (An erroring old build leaves
    residual internal margins e() state — content would be noise.)
  - Content (both rcs 0): shared objects compared at --rtol (default 1e-10, NO absolute
    floor). Stata missing codes in %21x (exponent-3ff patterns) decode to sentinel
    STRINGS ('.', '.a', ...) so a missing-code flip can never be explained by a bounded
    ledger row. New-only objects are informational ADDITIONS; old-only objects are
    REMOVALS = failures unless ledgered (old e() surface is a strict subset by audit).
  - Two canonicalisations happen before comparing: e(cmdline) has its uniform CORE-10
    `wooldid ` prefix stripped on load (so the REST of the command line is compared,
    not blanket-blessed), and an old MACRO whose name is posted as a SCALAR on new
    (batch 3, 605d6d7: joint_horiz*/joint_Choriz*) is float()ed and value-compared
    under kind `scalar' instead of being counted as a removal plus an addition.
  - Ledger rows explain differences: match by config glob + kind + object-name glob; a
    row with a numeric rtol only BOUNDS a finite numeric difference (never strings,
    removals, or NaN); empty rtol explains any difference (reserved for documented LARGE
    fixes). --verbose prints each bounded row's worst fired reldiff vs its bound.
  - STALE: a ledger row whose config glob matched a compared config but never fired
    fails the run (disable on subset runs with --no-stale-check). Rows whose glob
    matches NO compared config are reported UNSCOPED (fatal only when stale-checking).

Ledger CSV columns: id,config,kind,pattern,rtol,reason.
Kinds are the dump's own record types - scalar, macro, mat, matnames, esamp, varlist,
files (a `files` row's pattern is fnmatch'ed against the WHOLE new-side files| record,
so it can carry the expected file list; a bare `*` keeps the blanket) - plus two
synthetic kinds:
  - rc: the two ledgerable failures that correspond to no dump record at all - "rc 0/0 but
    expectations say X/Y (stale configs.csv expectation?)" and "rc 0 but the old dump has
    ZERO e() objects". No ledger row has ever used this kind, so this docstring is the only
    place it can be learned. (EXPLAINED[rc-by-design] in the output is NOT a ledger kind: it
    is an inline id for an rc pair matching configs.csv, and has no ledger row to find.)
  - mat_esblock, below:
  - mat_esblock: e(V) cells that are the PS-62 event-study covariance splice and nothing
    else. The comparer peels a cell into this kind only when the two stripes are identical
    AND the cell is off-diagonal AND it was exactly 0 on old AND it is finite and nonzero
    on new AND its row and column are two horizons of the SAME event-study family. Every
    other cell of the same matrix still goes through the ordinary 'mat' rows, so the
    bounded storage-noise tiers keep guarding the diagonal.
Exit 0 iff no UNEXPLAINED differences, no integrity/provenance failures, no STALE rows.
"""
import sys, os, glob, csv, fnmatch, re, math

TEMPVAR = re.compile(r"__0000[0-9A-Za-z]+")

def canon_tempvars(s):
    """Replace Stata tempvar tokens (__0000XX) with order-of-appearance placeholders.

    cfx matrices embed tempvar names in row/colnames; the two builds allocate different
    counters, so raw names always differ. Distinct tempvars stay distinguishable."""
    seen = {}
    def sub(m):
        t = m.group(0)
        if t not in seen:
            seen[t] = f"__TMP{len(seen) + 1}"
        return seen[t]
    return TEMPVAR.sub(sub, s)

def hex2f(s):
    """Parse Stata %21x to float. Stata MISSING codes ('.', '.a'..'.z') are stored as
    +1.<k in top 3 hex digits>X+3ff — decode them to sentinel strings so they compare
    by identity, never numerically. Other unparseables also stay strings."""
    s = s.strip()
    if s.startswith(".") and len(s) <= 2:
        return s
    try:
        t = s.replace(" ", "")
        neg = t.startswith("-")
        if t.startswith("+") or t.startswith("-"):
            t = t[1:]
        mant, exp = t.split("X")
        e = int(exp, 16)
        if e == 0x3ff and not neg and mant.startswith("1."):
            frac = int(mant[2:], 16)
            width = len(mant) - 2  # hex digits after the point (13 in %21x)
            top3, rest = divmod(frac, 16 ** (width - 3))
            if rest == 0 and top3 <= 26:
                return "." if top3 == 0 else "." + chr(ord("a") + top3 - 1)
        v = float.fromhex("0x" + mant + "p" + str(e))
        return -v if neg else v
    except Exception:
        return s

ES_ROW = re.compile(r"^(ES_.+?)_(pre\d+|post\d+|contemp)$")

def es_family(nm):
    """The event-study family a V row/column name belongs to, or None.

    ES_att_pre3 -> ES_att; ES_att_it_post0 -> ES_att_it; ES_ame_SG2_post1 -> ES_ame_SG2.
    The non-greedy stem plus the anchored horizon suffix is what separates ES_att from
    ES_att_it without a family list."""
    m = ES_ROW.match(nm)
    return m.group(1) if m else None

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
    assert hex2f("+1.0000000000000X+3ff") == "."          # system missing decodes
    assert hex2f("+1.0010000000000X+3ff") == ".a"         # extended missing decodes
    assert hex2f("+1.01a0000000000X+3ff") == ".z"
    assert isinstance(hex2f("+1.fffffffffffffX+3fe"), float)  # true maxdouble stays numeric
    assert isinstance(hex2f("+1.0010000000001X+3ff"), float)  # non-pattern 3ff stays numeric
    assert close(1.0, 1.0 + 1e-11, 1e-10) and not close(1.0, 1.0 + 1e-9, 1e-10)
    assert not close(1e-14, 2e-14, 1e-10), "no absolute floor"
    assert not close(0.0, 1e-300, 1e-10)
    assert not close(".", ".a", 1e-3) and close(".", ".", 0)

def load(path):
    d = {"scalar": {}, "esamp": {}, "macro": {}, "mat": {}, "matnames": {},
         "rc": None, "varlist": None, "prov": "", "end": False, "files": ""}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if not line:
            continue
        p = line.split("|")
        if p[0] == "rc": d["rc"] = p[1] if len(p) > 1 else ""
        elif p[0] == "prov": d["prov"] = "|".join(p[1:])
        elif p[0] == "end": d["end"] = True
        elif p[0] == "varlist": d["varlist"] = p[1] if len(p) > 1 else ""
        elif p[0] == "esamp": d["esamp"][p[1]] = hex2f(p[2])
        elif p[0] == "files": d["files"] = "|".join(p[1:])
        elif p[0] == "scalar": d["scalar"][p[1]] = hex2f(p[2])
        elif p[0] == "macro":
            v = "|".join(p[2:])
            # CORE-10: new e(cmdline) = "wooldid " + the v1.3 text (equal on all 399 rc-0/0
            # configs, 2026-09-12). Strip the uniform prefix so the REST of e(cmdline) is
            # compared instead of blanket-blessed by an all-config unbounded row.
            if p[1] == "cmdline" and v.startswith("wooldid "):
                v = v[len("wooldid "):]
            d["macro"][p[1]] = v
        elif p[0] == "matnames": d["matnames"][p[1]] = canon_tempvars("|".join(p[2:]))
        elif p[0] == "mat": d["mat"][(p[1], int(p[2]), int(p[3]))] = hex2f(p[4])
    return d

def load_ledger(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("id") or r["id"].startswith("#"):
                continue
            r["_rtol"] = float(r["rtol"]) if r.get("rtol", "").strip() else None
            r["_fired"] = False
            r["_scoped"] = False
            r["_worst"] = 0.0
            rows.append(r)
    return rows

def load_universe(paths, only):
    """Config universe from the CSVs: name -> (component, rc_old, rc_new)."""
    uni = {}
    for path in paths:
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                if not r.get("config") or r["config"].startswith("#"):
                    continue
                if r["config"] in uni:
                    sys.exit(f"duplicate config across CSVs: {r['config']}")
                if only and r["component"] not in only:
                    continue
                uni[r["config"]] = (r["component"],
                                    r.get("rc_old", "").strip() or "0",
                                    r.get("rc_new", "").strip() or "0")
    return uni

def main():
    _selftest()
    args = sys.argv[1:]
    def optval(flag, default=None):
        return args[args.index(flag) + 1] if flag in args else default
    def allvals(flag):
        return [args[i + 1] for i, a in enumerate(args) if a == flag]
    old_dir, new_dir = args[0], args[1]
    rtol = float(optval("--rtol", "1e-10"))
    ledger = load_ledger(optval("--ledger")) if "--ledger" in args else []
    only = set(optval("--only").split(",")) if "--only" in args else None
    universe = load_universe(allvals("--configs"), only)
    full_universe = load_universe(allvals("--configs"), None) if only else universe
    if not universe:
        sys.exit("empty config universe — pass --configs (and check --only)")
    verbose = "--verbose" in args
    stale_check = "--no-stale-check" not in args  # disable on SUBSET runs only

    def ledger_find(cfg, kind, name, reldiff):
        for r in ledger:
            if not fnmatch.fnmatch(cfg, r["config"]): continue
            if r["kind"] != kind: continue
            if not fnmatch.fnmatch(name, r["pattern"]): continue
            if r["_rtol"] is not None:
                # NOT first-match-wins, despite ledger.csv's header: a bounded row that
                # matched config+kind+pattern but whose diff EXCEEDED its bound falls
                # through to later rows, and because it never sets _fired it is reported at
                # the end as STALE - "never fired on any matching config". When such a row
                # is the only match for its cells, that STALE line means the opposite of
                # what it says: the divergence is still there and has grown past the bound.
                # The same cells also print as UNEXPLAINED in the per-config block hundreds
                # of lines earlier, which is where the real diagnosis is. Never repair a
                # STALE line by widening the bound until it fires - that is exactly how a
                # genuine numeric regression gets blessed by a row whose reason column
                # documents a different, smaller effect.
                if reldiff is None or math.isnan(reldiff) or reldiff > r["_rtol"]:
                    continue  # bounded rows only explain finite numeric diffs in-bound
            r["_fired"] = True
            if reldiff is not None and not math.isnan(reldiff):
                r["_worst"] = max(r["_worst"], reldiff)
            return r
        return None

    ocfg = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(old_dir, "*.csv"))}
    ncfg = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(new_dir, "*.csv"))}
    unexplained = explained = additions = 0
    verdicts = {}
    provs = {"old": set(), "new": set()}

    # integrity: stray dumps for configs outside the FULL universe (stale run dirs);
    # on --only runs, other components' dumps are expected and not strays
    for side, present in (("old", ocfg), ("new", ncfg)):
        for stray in sorted(present - set(full_universe)):
            print(f"[{stray}] STRAY dump on {side} side (not in the config universe — "
                  f"stale run dir?)")
            unexplained += 1

    for cfg in sorted(universe):
        comp, exp_old, exp_new = universe[cfg]
        for r in ledger:
            if fnmatch.fnmatch(cfg, r["config"]):
                # _scoped is set from the CONFIG GLOB ALONE - kind and pattern are not
                # consulted - which is what draws the STALE/UNSCOPED distinction. A row with
                # a typo in its kind or pattern column can therefore never report UNSCOPED:
                # its config glob still matches, so it reports STALE and points the reader
                # at the config axis while the fault is on one of the other two. And because
                # this runs BEFORE the MISSING continue below, a config whose shard crashed
                # still scopes its rows - one dead shard manufactures a block of STALE
                # failures on top of the real MISSING ones.
                r["_scoped"] = True
        if cfg not in ocfg or cfg not in ncfg:
            side = "BOTH" if (cfg not in ocfg and cfg not in ncfg) else \
                   ("new" if cfg in ocfg else "old")
            print(f"[{cfg}] MISSING dump ({side} side)")
            unexplained += 1
            verdicts[cfg] = "MISSING"
            continue
        O = load(os.path.join(old_dir, cfg + ".csv"))
        N = load(os.path.join(new_dir, cfg + ".csv"))
        provs["old"].add(O["prov"])
        provs["new"].add(N["prov"])
        if not O["end"] or not N["end"]:
            side = "old" if not O["end"] else "new"
            print(f"[{cfg}] MALFORMED dump ({side} side missing end| sentinel — "
                  f"truncated by a timeout kill?)")
            unexplained += 1
            verdicts[cfg] = "MALFORMED"
            continue
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

        rc_o, rc_n = O["rc"] or "", N["rc"] or ""
        rc_ok = (rc_o == exp_old and rc_n == exp_new)
        if rc_o != rc_n or rc_o != "0" or rc_n != "0":
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
            verdicts[cfg] = "RC_EXPLAINED" if rc_ok else "RC_UNEXPLAINED"
            print(f"[{cfg}] {verdicts[cfg]}")
            for _, l in lines:
                print("    " + l)
            continue
        if not rc_ok:
            classify("rc", "rc", f"rc 0/0 but expectations say {exp_old}/{exp_new} "
                                 f"(stale configs.csv expectation?)")
            # ONE-SIDED ON PURPOSE, and load-bearing rather than a redundant sanity check.
            # New-only objects are informational ADDITIONS and never enter `lines`, and the
            # verdict is MATCH iff `lines` is empty - so an OLD dump with no e() objects at
            # all (old build ran, rc 0, produced no estimation output) would classify every
            # new object as an addition, emit zero lines, and pass as
            # "MATCH (+N new-only objects)". This guard is the only thing standing between
            # that and a green run. The mirror case needs no guard because it is already
            # loud: an empty NEW dump turns every old object into a REMOVAL, which fails.
        if not O["scalar"] and not O["mat"] and not O["macro"]:
            classify("rc", "rc", "rc 0 but the old dump has ZERO e() objects — "
                                 "estimation output lost?")

        for kind in ("scalar", "esamp", "macro", "matnames"):
            for k in sorted(set(O[kind]) | set(N[kind])):
                if k not in O[kind]:
                    if not (kind == "scalar" and k in O["macro"]):
                        adds.append(f"{kind} {k}")
                elif k not in N[kind]:
                    if kind == "macro" and k in N["scalar"]:
                        # batch 3 (605d6d7): joint_horiz*/joint_Choriz* were re-posted as
                        # SCALARS under their old MACRO names. Kinds are keyed separately,
                        # so this used to be a ledgerable REMOVAL plus an informational
                        # ADDITION and the VALUE was never compared. Compare the pair
                        # numerically under kind `scalar' so the scalar rows govern it; an
                        # old macro that is not a number ('.', text) stays a string, which
                        # only an unbounded row can explain. Old macros carry the 12-15
                        # digits Stata printed, so rtol 1e-10 is safe; a macro printed at
                        # 10 digits or fewer would false-fire here and needs a floor of its
                        # own printing precision.
                        try:
                            a = float(O[kind][k])
                        except ValueError:
                            a = O[kind][k].strip()
                        b = N["scalar"][k]
                        if not close(a, b, rtol):
                            rd = (abs(a - b) / max(abs(a), abs(b))
                                  if not isinstance(a, str) and not isinstance(b, str)
                                  and math.isfinite(a) and math.isfinite(b)
                                  and max(abs(a), abs(b)) > 0 else None)
                            classify("scalar", k,
                                     f"macro->scalar {k}: [{O[kind][k]}] -> {b}", rd)
                        continue
                    classify(kind, k, f"{kind} {k}: REMOVED in new (was {O[kind][k]})")
                elif kind in ("scalar", "esamp"):
                    a, b = O[kind][k], N[kind][k]
                    if not close(a, b, rtol):
                        rd = (abs(a - b) / max(abs(a), abs(b))
                              if not isinstance(a, str) and not isinstance(b, str)
                              and max(abs(a), abs(b)) > 0 else None)
                        classify(kind, k, f"{kind} {k}: {a} -> {b}", rd)
                elif O[kind][k] != N[kind][k]:
                    rd = None
                    if kind == "macro":
                        try:
                            a, b = float(O[kind][k]), float(N[kind][k])
                            if (max(abs(a), abs(b)) > 0 and math.isfinite(a)
                                    and math.isfinite(b)):
                                rd = abs(a - b) / max(abs(a), abs(b))
                        except ValueError:
                            pass
                    classify(kind, k, f"{kind} {k}: [{O[kind][k]}] -> [{N[kind][k]}]", rd)
        if O["files"] != N["files"]:
            classify("files", N["files"] or "*", f"saveplots files: [{O['files']}] -> [{N['files']}]")
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
            # PS-62 (a366c46): the within-family event-study covariance blocks are spliced
            # into e(V), so hundreds of off-diagonals go 0 -> nonzero at reldiff 1. Peeled
            # off into their own kind (mat_esblock) under a STRICT predicate - off-diagonal,
            # exactly 0 on old, finite and nonzero on new, row and column two horizons of
            # the SAME family, stripes identical - so one documented row explains them while
            # every other e(V) cell still goes through the rows below and the bounded
            # storage-noise tiers keep guarding the diagonal. A blanket V row would have
            # masked all of e(V) on 360 configs.
            splice = []
            mn = O["matnames"].get(m)
            if m == "V" and mn and mn == N["matnames"].get(m) and len(mn.split("|")) >= 3:
                rn = mn.split("|")[1].split()
                cn = mn.split("|")[2].split()
                keep = []
                for kk in ks:
                    a, b = O["mat"].get(kk), N["mat"].get(kk)
                    r, c = kk[1], kk[2]
                    fr = es_family(rn[r - 1]) if 1 <= r <= len(rn) else None
                    fc = es_family(cn[c - 1]) if 1 <= c <= len(cn) else None
                    if (r != c and fr is not None and fr == fc
                            and isinstance(a, float) and isinstance(b, float)
                            and a == 0.0 and b == b and b != 0.0):
                        splice.append(kk)
                    else:
                        keep.append(kk)
                ks = keep
            if splice:
                classify("mat_esblock", m,
                         f"matrix {m}: {len(splice)} event-study within-family "
                         f"off-diagonal cell(s) 0 -> nonzero", None)
            if not ks:
                continue
            worst = 0.0
            # has_str is TRUE FOR ABSENT CELLS as well as string cells - the name and the
            # printed "n/a(str)" both undersell it - and that is the mechanism that makes a
            # matrix SHAPE CHANGE impossible for any bounded ledger row to absorb: an absent
            # cell yields reldiff None, and ledger_find refuses None for every row carrying a
            # numeric rtol, so only an unbounded row (reserved for documented large fixes)
            # can ever explain it. Deliberate safety property. It also means that
            # "matrix V: 12 differing cell(s), worst reldiff n/a(str)" on a matrix containing
            # no strings is the shape-change message.
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
        verdicts[cfg] = ("MATCH" if not lines else
                         "EXPLAINED" if n_un == 0 else "DIFFERS")
        additions += len(adds)
        print(f"[{cfg}] {verdicts[cfg]}" + (f" (+{len(adds)} new-only objects)" if adds else ""))
        for t, l in lines:
            if t == "U" or verbose:
                print("    " + l)
        if verbose and adds:
            print("    ADDITIONS: " + ", ".join(adds))

    # provenance uniformity per side
    prov_fail = 0
    for side in ("old", "new"):
        ps = provs[side] - {""}
        if len(ps) > 1:
            print(f"PROVENANCE FAILURE: {side}-side dumps carry {len(ps)} distinct "
                  f"prov stamps: {sorted(ps)}")
            prov_fail += 1
        elif ps:
            p = next(iter(ps))
            if not p.startswith(side + "|"):
                print(f"PROVENANCE FAILURE: {side}-side dumps labeled '{p}'")
                prov_fail += 1
            else:
                print(f"prov[{side}] = {p}")

    stale = [r for r in ledger if r["_scoped"] and not r["_fired"]] if stale_check else []
    for r in stale:
        print(f"STALE ledger row {r['id']} ({r['config']} {r['kind']} {r['pattern']}): "
              f"never fired on any matching config")
    unscoped = [r for r in ledger if not r["_scoped"]]
    for r in unscoped:
        print(f"UNSCOPED ledger row {r['id']} ({r['config']}): glob matched no compared "
              f"config" + (" [fatal]" if stale_check else " [subset run — informational]"))
    if verbose:
        for r in ledger:
            if r["_fired"] and r["_rtol"] is not None:
                print(f"ledger {r['id']}: worst fired reldiff {r['_worst']:.2e} "
                      f"(bound {r['_rtol']:.0e})")
    counts = {}
    for v in verdicts.values():
        counts[v] = counts.get(v, 0) + 1
    print("VERDICTS " + " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"COMPARED {len(verdicts)} of {len(universe)} configured")
    print(f"SUMMARY unexplained={unexplained} explained={explained} "
          f"additions={additions} stale_ledger={len(stale)} "
          f"unscoped={len(unscoped)} prov_fail={prov_fail}")
    bad_exit = (unexplained > 0 or prov_fail > 0 or bool(stale)
                or (stale_check and bool(unscoped)) or len(verdicts) == 0)
    sys.exit(1 if bad_exit else 0)

main()
