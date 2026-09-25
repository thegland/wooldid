#!/usr/bin/env python3
"""degeneracy_scan.py - the `degeneracy` suite component: hidden-coverage assertion over the
oldnew dumps (no Stata; ~10 s).  Adopted from the Stage A suite audit of 2026-09-12/13
(degeneracy_scan F12, critic amendments A2/A3/A4, F9 digest diff, F11 seed rows).

What it answers: "does every clean-rc config on well-behaved/stress data still exercise the
surface it was written to test, or has that surface gone degenerate (b-only, -1e100 SEs,
blanked pwsgs, empty ES window, ...) on BOTH builds, where oldnew cannot see it?"

Usage (the run_full.sh hook passes exactly these; every flag takes one value):
  degeneracy_scan.py --dumps OUT/run_new --configs configs.csv [--configs configs_ext.csv]
                     --data-classes data_classes.csv --ado wooldid.ado --out OUT/degeneracy
                     --expect degeneracy_expected.csv
                     [--ledger ledger.csv]        default: ledger.csv beside the first --configs
                     [--logs-new OUT/wd_new]      shard logs -> per-config caution digest
                                                  (--logs is accepted as a synonym)
                     [--logs-old OUT/wd_old]      old-side digest + caution_diff_old_new.csv
                     [--dumps-old OUT/run_old]    old-side rc/properties annotation column
  Any other flag is a usage error (exit 2).  Without --expect it only scans (exit 3).

Per config it reads the %21x dump and classifies every result surface as LIVE / INFDEG
(estimate finite, inference cells missing or -1e100) / VOID (estimate itself missing), reads
the feature flags the ado posts (perm_*/boot_*/es_ndegen_*/wdid_sgdegen/...), checks the
realized event-study window against the requested one, counts exact-zero b cells and
non-positive e(V) diagonals, joins the config's options (canonicalised with the same
abbreviation rule coverage_check.py uses), its data class, the LARGE ledger rows that
degeneracy_expected.csv cites as explaining a degeneracy, and the column-0 `Caution:` lines of
its shard log.

Classes (per clean-rc config; see DEGENERACY.md):
  CLEAN            no flag
  BY_CONSTRUCTION  only ADV_IMPOSSIBLE_WINDOW / ADV_POOLED_NOP / JT_WINDOW_SHORT / SG_ONE_LEVEL
                   (the config's own options make the surface impossible)
  INFO             only NOSE_BY_DESIGN / NCLUST_TINY / BES_ABSENT(impvalues|permute) /
                   SGDEGEN(:report) - the values are present or the mode computes no variance
  RULED            every substantive flag's e() object is named by a LARGE ledger row that
                   degeneracy_expected.csv CITES for the config (the citation column gates the
                   join - see ledger_known; before 2026-09-14 any row on the same matrix
                   counted, which credited rows written about something else)
  PART_RULED       the cited rows name some of them
  HIDDEN           none of them - nobody has ruled this degeneracy
  (degeneracy_expected.csv may additionally class a row DATASET_DESIGN / DATA_EDGE: the data was
  built to lack the cells; those are HIDDEN or PART_RULED to the scanner.)

Assertion (--expect degeneracy_expected.csv; columns config,class,stems,citation[,ncaution]):
  DEGEN_FAIL   <cfg>: NO DUMP ...            a config of the universe has no dump (A2: a timed-out
                                             shard or an ONLY= subset run; never silently skipped)
  DEGEN_FAIL   <cfg>: NEW degeneracy [...]   a flag stem fires that the row does not allow
  DEGEN_FAIL   <cfg>: rc ... / not classified   the surface cannot be assessed
  DEGEN_STALE  <cfg>: STALE expectation [..] an allowed stem no longer fires -> retire it
  DEGEN_STALE  <cfg>: ... no longer exists   the row's config is not in configs.csv/configs_ext.csv
                                             (A3: mirrors the ledger's STALE rule)
  DEGEN_BADROW ...                           malformed row (duplicate config, unknown class, no stems)
  DEGEN_WARN   ...                           never fails: degen-only options / option pairs, class
                                             moved, dangling ledger citation, caution count moved,
                                             UNRULED rows pending a ruling
  DEGEN_NOTE   ...                           information (old-vs-new caution diff, skipped digest)
  DEGEN_SUMMARY ...                          one line of counts, then DEGEN_CHECK_OK (last line,
                                             exit 0) or DEGEN_CHECK_FAIL (exit 1)
Only clean-rc configs on wellbehaved/stress data are asserted for stems (pathological data answers
"correct or observable failure?" only - README "Role of pathological data"); a missing dump fails
for every config regardless of class.

Outputs in --out:
  degeneracy_table.csv     one row per config, all signals + CLASS + FLAGS
  option_inverted.txt      coverage_check-style per-option view with the extra column "clean on
                           wellbehaved+stress AND live on the option's surface", DEGEN-ONLY
                           options and DEGEN-ONLY option PAIRS
  cautions.csv             per-config caution digest of the new shard logs (empty without --logs-new)
  cautions_old.csv, caution_diff_old_new.csv   with --logs-old
  summary.txt              flag histogram, class x data-class table, per-class config lists
  expectations.csv         GENERATED seed in the degeneracy_expected.csv format (never asserted
                           against itself; copy rows from it when curating)

The e() names it keys on are the ones the ado posts as of HEAD 08ec0b4; the SURFACE table
below is the contract "which e() objects prove option X actually ran".  A rename in the ado
needs a scanner update - the same maintenance class as coverage_check's syntax-line parser.
"""
import sys, os, csv, re, glob, fnmatch, collections

# ---------------------------------------------------------------------------------------
# %21x parser: copied from test/suite/compare_oldnew.py:hex2f (decodes Stata missing codes
# to sentinel strings). test/compare.py's hex2f does NOT decode them (parses "." as 2^1023).
# ---------------------------------------------------------------------------------------
def hex2f(s):
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
            width = len(mant) - 2
            top3, rest = divmod(frac, 16 ** (width - 3))
            if rest == 0 and top3 <= 26:
                return "." if top3 == 0 else "." + chr(ord("a") + top3 - 1)
        v = float.fromhex("0x" + mant + "p" + str(e))
        return -v if neg else v
    except Exception:
        return s

SENTINEL = -1e99   # anything <= this is the package's -1e100 "no variance" sentinel

def is_missing(v):
    return isinstance(v, str)

def is_dead(v):
    """missing, -1e100 sentinel, or non-finite"""
    if isinstance(v, str):
        return True
    return v <= SENTINEL or v != v or v in (float("inf"), float("-inf"))

# ---------------------------------------------------------------------------------------
# option canonicalisation: copied from test/suite/coverage_check.py (same abbreviation rule)
# ---------------------------------------------------------------------------------------
def parse_syntax_options(ado_path):
    txt = open(ado_path, encoding="utf-8", errors="replace").read()
    m = re.search(r"^\s*syntax\s+\[varlist.*?,\s*\[(.*?)\]\s*$", txt, re.M | re.S)
    if not m:
        sys.exit("could not locate the wooldid syntax line")
    opts = {}
    for tok in re.finditer(r"([A-Za-z_][A-Za-z_0-9]*)(?:\([^)]*\))?", m.group(1)):
        name = tok.group(1)
        upper = re.match(r"[A-Z_0-9]+", name)
        abbrev = upper.group(0).lower() if upper and name[0].isupper() else name.lower()
        opts[name.lower()] = abbrev
    return opts

def match_option(token, opts):
    t = token.lower()
    for full, ab in opts.items():
        if t.startswith(ab) and full.startswith(t):
            return full
    return None

def config_tokens_with_args(optstring):
    """[(token, argstring)] split on spaces at paren depth 0."""
    toks, depth, cur = [], 0, ""
    for ch in optstring:
        if ch == "(": depth += 1
        elif ch == ")": depth -= 1
        if ch == " " and depth == 0:
            if cur: toks.append(cur)
            cur = ""
        else:
            cur += ch
    if cur: toks.append(cur)
    out = []
    for t in toks:
        if "(" in t:
            name, arg = t.split("(", 1)
            out.append((name, arg.rstrip(")")))
        else:
            out.append((t, ""))
    return out

def canon_options(optstring, opts):
    """{canonical: arg} for every recognised token; unknown tokens under '?<tok>'."""
    d = {}
    for name, arg in config_tokens_with_args(optstring):
        full = match_option(name, opts)
        d[full if full else "?" + name] = arg
    return d

# ---------------------------------------------------------------------------------------
# dump loader
# ---------------------------------------------------------------------------------------
def load_dump(path):
    d = {"scalar": {}, "macro": {}, "mat": collections.defaultdict(dict), "matnames": {},
         "rc": None, "end": False, "esampN": None, "files": ""}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if not line: continue
        p = line.split("|")
        k = p[0]
        if k == "rc": d["rc"] = p[1]
        elif k == "esamp" and p[1] == "__N": d["esampN"] = hex2f(p[2])
        elif k == "scalar": d["scalar"][p[1]] = hex2f(p[2])
        elif k == "macro": d["macro"][p[1]] = "|".join(p[2:])
        elif k == "matnames":
            R, C = p[2].split("x")
            d["matnames"][p[1]] = {"R": int(R), "C": int(C),
                                   "rows": p[3].split() if len(p) > 3 else [],
                                   "cols": p[4].split() if len(p) > 4 else []}
        elif k == "mat": d["mat"][p[1]][(int(p[2]), int(p[3]))] = hex2f(p[4])
        elif k == "files": d["files"] = "|".join(p[1:])
        elif k == "end": d["end"] = True
    return d

def load_configs(files):
    rows = []
    for f in files:
        for r in csv.DictReader(open(f, newline="")):
            if not r.get("config") or r["config"].startswith("#"): continue
            r["_src"] = os.path.basename(f)
            rows.append(r)
    return rows

def load_classes(path):
    classes = {}
    for r in csv.DictReader(open(path, newline="")):
        name = (r.get("name") or "").strip()
        if not name or name.startswith("#"): continue
        classes[name] = (r.get("class") or "").strip()
    return classes

def config_class(r, classes):
    over = "config:" + r["config"].strip()
    ds = (r.get("dataset") or "").strip()
    if over in classes: return classes[over]
    if ds in classes: return classes[ds]
    return "UNCLASSIFIED"

def load_ledger(path):
    rows = []
    for r in csv.DictReader(open(path, newline="")):
        if not r.get("id") or r["id"].startswith("#"): continue
        rows.append(r)
    return rows

# ---------------------------------------------------------------------------------------
# results-matrix row classification
# ---------------------------------------------------------------------------------------
RESULT_MATS = re.compile(r"^(wdidmainresults|wdidsgresults_\d+|pwsgs|gteffects|cfx.*|"
                         r"boot_(pctl|bc|bca|wald|tpctl|percentile)_(mainresults|sgresults_\d+|pwsgs|gteffects)|"
                         r"perm_analytic_(mainresults|cfx.*))$")

def family_of(rowlabel, matname):
    """estimand family of a results row: att / ame / att_it / ame_it / att_i / ..."""
    if matname == "gteffects" or matname.endswith("_gteffects"):
        return "gt"
    if matname.startswith("cfx") or "cfx" in matname:
        m = re.match(r"(?:perm_analytic_)?cfx(att|ame)(_it|_i|_itime)?_", matname)
        return (m.group(1) + (m.group(2) or "")) if m else "cfx"
    # wdidmainresults / sgresults rows: 'main:att', 'ES_ame:pre2', 'pretreatment:att_it'
    # pwsgs rows: 'att:2vs1bn.region', 'pre_att:3vs1.region', 'ame_it:..'
    lab = rowlabel
    if ":" in lab:
        a, b = lab.split(":", 1)
        if a in ("main", "pretreatment"): return b
        if a.startswith("ES_"): return a[3:]
        if a.startswith("pre_"): return a[4:]
        return a
    return lab

def classify_results(d, want_analytic_se=True):
    """per results matrix: {(mat, family): {'live':n,'infdeg':n,'void':n, 'rows':n}} plus
    a list of (mat, rowlabel, kind) for the dead rows. Under permute() fast mode the
    se_analytic column is missing BY DESIGN (permanalyticse restores it), so it is only
    read when want_analytic_se."""
    out = {}
    dead = []
    for m, meta in d["matnames"].items():
        if not RESULT_MATS.match(m): continue
        cols = meta["cols"]; rows = meta["rows"]
        cidx = {c: i + 1 for i, c in enumerate(cols)}
        est = cidx.get("estimate", cidx.get("b"))
        se = cidx.get("se")
        if se is None and "se_analytic" in cidx and want_analytic_se:
            se = cidx["se_analytic"]
        pcol = cidx.get("p", cidx.get("pvalue", cidx.get("p_perm")))
        lb = next((cidx[c] for c in cols if c.startswith("lb_") or c == "ll"), None)
        ub = next((cidx[c] for c in cols if c.startswith("ub_") or c == "ul"), None)
        inwin = cidx.get("inwindow")
        cells = d["mat"][m]
        if est is None: continue
        nose = se is None and pcol is None  # b-only matrices (cfxplotnose, cfx*_i)
        for r in range(1, meta["R"] + 1):
            lab = rows[r - 1] if r - 1 < len(rows) else f"r{r}"
            fam = family_of(lab, m)
            key = (m, fam)
            rec = out.setdefault(key, {"live": 0, "infdeg": 0, "void": 0, "rows": 0, "nose": nose})
            if inwin is not None and not is_missing(cells.get((r, inwin), ".")) and cells[(r, inwin)] != 1:
                continue  # gteffects: only in-window rows carry inference by contract
            rec["rows"] += 1
            e = cells.get((r, est), ".")
            if is_dead(e):
                rec["void"] += 1; dead.append((m, lab, "void")); continue
            if nose:
                rec["live"] += 1; continue
            infcells = [cells.get((r, c), ".") for c in (se, pcol, lb, ub) if c is not None]
            if any(is_dead(v) for v in infcells):
                rec["infdeg"] += 1; dead.append((m, lab, "infdeg"))
            else:
                rec["live"] += 1
    return out, dead

def vdiag(d, m):
    """counts over the diagonal of matrix m: positive / zero / sentinel / missing"""
    meta = d["matnames"].get(m)
    if not meta: return None
    cells = d["mat"][m]
    c = {"n": 0, "pos": 0, "zero": 0, "sentinel": 0, "missing": 0}
    for i in range(1, min(meta["R"], meta["C"]) + 1):
        v = cells.get((i, i), ".")
        c["n"] += 1
        if is_missing(v): c["missing"] += 1
        elif v <= SENTINEL: c["sentinel"] += 1
        elif v == 0: c["zero"] += 1
        elif v > 0: c["pos"] += 1
        else: c["missing"] += 1
    return c

def matrix_dead_cells(d, m):
    meta = d["matnames"].get(m)
    if not meta: return None
    cells = d["mat"][m]
    n = meta["R"] * meta["C"]
    dead = sum(1 for v in cells.values() if is_dead(v))
    return n, dead

# ---------------------------------------------------------------------------------------
# per-config analysis
# ---------------------------------------------------------------------------------------
def sc(d, name):
    return d["scalar"].get(name)

def num(v):
    return v if isinstance(v, float) else None

def analyze(cfg, d, opts):
    o = canon_options(cfg["options"], opts)
    S = {"config": cfg["config"], "component": cfg["component"],
         "dataset": (cfg.get("dataset") or "").strip() or ("prep:" + (cfg.get("prep") or "").strip()),
         "cls": cfg["_cls"], "rc_expected": (cfg.get("rc_new") or "").strip() or "0",
         "rc": d["rc"], "options": cfg["options"]}
    S["opts"] = o
    flags = []
    F = flags.append
    S["rc_ok"] = (d["rc"] == "0")
    if d["rc"] != S["rc_expected"]:
        F(f"RC_MISMATCH({S['rc_expected']}->{d['rc']})")
    if d["rc"] != "0":
        S["flags"] = flags
        S["dead"] = []
        return S
    # ---- generic e() state
    S["properties"] = d["macro"].get("properties", "")
    S["vcvstate"] = d["macro"].get("wdid_vcvstate", "")
    S["has_V"] = "V" in d["matnames"]
    S["N"] = num(sc(d, "N"))
    S["N_clust"] = num(sc(d, "N_clust"))
    S["estimands"] = d["macro"].get("wdid_estimands", "").split()
    S["sglist"] = d["macro"].get("wdid_sglist", "").split()
    S["suppressjt"] = d["macro"].get("wdid_suppressjt", "")
    S["permest"] = d["macro"].get("wdid_permest", "").split()
    S["fastnose"] = d["macro"].get("wdid_fastnose", "")
    # nose modes: fastnose, and permute() without permanalyticse (fast-mode default, author
    # ruling 2026-09-05) compute no analytic variance BY DESIGN
    nose = bool(S["fastnose"]) or ("permute" in o and "permanalyticse" not in o)
    S["nose"] = nose
    if nose: F("NOSE_BY_DESIGN(" + ("fastnose" if S["fastnose"] else "permute-fast") + ")")
    if not S["has_V"] and not nose:
        F(f"V_ABSENT({S['vcvstate'] or 'nostate'})")
    # ---- results matrices
    res, dead = classify_results(d, want_analytic_se=("permanalyticse" in o or "permute" not in o))
    S["res"] = res; S["dead"] = dead
    def agg(prefix_pred):
        t = {"live": 0, "infdeg": 0, "void": 0, "rows": 0}
        fams = collections.defaultdict(lambda: {"live": 0, "infdeg": 0, "void": 0})
        for (m, fam), r in res.items():
            if prefix_pred(m):
                for k in t: t[k] += r[k]
                for k in ("live", "infdeg", "void"): fams[fam][k] += r[k]
        return t, fams
    main, mainf = agg(lambda m: m == "wdidmainresults")
    S["main"] = main; S["mainf"] = dict(mainf)
    if "wdidmainresults" not in d["matnames"]:
        F("MAIN_ABSENT")
    else:
        if main["void"]: F(f"MAIN_VOID({main['void']}/{main['rows']})")
        if main["infdeg"] and not S["fastnose"]:
            det = ",".join(f"{f}:{v['infdeg']}" for f, v in sorted(mainf.items()) if v["infdeg"])
            F(f"MAIN_INFDEG({main['infdeg']}/{main['rows']};{det})")
    sg, sgf = agg(lambda m: m.startswith("wdidsgresults_"))
    S["sg"] = sg
    if sg["rows"] and (sg["void"] or sg["infdeg"]):
        F(f"SG_DEAD({sg['void']}v+{sg['infdeg']}i/{sg['rows']})")
    pw, pwf = agg(lambda m: m == "pwsgs")
    S["pw"] = pw
    if pw["rows"] and (pw["void"] or pw["infdeg"]):
        F(f"PWSGS_DEAD({pw['void']}v+{pw['infdeg']}i/{pw['rows']})")
    gt, _ = agg(lambda m: m == "gteffects")
    S["gt"] = gt
    if gt["rows"] and (gt["void"] or gt["infdeg"]):
        F(f"GT_DEAD({gt['void']}v+{gt['infdeg']}i/{gt['rows']}inwin)")
    cfx, _ = agg(lambda m: m.startswith("cfx"))
    S["cfx"] = cfx
    # cfxplotnose asks for point-only plot cells, so their inference columns are missing by design
    if cfx["rows"] and (cfx["void"] or (cfx["infdeg"] and "cfxplotnose" not in o)):
        F(f"CFX_DEAD({cfx['void']}v+{cfx['infdeg']}i/{cfx['rows']})")
    bres, _ = agg(lambda m: m.startswith("boot_") and RESULT_MATS.match(m))
    S["bres"] = bres
    if bres["rows"] and (bres["void"] or bres["infdeg"]):
        F(f"BOOTRES_DEAD({bres['void']}v+{bres['infdeg']}i/{bres['rows']})")
    # ---- e(V) diagonals
    for vm in ("V", "V_boot", "V_analytic"):
        c = vdiag(d, vm)
        if c:
            S[f"diag_{vm}"] = c
            bad = c["zero"] + c["sentinel"] + c["missing"]
            if bad: F(f"{vm}_DIAG_DEAD({c['sentinel']}s+{c['zero']}z+{c['missing']}m/{c['n']})")
    ves = {}
    for m in d["matnames"]:
        if m.startswith("V_es_"):
            c = vdiag(d, m); ves[m] = c
            bad = c["zero"] + c["sentinel"] + c["missing"]
            if bad: F(f"{m}_DIAG_DEAD({bad}/{c['n']})")
    S["ves"] = ves
    nd = {k: v for k, v in d["scalar"].items() if k.startswith("es_ndegen_")}
    S["es_ndegen"] = nd
    if nd: F("ES_NDEGEN(" + ",".join(f"{k[10:]}={int(v)}" for k, v in nd.items()) + ")")
    # b_es present but V_es absent for the family = "no identified variance" caution path
    if not nose:
        for m in d["matnames"]:
            if m.startswith("b_es_") and ("V_es_" + m[5:]) not in d["matnames"]:
                F(f"VES_ABSENT({m[5:]})")
    # ---- exact-zero b cells
    bz = []
    if "b" in d["matnames"]:
        cols = d["matnames"]["b"]["cols"]
        for (r, c), v in d["mat"]["b"].items():
            if isinstance(v, float) and v == 0.0:
                bz.append(cols[c - 1] if c - 1 < len(cols) else str(c))
    S["b_zero"] = bz
    if bz: F(f"B_ZERO({len(bz)}:{','.join(bz[:4])})")
    # ---- event-study window: requested vs realized
    def intarg(k, default=0):
        a = o.get(k, "")
        try: return int(a)
        except Exception: return default
    pre_req = intarg("esprelength", 0)
    post_req = intarg("espostlength", 0)
    posttrip = post_req == -999
    if posttrip: post_req = 0
    if "subgroups" in o and pre_req > 0 and post_req == 0 and not posttrip:
        post_req = 1  # the ado's documented rewrite (Caution: espostlength set to 1 ...)
    S["pre_req"] = pre_req; S["post_req"] = post_req
    win = {}
    for m, meta in d["matnames"].items():
        if not m.startswith("b_es_"): continue
        fam = m[5:]
        npre = sum(1 for c in meta["cols"] if c.startswith("pre"))
        npost = sum(1 for c in meta["cols"] if c.startswith("post"))
        win[fam] = (npre, npost)
        short = []
        if npre < pre_req: short.append(f"pre{npre}/{pre_req}")
        if npost < post_req: short.append(f"post{npost}/{post_req}")
        if short: F(f"WINDOW_SHORT({fam}:{'+'.join(short)})")
    S["win"] = win
    # realized pre-horizon count for the att family, in the ado's `espre' sense (the
    # omitted fixed reference counts when it sits inside the requested window, ado ~4979)
    esf_bump = ("esfixedbaseperiod" in o and abs(intarg("esrelativeto", -1)) <= pre_req)
    S["espre_realized"] = (win.get("att", win.get("ame", (0, 0)))[0] + (1 if esf_bump else 0)) if win else (0 if pre_req == 0 else None)
    if (pre_req + post_req) > 0 and not win and "b" in d["matnames"] and not S["fastnose"]:
        # ES requested but no b_es_* posted at all: under impvalues and permute fast mode the
        # F2 returns are not built (mode-specific), elsewhere it is a finding
        mode = "impvalues" if "impvalues" in o else ("permute" if "permute" in o else "")
        F(f"BES_ABSENT({mode})" if mode else "BES_ABSENT")
    # ---- subgroups
    if "subgroups" in o:
        nsg = len(S["sglist"])
        S["nsg"] = nsg
        if nsg < 2: F(f"SG_ONE_LEVEL({nsg})")
        if nsg >= 2 and "pwsgs" not in d["matnames"] and "suppresssgprimaryeffects" not in o:
            F("PWSGS_ABSENT")
        if "wdid_sgdegen" in d["scalar"]:
            F(f"SGDEGEN({int(sc(d,'wdid_sgdegen'))}:{d['macro'].get('wdid_sgdegenact','')})")
    # ---- joint tests
    if "jointtests" in o:
        js = {k: v for k, v in d["scalar"].items() if k.startswith("joint_")}
        S["joint"] = js
        if S["suppressjt"] == "suppressjointtests":
            F("JT_SUPPRESSED")
        elif not js:
            # the joint pre-trend test needs a realized `espre' > 1 (ado ~2164): a window with
            # one pre horizon can never carry it, and the ado says nothing when that happens
            er = S.get("espre_realized")
            if er is not None and er <= 1: F(f"JT_WINDOW_SHORT(espre={er})")
            else: F("JT_NONE_POSTED")
        else:
            miss = [k for k, v in js.items() if is_dead(v)]
            if miss: F(f"JT_MISSING({len(miss)}/{len(js)}:{','.join(miss[:3])})")
            # across-subgroup equality tests expected with >= 2 subgroups
            if len(S["sglist"]) >= 2 and not any(k.endswith("_allsgssame") for k in js):
                F("JT_SGEQ_ABSENT")
    # ---- adversarial
    if "adversarial" in o:
        adv = {k: v for k, v in d["scalar"].items() if k.startswith("adversarial_pval_")}
        S["adv"] = adv
        if pre_req <= 0 or post_req <= 0 or (posttrip):
            F("ADV_IMPOSSIBLE_WINDOW")
        elif "esfixedbaseperiod" not in o:
            # BS-12 (2026-09-08): the adversarial null is a line through the single omitted
            # reference period, so a pooled-reference run posts no p (cautioned)
            F("ADV_POOLED_NOP")
        elif not adv:
            F("ADV_NONE_POSTED")
        else:
            miss = [k for k, v in adv.items() if is_dead(v)]
            if miss: F(f"ADV_MISSING({len(miss)}/{len(adv)})")
            # every event-study estimand family (att / att_it, ame / ame_it under contreat)
            # and every subgroup should carry one; att_i / att_itime have no horizons
            fams = [f for f in S["estimands"] if f in ("att", "att_it")]
            if "contreat" in o: fams += ["ame" + f[3:] for f in fams]
            want = set()
            for f in fams:
                want.add(f"adversarial_pval_{f}")
                for s in S["sglist"]: want.add(f"adversarial_pval_{f}_s{s}")
            lack = sorted(want - set(adv))
            if lack: F(f"ADV_PARTIAL({len(adv)}/{len(want)})")
    # ---- permute
    if "permute" in o:
        P = {k: v for k, v in d["scalar"].items() if k.startswith("perm_")}
        S["perm"] = P
        reps, comp, fail = num(sc(d, "perm_reps")), num(sc(d, "perm_completed")), num(sc(d, "perm_failed"))
        if reps is None: F("PERM_NONE_POSTED")
        else:
            if fail: F(f"PERM_FAILED({int(fail)})")
            if comp is not None and reps and comp < reps: F(f"PERM_INCOMPLETE({int(comp)}/{int(reps)})")
            if num(sc(d, "perm_blank")): F(f"PERM_BLANK({int(sc(d,'perm_blank'))})")
            for fl in ("perm_gt_blank", "perm_cfx_blank"):
                if sc(d, fl): F(fl.upper())
            for fl in ("perm_offcentre", "perm_incomplete", "perm_addone"):
                if d["macro"].get(fl): F(f"{fl.upper()}[{d['macro'][fl][:30]}]")
            for est in S["permest"]:
                if f"perm_p_{est}" not in d["scalar"]: F(f"PERM_P_ABSENT({est})")
                elif is_dead(sc(d, f"perm_p_{est}")): F(f"PERM_P_MISSING({est})")
            for pm in ("perm_p", "perm_ci"):
                nd_ = matrix_dead_cells(d, pm)
                if nd_ and nd_[1]: F(f"{pm.upper()}_DEAD({nd_[1]}/{nd_[0]})")
            if "permexact" in o and not num(sc(d, "perm_exact")): F("PERMEXACT_NOT_EXACT")
    # ---- bootstrap
    if "bootstrap" in o:
        B = {k: v for k, v in d["scalar"].items() if k.startswith("boot_")}
        S["boot"] = B
        reps, comp, fail = num(sc(d, "boot_reps")), num(sc(d, "boot_completed")), num(sc(d, "boot_failed"))
        if reps is None: F("BOOT_NONE_POSTED")
        else:
            if fail: F(f"BOOT_FAILED({int(fail)}/{int(reps)})")
            if comp is not None and comp < reps: F(f"BOOT_INCOMPLETE({int(comp)}/{int(reps)})")
            if num(sc(d, "boot_v_posted")) != 1: F("BOOT_V_NOT_POSTED")
            for fl in ("boot_adv_absent", "boot_jt_absent", "boot_sgjt_absent", "boot_ves_failed"):
                if d["macro"].get(fl): F(f"{fl.upper()}[{d['macro'][fl][:30]}]")
            for bm in ("boot_ci_percentile", "boot_ci_bc", "boot_ci_bca", "boot_p_percentile", "boot_p_bc"):
                nd_ = matrix_dead_cells(d, bm)
                if nd_ and nd_[1]: F(f"{bm.upper()}_DEAD({nd_[1]}/{nd_[0]})")
            if "jointtests" in o and "boot_jointW" not in d["matnames"]: F("BOOT_JOINTW_ABSENT")
            if "adversarial" in o and "boot_advW" not in d["matnames"]: F("BOOT_ADVW_ABSENT")
            if len(S["sglist"]) >= 2 and "jointtests" in o and "boot_sgW" not in d["matnames"]: F("BOOT_SGW_ABSENT")
    # ---- simultaneous bands
    if "simultaneousbands" in o:
        eb = [m for m in d["matnames"] if m.startswith("esband_")]
        crit = [k for k in d["scalar"] if k.startswith("esband_crit_")]
        if not eb or not crit: F("ESBAND_ABSENT")
        else:
            for m in eb:
                nd_ = matrix_dead_cells(d, m)
                if nd_[1]: F(f"{m.upper()}_DEAD({nd_[1]}/{nd_[0]})")
    # ---- pretrendpower
    if "pretrendpower" in o:
        if num(sc(d, "roth_ok")) != 1: F(f"ROTH_NOT_OK({sc(d,'roth_ok')})")
        if not any(m.startswith("roth_power_") for m in d["matnames"]): F("ROTH_POWER_ABSENT")
    # ---- impvalues
    if "impvalues" in o:
        if "imdf_min" not in d["scalar"]: F("IMDF_ABSENT")
    # ---- gteffects
    if "gteffects" in o and "gteffects" not in d["matnames"]: F("GTEFFECTS_ABSENT")
    # ---- cfx plots
    if "makecfxplots" in o:
        cm = [m for m in d["matnames"] if m.startswith("cfx")]
        if not cm: F("CFX_ABSENT")
        elif "cfxplotnose" not in o and all(res[k].get("nose") for k in res if k[0].startswith("cfx")):
            F("CFX_NOSE_UNREQUESTED")
    # ---- contreat: ame rows should carry inference
    if "contreat" in o and "wdidmainresults" in d["matnames"]:
        amef = {f: v for f, v in mainf.items() if f.startswith("ame")}
        tot = sum(v["live"] + v["infdeg"] + v["void"] for v in amef.values())
        livec = sum(v["live"] for v in amef.values())
        if tot and livec == 0: F("AME_NO_LIVE_INFERENCE")
        if not amef: F("AME_ROWS_ABSENT")
    # ---- small-cluster / tiny N
    if S["N_clust"] is not None and S["N_clust"] < 10: F(f"NCLUST_TINY({int(S['N_clust'])})")
    S["flags"] = flags
    return S

# ---------------------------------------------------------------------------------------
# surface predicates: option -> (live?, detail). "live" = the dump proves the option's
# own surface produced non-degenerate output.
# ---------------------------------------------------------------------------------------
def main_live(S, fams=None):
    mf = S.get("mainf", {})
    if fams is None:
        fams = [f for f in mf]
    tot = sum(mf[f]["live"] + mf[f]["infdeg"] + mf[f]["void"] for f in fams if f in mf)
    live = sum(mf[f]["live"] for f in fams if f in mf)
    return (tot > 0 and live > 0), f"main {live}/{tot} live"

def _noflag(S, *prefixes):
    return not any(f.startswith(p) for p in prefixes for f in S["flags"])

def surface(opt, S):
    o = S["opts"]; fl = S["flags"]
    if opt == "subgroups":
        # suppresssgprimaryeffects withholds the per-subgroup rows by design; the pairwise table is then the surface
        ok = S.get("nsg", 0) >= 2 and (S["sg"]["live"] > 0 or "suppresssgprimaryeffects" in o) and S["pw"]["live"] > 0
        return ok, f"nsg={S.get('nsg')} sg_live={S['sg']['live']} pw_live={S['pw']['live']}"
    if opt == "jointtests":
        js = S.get("joint", {})
        ok = bool(js) and _noflag(S, "JT_")
        return ok, f"joint={len(js)} suppressed={S.get('suppressjt')!r}"
    if opt == "adversarial" or opt == "altmethod_adversarial":
        adv = S.get("adv", {})
        ok = bool(adv) and _noflag(S, "ADV_")
        return ok, f"adv={len(adv)}"
    if opt in ("permute", "permreps", "permexact", "permstrata", "permtstat", "permtci", "permside",
               "permstepdown", "permabsp", "permanalyticse", "permsave"):
        ok = "perm" in S and _noflag(S, "PERM_")
        if opt == "permstepdown": ok = ok and any(m == "perm_stepdown" for m in S.get("mats", []))
        if opt == "permanalyticse": ok = ok and "perm_analytic_mainresults" in S.get("mats", [])
        if opt == "permtstat": ok = ok and ("perm_stat" in S.get("macros", {}) and S["macros"]["perm_stat"] == "t")
        return ok, ",".join(f for f in fl if f.startswith("PERM"))
    if opt in ("bootstrap", "reps", "rseed", "boottstat", "bootdisplay", "bootdf",
               "bootsave"):
        ok = "boot" in S and _noflag(S, "BOOT_", "BOOTRES_", "V_boot_DIAG")
        if opt == "boottstat": ok = ok and "boot_tstud" in S.get("mats", [])
        if opt == "rseed" and "boot" not in S:  # rseed also used by permute configs
            ok = "perm" in S and _noflag(S, "PERM_")
        return ok, ",".join(f for f in fl if f.startswith("BOOT"))
    if opt == "gteffects":
        ok = S["gt"]["rows"] > 0 and S["gt"]["live"] == S["gt"]["rows"]
        return ok, f"gt {S['gt']['live']}/{S['gt']['rows']} inwin live"
    if opt in ("makecfxplots", "cfxplottypes", "cfxlattice", "makecfxplotsbysg", "cfxplotnose"):
        if "cfxplotnose" in o or opt == "cfxplotnose":
            ok = S["cfx"]["rows"] > 0 and S["cfx"]["void"] == 0
        else:
            ok = S["cfx"]["rows"] > 0 and S["cfx"]["live"] == S["cfx"]["rows"]
        return ok, f"cfx {S['cfx']['live']}/{S['cfx']['rows']}"
    if opt == "simultaneousbands":
        return _noflag(S, "ESBAND"), ""
    if opt == "pretrendpower":
        return _noflag(S, "ROTH"), ""
    if opt == "impvalues":
        ok = _noflag(S, "IMDF") and main_live(S)[0]
        return ok, ""
    if opt in ("contreat", "contreatcontrols", "contreatpoly", "contreatwithin", "contreatelasticitytype",
               "lattice", "latticeignoreweights", "ccnoabsorb"):
        ok, det = main_live(S, [f for f in S.get("mainf", {}) if f.startswith("ame")])
        return ok, "ame " + det
    if opt == "timetrends":
        return main_live(S)[0] and S.get("has_V", False), f"has_V={S.get('has_V')}"
    if opt == "histogramcohorteffects":
        return "histogramestimates" in S.get("mats", []), ""
    if opt == "summarystats":
        return "summarystats" in S.get("mats", []), ""
    if opt == "esprelength":
        return not any(f.startswith("WINDOW_SHORT") and "pre" in f for f in fl) and S.get("pre_req", 0) > 0, ""
    if opt == "espostlength":
        return not any(f.startswith("WINDOW_SHORT") and "post" in f for f in fl) and S.get("post_req", 0) > 0, ""
    if opt in ("esfixedbaseperiod", "esrelativeto"):
        return not any(f.startswith("WINDOW_SHORT") for f in fl) and main_live(S)[0], ""
    if opt in ("att_it", "att_i", "att_itime"):
        return main_live(S, [opt, "ame" + opt[3:]])
    if opt == "fastnose":
        return S.get("main", {}).get("rows", 0) > 0 and S["main"]["void"] == 0, "b-only by design"
    if opt == "unconditionalse" or opt == "customvce" or opt == "cluster" or opt == "robust":
        if S.get("nose") and not S.get("has_V", False):
            return None, "n/a: b-only mode (fastnose / permute fast) computes no VCE"   # tri-state: not counted
        return main_live(S)[0] and S.get("has_V", False), ""
    if opt == "level":
        return main_live(S)[0] or ("boot" in S and _noflag(S, "BOOT_")), ""
    # default: the run produced live inference on its main rows (or is a b-only mode by design)
    if S.get("nose"):
        return S.get("main", {}).get("rows", 0) > 0 and S["main"]["void"] == 0, "nose-mode default"
    return main_live(S)

# ---------------------------------------------------------------------------------------
# classification vocabulary (explicit, so the expected CSV and the classifier agree)
# ---------------------------------------------------------------------------------------
# flags whose only cause is the config's own options (surface impossible by construction)
BY_CONSTRUCTION_STEMS = ("ADV_IMPOSSIBLE_WINDOW", "ADV_POOLED_NOP", "JT_WINDOW_SHORT", "SG_ONE_LEVEL")
# flags that carry information but no degeneracy: the values are present (sgdegen 'report'),
# the mode computes no variance by design, or a smallness note. Never asserted, never a stem.
INFO_FLAG_PREFIXES = ("NOSE_BY_DESIGN", "BES_ABSENT(impvalues)", "BES_ABSENT(permute)", "NCLUST_TINY")
GOOD_DATA = ("wellbehaved", "stress")
# classes the scanner computes
SCANNER_CLASSES = ("CLEAN", "BY_CONSTRUCTION", "INFO", "RULED", "PART_RULED", "HIDDEN")
# classes degeneracy_expected.csv may carry: the scanner's, plus the curated ones that record
# "the data was built to lack the cells" (scanner sees HIDDEN or PART_RULED for those)
CURATED_CLASSES = ("DATASET_DESIGN", "DATA_EDGE")
EXPECTED_CLASSES = SCANNER_CLASSES + CURATED_CLASSES

def is_info_flag(f):
    return f.startswith(INFO_FLAG_PREFIXES) or (f.startswith("SGDEGEN") and f.endswith(":report)"))

def stems(S):
    """the flag STEMS a config shows, minus INFO flags: this is what the expected CSV lists."""
    return sorted({re.sub(r"\(.*", "", f) for f in S["flags"] if not is_info_flag(f)})

def class_compatible(expected, computed):
    if expected == computed: return True
    return expected in CURATED_CLASSES and computed in ("HIDDEN", "PART_RULED")

# ---------------------------------------------------------------------------------------
# ledger join
# ---------------------------------------------------------------------------------------
def ledger_known(cfgname, S, ledger, candidates=None):
    """ledger row ids whose (config glob, kind, pattern) covers one of this config's dead
    objects; returns (ids, uncovered_objects).
    CANDIDATES (author ruling 2026-09-14): when not None, only ledger ids in this set may
    credit coverage. It is the CITATION column of this config's degeneracy_expected.csv row -
    i.e. the ids a human has actually tied to this config's surfaces. Before this, ANY ledger
    row matching the same MATRIX credited coverage, so a row written to explain one thing was
    credited with ruling another: `sg_degpre`'s WINDOW_SHORT(att_s2:pre0/2) was promoted to
    RULED by F1SG-degpre-m-b, a row about a positional insertion in `b`, and
    `ext_cov10_rcs_pweight_ddd_c1` by DDG-pw-1 on the same shape. Nothing in `kind`, `pattern`
    or the matrix name separates those from the honest cases (BS15-pw-esptrip really is about
    esptrip's pwsgs deadness) - the difference lives only in the prose, so the join cannot
    recover it and a human has to assert it. The citation column already IS that assertion.
    NOT implemented by qualifying `pattern` with the dead stem, which was the first design:
    `pattern` is also compare_oldnew.py's cell-match key (compare_oldnew.py:216), so narrowing
    it there would have stopped those rows explaining the cells they exist to explain and
    turned the comparison red.
    A4 (critique 2026-09-13): flags that name no e() matrix of their own - WINDOW_SHORT, B_ZERO,
    ES_NDEGEN, VES_ABSENT, BES_ABSENT - are a degeneracy of the estimate vector, so a LARGE
    ledger row on `mat b` is required to rule them. Before this they contributed no object, and
    a config carrying one of them plus any ruled object came out RULED (ext_b4_h55_p).
    INFO flags (BES_ABSENT(impvalues|permute), NOSE_BY_DESIGN, ...) contribute no object: they are
    not degeneracies, and the literal prefix match would otherwise un-rule trig_imdegen."""
    objs = set()
    for m, lab, kind in S.get("dead", []):
        objs.add(("mat", m))
    for f in S["flags"]:
        if is_info_flag(f): continue
        if f.startswith("V_DIAG") or f.startswith("V_ABSENT"): objs.add(("mat", "V")); objs.add(("macro", "properties"))
        if f.startswith("V_es_"): objs.add(("mat", f.split("_DIAG")[0]))
        if f.startswith("ADV_"): objs.add(("scalar", "adversarial_pval_*"))
        if f.startswith("JT_"): objs.add(("scalar", "joint_*"))
        if f.startswith("PWSGS"): objs.add(("mat", "pwsgs"))
        if f.startswith(("WINDOW_SHORT", "B_ZERO", "ES_NDEGEN", "VES_ABSENT", "BES_ABSENT")): objs.add(("mat", "b"))
    ids, covered = [], set()
    for r in ledger:
        if candidates is not None and r["id"] not in candidates: continue
        if (r.get("rtol") or "").strip(): continue  # bounded rows explain noise, never a removal/missing
        if not fnmatch.fnmatch(cfgname, r["config"]): continue
        for kind, name in objs:
            if r["kind"] == kind and (fnmatch.fnmatch(name, r["pattern"]) or fnmatch.fnmatch(r["pattern"], name)):
                ids.append(r["id"]); covered.add((kind, name))
    return sorted(set(ids)), sorted(o for o in objs if o not in covered)

CLASS_RANK = {"HIDDEN": 0, "PART_RULED": 1, "RULED": 2}

def classify(S):
    """CLASS of a clean-rc config from its flags and ledger join (see module docstring)."""
    subst = [f for f in S["flags"] if not f.startswith(BY_CONSTRUCTION_STEMS) and not is_info_flag(f)]
    if not S["flags"]: return "CLEAN"
    if not subst:
        return "BY_CONSTRUCTION" if any(f.startswith(BY_CONSTRUCTION_STEMS) for f in S["flags"]) else "INFO"
    if S["ids"] and not S["uncovered"]: return "RULED"
    if S["ids"]: return "PART_RULED"
    return "HIDDEN"

# ---------------------------------------------------------------------------------------
# caution digest from the shard logs (F9). Every shard log opens with the echoed ado source
# (`do <ado>`), so a raw grep for "Caution" is useless (11,417 hits vs 388 real on the
# 2026-09-11 run); only column-0 message lines AFTER the config marker mkrun.py writes
# (`. * ---- <cfg> (<comp>) ----`) count, with the 80-column `> ` wrap re-joined.
# ---------------------------------------------------------------------------------------
MARK = re.compile(r"^\. \* ---- (\S+) \((\S+)\) ----")

def caution_digest(logdir):
    """{config: [message lines]} over every <logdir>/*/*.log"""
    out = {}
    for lg in sorted(glob.glob(os.path.join(logdir, "*", "*.log"))):
        cur = None
        buf = []
        def flush():
            if cur is not None:
                out[cur] = buf[:]
        cont = False   # True only while the previous physical line was a buffered message
        for line in open(lg, encoding="utf-8", errors="replace"):
            line = line.rstrip("\n")
            m = MARK.match(line)
            if m:
                flush(); cur = m.group(1); buf = []; cont = False
                continue
            if cur is None: continue
            if line.startswith("Caution:") or line.startswith("Note:") or line.startswith("Error:") \
                    or line.startswith("Warning:"):
                buf.append(line); cont = True
            elif line.startswith("> ") and cont:
                # batch logs hard-wrap at 80 cols; a "> " line continues the line just above
                buf[-1] += line[2:]
            else:
                cont = False
        flush()
    return out

def n_caution(lines):
    return sum(1 for x in lines if x.startswith("Caution:"))

def write_cautions_csv(path, configs, digest):
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["config", "n_caution", "n_note", "n_error", "lines"])
        for cfg in configs:
            cl = digest.get(cfg, [])
            w.writerow([cfg, n_caution(cl), sum(1 for x in cl if x.startswith("Note:")),
                        sum(1 for x in cl if x.startswith("Error:")), " || ".join(x[:160] for x in cl)])

# ---------------------------------------------------------------------------------------
# expected CSV (the curated contract)
# ---------------------------------------------------------------------------------------
def load_expected(path):
    """rows of degeneracy_expected.csv -> ({config: row}, [DEGEN_BADROW messages]).
    '#' comment rows are skipped; columns config,class,stems,citation are required, ncaution is
    optional (blank = not pinned)."""
    exp, bad = {}, []
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        need = {"config", "class", "stems", "citation"}
        if not need.issubset(set(rd.fieldnames or [])):
            bad.append(f"{path}: header must contain {sorted(need)}; got {rd.fieldnames}")
            return exp, bad
        for i, r in enumerate(rd, start=2):
            cfg = (r.get("config") or "").strip()
            if not cfg or cfg.startswith("#"): continue
            if cfg in exp:
                bad.append(f"{cfg}: duplicate row (line {i})"); continue
            klass = (r.get("class") or "").strip()
            if klass not in EXPECTED_CLASSES:
                bad.append(f"{cfg}: class '{klass}' not in {list(EXPECTED_CLASSES)} (line {i})"); continue
            st = set((r.get("stems") or "").split())
            if not st:
                bad.append(f"{cfg}: no stems listed (line {i}); a row must allow at least one stem"); continue
            nc = (r.get("ncaution") or "").strip()
            if nc and not nc.isdigit():
                bad.append(f"{cfg}: ncaution '{nc}' is not an integer (line {i})"); continue
            exp[cfg] = {"class": klass, "stems": st, "citation": (r.get("citation") or "").strip(),
                        "ncaution": int(nc) if nc else None, "line": i}
    return exp, bad

# ---------------------------------------------------------------------------------------
KNOWN_FLAGS = ("--dumps", "--dumps-old", "--configs", "--data-classes", "--ado", "--ledger",
               "--logs", "--logs-new", "--logs-old", "--out", "--expect")

def usage_error(msg):
    print("DEGEN_USAGE_ERROR: " + msg)
    print("usage: see the module docstring (python3 degeneracy_scan.py --help)")
    sys.exit(2)

def parse_args(argv):
    if not argv or "--help" in argv or "-h" in argv:
        print(__doc__); sys.exit(0 if argv else 2)
    args = collections.defaultdict(list)
    i = 0
    while i < len(argv):
        a = argv[i]
        if a not in KNOWN_FLAGS: usage_error(f"unknown argument {a!r}; known: {' '.join(KNOWN_FLAGS)}")
        if i + 1 >= len(argv) or argv[i + 1].startswith("--"): usage_error(f"{a} needs a value")
        args[a].append(argv[i + 1]); i += 2
    for req in ("--dumps", "--configs", "--data-classes", "--ado"):
        if not args[req]: usage_error(f"{req} is required")
    for one in ("--dumps", "--dumps-old", "--data-classes", "--ado", "--ledger", "--logs", "--logs-new",
                "--logs-old", "--out", "--expect"):
        if len(args[one]) > 1: usage_error(f"{one} given more than once")
    if args["--logs"] and args["--logs-new"]: usage_error("--logs and --logs-new are synonyms; give one")
    return args

def main():
    args = parse_args(sys.argv[1:])
    one = lambda k, d=None: (args[k][0] if args[k] else d)
    dumps, dumps_old, cfgs = one("--dumps"), one("--dumps-old"), args["--configs"]
    for p in cfgs + [one("--data-classes"), one("--ado")]:
        if not os.path.isfile(p): usage_error(f"file not found: {p}")
    ledger_path = one("--ledger") or os.path.join(os.path.dirname(os.path.abspath(cfgs[0])), "ledger.csv")
    logs_new = one("--logs-new") or one("--logs")
    logs_old = one("--logs-old")
    outdir = one("--out", "."); os.makedirs(outdir, exist_ok=True)
    expect = one("--expect")
    notes, warns, classmoved = [], [], []
    classes = load_classes(one("--data-classes"))
    opts = parse_syntax_options(one("--ado"))
    if os.path.isfile(ledger_path):
        ledger = load_ledger(ledger_path)
    else:
        ledger = []; notes.append(f"no ledger at {ledger_path}: RULED/PART_RULED cannot be derived (all substantive flags read HIDDEN)")
    # --expect is read HERE, not only at the assertion below: from 2026-09-14 the citation
    # column gates which ledger rows may credit degeneracy coverage (see ledger_known).
    exp, bad = (load_expected(expect) if expect else ({}, []))
    rows = load_configs(cfgs)
    if not rows: usage_error("no config rows read from " + " ".join(cfgs))
    if not os.path.isdir(dumps):
        # the component runs after oldnew; no dumps dir means oldnew did not run at all
        print(f"DEGEN_FAIL *: dumps directory {dumps} does not exist (oldnew has not run)")
        print(f"DEGEN_SUMMARY configs={len(rows)} dumps_dir=missing")
        print("DEGEN_CHECK_FAIL: no dumps"); sys.exit(1)
    # ---- caution digests (F9): new side always when logs are given, old side optional
    digest_new = None
    if logs_new:
        if os.path.isdir(logs_new): digest_new = caution_digest(logs_new)
        else: warns.append(f"shard-log directory {logs_new} not found: caution digest skipped")
    else:
        notes.append("no --logs-new: caution digest skipped (ncaution blank, caution pins not checked)")
    digest_old = None
    if logs_old:
        if os.path.isdir(logs_old): digest_old = caution_digest(logs_old)
        else: warns.append(f"old shard-log directory {logs_old} not found: caution diff skipped")
    # ---- per-config analysis
    table = []
    for r in rows:
        r["_cls"] = config_class(r, classes)
        p = os.path.join(dumps, r["config"] + ".csv")
        if not os.path.exists(p):
            table.append({"config": r["config"], "component": r["component"], "cls": r["_cls"],
                          "rc_expected": (r.get("rc_new") or "").strip() or "0",
                          "rc": "NODUMP", "flags": ["NO_DUMP"], "opts": canon_options(r["options"], opts),
                          "dead": [], "options": r["options"], "dataset": r.get("dataset", ""),
                          "ncaution": None})
            continue
        d = load_dump(p)
        S = analyze(r, d, opts)
        S["mats"] = list(d["matnames"]); S["macros"] = d["macro"]; S["scalars_seen"] = list(d["scalar"])
        # all-ledger pass: DISCOVERY only. It no longer classifies (it is what promoted
        # sg_degpre on a row about something else); it is reported as a DEGEN_NOTE so
        # "the ledger now covers this, consider citing it" is still visible.
        S["ids_all"], S["uncovered_all"] = ledger_known(r["config"], S, ledger)
        if expect:
            _e = exp.get(r["config"])
            _cand = {c.strip() for c in ((_e["citation"] or "").split(";") if _e else []) if c.strip()}
            S["ids"], S["uncovered"] = ledger_known(r["config"], S, ledger, _cand)
        else:
            # scan-only mode seeds expectations.csv and has no citations to gate on
            S["ids"], S["uncovered"] = S["ids_all"], S["uncovered_all"]
        cl = digest_new.get(r["config"], []) if digest_new is not None else []
        S["ncaution"] = n_caution(cl) if digest_new is not None else None
        S["cautions"] = cl
        # old-side annotation: was the old build b-only / erroring on the same config?
        S["old"] = ""
        if dumps_old:
            po = os.path.join(dumps_old, r["config"] + ".csv")
            if os.path.exists(po):
                do = load_dump(po)
                S["old"] = f"rc{do['rc']}:{do['macro'].get('properties', '-').replace(' ', '')}"
        S["class"] = classify(S) if S["rc"] == "0" else ""
        S["class_all"] = classify(dict(S, ids=S["ids_all"], uncovered=S["uncovered_all"])) \
            if S["rc"] == "0" else ""
        table.append(S)
    names = [S["config"] for S in table]
    # DISCOVERY signal, deliberately not a classification: configs the whole ledger would
    # cover further than their cited ids do. Each is a candidate for citing the row in
    # degeneracy_expected.csv - after a human checks the row is really about that surface.
    promo = [S["config"] for S in table
             if CLASS_RANK.get(S.get("class_all"), -1) > CLASS_RANK.get(S.get("class"), -1)]
    if promo:
        notes.append(f"{len(promo)} config(s) where the whole ledger would credit more coverage than the "
                     f"cited ids do - NOT applied, cite the row in degeneracy_expected.csv if it really "
                     f"rules that surface: " + " ".join(sorted(promo)[:12]) + (" ..." if len(promo) > 12 else ""))
    # ---- per-config table
    cols = ["config", "component", "dataset", "cls", "rc_expected", "rc", "vcvstate", "properties",
            "N", "N_clust", "nsg", "main_live", "main_infdeg", "main_void", "sg_live", "sg_dead",
            "pw_live", "pw_dead", "gt_live", "gt_dead", "cfx_live", "cfx_dead", "pre_req", "post_req",
            "win", "b_zero", "ncaution", "old_side", "CLASS", "ledger_ids", "unledgered", "FLAGS", "options"]
    with open(os.path.join(outdir, "degeneracy_table.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for S in table:
            g = lambda k: S.get(k, "")
            def cnt(k, a): return S[k][a] if k in S else ""
            w.writerow([g("config"), g("component"), g("dataset"), g("cls"), g("rc_expected"), g("rc"),
                        g("vcvstate"), g("properties"), g("N") and int(S["N"]), g("N_clust") and int(S["N_clust"]),
                        g("nsg"), cnt("main", "live"), cnt("main", "infdeg"), cnt("main", "void"),
                        cnt("sg", "live"), (S["sg"]["infdeg"] + S["sg"]["void"]) if "sg" in S else "",
                        cnt("pw", "live"), (S["pw"]["infdeg"] + S["pw"]["void"]) if "pw" in S else "",
                        cnt("gt", "live"), (S["gt"]["infdeg"] + S["gt"]["void"]) if "gt" in S else "",
                        cnt("cfx", "live"), (S["cfx"]["infdeg"] + S["cfx"]["void"]) if "cfx" in S else "",
                        g("pre_req"), g("post_req"),
                        ";".join(f"{k}={v[0]}p{v[1]}" for k, v in sorted(S.get("win", {}).items())),
                        len(S.get("b_zero", [])), "" if S.get("ncaution") is None else S["ncaution"],
                        g("old"), g("class"), ";".join(S.get("ids", [])),
                        ";".join(f"{k}:{n}" for k, n in S.get("uncovered", [])),
                        " ".join(S["flags"]), g("options")])
    # ---- cautions.csv (+ old side and diff, F9)
    # header-only when no shard logs were given: a 0 there must never read as "no cautions"
    write_cautions_csv(os.path.join(outdir, "cautions.csv"), names if digest_new is not None else [], digest_new or {})
    n_changed = None
    if digest_old is not None and digest_new is not None:
        write_cautions_csv(os.path.join(outdir, "cautions_old.csv"), names, digest_old)
        up = down = 0
        with open(os.path.join(outdir, "caution_diff_old_new.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["config", "n_old", "n_new", "delta"])
            for cfg in names:
                no, nn = n_caution(digest_old.get(cfg, [])), n_caution(digest_new.get(cfg, []))
                if no != nn:
                    w.writerow([cfg, no, nn, nn - no])
                    if nn > no: up += 1
                    else: down += 1
        n_changed = up + down
        # a permanent block when old = main = v1.3 (critique F9): one NOTE line, the CSV has the rows
        notes.append(f"caution counts old-vs-new: {n_changed} configs changed (up {up}, down {down}); "
                     f"total old {sum(n_caution(v) for v in digest_old.values())}, new "
                     f"{sum(n_caution(v) for v in digest_new.values())}; see {outdir}/caution_diff_old_new.csv")
    # ---- inverted per-option view
    counts = {k: [0, 0, 0, 0, 0] for k in opts}  # total / clean / clean good-data (= coverage_check's column) / clean good-data AND live / n-a (surface not applicable in that mode)
    live_examples = collections.defaultdict(list); dead_examples = collections.defaultdict(list)
    for S in table:
        clean = (S.get("rc") == "0")
        good = S.get("cls") in GOOD_DATA
        for opt in S["opts"]:
            if opt.startswith("?") or opt not in counts: continue
            counts[opt][0] += 1
            if clean:
                counts[opt][1] += 1
                if good:
                    counts[opt][2] += 1
                    ok, det = surface(opt, S)
                    if ok is None:
                        counts[opt][4] += 1   # surface not applicable in this mode: neither live nor dead
                    elif ok:
                        counts[opt][3] += 1; live_examples[opt].append(S["config"])
                    else:
                        dead_examples[opt].append(f"{S['config']}[{det or ' '.join(S['flags'])[:60]}]")
    with open(os.path.join(outdir, "option_inverted.txt"), "w") as f:
        f.write("option coverage: total / clean-rc / clean on wellbehaved+stress (= coverage_check's 3rd column) / ... AND live on the option's own surface [/ n-a: b-only mode, surface not applicable]\n")
        degen_only = []
        for k in sorted(counts, key=lambda x: (counts[x][3], counts[x][2], x)):
            t, c, g, l, na = counts[k]
            mark = ""
            if t == 0: mark = "  (zero)"
            elif c == 0: mark = "  <-- ERR-ONLY"
            elif g == 0: mark = "  <-- PATHO-ONLY"
            elif l == 0: mark = "  <-- DEGEN-ONLY (every clean good-data config is degenerate on this option's surface)"; degen_only.append(k)
            elif l + na < g: mark = f"  ({g - l - na} of {g - na} degenerate on surface)"
            f.write(f"  {k:28s} {t:4d} /{c:4d} /{g:4d} /{l:4d}" + (f" / n-a {na}" if na else "") + f"{mark}\n")
        f.write("\nDEGEN-ONLY: " + (", ".join(degen_only) if degen_only else "none") + "\n")
        f.write("\nper-option degenerate configs (clean rc, good data, surface not live):\n")
        for k in sorted(dead_examples):
            f.write(f"  {k} ({len(dead_examples[k])}): " + "; ".join(dead_examples[k][:12]) +
                    (" ..." if len(dead_examples[k]) > 12 else "") + "\n")
        # ---- pairwise inversion: an option PAIR that co-occurs in clean good-data configs but is
        # never live on BOTH surfaces at once is a hidden gap the per-option view cannot see
        # (e.g. timetrends x cluster: each option is live somewhere, never together).
        pair = collections.defaultdict(lambda: [0, 0]); pair_dead = collections.defaultdict(list)
        struct = {"att", "safety", "esprelength", "espostlength", "makeplots", "dots", "verbose", "replace", "saveplots"}
        for S in table:
            if S.get("rc") != "0" or S.get("cls") not in GOOD_DATA: continue
            ol = sorted(o for o in S["opts"] if not o.startswith("?") and o in counts and o not in struct)
            ok = {o: surface(o, S)[0] for o in ol}
            ol = [o for o in ol if ok[o] is not None]
            for i in range(len(ol)):
                for j in range(i + 1, len(ol)):
                    k = (ol[i], ol[j]); pair[k][0] += 1
                    if ok[ol[i]] and ok[ol[j]]: pair[k][1] += 1
                    else: pair_dead[k].append(S["config"])
        dead_pairs = sorted((k for k in pair if pair[k][1] == 0), key=lambda k: (-pair[k][0], k))
        f.write("\nDEGEN-ONLY PAIRS (option pairs co-occurring in clean good-data configs, never live on both surfaces at once;"
                " structural options att/safety/es*length/makeplots/dots/verbose/replace/saveplots excluded):\n")
        for k in dead_pairs:
            f.write(f"  {k[0]} x {k[1]:24s} {pair[k][0]:3d} configs, 0 live: " + "; ".join(pair_dead[k][:6]) +
                    (" ..." if len(pair_dead[k]) > 6 else "") + "\n")
        if not dead_pairs: f.write("  none\n")
    # ---- summary.txt
    n = len(table); nclean = sum(1 for S in table if S.get("rc") == "0")
    nodump = [S["config"] for S in table if S.get("rc") == "NODUMP"]
    cc_good = collections.Counter(S["class"] for S in table if S.get("rc") == "0" and S.get("cls") in GOOD_DATA)
    with open(os.path.join(outdir, "summary.txt"), "w") as f:
        fc = collections.Counter()
        for S in table:
            for fl in S["flags"]: fc[re.sub(r"\(.*", "", fl)] += 1
        f.write(f"configs {n}, clean rc {nclean}, no dump {len(nodump)}\n")
        f.write("flag histogram (all configs):\n")
        for k, v in fc.most_common(): f.write(f"  {k:32s} {v}\n")
        cc = collections.Counter((S.get("class", "?"), S.get("cls")) for S in table if S.get("rc") == "0")
        f.write("\nclass histogram over clean-rc configs (degeneracy class x data class):\n")
        for k, v in sorted(cc.items()): f.write(f"  {k[0]:16s} {k[1]:13s} {v}\n")
        for klass in ("HIDDEN", "PART_RULED", "RULED", "BY_CONSTRUCTION", "INFO"):
            f.write(f"\n== {klass}: clean-rc configs on wellbehaved/stress data ==\n")
            for S in table:
                if S.get("rc") == "0" and S.get("cls") in GOOD_DATA and S.get("class") == klass:
                    f.write(f"  {S['config']:36s} {S['cls']:12s} old={S.get('old',''):10s} ledger={','.join(S.get('ids', [])) or '-':28s} {' '.join(S['flags'])}\n")
        f.write("\n== pathological-class configs (info only) ==\n")
        for S in table:
            if S.get("rc") == "0" and S.get("cls") == "pathological":
                f.write(f"  {S['config']:36s} old={S.get('old',''):10s} {' '.join(S['flags']) or 'CLEAN'}\n")
        if nodump:
            f.write("\n== configs with NO dump ==\n  " + " ".join(nodump) + "\n")
        # per-dataset view: which surfaces are ever live on each dataset
        f.write("\n== per-dataset live surfaces (clean-rc configs) ==\n")
        ds = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
        for S in table:
            if S.get("rc") != "0": continue
            for opt in S["opts"]:
                if opt.startswith("?") or opt not in opts: continue
                ok, _ = surface(opt, S)
                ds[S["dataset"]][opt][0] += 1
                if ok: ds[S["dataset"]][opt][1] += 1
        for dsn in sorted(ds):
            dead = [f"{o}({v[1]}/{v[0]})" for o, v in sorted(ds[dsn].items()) if v[1] < v[0]]
            never = [o for o, v in sorted(ds[dsn].items()) if v[1] == 0]
            f.write(f"  {dsn}: never-live={never}; partially-live={[x for x in dead if x.split('(')[0] not in never]}\n")
        # ado flags never posted by any dump (F10 stays visible)
        try:
            ado_src = open(one("--ado"), encoding="utf-8", errors="replace").read()
            posted = set()
            for S in table:
                posted.update(S.get("mats", [])); posted.update(S.get("macros", {}).keys())
                posted.update(S.get("scalars_seen", []))
            ado_flags = sorted(set(re.findall(r"ereturn (?:hidden )?(?:scalar|local|matrix) (perm_[A-Za-z0-9_]+|boot_[A-Za-z0-9_]+)", ado_src)))
            never = [x for x in ado_flags if x not in posted and not any(m.startswith(x) for m in posted)]
            f.write("\n== ado perm_*/boot_* e() names never posted by any dump (matrix/macro names only; scalars are not tracked here) ==\n  "
                    + (" ".join(never) if never else "none") + "\n")
        except OSError:
            pass
    # ---- expectations.csv: GENERATED seed in the expected-CSV format (never asserted)
    with open(os.path.join(outdir, "expectations.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["config", "class", "stems", "citation", "ncaution"])
        for S in table:
            if S.get("rc") != "0" or S.get("cls") not in GOOD_DATA: continue
            st = stems(S)
            if not st: continue
            cite = ";".join(S.get("ids", [])) or ("BY_CONSTRUCTION" if S["class"] == "BY_CONSTRUCTION" else "UNRULED")
            if S["class"] == "PART_RULED": cite += ";UNRULED: " + " ".join(f"{k}:{v}" for k, v in S.get("uncovered", []))
            w.writerow([S["config"], S["class"], " ".join(st), cite, "" if S.get("ncaution") is None else S["ncaution"]])
    wrote = f"wrote {outdir}/degeneracy_table.csv, option_inverted.txt, cautions.csv" + \
            (", cautions_old.csv, caution_diff_old_new.csv" if n_changed is not None else "") + ", summary.txt, expectations.csv"
    def summary_line(extra=""):
        cl = " ".join(f"{k}={cc_good.get(k, 0)}" for k in SCANNER_CLASSES)
        return (f"DEGEN_SUMMARY configs={n} clean_rc={nclean} nodump={len(nodump)} good_data_clean={sum(cc_good.values())} "
                f"[{cl}]" + (f" cautions_new={sum(n_caution(v) for v in digest_new.values())}" if digest_new is not None else "") + extra)
    if not expect:
        print(wrote)
        for x in notes: print("DEGEN_NOTE " + x)
        for x in warns: print("DEGEN_WARN " + x)
        print(summary_line())
        print("DEGEN_SCAN_ONLY: no --expect given; nothing asserted (expectations.csv written as a seed)")
        sys.exit(3)
    # ---- assertion --------------------------------------------------------------------
    # exp/bad were loaded before the per-config loop (the citations gate the ledger join)
    ledger_ids = {r["id"] for r in ledger}
    fails, stale = [], []
    by_name = {S["config"]: S for S in table}
    for S in table:
        cfg = S["config"]
        if S["rc"] == "NODUMP":
            # A2: never skip silently - a timed-out shard or an ONLY= subset run leaves no dump
            fails.append(f"{cfg}: NO DUMP in {dumps} (shard did not finish, or oldnew ran an ONLY= subset)"); continue
        good = S["cls"] in GOOD_DATA
        if S["rc"] == "0" and S["cls"] not in GOOD_DATA + ("pathological",):
            fails.append(f"{cfg}: data class '{S['cls']}' (not in data_classes.csv) - surface cannot be asserted"); continue
        if not good: continue          # pathological data: "correct or observable failure?" only
        if S["rc"] != S["rc_expected"]:
            fails.append(f"{cfg}: rc {S['rc']} but configs expect rc {S['rc_expected']} - surface not assessable"); continue
        if S["rc"] != "0": continue    # an expected error has no surface
        got = set(stems(S)); e = exp.get(cfg); want = e["stems"] if e else set()
        new = got - want; gone = want - got
        if new: fails.append(f"{cfg}: NEW degeneracy {sorted(new)} (allowed {sorted(want)}; flags {' '.join(S['flags'])})")
        if gone: stale.append(f"{cfg}: STALE expectation {sorted(gone)} (now shows {sorted(got) or 'nothing'})")
        if e:
            if not class_compatible(e["class"], S["class"]):
                # RF-12 (2026-09-13): a curated class moving is a FAILURE, not a warning. It
                # was a warning until a real one fired on real_sg during the Stage B triage and
                # the component still exited 0 - the contract this file asserts is that a
                # curated expectation still describes the run, and a class that has moved does
                # not. Kept in `warns` too so the DEGEN_WARN line still names it.
                warns.append(f"class moved: {cfg} expected={e['class']} got={S['class']} (update degeneracy_expected.csv)")
                classmoved.append(cfg)
            if e["class"] in ("RULED", "PART_RULED"):
                for cid in e["citation"].split(";"):
                    cid = cid.strip()
                    if cid and not cid.startswith("UNRULED") and not cid.startswith("PINNED_DEGENERATE") and cid not in ledger_ids:
                        # RF-12: a citation that resolves to nothing is a broken contract, same tier.
                        warns.append(f"dangling ledger citation: {cfg} cites '{cid}' which is not a ledger.csv id")
                        classmoved.append(cfg)
            if e["ncaution"] is not None and S.get("ncaution") is not None and e["ncaution"] != S["ncaution"]:
                warns.append(f"caution count moved: {cfg} expected={e['ncaution']} got={S['ncaution']}")
    # A3: rows whose config is gone, or which can never be asserted, are STALE (ledger discipline)
    for cfg in sorted(exp):
        if cfg not in by_name:
            stale.append(f"{cfg}: expectation for a config that no longer exists in " + "/".join(os.path.basename(c) for c in cfgs))
        else:
            S = by_name[cfg]
            if S["cls"] not in GOOD_DATA or S["rc_expected"] != "0":
                stale.append(f"{cfg}: expectation on a config outside the assertion universe (class={S['cls']}, rc_new={S['rc_expected']})")
    n_unruled = sum(1 for e in exp.values() if "UNRULED" in e["citation"])
    if n_unruled: warns.append(f"{n_unruled} degeneracy_expected.csv rows carry UNRULED citations (pending an author ruling)")
    # Ruled 2026-09-13: keep these as pins of their degenerate surface rather than re-hosting
    # or silently accepting them. They stay WARNed on purpose - the author's stated reason for
    # keeping them was that the warning keeps them visible rather than forgotten - but they are
    # no longer PENDING, so the wording above would be false for them.
    n_pinned = sum(1 for e in exp.values() if "PINNED_DEGENERATE" in e["citation"])
    if n_pinned: warns.append(f"{n_pinned} degeneracy_expected.csv rows are PINNED degenerate surfaces (ruled 2026-09-13: kept as pins, kept visible deliberately)")
    # RF-13 (2026-09-13): DEGENERACY.md states a census of this CSV in prose, and prose goes
    # stale in silence - it said 165 rows and 10 UNRULED within an hour of nine rows being
    # added, and nothing in the suite noticed because no component checks a document's
    # arithmetic against the data it documents. The census sentence is machine-read here and
    # a disagreement is a FAILURE, on the same footing as a stale expectation row.
    docpath = os.path.join(os.path.dirname(os.path.abspath(expect)), "DEGENERACY.md")
    if not os.path.exists(docpath):
        warns.append(f"DEGENERACY.md not found beside {os.path.basename(expect)}: its census is unchecked")
    else:
        with open(docpath, encoding="utf-8") as fh: doc = " ".join(fh.read().split())
        mdoc = re.search(r"holds\s+(\d+)\s+rows:\s*([^.]*)\.", doc)
        if not mdoc:
            warns.append("DEGENERACY.md has no 'holds N rows: ...' census sentence, so its arithmetic is unchecked")
        else:
            cur = collections.Counter(e["class"] for e in exp.values())
            said = {k: int(n) for n, k in re.findall(r"(\d+)\s+(PART_RULED|BY_CONSTRUCTION|DATASET_DESIGN|DATA_EDGE|HIDDEN|RULED)\b", mdoc.group(2))}
            munr = re.search(r"(\d+)\s+of those rows carry an", doc)
            docbad = []
            if int(mdoc.group(1)) != len(exp): docbad.append(f"total says {mdoc.group(1)}, file has {len(exp)}")
            for k in sorted(set(said) | set(cur)):
                if said.get(k, 0) != cur.get(k, 0): docbad.append(f"{k} says {said.get(k, 0)}, file has {cur.get(k, 0)}")
            if munr and int(munr.group(1)) != n_unruled: docbad.append(f"UNRULED says {munr.group(1)}, file has {n_unruled}")
            if docbad:
                fails.append("DEGENERACY.md: its census disagrees with " + os.path.basename(expect) + " (" + "; ".join(docbad) + ") - recompute the note from the CSV")
    do = [k for k in counts if counts[k][2] > 0 and counts[k][3] == 0]
    if do: warns.append("degen-only options (clean good-data coverage is degenerate on their own surface): " + ", ".join(sorted(do)))
    if dead_pairs: warns.append("degen-only option PAIRS (never live on both surfaces at once): " + ", ".join(f"{a}x{b}" for a, b in dead_pairs))
    for x in bad: print("DEGEN_BADROW " + x)
    for x in fails: print("DEGEN_FAIL " + x)
    for x in stale: print("DEGEN_STALE " + x)
    for x in notes: print("DEGEN_NOTE " + x)
    for x in warns: print("DEGEN_WARN " + x)
    print(wrote)
    # RF-11 (2026-09-13): the class bracket above prints SCANNER classes only, so curated
    # DATASET_DESIGN rows are counted under HIDDEN and the line cannot be reconciled
    # against degeneracy_expected.csv. curated_* below reports the file's own tally.
    print(summary_line(f" expected_rows={len(exp)}"
        f" curated=[" + " ".join(f"{k}={v}" for k,v in sorted(__import__("collections").Counter(e['class'] for e in exp.values()).items())) + "]"
        f" fail={len(fails)} stale={len(stale)} badrow={len(bad)} classmoved={len(classmoved)} warn={len(warns)}"))
    if fails or stale or bad or classmoved:
        print(f"DEGEN_CHECK_FAIL: {len(fails)} fail, {len(stale)} stale, {len(bad)} bad rows, {len(classmoved)} class-moved/dangling"); sys.exit(1)
    print("DEGEN_CHECK_OK")

if __name__ == "__main__":
    main()
