#!/usr/bin/env python3
"""Assert option-surface coverage of the suite configs.

Usage: coverage_check.py --ado ../wooldidinstall/wooldid.ado --configs configs.csv
                         [--configs configs_ext.csv ...] [--allow-zero fastnose,...]
                         [--data-classes data_classes.csv]

Parses the wooldid syntax line (the worker's), derives each option's minimum abbreviation
from its capitalization, then scans every config's option string (plus weight/if columns)
and reports:
  - per-option config counts (canonical name),
  - OPTIONS WITH ZERO COVERAGE (failures unless --allow-zero),
  - weight-class and if-condition coverage.

With --data-classes (the dataset classification table; see the README's "Role of
pathological data"): joins each config to its dataset's class and additionally reports
  - PATHO-ONLY options (warn): clean-rc coverage exists ONLY on pathological-class data,
    i.e. the feature has never been shown working where a trustworthy answer is expected;
  - options with no pathological-data exposure (informational — NOT a defect; the
    effort-allocation principle says not every feature needs a pathological probe);
  - sync failures: a config whose dataset (or inline-prep config name) is missing from
    the table, or a table row no config references, both fail the run so the marking
    cannot silently go stale.
Exit 0 iff every option not in --allow-zero has at least one config and (when given)
the data-class table is in sync.
"""
import sys, csv, re

def parse_syntax_options(ado_path):
    """Return {canonical_lowercase_name: min_abbrev_lowercase} from the first full syntax line."""
    txt = open(ado_path, encoding="utf-8", errors="replace").read()
    m = re.search(r"^\s*syntax\s+\[varlist.*?,\s*\[(.*?)\]\s*$",
                  txt, re.M | re.S)
    if not m:
        sys.exit("could not locate the wooldid syntax line")
    body = m.group(1)
    opts = {}
    for tok in re.finditer(r"([A-Za-z_][A-Za-z_0-9]*)(?:\([^)]*\))?", body):
        name = tok.group(1)
        upper = re.match(r"[A-Z_0-9]+", name)
        abbrev = upper.group(0).lower() if upper and name[0].isupper() else name.lower()
        opts[name.lower()] = abbrev
    return opts

def match_option(token, opts):
    """Stata abbreviation rule: token must extend the min abbrev and prefix the full name."""
    t = token.lower()
    for full, ab in opts.items():
        if t.startswith(ab) and full.startswith(t):
            return full
    return None

def config_tokens(optstring):
    """Split an option string on spaces at paren depth 0; strip any (args)."""
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
    return [t.split("(")[0] for t in toks if t]

def load_data_classes(path):
    """Return {name: class} from the classification table. Names are dataset names as
    spelled in the configs' dataset column, or config:<name> for inline-prep configs
    (and for deliberate per-config overrides, which win over the dataset lookup)."""
    classes = {}
    for r in csv.DictReader(open(path, newline="")):
        name = (r.get("name") or "").strip()
        if not name or name.startswith("#"):
            continue
        cls = (r.get("class") or "").strip()
        if cls not in ("wellbehaved", "stress", "pathological"):
            sys.exit(f"data_classes: bad class '{cls}' for '{name}'")
        if name in classes:
            sys.exit(f"data_classes: duplicate row for '{name}'")
        classes[name] = cls
    return classes

def main():
    args = sys.argv[1:]
    ado = args[args.index("--ado") + 1]
    cfg_files = [args[i + 1] for i, a in enumerate(args) if a == "--configs"]
    allow = set()
    if "--allow-zero" in args:
        allow = set(args[args.index("--allow-zero") + 1].split(","))
    classes = None
    if "--data-classes" in args:
        classes = load_data_classes(args[args.index("--data-classes") + 1])
    opts = parse_syntax_options(ado)
    counts = {k: 0 for k in opts}       # any config
    clean = {k: 0 for k in opts}        # configs expected to run rc 0 on the NEW build
    goodclean = {k: 0 for k in opts}    # clean-rc configs on wellbehaved/stress data
    pathoexp = {k: 0 for k in opts}     # configs on pathological-class data (any rc)
    wclasses = {"aw": 0, "pw": 0, "none": 0}
    class_hist = {"wellbehaved": 0, "stress": 0, "pathological": 0}
    if_count = 0
    nconfigs = 0
    unknown = set()
    unclassified = []
    referenced = set()
    for f in cfg_files:
        for r in csv.DictReader(open(f, newline="")):
            if not r.get("config") or r["config"].startswith("#"):
                continue
            nconfigs += 1
            is_clean = (r.get("rc_new", "").strip() or "0") == "0"
            cls = None
            if classes is not None:
                over = "config:" + r["config"].strip()
                ds = (r.get("dataset") or "").strip()
                if over in classes:
                    cls, key = classes[over], over
                elif ds and ds in classes:
                    cls, key = classes[ds], ds
                else:
                    unclassified.append(r["config"].strip())
                    key = None
                if key:
                    referenced.add(key)
                    class_hist[cls] += 1
            for tok in config_tokens(r["options"]):
                full = match_option(tok, opts)
                if full:
                    counts[full] += 1
                    if is_clean:
                        clean[full] += 1
                        if cls in ("wellbehaved", "stress"):
                            goodclean[full] += 1
                    if cls == "pathological":
                        pathoexp[full] += 1
                else:
                    unknown.add(tok)
            w = (r.get("weight") or "").strip()
            wclasses["aw" if w.startswith("aw") else "pw" if w.startswith("pw") else "none"] += 1
            if (r.get("ifcond") or "").strip():
                if_count += 1
    zero = sorted(k for k, v in counts.items() if v == 0 and k not in allow)
    dead_allow = sorted(k for k in allow if k in counts and counts[k] > 0)
    err_only = sorted(k for k in counts
                      if counts[k] > 0 and clean[k] == 0 and k not in allow)
    patho_only = sorted(k for k in counts
                        if clean[k] > 0 and goodclean[k] == 0
                        and classes is not None and k not in allow)
    print(f"configs scanned: {nconfigs}")
    print(f"weights: aw={wclasses['aw']} pw={wclasses['pw']} unweighted={wclasses['none']}; "
          f"if-restricted configs: {if_count}")
    if classes is not None:
        print(f"data classes: wellbehaved={class_hist['wellbehaved']} "
              f"stress={class_hist['stress']} pathological={class_hist['pathological']}")
    print("option coverage (canonical name = total / clean-rc configs"
          + (" / clean on wellbehaved+stress data" if classes is not None else "") + "):")
    for k in sorted(counts, key=lambda x: (clean[x], counts[x], x)):
        mark = ("  <-- ZERO" if k in zero else
                "  (allowed zero)" if counts[k] == 0 else
                "  <-- ERR-ONLY (never exercised on a clean run)" if k in err_only else
                "  <-- PATHO-ONLY (clean runs only on pathological data)" if k in patho_only
                else "")
        gc = f" /{goodclean[k]:4d}" if classes is not None else ""
        print(f"  {k:28s} {counts[k]:4d} /{clean[k]:4d}{gc}{mark}")
    if unknown:
        print("UNRECOGNIZED option tokens (check spelling/abbreviation): " + ", ".join(sorted(unknown)))
    fail = False
    if zero:
        print(f"COVERAGE_FAIL: {len(zero)} option(s) with zero configs: " + ", ".join(zero))
        fail = True
    if dead_allow:
        print(f"COVERAGE_FAIL: --allow-zero entries that actually HAVE coverage "
              f"(remove them so future coverage loss is caught): " + ", ".join(dead_allow))
        fail = True
    if classes is not None:
        if unclassified:
            print(f"COVERAGE_FAIL: {len(unclassified)} config(s) whose dataset/prep has no "
                  f"data_classes row: " + ", ".join(sorted(set(unclassified))))
            fail = True
        stale = sorted(set(classes) - referenced)
        if stale:
            print(f"COVERAGE_FAIL: data_classes row(s) no config references (stale): "
                  + ", ".join(stale))
            fail = True
    if err_only:
        print("COVERAGE_WARN err-only options (exercised, but never on a clean run): "
              + ", ".join(err_only))
    if patho_only:
        print("COVERAGE_WARN patho-only options (clean coverage exists ONLY on "
              "pathological-class data — feature never shown working where a good "
              "answer is expected): " + ", ".join(patho_only))
    if classes is not None:
        nopatho = sorted(k for k in counts if counts[k] > 0 and pathoexp[k] == 0)
        print("info: options with no pathological-data exposure (not a defect; see the "
              "effort-allocation principle): " + (", ".join(nopatho) if nopatho else "none"))
    if fail:
        sys.exit(1)
    print("COVERAGE_OK")

main()
