#!/usr/bin/env python3
"""
# REGENERATION POLICY: this miner writes configs_ext.GEN.csv and preps_ext.GEN.json —
# never the curated configs_ext.csv / preps.json, which carry hand-verified rc pins and
# hand-authored prep keys. To adopt regenerated output, DIFF it against the curated
# files and merge deliberately (the curated rc pins win unless re-verified empirically).
External-config MINER for the wooldid dual-build test suite.

Reads the external colleagues' bundle
  wooldid_external/wooldid_test_comparison_code
(manifests + case .do files) and emits runner configs:

  configs_ext.gen.csv  columns:
    config,component,dataset,weight,ifcond,prep,varlist,options,rc_old,rc_new,msg_old,slow
  preps_ext.gen.json       prep key -> list of Stata lines (run after loading dataset,
                   before the wooldid call)

Design notes
------------
* Battery (216) and lattice (24) configs are constructed by replicating the
  template logic of tests/cases/binary_batteries.do, continuous_batteries.do
  and continuous_extensions.do::wdbt_parse_lattice_case/wdbt_case_lattice,
  enumerated from tests/manifest/battery-cases.csv / lattice-cases.csv.
* Authored cases are STATICALLY extracted from the case .do files: the miner
  joins /// continuations, finds each case's if/else-if block inside the
  wdbt_case_* dispatch programs, resolves trivial `local X "Y"` macros, and
  pulls out literal wooldid calls plus preceding data-mutation lines (preps).
  Cases whose call cannot be resolved statically are SKIPPED and listed.
* Multi-call cases are split into <config>_c1, <config>_c2 ... (one config per
  wooldid call; preps are cumulative). Only rc_old==0 cases may be split.
* SKIPPED case families (by design): harness (serializer self-tests),
  pilot (redundant smoke calls), oracle (known-value/metamorphic assertion
  library -- port as future work).
* rc_old = the bundle's expected_rc (the OLD-build v1.3 contract, verified by
  their b003 double-run). rc_new = rc_old unless the intentional-divergence
  ledger (/tmp/wooldid_suite_survey/divergence_ledger_draft.md) documents an
  rc change; where the ledger is silent on a contract-deviation site, rc_new
  is left BLANK = assert later.
* dataset/setup mapping is transcribed from tests/cases/registry.do (dispatch)
  and coverage_extensions.do::wdbt_dispatch_coverage10; the `use ...dta` line
  and any setup-time mutations are parsed from the wdbt_setup_* programs
  themselves. The miner asserts every included case got a dataset.

Stdlib-only; re-runnable; degrades gracefully (note + exit 0) if the bundle
is absent.
"""

import csv
import json
import re
import sys
from pathlib import Path

BUNDLE = Path(
    "/Users/hegland/Documents/Research/wooldid_dev/wooldid_external/"
    "wooldid_test_comparison_code"
)
OUTDIR = Path(__file__).resolve().parent

# --------------------------------------------------------------------------
# policy tables
# --------------------------------------------------------------------------

SKIP_FAMILIES = {"harness", "pilot", "oracle"}
SKIP_REASON = {
    "harness": "serializer self-test, not a wooldid case",
    "pilot": "smoke pilot, redundant with batteries",
    "oracle": "oracle assertion library -- future work",
}

# rc_new: divergence-ledger-driven overrides (ledger row ids in comments)
RC_NEW_OVERRIDE = {
    "contract_fractional_subgroup": "411",   # B-B20 7f7f73c: cryptic r198 -> clear r411
    "contract10_inconsistent_ttre": "459",   # B-F13 c4d7e62: ttre varies within i, safety(on) -> r459
    "contract10_if_outside_negative": "0",   # B-B23 8d7bd01: bad values excluded via if -> rc 0
    "val10_fractional_ttre": "459",          # B-S5scr 671dba1: in-window fractional onset -> r459
}
# contract-deviation sites the ledger does not resolve: leave rc_new blank
RC_NEW_BLANK = {
    "contract_bad_elasticity",
    "contract_fractional_lattice",
    "contract10_nonpositive_time",
    "contract10_fractional_time",
    "contract10_poly_zero",
    "contract10_bad_safety",
    "contract10_poisexp_without_poisson",
    "contract10_rcs_factor_fe",
    "path10_ppml_separation",              # A7-B10 may turn separation crash into success
}

# case_id -> setup program (transcribed from registry.do + wdbt_dispatch_coverage10)
SETUP_BY_CASE = {
    # registry.do lines 80-86 (binary extensions; default binary_balanced)
    "bin_ext_fine": "wdbt_setup_binary_fine",
    "bin_ext_mixed_imbalanced": "wdbt_setup_binary_imbalanced",
    "bin_ext_cleanmatrices": "wdbt_setup_binary_sentinel",
    # registry.do lines 112-122 (continuous extensions; default continuous_control)
    "cont_ext_poly2": "wdbt_setup_continuous_dense",
    "cont_ext_poly3": "wdbt_setup_continuous_dense",
    "cont_ext_cfx": "wdbt_setup_continuous_dense",
    "cont_ext_cfx_subgroups": "wdbt_setup_continuous_dense",
    "cont_ext_interactive": "wdbt_setup_continuous_dense",
    "cont_ext_custom_lattice": "wdbt_setup_continuous_dense",
    "cont_ext_weighted_lattice": "wdbt_setup_continuous_dense",
    "cont_ext_untreatedzero": "wdbt_setup_cont_untreated0_dense",
    "cont_ext_imbalanced": "wdbt_setup_cont_imbalanced_dense",
    # registry.do line 144 (sanitation)
    "san_input_missing": "wdbt_setup_sanitation",
    # coverage_extensions.do wdbt_dispatch_coverage10 lines 755-777
    # (default wdbt_setup_cov_binary)
    "cov10_rcs_pweight_ddd": "wdbt_setup_cov_survey",
    "cov10_rcs_twoway_cluster": "wdbt_setup_cov_survey",
    "contract10_rcs_factor_fe": "wdbt_setup_cov_survey",
    "oracle10_pweight_target": "wdbt_setup_cov_survey",       # skipped (oracle)
    "cov10_binary_mixed_grain": "wdbt_setup_cov_mixed",
    "oracle10_aggregation_exact": "wdbt_setup_cov_aggregation",  # skipped
    "oracle10_im_unequal_clusters": "wdbt_setup_cov_im_unequal", # skipped
    "cov10_cont_within_vary": "wdbt_setup_cov_dense_within",
    "cov10_cont_single_sg": "wdbt_setup_cov_dense",
    "cov10_cont_panel_fe": "wdbt_setup_cov_dense",
    "cov10_cont_impvalues": "wdbt_setup_cov_dense",
    "cov10_cont_elasticity_eydx": "wdbt_setup_cov_dense",
    "cov10_cont_elasticity_dyex": "wdbt_setup_cov_dense",
    "cov10_cfx_inference": "wdbt_setup_cov_dense",
    "cov10_rare_binary": "wdbt_setup_cov_dense",
    "cov10_skewed_cost": "wdbt_setup_cov_dense",
    "oracle10_ppml_xb": "wdbt_setup_cov_dense",                # skipped
    "oracle10_ppml_exp": "wdbt_setup_cov_dense",               # skipped
    "contract10_ppml_d_collision": "wdbt_setup_cov_dense",
    "contract10_poisson_converged": "wdbt_setup_cov_dense",
    "path10_ppml_separation": "wdbt_setup_cov_dense",
    "path10_cont_constant_dose": "wdbt_setup_cov_dense",
    "contract10_poly_zero": "wdbt_setup_cov_dense",
    "val10_fractional_lattice_count": "wdbt_setup_cov_dense",
}


def default_setup(case_id: str) -> str:
    """Defaults per registry.do dispatch branches."""
    if case_id.startswith("bin_ext_"):
        return "wdbt_setup_binary_balanced"
    if case_id.startswith("cont_ext_"):
        return "wdbt_setup_continuous_control"
    if re.match(r"^(cov10|path10|contract10|val10|oracle10)_", case_id):
        return "wdbt_setup_cov_binary"
    # val_*, char_*, contract_* (phase 7), sanitation fallback
    return "wdbt_setup_binary_balanced"


# case_id -> the wdbt_case_* program that contains its block
def case_program(case_id: str) -> str:
    if case_id == "san_input_missing":
        return "wdbt_case_sanitation"
    if case_id.startswith("bin_ext_"):
        return "wdbt_case_binary_extension"
    if case_id.startswith("cont_ext_"):
        return "wdbt_case_continuous_extension"
    if re.match(r"^(cov10|path10|contract10|val10|oracle10)_", case_id):
        return "wdbt_case_coverage10"
    if case_id.startswith("val_"):
        return "wdbt_case_validation"
    # char_* and phase-7 contract_* live in pathologies.do
    return "wdbt_case_characterization"


CASE_FILES = [
    "binary_extensions.do",
    "continuous_extensions.do",
    "coverage_extensions.do",
    "validation.do",
    "pathologies.do",
    "pilots.do",  # setup programs only
]

MUTATION_PREFIXES = (
    "replace ", "generate ", "gen ", "recast ", "expand ",
    "keep if ", "keep in ", "drop if ", "drop in ",
    "set emptycells", "twoway ", "egen ",
)
IGNORABLE_PREFIXES = (
    "matrix ", "scalar ", "local ", "quietly ", "assert ", "confirm ",
    "count ", "summarize ", "capture ", "display ", "estimates ",
    "graph ", "mata", "isid ", "bysort ", "preserve", "restore",
    "tempvar ", "version ", "wdbt_assert_close",
)

# --------------------------------------------------------------------------
# .do file parsing
# --------------------------------------------------------------------------


def read_statements(path: Path):
    """Join /// continuations; strip blank/comment lines; return statements."""
    raw = path.read_text(encoding="utf-8", errors="replace").splitlines()
    stmts, buf = [], ""
    for line in raw:
        s = line.rstrip()
        joined = False
        if s.rstrip().endswith("///"):
            buf += s.rstrip()[:-3].rstrip() + " "
            joined = True
        if not joined:
            stmt = (buf + s).strip()
            buf = ""
            if not stmt or stmt.startswith("//") or stmt.startswith("*"):
                continue
            stmts.append(stmt)
    return stmts


def program_bodies(stmts):
    """Map program name -> list of statements (between define and end)."""
    bodies, name, body = {}, None, []
    for st in stmts:
        m = re.match(r"^program define (\S+)", st)
        if m:
            name, body = m.group(1).rstrip(","), []
            continue
        if st == "end" and name:
            bodies[name] = body
            name, body = None, []
            continue
        if name:
            body.append(st)
    return bodies


def find_case_block(body, case_id):
    """Return the statements of case_id's if/else-if block inside a dispatch
    program body, or None. Single-branch programs (no case_id dispatch)
    return the whole body."""
    if not any(re.match(r"^(if|else if) ", st) and "case_id" in st for st in body):
        return list(body)  # e.g. wdbt_case_sanitation
    i = 0
    while i < len(body):
        st = body[i]
        m = re.match(r"^(?:if|else if) (.*)\{$", st)
        if m and "case_id" in m.group(1):
            ids = re.findall(r'"([A-Za-z0-9_]+)"', m.group(1))
            depth, j, block = 1, i + 1, []
            while j < len(body) and depth > 0:
                depth += body[j].count("{") - body[j].count("}")
                if depth > 0:
                    block.append(body[j])
                j += 1
            if case_id in ids:
                return block
            i = j
            continue
        i += 1
    return None


class DynamicCase(Exception):
    pass


def extract_calls(block, case_id):
    """Walk a case block; return (calls, note) where calls is a list of dicts
    {varlist, weight, ifcond, preps}. Raises DynamicCase if not static."""
    locals_map = {}
    preps = []       # cumulative mutations before the *next* call
    calls = []
    depth = 0
    for st in block:
        opens, closes = st.count("{"), st.count("}")
        if depth > 0:
            # inside a nested braced sub-block: must be pure bookkeeping
            if st.startswith("wooldid ") or st.startswith(MUTATION_PREFIXES):
                raise DynamicCase(f"mutation/call inside nested block: {st[:60]}")
            depth += opens - closes
            continue
        if opens > closes:
            # nested block header (foreach/forvalues/if) -- bookkeeping only
            depth = opens - closes
            continue

        # one-line `if <cond> <stmt>`: resolve case_id conditions
        m = re.match(r'^if\s+(.+?)\s+((?:local|assert|display)\s.+)$', st)
        if m and "case_id" in m.group(1):
            ids = re.findall(r'"([A-Za-z0-9_]+)"', m.group(1))
            if case_id not in ids:
                continue
            st = m.group(2)
        elif re.match(r"^if\s+", st):
            # non-case-id one-liner: fine if the inner statement is bookkeeping
            if re.search(r"\s(assert|display|local)\s", st):
                continue
            raise DynamicCase(f"unresolvable one-line if: {st[:60]}")

        m = re.match(r'^local\s+(\w+)\s+"([^"]*)"\s*$', st)
        if m:
            locals_map[m.group(1)] = m.group(2)
            continue

        if st.startswith("wooldid ") or st.startswith("wooldid,"):
            line = st
            for k, v in locals_map.items():
                line = line.replace(f"`{k}'", v)
            if "`" in line or "$" in line:
                raise DynamicCase(f"unresolved macro in call: {line[:80]}")
            calls.append({"line": re.sub(r"\s+", " ", line).strip(),
                          "preps": list(preps)})
            continue

        if st.startswith(MUTATION_PREFIXES):
            preps.append(re.sub(r"\s+", " ", st))
            continue
        if st.startswith(IGNORABLE_PREFIXES) or st in ("preserve", "restore"):
            continue
        raise DynamicCase(f"unclassified statement: {st[:60]}")
    if not calls:
        raise DynamicCase("no wooldid call found")
    return calls


def split_call(line):
    """'wooldid <varlist> [w=x] if <cond>, <opts>' -> (varlist, weight, if, opts)"""
    body = line[len("wooldid"):].strip()
    # find first comma outside () and []
    depth = 0
    cut = None
    for i, ch in enumerate(body):
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif ch == "," and depth == 0:
            cut = i
            break
    left, opts = (body, "") if cut is None else (body[:cut], body[cut + 1:])
    left, opts = left.strip(), re.sub(r"\s+", " ", opts).strip()
    weight = ""
    m = re.search(r"\[(\w+)\s*=\s*([^\]]+)\]", left)
    if m:
        weight = f"{m.group(1)}={m.group(2).strip()}"
        left = (left[:m.start()] + left[m.end():]).strip()
    ifcond = ""
    m = re.search(r"\bif\s+(.+)$", left)
    if m:
        ifcond = m.group(1).strip()
        left = left[:m.start()].strip()
    varlist = re.sub(r"\s+", " ", left)
    return varlist, weight, ifcond, opts


def parse_setups(bodies):
    """wdbt_setup_* -> (dataset or None, setup prep lines)."""
    setups = {}
    for name, body in bodies.items():
        if not name.startswith("wdbt_setup_"):
            continue
        dataset, preps = None, []
        for st in body:
            m = re.search(r"generated/([\w]+\.dta)", st)
            if m and st.lstrip().startswith("use"):
                dataset = m.group(1)
                continue
            if st.startswith(("version ", "clear", "quietly ", "scalar ",
                              "isid ", "assert ")):
                continue
            if st.startswith("bysort ") and ": assert" in st:
                continue
            if st.startswith(("set obs", "generate ", "gen ", "expand ",
                              "collapse ", "egen ", "preserve", "restore",
                              "keep ", "drop ")) and dataset is None:
                # in-memory generator (oracle setups) -- flagged by dataset None
                preps.append(st)
                continue
            preps.append(st)
        setups[name] = (dataset, preps)
    return setups


# --------------------------------------------------------------------------
# battery / lattice templates (replicated from the bundle case drivers)
# --------------------------------------------------------------------------

ALL_ESTIMANDS = "att att_i att_it att_itime"


def battery_config(case_id):
    """Replicates wdbt_parse_battery_case + wdbt_case_binary_battery /
    wdbt_case_continuous_battery exactly."""
    m = re.match(r"^(b[1-4]|c[1-2])_h([0-5])([0-5])_(p|f1|f2|f4)$", case_id)
    if not m:
        raise ValueError(f"bad battery id {case_id}")
    battery, pre, post, ref = m.groups()
    horizon = pre + post
    if horizon not in {"00", "10", "20", "30", "03", "11", "22", "33", "55"}:
        raise ValueError(f"unregistered horizon {case_id}")

    estimands = "att"
    if horizon == "55" and ref == "f4":
        estimands = "att att_it"
    if battery in ("b1", "b3", "c1") and horizon == "22" and ref == "p":
        estimands = ALL_ESTIMANDS
    if battery in ("b2", "b4", "c2") and horizon == "22" and ref == "f2":
        estimands = ALL_ESTIMANDS

    refopts = ""
    if ref != "p":
        refopts = f"esfixedbaseperiod esrelativeto(-{ref[1:]})"

    if battery == "b1":
        varlist = "y_single cohort t ttre_single"
    else:
        varlist = "y cohort t ttre"

    overlay = {
        "b1": "",
        "b2": "subgroups(sg_one)",
        "b3": "subgroups(sg_three) timetrends",
        "b4": "subgroups(sg_ttre)",
        "c1": "contreat(dose) contreatcontrols(basic) contreatwithin",
        "c2": "contreat(dose) contreatcontrols(basic) subgroups(sg_three) timetrends",
    }[battery]

    parts = [estimands, f"esprelength({int(pre)})", f"espostlength({int(post)})",
             refopts, overlay, "jointtests adversarial makeplots safety(off)"]
    options = " ".join(p for p in parts if p)

    dataset = {"b": "binary_balanced.dta"}.get(battery[0])
    if battery == "c1":
        dataset = "continuous_within_compact.dta"   # wdbt_setup_continuous_within
    elif battery == "c2":
        dataset = "continuous_control_compact.dta"  # wdbt_setup_continuous_control

    has_sg = "subgroups(" in options
    slow = 1 if has_sg and (int(pre) + int(post) >= 6 or estimands == ALL_ESTIMANDS) else 0
    return varlist, options, dataset, slow


LATTICE_TYPES = {
    "all": "all", "btr": "baseTreated", "bpost": "baseTreatedPost",
    "bpre": "baseTreatedPre", "bref": "baseTreatedRef", "cohort": "byCohort",
    "time": "byTime", "cohorttime": "byCohortByTime", "treated": "byTreated",
    "treatedtime": "byTreatedByTime", "ttre": "byTtre", "ttretime": "byTtreByTime",
}


def lattice_config(token, bins):
    ltype = LATTICE_TYPES[token]
    options = (f"att contreat(dose_common) contreatcontrols(basic) "
               f"contreatpoly(1) lattice({ltype} {bins}) "
               f"esprelength(1) espostlength(1) robust safety(off)")
    return "y_cont_common cohort t ttre", options, "continuous_control_dense.dta"


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def main():
    if not BUNDLE.is_dir():
        print(f"NOTE: external bundle not found at {BUNDLE}; nothing to mine. "
              f"configs_ext.gen.csv left untouched.")
        return 0

    manifest = BUNDLE / "tests" / "manifest"
    cases = list(csv.DictReader(open(manifest / "cases.csv", encoding="utf-8-sig")))
    batteries = list(csv.DictReader(open(manifest / "battery-cases.csv", encoding="utf-8-sig")))
    lattices = list(csv.DictReader(open(manifest / "lattice-cases.csv", encoding="utf-8-sig")))
    # messages.csv contains unquoted commas inside expected_text (e.g.
    # "must consist of positive, non-zero integers") -- rejoin extra fields.
    messages = {}
    with open(manifest / "messages.csv", encoding="utf-8-sig") as fh:
        rd = csv.reader(fh)
        next(rd)  # header
        for rec in rd:
            if rec:
                messages[rec[0]] = ",".join(rec[1:])

    # parse all case files
    bodies = {}
    for f in CASE_FILES:
        bodies.update(program_bodies(read_statements(BUNDLE / "tests" / "cases" / f)))
    setups = parse_setups(bodies)

    rows, preps_out, skipped = [], {}, []

    def poisson_slow(options, rc_old, rc_new):
        if re.search(r"\bpoisson\b", options):
            # parse-time failures on both builds never reach estimation
            if rc_old != "0" and rc_new not in ("0", ""):
                return 0
            return 1
        return 0

    # ---- extgrid ----------------------------------------------------------
    for r in batteries:
        cid = r["case_id"]
        varlist, options, dataset, slow = battery_config(cid)
        # cross-check the manifest's own decomposition
        assert r["battery"] == cid[:2] and r["horizon"] == cid[4:6] \
            and r["reference"] == cid.split("_")[-1], f"manifest mismatch {cid}"
        assert r["expected_rc"] == "0"
        rows.append(dict(config=f"ext_{cid}", component="extgrid", dataset=dataset,
                         weight="", ifcond="", prep="", varlist=varlist,
                         options=options, rc_old="0", rc_new="0", msg_old="",
                         slow=slow))

    # ---- extlattice --------------------------------------------------------
    for r in lattices:
        cid = r["case_id"]
        varlist, options, dataset = lattice_config(r["lattice_token"], r["bins"])
        assert r["expected_rc"] == "0"
        rows.append(dict(config=f"ext_{cid}", component="extlattice", dataset=dataset,
                         weight="", ifcond="", prep="", varlist=varlist,
                         options=options, rc_old="0", rc_new="0", msg_old="",
                         slow=0))

    # ---- authored cases ----------------------------------------------------
    for r in cases:
        cid, family, rc_old = r["case_id"], r["family"], r["expected_rc"]
        if family in SKIP_FAMILIES:
            skipped.append((cid, SKIP_REASON[family]))
            continue
        component = "extval" if family == "validation" else "extcase"

        prog = case_program(cid)
        body = bodies.get(prog)
        if body is None:
            skipped.append((cid, f"case program {prog} not found"))
            continue
        block = find_case_block(body, cid)
        if block is None:
            skipped.append((cid, f"no static block for {cid} in {prog}"))
            continue
        try:
            calls = extract_calls(block, cid)
        except DynamicCase as e:
            skipped.append((cid, f"dynamic: {e}"))
            continue

        setup_name = SETUP_BY_CASE.get(cid, default_setup(cid))
        if setup_name not in setups:
            skipped.append((cid, f"setup {setup_name} not found"))
            continue
        dataset, setup_preps = setups[setup_name]
        if dataset is None:
            skipped.append((cid, f"in-memory setup {setup_name} (no fixture)"))
            continue

        if len(calls) > 1 and rc_old != "0":
            skipped.append((cid, "multi-call with nonzero rc -- ambiguous"))
            continue

        rc_new = RC_NEW_OVERRIDE.get(cid, "" if cid in RC_NEW_BLANK else rc_old)
        msg_old = messages.get(cid, "")

        for k, call in enumerate(calls, start=1):
            name = f"ext_{cid}" if len(calls) == 1 else f"ext_{cid}_c{k}"
            varlist, weight, ifcond, options = split_call(call["line"])
            prep_lines = list(setup_preps) + call["preps"]
            prep_key = ""
            if prep_lines:
                prep_key = name
                preps_out[prep_key] = prep_lines
            slow = poisson_slow(options, rc_old, rc_new)
            if cid == "cont_ext_cfx_subgroups":
                slow = 1   # subgroup+cfx-heavy
            rows.append(dict(config=name, component=component, dataset=dataset,
                             weight=weight, ifcond=ifcond, prep=prep_key,
                             varlist=varlist, options=options,
                             rc_old=rc_old, rc_new=rc_new, msg_old=msg_old,
                             slow=slow))

    # ---- checks ------------------------------------------------------------
    for cid, why in skipped:
        print(f"SKIP {cid}: {why}")
    n = {"extgrid": 0, "extlattice": 0, "extval": 0, "extcase": 0}
    for row in rows:
        n[row["component"]] += 1
        assert row["dataset"], f"missing dataset for {row['config']}"
    total_cases = len(batteries) + len(lattices) + len(cases)
    mined_cases = len(batteries) + len(lattices) + (len(cases) - len(skipped))
    print(f"bundle cases: {total_cases}  mined cases: {mined_cases}  "
          f"skipped: {len(skipped)}  configs: {len(rows)}")
    print(f"components: {n}")
    assert total_cases == 382, total_cases
    assert n["extgrid"] == 216 and n["extlattice"] == 24, n
    assert n["extval"] == 35, n
    assert n["extcase"] == 83, n
    assert len(rows) == 358, len(rows)

    # ---- emit ---------------------------------------------------------------
    cols = ["config", "component", "dataset", "weight", "ifcond", "prep",
            "varlist", "options", "rc_old", "rc_new", "msg_old", "slow"]
    with open(OUTDIR / "configs_ext.gen.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    with open(OUTDIR / "preps_ext.gen.json", "w", encoding="utf-8") as fh:
        json.dump(preps_out, fh, indent=1, sort_keys=True)
    with open(OUTDIR / "skipped_cases.txt", "w", encoding="utf-8") as fh:
        for cid, why in skipped:
            fh.write(f"{cid}\t{why}\n")
    print(f"wrote {OUTDIR/'configs_ext.gen.csv'} ({len(rows)} configs), "
          f"preps_ext.gen.json ({len(preps_out)} keys), skipped_cases.txt ({len(skipped)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
