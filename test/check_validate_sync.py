#!/usr/bin/env python3
"""S5 guards for the wooldid_worker / wooldid_validate split. Run from anywhere:

    python3 test/check_validate_sync.py

Two checks, both static:

1. SYNTAX SYNC - the two programs re-parse the same `0', so their syntax lines
   must be BYTE-IDENTICAL. A drift (an option added to one but not the other)
   makes the validator reject every call using the new option - loud - or, if a
   type/default drifts, silently changes parsing.

2. C_LOCAL CLOSURE - every local ASSIGNED inside wooldid_validate must be either
   (a) in the blanket c_local return list, or (b) never read by the worker after
   the wooldid_validate call, or (c) assigned by the worker itself before its
   first read there. This is the failure class caught during implementation
   (suppressjointtests was mutated in the validator but not returned, silently
   enabling joint tests): a future edit that mutates one more local in the
   validator without adding it to the list would otherwise strand the worker on
   the stale syntax-parsed value with no error anywhere.
"""
import os
import re
import sys

ado = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "wooldidinstall", "wooldid.ado")
src = open(ado, encoding="utf-8").read().split("\n")

fail = False

# ---- check 1: syntax lines byte-identical -----------------------------------
syn = [l for l in src if l.startswith("syntax [varlist(numeric default=none)]")]
if len(syn) != 2:
    print(f"FAIL: expected exactly 2 wooldid syntax lines, found {len(syn)}")
    fail = True
elif syn[0] != syn[1]:
    print("FAIL: the wooldid_worker and wooldid_validate syntax lines differ:")
    print("  a: " + syn[0][:120] + " ...")
    print("  b: " + syn[1][:120] + " ...")
    fail = True
else:
    print("OK: worker and validate syntax lines are byte-identical")

# ---- check 2: c_local closure ------------------------------------------------
# TWO PRECONDITIONS OF THIS CHECK, both properties of today's ado rather than of the
# checker, so they need re-verifying whenever the validator changes:
#  (1) span() ends a program at the FIRST line whose strip() is "end", so a nested end - a
#      mata: block, an inner program define - would silently truncate the body and shrink
#      the set of locals checked. (It resolves correctly today, including the worker's
#      INDENTED "  end", which only the .strip() makes findable.)
#  (2) the `assigned' scan below recognizes ONLY the `local' keyword, so a local created by
#      gettoken, by Mata's st_local(), or by c_local is invisible to the closure check. The
#      current wooldid_validate body contains none of those three; add one and this script
#      still prints "OK: c_local closure holds" while checking nothing about it.
def span(name):
    s = e = None
    for i, l in enumerate(src):
        if l.startswith(f"prog def {name}"):
            s = i
        elif s is not None and l.strip() == "end" and i > s:
            e = i
            break
    return s, e

ws, we = span("wooldid_worker")
vs, ve = span("wooldid_validate")
if None in (ws, we, vs, ve):
    print("FAIL: could not locate wooldid_worker / wooldid_validate program spans")
    sys.exit(1)
vbody = src[vs:ve]

clist = None
for l in vbody:
    if l.strip().startswith("foreach wdv_ret in"):
        clist = set(l.split(" in ", 1)[1].split(" {")[0].split())
if not clist:
    print("FAIL: could not find the c_local return list (foreach wdv_ret in ...)")
    sys.exit(1)

assigned = set()
for l in vbody:
    # mid-line mutations ("if COND local x ..." / "local ++x") are this file's dominant
    # idiom - a line-anchored match missed them (verifier s5_v1, 2026-07-26). skip
    # comment lines so prose after "* local ..." is not picked up.
    if l.lstrip().startswith("*"):
        continue
    m = re.search(r"(?:^|\s)(?:cap\s+)?local\s+(?:\+\+)?([A-Za-z_][A-Za-z0-9_]*)", l)
    if m:
        assigned.add(m.group(1))
assigned.discard("wdv_ret")

callidx = None
for i in range(ws, we):
    if "wooldid_validate `0'" in src[i]:
        callidx = i
        break
if callidx is None:
    print("FAIL: worker does not call wooldid_validate `0'")
    sys.exit(1)

worker_after = src[callidx + 1:we]
missing = []
for name in sorted(assigned - clist):
    first_read = first_assign = None
    pat_read = re.compile(r"`" + re.escape(name) + r"'")
    # THE NAME MUST BE THE THING BEING ASSIGNED, not merely present on the line (author
    # ruling 2026-09-11, documentation review D-14). This used to be
    # `^...(?:local|tempvar|tempname)\s.*\bNAME\b', which counted a line that only READS
    # the local as assigning it - so an ordinary `local N = `nsg'' registered as an
    # assignment of nsg and suppressed the finding, in the script whose whole job is to
    # catch that class. Two forms are recognized now: `local NAME ...' (the macro is the
    # first token, optionally ++-prefixed) and `tempvar/tempname ... NAME ...' (any token
    # in the declaration list, but only a BARE one - a `NAME' dereference does not count,
    # which is what used to suppress `item' via `tempvar `item''.)
    # foreach/forvalues are recognized too, and had to be: tightening the regex immediately
    # surfaced `item', which the worker assigns fourteen times as a `foreach item in ...'
    # loop variable and the old pattern had been suppressing for the wrong reason. A loop
    # variable is a real assignment of the macro, so it belongs here; anything the worker
    # sets by some OTHER idiom (gettoken, mata st_local, c_local) is still invisible to this
    # check and would show up as a false positive here rather than passing silently.
    pat_assign = re.compile(
        r"^\s*(?:qui(?:etly)?\s+|cap(?:ture)?\s+|noi(?:sily)?\s+)*"
        r"(?:local\s+(?:\+\+)?" + re.escape(name) + r"(?![A-Za-z0-9_`])"
        r"|(?:tempvar|tempname)\s+(?:[A-Za-z_][A-Za-z0-9_]*\s+)*"
        + re.escape(name) + r"(?![A-Za-z0-9_`])"
        r"|(?:foreach|forv(?:al(?:ues)?)?)\s+" + re.escape(name) + r"(?![A-Za-z0-9_`]))"
    )
    for j, l in enumerate(worker_after):
        if first_read is None and pat_read.search(l):
            first_read = j
        if first_assign is None and pat_assign.match(l):
            first_assign = j
        if first_read is not None and first_assign is not None:
            break
    if first_read is not None and (first_assign is None or first_assign > first_read):
        missing.append(name)

# typo check: every c_local name must trace to the syntax line or an assignment
synline = syn[0].lower() if syn else ""
typos = [n for n in sorted(clist) if n not in assigned and n.lower() not in synline]

# BLANKET COVERAGE: every option the syntax line declares must appear in the c_local
# list. The blanket's own comment claims this is true "by construction", but nothing
# enforced it, and a new option that is merely LEFT OUT is invisible to the closure
# check above (the worker parses the syntax itself, so it reads a correct value - only
# the SHELL, which relies on the c_local, silently sees an empty one). That is how
# fastcollapse initially failed to reach the permute/bootstrap replication option
# strings while every direct run behaved perfectly.
opt_names = []
for tok in re.findall(r"[A-Za-z_]+", syn[0].split("]", 1)[-1] if syn else ""):
    opt_names.append(tok)
# recover each option's full name from the syntax line: NAMEabbrev(...) or NAMEabbrev
declared = set()
for m in re.finditer(r"(?<![A-Za-z_])([A-Za-z][A-Za-z_0-9]*)(?=\s*\(|\s|$)", syn[0] if syn else ""):
    declared.add(m.group(1).lower())
# syntax-line noise that is not an option
NOISE = {"syntax", "varlist", "numeric", "default", "none", "if", "aw", "fw", "pw", "string",
         "varname", "integer", "int", "numlist", "real", "name", "namelist", "passthru"}
uncovered = []
for opt in sorted(declared - NOISE):
    if opt in ("y", "i", "t"):
        continue
    if opt in {c.lower() for c in clist}:
        continue
    # options the blanket documents as deliberately excluded
    if opt in {"varlist", "weight", "exp"}:
        continue
    uncovered.append(opt)

if missing:
    print("FAIL: locals assigned in wooldid_validate, read by the worker after the "
          "call, but not in the c_local list (the suppressjointtests defect class): "
          + " ".join(missing))
    fail = True
else:
    print(f"OK: c_local closure holds ({len(assigned)} assigned, {len(clist)} returned)")
if typos:
    print("FAIL: c_local names that are neither syntax options nor assigned in the "
          "validator (typo?): " + " ".join(typos))
    fail = True
if uncovered:
    print("FAIL: syntax options not in the c_local blanket (the shell reads these as "
          "empty, so rep option strings and shell-side gates silently lose them): "
          + " ".join(uncovered))
    fail = True
else:
    print("OK: every syntax option is covered by the c_local blanket")

# ---- the permutation draw-matrix cap, written three times (added 2026-09-22) ----
# When the author raised the permexact enumeration ceiling to the Stata matrix row limit, an
# exact run between 10,001 and 65,534 assignments became legal for the first time, and such a
# run silently loses its draw matrices to a SEPARATE cap of 10,000 completed replications. That
# cap now appears three times: the shell's e(permdist) test, wooldid_permwire's `done' test, and
# the caution the shell prints when it bites. They share no storage. If they drift, the run
# either announces a loss that did not happen or loses the matrices in silence again - which is
# exactly the defect the caution was added to close. This is not testable at run time without
# actually completing 10,001 replications, so it is pinned here instead.
capvals = {}
joined = "\n".join(src)
# two of the three are unique one-liners
for label, pat in [
    ("shell e(permdist) post", r"if\s+`wdp_r_done'\s*<=\s*(\d+)\s+cap ereturn matrix permdist"),
    ("shell draw-matrix caution", r"if\s+`wdp_r_done'\s*>\s*(\d+)\s+dis as text"),
]:
    hits = re.findall(pat, joined)
    if len(hits) != 1:
        print(f"FAIL: expected exactly one '{label}' site for the permutation draw-matrix cap, "
              f"found {len(hits)}")
        fail = True
    else:
        capvals[label] = hits[0]
# the third shares its shape with wooldid_bootwire's e(bootdist) twin, so it is identified by
# the matrices posted inside the block rather than by the test itself
permwire_hits = []
for idx, line in enumerate(src):
    m = re.search(r"if\s+`done'\s*<=\s*(\d+)\s*\{", line)
    if m and any("permdist" in l for l in src[idx:idx + 14]):
        permwire_hits.append(m.group(1))
if len(permwire_hits) != 1:
    print(f"FAIL: expected exactly one wooldid_permwire draw-matrix cap site, found "
          f"{len(permwire_hits)}")
    fail = True
else:
    capvals["permwire draw-matrix post"] = permwire_hits[0]
if len(capvals) == 3:
    if len(set(capvals.values())) != 1:
        print("FAIL: the permutation draw-matrix cap disagrees across its three sites: "
              + ", ".join(f"{k}={v}" for k, v in capvals.items()))
        fail = True
    else:
        print(f"OK: the permutation draw-matrix cap is {next(iter(capvals.values()))} at all three sites")

# ---- the bootstrap draw-matrix cap and its caution ----
# wooldid_bootwire posts e(bootdist)/e(bootdist_se) only up to the same kind of cap, and the
# shell cautions when the matrices are withheld, keyed on the same completed-replication count
# (bootwire's `done' is the shell's `wdb_r_done'). As in the permutation block the cap is then
# written three times - bootwire's test, the shell caution's test and the number in its text -
# and they share no storage, so they can drift. Pinned here for the permutation twin's reason:
# exercising it needs 10,001 completed replications. The bootwire site is the `done' test with
# bootdist (not permdist) posted inside it. Also pinned: the caution prints once, from the
# shell's bootstrap block below its part(boot) display - below the result tables - and before
# that block clears the replay marker.
bootwire_sites = []
for idx, line in enumerate(src):
    m = re.search(r"if\s+`done'\s*<=\s*(\d+)\s*\{", line)
    if m and any("ereturn matrix bootdist" in l for l in src[idx:idx + 20]) \
            and not any("permdist" in l for l in src[idx:idx + 20]):
        bootwire_sites.append((idx, m.group(1)))
bshell_sites = []
for idx, line in enumerate(src):
    m = re.match(r"\s*if\s+`wdb_r_done'\s*>\s*(\d+)\s+dis as text \"Caution: `wdb_r_done' completed "
                 r"replications exceeds the ([\d,]+) cap on posting draw matrices, so e\(bootdist\)",
                 line)
    if m:
        bshell_sites.append((idx, m.group(1), m.group(2)))
bcaut_n = sum(1 for l in src if "so e(bootdist) and e(bootdist_se) are not posted" in l)
if len(bootwire_sites) != 1:
    print(f"FAIL: expected exactly one wooldid_bootwire draw-matrix cap site, found "
          f"{len(bootwire_sites)}")
    fail = True
elif len(bshell_sites) != 1 or bcaut_n != 1:
    print(f"FAIL: expected the bootstrap draw-matrix caution exactly once, in the shell keyed on "
          f"`wdb_r_done'; found {len(bshell_sites)} shell site(s) and {bcaut_n} line(s) with its text")
    fail = True
else:
    widx, bcap = bootwire_sites[0]
    sidx, scap, tcap = bshell_sites[0]
    # the nearest display call above the caution must be part(boot), with no marker clear between
    disp = None
    for j in range(sidx - 1, -1, -1):
        if re.match(r"\s*wooldid_display\s*,\s*part\(", src[j]):
            disp = j
            break
    cleared = disp is not None and any(re.match(r"\s*global WOOLDID_LASTRC\s*$", src[j])
                                       for j in range(disp, sidx))
    if disp is None or not re.match(r"\s*wooldid_display\s*,\s*part\(boot\)", src[disp]) or cleared:
        print("FAIL: the bootstrap draw-matrix caution does not print below the part(boot) display "
              "inside the shell's bootstrap block")
        fail = True
    elif len({bcap, scap, tcap.replace(",", "")}) != 1:
        print(f"FAIL: the bootstrap draw-matrix cap disagrees across its three sites: bootwire "
              f"test {bcap}, shell caution test {scap}, caution text {tcap}")
        fail = True
    else:
        print(f"OK: the bootstrap draw-matrix cap is {bcap} in the e(bootdist) test, the shell "
              f"caution's test and its text, and the caution prints below part(boot)")

sys.exit(1 if fail else 0)
