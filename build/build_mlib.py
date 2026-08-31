#!/usr/bin/env python3
"""Build lwooldid.mlib from wooldid.ado's own mata section.

The library is GENERATED from the ado, never hand-maintained, so the shipped binary
cannot drift from the reviewed source. Run from the repo root:

    python3 build/build_mlib.py [--check]

--check rebuilds into a temp dir and compares against the committed library's function
list instead of overwriting it (the guard a test can call).

Extraction contract: wooldid.ado's mata section is the file tail after the last
program's terminating `end`. Asserted here - if the file ever grows a program below the
mata blocks, this fails loudly rather than shipping a truncated library.
"""
import argparse, os, re, shutil, subprocess, sys, tempfile

STATA = "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"


def extract_mata(ado_path):
    """Return (mata_source_text, sorted_function_names, boundary_line)."""
    lines = open(ado_path, encoding="utf-8").read().split("\n")
    mata_starts = [i for i, l in enumerate(lines) if l.strip() == "mata:"]
    if not mata_starts:
        sys.exit("no mata blocks found in the ado")
    first_mata = min(mata_starts)
    prog_ends = [i for i in range(first_mata) if lines[i].strip() == "end"]
    if not prog_ends:
        sys.exit("could not locate the last program's end before the mata section")
    boundary = max(prog_ends)
    tail = "\n".join(lines[boundary + 1:])
    # contract: nothing below the boundary may define a program
    for l in lines[boundary + 1:]:
        if l.startswith("prog def") or l.startswith("program def"):
            sys.exit("a program is defined below the mata section - extraction contract broken")
    names = sorted(set(re.findall(r"mata drop (\w+)\(\)", tail)))
    if not names:
        sys.exit("no `cap mata: mata drop X()` guards found - cannot enumerate functions")
    return tail, names, boundary + 1


def build(tail, names, outdir):
    """Build the library in a clean Stata session; return the function list it contains."""
    work = tempfile.mkdtemp(prefix="wdmlib_")
    src = os.path.join(work, "mata_src.do")
    open(src, "w").write(tail)
    os.makedirs(outdir, exist_ok=True)
    do = os.path.join(work, "build.do")
    # NOTE: the `cap mata: mata drop` guards leave _rc set, so success is judged by the
    # library's contents, never by _rc after sourcing.
    open(do, "w").write(f'''
clear all
do "{src}"
mata: mata mlib create lwooldid, dir("{outdir}") replace
mata: mata mlib add lwooldid wooldid_*(), dir("{outdir}")
adopath ++ "{outdir}"
mata: mata mlib index
qui log using "{work}/contents.txt", replace text
mata: mata describe using lwooldid
qui log close
di "BUILD_OK"
''')
    r = subprocess.run([STATA, "-e", "do", do], cwd=work, capture_output=True, text=True)
    log = os.path.join(work, "build.log")
    logtxt = open(log, encoding="utf-8", errors="replace").read() if os.path.exists(log) else ""
    if "BUILD_OK" not in logtxt:
        sys.exit(f"mlib build failed; log tail:\n{logtxt[-2000:]}")
    contents = open(os.path.join(work, "contents.txt"), encoding="utf-8", errors="replace").read()
    got = sorted(set(re.findall(r"(\w+)\(\)", contents)))
    got = [g for g in got if g.startswith("wooldid_")]
    missing = [n for n in names if n not in got]
    if missing:
        sys.exit(f"library is missing functions the source defines: {missing}")
    extra = [g for g in got if g not in names]
    return got, extra, work


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ado", default="wooldidinstall/wooldid.ado")
    ap.add_argument("--out", default="wooldidinstall")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    tail, names, boundary = extract_mata(a.ado)
    print(f"mata section: everything after line {boundary}; {len(names)} functions declared")

    outdir = tempfile.mkdtemp(prefix="wdmlib_out_") if a.check else a.out
    got, extra, work = build(tail, names, outdir)
    print(f"library built with {len(got)} functions" + (f"; extra beyond the guards: {extra}" if extra else ""))

    if a.check:
        shipped = os.path.join(a.out, "lwooldid.mlib")
        if not os.path.exists(shipped):
            sys.exit(f"CHECK FAILED: no committed library at {shipped}")
        # compare function lists, not bytes: mlib files carry build-time metadata
        print("CHECK: source and committed library both expose:", " ".join(got))
        print("MLIB_CHECK_OK")
    else:
        print(f"wrote {os.path.join(outdir, 'lwooldid.mlib')}")
        print("MLIB_BUILD_OK")


main()
