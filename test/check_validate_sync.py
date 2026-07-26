#!/usr/bin/env python3
"""S5 guard: wooldid_worker and wooldid_validate re-parse the same `0', so their
syntax lines must be BYTE-IDENTICAL. A drift (an option added to one but not the
other) makes the validator reject every call using the new option - loud - or, if
a type/default drifts, silently changes parsing. Run from anywhere:

    python3 test/check_validate_sync.py
"""
import os
import sys

ado = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "wooldidinstall", "wooldid.ado")
lines = [l for l in open(ado, encoding="utf-8").read().split("\n")
         if l.startswith("syntax [varlist(numeric default=none)]")]
if len(lines) != 2:
    print(f"FAIL: expected exactly 2 wooldid syntax lines, found {len(lines)}")
    sys.exit(1)
if lines[0] != lines[1]:
    print("FAIL: the wooldid_worker and wooldid_validate syntax lines differ:")
    print("  worker:   " + lines[0][:120] + " ...")
    print("  validate: " + lines[1][:120] + " ...")
    sys.exit(1)
print("OK: worker and validate syntax lines are byte-identical")
