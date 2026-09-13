#!/usr/bin/env python3
"""Demote the analyzer's own heading levels so its tables nest under section 5 of the
record instead of clashing with the prose's section numbering."""
import re, sys
MAP = {"0": "5.1", "1": "5.2", "2": "5.3", "3": "5.4", "4": "5.5", "5": "5.6", "6": "5.7"}
for line in sys.stdin:
    if line.startswith("### "):
        line = "#" + line
    elif line.startswith("## "):
        m = re.match(r"## (\d+)\.\s*(.*)$", line)
        if m and m.group(1) in MAP:
            line = f"### {MAP[m.group(1)]} {m.group(2)}\n"
        else:
            line = "#" + line
    sys.stdout.write(line)
