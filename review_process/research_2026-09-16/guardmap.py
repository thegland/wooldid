import re, sys
src = open('wooldidinstall/wooldid.ado', encoding='utf-8', errors='replace').read().split('\n')

EXP = re.compile(r'\b(gegen|gstats|gduplicates|gunique|glevelsof|gcollapse|gsort|gisid|greshape)\b')
prog = None
stack = []          # list of (depth, condition-text, lineno)
depth = 0
rows = []
progline = {}

for n, raw in enumerate(src, 1):
    line = raw.strip()
    m = re.match(r'^(?:cap\w*\s+)?prog(?:ram)?\s+def(?:ine)?\s+(\w+)', line)
    if m:
        prog = m.group(1); depth = 0; stack = []; progline[prog] = n
        continue
    if re.match(r'^end\b', line):
        prog = None; continue
    if line.startswith('*') or not line:
        continue
    code = line.split('//')[0]
    opens = code.count('{'); closes = code.count('}')
    # record the guard BEFORE adjusting depth for this line's own braces
    if EXP.search(code) and prog:
        conds = [c for (_d, c, _l) in stack]
        rows.append((prog, n, conds[-3:], code[:90]))
    if opens:
        cm = re.match(r'^(?:qui\w*\s+)?(if|else if|else|foreach|forv\w*|while)\b(.*)$', code)
        cond = (cm.group(1) + ' ' + cm.group(2).replace('{','').strip())[:78] if cm else code.replace('{','').strip()[:78]
        for _ in range(opens):
            stack.append((depth, cond, n)); depth += 1
    for _ in range(closes):
        if stack: stack.pop()
        depth = max(0, depth-1)

import collections
bycond = collections.Counter()
for prog_, n, conds, code in rows:
    inner = conds[-1] if conds else '<TOP OF PROGRAM - runs on every call into it>'
    bycond[(prog_, inner)] += 1

print("=== expensive-op sites by (program, innermost enclosing condition) ===")
for (p, c), k in sorted(bycond.items(), key=lambda x: -x[1]):
    if k >= 1:
        print("%3d  %-22s %s" % (k, p, c))

print()
print("=== ALWAYS-ON-PATH sites, with code ===")
TARGETS = {('wooldid_worker','<TOP OF PROGRAM - runs on every call into it>'),
           ('wooldid_worker','if "`safetyinput\'" != "off"'),
           ('wooldid_worker','if "`safety\'" != "off"'),
           ('wooldid_reltime','<TOP OF PROGRAM - runs on every call into it>'),
           ('wooldid_fcol','<TOP OF PROGRAM - runs on every call into it>'),
           ('wooldid_validate','<TOP OF PROGRAM - runs on every call into it>')}
for prog_, n, conds, code in rows:
    inner = conds[-1] if conds else '<TOP OF PROGRAM - runs on every call into it>'
    if (prog_, inner) in TARGETS:
        print("%-18s L%-6d [%s]" % (prog_, n, inner[:40]))
        print("        %s" % code)

print()
print("=== sites in the replication engines, with FULL enclosing stack ===")
for prog_, n, conds, code in rows:
    if prog_ in ('wooldid_permrun','wooldid_bootrun','wooldid_permctl','wooldid_bootctl'):
        print("%-16s L%-6d %s" % (prog_, n, ' >> '.join(c[:46] for c in conds) or '<TOP>'))
        print("      %s" % code[:86])
