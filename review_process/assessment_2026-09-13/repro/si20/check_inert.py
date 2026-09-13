"""Inertness check: my instrumented-build dumps vs the suite's 2026-09-11 dumps, exact %21x string
equality on every e(b)/e(V)/e(wdidmainresults) cell and every scalar (rtol 0)."""
import sys
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
suite = _REPO + '/test/suite/out'
def load(path):
    d = {}
    for line in open(path, encoding='utf-8', errors='replace'):
        p = line.rstrip('\n').split('|')
        if p[0] == 'scalar': d[('scalar', p[1])] = p[2]
        elif p[0] == 'mat': d[('mat', p[1], int(p[2]), int(p[3]))] = p[4]
    return d
W = sys.argv[1]; build = sys.argv[2]
for cfg in sys.argv[3:]:
    a = load(f'{W}/run_{build}/{cfg}_{build}.csv'); b = load(f'{suite}/run_{build}/{cfg}.csv')
    keys = [k for k in a if k[0] == 'scalar' or (k[0] == 'mat' and k[1] in ('b', 'V', 'wdidmainresults'))]
    same = sum(1 for k in keys if b.get(k) == a[k]); missing = [k for k in keys if k not in b]
    diff = [k for k in keys if k in b and b[k] != a[k]]
    print(f'INERT|{build}|{cfg}|cells_checked={len(keys)}|bit_identical={same}|differ={len(diff)}|absent_in_suite={len(missing)}')
    for k in diff[:8]: print('   differs:', k, a[k], b[k])
