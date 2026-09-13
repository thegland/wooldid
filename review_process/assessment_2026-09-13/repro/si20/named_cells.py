import sys
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
base = _REPO + '/test/suite/out'
def hex2f(s):
    s = s.strip()
    if s.startswith('.') and len(s) <= 2: return s
    try:
        t = s.replace(' ', ''); neg = t.startswith('-')
        if t[0] in '+-': t = t[1:]
        mant, exp = t.split('X'); v = float.fromhex('0x' + mant + 'p' + str(int(exp, 16)))
        return -v if neg else v
    except Exception: return s
def load(path):
    d = {'scalar': {}, 'mat': {}, 'names': {}, 'rc': None, 'macro': {}}
    for line in open(path, encoding='utf-8', errors='replace'):
        p = line.rstrip('\n').split('|')
        if p[0] == 'rc': d['rc'] = p[1]
        elif p[0] == 'scalar': d['scalar'][p[1]] = hex2f(p[2])
        elif p[0] == 'macro': d['macro'][p[1]] = p[2] if len(p) > 2 else ''
        elif p[0] == 'mat': d['mat'][(p[1], int(p[2]), int(p[3]))] = hex2f(p[4])
        elif p[0] == 'matnames': d['names'][p[1]] = (p[2], p[3].split(), p[4].split())
    return d
scal = ['N','df_r','df_a','F','joint_horiz','joint_Choriz','ll','rss','r2','r2_within','r2_withinadj','r2adj','ll_0','rmse']
for cfg in sys.argv[1:]:
    do, dn = load(f'{base}/run_old/{cfg}.csv'), load(f'{base}/run_new/{cfg}.csv')
    print(f'##### {cfg}   rc old={do["rc"]} new={dn["rc"]}  props old=[{do["macro"].get("properties")}] new=[{dn["macro"].get("properties")}]')
    for sc in scal:
        if sc in do['scalar'] or sc in dn['scalar']:
            vo, vn = do['scalar'].get(sc), dn['scalar'].get(sc)
            flag = ''
            if isinstance(vo, float) and isinstance(vn, float) and abs(vo-vn) > 1e-10*max(abs(vo),abs(vn)): flag = ' <<<'
            print(f'  scalar {sc:13s} old={vo!r:>22} new={vn!r:>22}{flag}')
    for mat in ['b', 'wdidmainresults']:
        if mat not in do['names']: continue
        dim, rn, cn = do['names'][mat]
        print(f'  matrix {mat} {dim}')
        keys = sorted(k for k in do['mat'] if k[0] == mat)
        for (m, r, c) in keys:
            vo, vn = do['mat'][(m,r,c)], dn['mat'].get((m,r,c))
            rname = cn[c-1] if mat == 'b' else rn[r-1]
            cname = '' if mat == 'b' else cn[c-1]
            if isinstance(vo, float) and isinstance(vn, float):
                rd = abs(vo-vn)/max(abs(vo),abs(vn),1e-300); flag = ' <<<' if rd > 1e-10 else ''; rds = f'{rd:.2e}'
            else:
                rds = 'NA'; flag = ' <<<' if vo != vn else ''
            if mat == 'b' or flag:
                print(f'    [{r},{c}] {rname:18s} {cname:8s} old={vo!r:>24} new={vn!r:>24} rel={rds}{flag}')
