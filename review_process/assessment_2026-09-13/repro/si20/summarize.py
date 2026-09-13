"""Per-config 'This run' summary from the suite's 2026-09-11 dumps (run_old = v1.3 301b706, run_new = 14c30d3)."""
import sys, math
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
base = _REPO + '/test/suite/out'
def hex2f(s):
    s = s.strip()
    if s.startswith('.') and len(s) <= 2: return None
    try:
        t = s.replace(' ', ''); neg = t.startswith('-')
        if t[0] in '+-': t = t[1:]
        mant, exp = t.split('X'); v = float.fromhex('0x' + mant + 'p' + str(int(exp, 16)))
        return -v if neg else v
    except Exception: return None
def load(path):
    d = {'scalar': {}, 'b': {}, 'bn': [], 'rc': None, 'props': ''}
    for line in open(path, encoding='utf-8', errors='replace'):
        p = line.rstrip('\n').split('|')
        if p[0] == 'rc': d['rc'] = p[1]
        elif p[0] == 'scalar': d['scalar'][p[1]] = hex2f(p[2])
        elif p[0] == 'macro' and p[1] == 'properties': d['props'] = p[2]
        elif p[0] == 'matnames' and p[1] == 'b': d['bn'] = p[4].split()
        elif p[0] == 'mat' and p[1] == 'b': d['b'][int(p[3])] = hex2f(p[4])
    d['bb'] = {d['bn'][c-1]: v for c, v in d['b'].items() if c-1 < len(d['bn'])}
    return d
def rd(a, b):
    if a is None or b is None: return float('nan')
    return abs(a-b)/max(abs(a), abs(b), 1e-300)
mode = sys.argv[1]
print(f'{"config":30s} {"df_r o->n":>11s} {"main_ame old->new":>36s} {"ame shift(old-new) min..max":>30s} {"worst att b rel":>15s} {"joint_horiz old->new":>28s} {"joint_Choriz old->new":>28s}')
for cfg in sys.argv[2:]:
    try:
        o, n = load(f'{base}/run_old/{cfg}.csv'), load(f'{base}/run_new/{cfg}.csv')
    except FileNotFoundError as e:
        print(cfg, 'MISSING', e); continue
    ame = [k for k in o['bb'] if 'ame' in k and k in n['bb']]
    att = [k for k in o['bb'] if 'att' in k and k in n['bb']]
    shifts = [o['bb'][k]-n['bb'][k] for k in ame if o['bb'][k] is not None and n['bb'][k] is not None]
    watt = max([rd(o['bb'][k], n['bb'][k]) for k in att], default=float('nan'))
    ma = f"{o['bb'].get('main_ame')!s:>17.17} -> {n['bb'].get('main_ame')!s:<17.17}" if 'main_ame' in o['bb'] else '(no ame rows)'
    sh = f"{min(shifts):.6f}..{max(shifts):.6f}" if shifts else '-'
    jh = f"{o['scalar'].get('joint_horiz')!s:>12.12} -> {n['scalar'].get('joint_horiz')!s:<12.12}" if 'joint_horiz' in o['scalar'] else '-'
    jc = f"{o['scalar'].get('joint_Choriz')!s:>12.12} -> {n['scalar'].get('joint_Choriz')!s:<12.12}" if 'joint_Choriz' in o['scalar'] else '-'
    dfr = f"{o['scalar'].get('df_r')}->{n['scalar'].get('df_r')}"
    print(f'{cfg:30s} {dfr:>11s} {ma:>36s} {sh:>30s} {watt:15.2e} {jh:>28s} {jc:>28s}   props {o["props"]}->{n["props"]}')
