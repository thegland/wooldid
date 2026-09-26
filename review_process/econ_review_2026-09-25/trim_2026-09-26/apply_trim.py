#!/usr/bin/env python3
"""Apply the verified trim proposals to wooldid.sthlp and run the gates.
Usage: apply_trim.py <repo> <tmpdir> <out.sthlp> [--groups G1,G2,...]
Reads <tmpdir>/blocks_f556138.json, proposals_<G>.json and verdicts_<G>.json.
accept -> proposal.new_smcl (or delete); fixed_replace -> verdict.new_smcl; reject -> unchanged.
"""
import json, os, re, sys
repo, T, out = sys.argv[1], sys.argv[2], sys.argv[3]
groups = sys.argv[sys.argv.index('--groups')+1].split(',') if '--groups' in sys.argv else ['G1','G2','G3','G4']
sys.path.insert(0, os.path.join(repo, 'review_process/tools/help_editor'))
from parse_sthlp import parse
from apply_edits import rewrap
src = os.path.join(repo, 'wooldidinstall/wooldid.sthlp')
text = open(src, encoding='utf-8').read()
blocks = parse(text)
ref = json.load(open(os.path.join(T, 'blocks_f556138.json')))
assert [b['hash'] for b in blocks] == [b['hash'] for b in ref], 'the help changed since the blocks were cut'
byidx = {b['idx']: b for b in blocks}
new = {b['idx']: b['smcl'] for b in blocks}
deleted = set()
stats = {'accept': 0, 'fixed_replace': 0, 'reject': 0, 'skipped': 0}
problems = []
for g in groups:
    pp = os.path.join(T, f'proposals_{g}.json'); vp = os.path.join(T, f'verdicts_{g}.json')
    if not (os.path.exists(pp) and os.path.exists(vp)):
        problems.append(f'{g}: missing proposals or verdicts'); continue
    props = {p['idx']: p for p in json.load(open(pp))['proposals']}
    verds = json.load(open(vp))['verdicts']
    seen = set()
    for v in verds:
        idx = v['idx']; seen.add(idx)
        p = props.get(idx)
        if p is None:
            problems.append(f'{g}: verdict for idx {idx} has no proposal'); continue
        b = byidx.get(idx)
        if b is None or b['hash'] != p['hash']:
            problems.append(f'{g}: hash mismatch at idx {idx}'); stats['skipped'] += 1; continue
        if b['kind'] in ('title', 'dlgtab', 'meta', 'raw', 'blank'):
            problems.append(f'{g}: idx {idx} is kind {b["kind"]}, not editable'); stats['skipped'] += 1; continue
        verdict = v['verdict']
        if verdict == 'reject':
            stats['reject'] += 1; continue
        if verdict == 'accept':
            smcl = '' if p['action'] == 'delete' else p['new_smcl']
        elif verdict == 'fixed_replace':
            smcl = v.get('new_smcl') or ''
            if not smcl:
                problems.append(f'{g}: fixed_replace at idx {idx} without new_smcl'); stats['skipped'] += 1; continue
        else:
            problems.append(f'{g}: unknown verdict {verdict} at idx {idx}'); continue
        if smcl == '':
            if b['kind'] in ('synopt', 'p2col'):
                problems.append(f'{g}: delete of a {b["kind"]} row at idx {idx} refused'); stats['skipped'] += 1; continue
            deleted.add(idx); new[idx] = ''
        else:
            if not smcl.endswith('\n'): smcl += '\n'
            # directive and p_end checks
            if not re.match(r'^\{(pstd|phang2|phang|pmore2|pmore|p\s+\d|synopt\s*:|p2col\s)', smcl):
                problems.append(f'{g}: idx {idx} replacement lacks its opening directive'); stats['skipped'] += 1; continue
            if '{p_end}' in b['smcl'] and '{p_end}' not in smcl:
                problems.append(f'{g}: idx {idx} replacement lost {{p_end}}'); stats['skipped'] += 1; continue
            if smcl.count('{') != smcl.count('}'):
                problems.append(f'{g}: idx {idx} replacement has unbalanced braces'); stats['skipped'] += 1; continue
            new[idx] = rewrap(smcl)
        stats[verdict] += 1
    for idx in props:
        if idx not in seen: problems.append(f'{g}: proposal idx {idx} has no verdict (left unchanged)')
# hand fixes by the coordinator, applied last (keyed by original block idx; text is the full final block)
mf = os.path.join(T, 'manual_fixes.json')
if os.path.exists(mf):
    for k, smcl in json.load(open(mf)).items():
        i = int(k)
        if smcl.count('{') != smcl.count('}') or '{p_end}' not in smcl:
            problems.append(f'manual fix at idx {i} fails the SMCL checks'); continue
        new[i] = rewrap(smcl if smcl.endswith('\n') else smcl + '\n'); stats['manual'] = stats.get('manual', 0) + 1
# collapse the blank a deleted paragraph leaves behind (blank, deleted, blank -> one blank)
for idx in sorted(deleted):
    prev, nxt = idx - 1, idx + 1
    if prev in new and nxt in new and byidx[prev]['kind'] == 'blank' and byidx[nxt]['kind'] == 'blank' and new[prev] == '\n' and new[nxt] == '\n':
        new[nxt] = ''
newtext = ''.join(new[i] for i in sorted(new))
open(out, 'w', encoding='utf-8').write(newtext)
print('stats', stats, 'deleted blocks', len(deleted))
for p in problems: print('PROBLEM:', p)
# gates
def inv(t):
    return dict(opt=set(re.findall(r'\{opt\s+([A-Za-z_0-9]+)', t)), markers=set(re.findall(r'\{marker\s+(\w+)\}', t)),
                links=set(re.findall(r'\{help wooldid##(\w+)', t)), eret=set(re.findall(r'e\(([A-Za-z_0-9*]+)\)', t)),
                opth=set(re.findall(r'\{opth?\s+([A-Za-z_0-9]+)', t)), cmd=set(re.findall(r'\{cmd:([^}]+)\}', t)))
a, b2 = inv(text), inv(newtext)
ok = True
for k in a:
    lost = a[k] - b2[k]; added = b2[k] - a[k]
    if lost or added:
        ok = False; print(f'COVERAGE {k}: lost {sorted(lost)} added {sorted(added)}')
print('coverage gate', 'OK' if ok else 'FAILED')
nb = parse(newtext)
para = [x for x in nb if x['kind'] in ('pstd','phang','phang2','pmore','pmore2','synopt','p2col')]
print('paragraph blocks without p_end:', [(x['line'], x['kind']) for x in para if '{p_end}' not in x['smcl']])
print('brace balance:', newtext.count('{') == newtext.count('}'), 'dangling links:', sorted(b2['links'] - b2['markers']))
def metrics(t):
    body = re.sub(r'\{[^}]*\}', '', t); w = body.split(); n = len(w)
    return dict(words=n, wc=len(t.split()), dash=round(t.count(' - ')*1000/n, 2), paren=round(t.count('(')*1000/n, 1), semis=t.count(';'), ul=t.count('{ul:'), emdash=t.count('—'))
print('before', metrics(text)); print('after ', metrics(newtext))
# per-section words
def secwords(t):
    sec = '(pre)'; c = {}
    for l in t.split('\n'):
        m = re.match(r'\{marker (\w+)\}', l)
        if m: sec = m.group(1)
        c[sec] = c.get(sec, 0) + len(re.sub(r'\{[^}]*\}', '', l).split())
    return c
sa, sb = secwords(text), secwords(newtext)
for s in sa: 
    if sa[s] != sb.get(s, 0): print(f'  {s:18s} {sa[s]:6d} -> {sb.get(s,0):6d}  ({sb.get(s,0)-sa[s]:+d})')
