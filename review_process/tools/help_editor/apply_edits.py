#!/usr/bin/env python3
"""Apply the help editor's saved edits to wooldid.sthlp.

Usage: apply_edits.py <edits_dir> <sthlp> [--ver TAG] [--write]

<edits_dir> holds one JSON file per edited block, as ArtifactData's out_dir writes them
(collection edits, documents b<idx>): fields ver, idx, hash, smcl_original, smcl.
The current sthlp is re-parsed with parse_sthlp.parse; each edit must match a block by
hash (its position may have moved). Edited blocks are re-wrapped at WIDTH columns at
spaces outside braces. Without --write it only prints the diff.
"""
import difflib, glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_sthlp import parse

WIDTH = 170

def rewrap(smcl):
    if not smcl.endswith('\n'):
        smcl += '\n'
    body = smcl[:-1]
    if '\n' in body or len(body) <= WIDTH:
        return smcl
    out, line, depth, i = [], '', 0, 0
    last_space = -1
    while i < len(body):
        ch = body[i]
        if ch == '{': depth += 1
        elif ch == '}': depth -= 1
        line += ch
        if ch == ' ' and depth == 0:
            last_space = len(line) - 1
        if len(line) >= WIDTH and last_space > 0:
            out.append(line[:last_space].rstrip())
            line = line[last_space + 1:]
            last_space = -1
        i += 1
    if line:
        out.append(line)
    return '\n'.join(out) + '\n'

def main():
    edits_dir, sthlp = sys.argv[1], sys.argv[2]
    ver = sys.argv[sys.argv.index('--ver') + 1] if '--ver' in sys.argv else None
    write = '--write' in sys.argv
    text = open(sthlp, encoding='utf-8').read()
    blocks = parse(text)
    by_hash = {}
    for b in blocks:
        by_hash.setdefault(b['hash'], []).append(b)
    edits = []
    for f in sorted(glob.glob(os.path.join(edits_dir, '*.json'))):
        d = json.load(open(f, encoding='utf-8'))
        if ver and d.get('ver') != ver:
            print('SKIP (other version %s):' % d.get('ver'), f); continue
        edits.append(d)
    applied, problems = [], []
    new_blocks = [b['smcl'] for b in blocks]
    for d in edits:
        cands = by_hash.get(d['hash'], [])
        cands = [b for b in cands if b['smcl'] == d['smcl_original']]
        if len(cands) != 1:
            problems.append('block for edit b%s (line %s) not found once in the current file (%d matches)' % (d['idx'], d.get('line'), len(cands)))
            continue
        b = cands[0]
        if d['smcl'] == b['smcl']:
            continue
        new_blocks[b['idx']] = rewrap(d['smcl'])
        applied.append((b['idx'], b['line']))
    new_text = ''.join(new_blocks)
    diff = difflib.unified_diff(text.splitlines(True), new_text.splitlines(True), 'a/' + os.path.basename(sthlp), 'b/' + os.path.basename(sthlp), n=1)
    sys.stdout.write(''.join(diff))
    print('\napplied %d edit(s) at blocks %s' % (len(applied), applied))
    for p in problems:
        print('PROBLEM:', p)
    if write and not problems:
        open(sthlp, 'w', encoding='utf-8').write(new_text)
        print('written', sthlp)
    elif write:
        print('NOT written: fix the problems first')

if __name__ == '__main__':
    main()
