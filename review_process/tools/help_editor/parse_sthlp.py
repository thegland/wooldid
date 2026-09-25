#!/usr/bin/env python3
"""Split a Stata .sthlp (SMCL) file into editable blocks that round-trip byte-exact.

Usage: parse_sthlp.py <in.sthlp> <out.json> [--ver TAG]

A block is one logical unit: a paragraph ({pstd} ... {p_end}, possibly over several
physical lines), a synopt table row, a title or dlgtab line, a one-line meta directive
({marker}, {viewerjumpto}, {synoptset}, ...), a blank line, or a raw line the parser does
not recognise. Concatenating every block's `smcl` reproduces the file exactly.
"""
import hashlib, json, re, sys

PARA_START = re.compile(r"^\{(pstd|phang2|phang|pmore2|pmore|p\s+\d|synopt\s*:|p2col\s)")
SINGLE = re.compile(r'^\{(title:|dlgtab:|marker\s|synoptset|synopthdr|p2colreset|hline|smcl\}|\*|vieweralsosee|viewerjumpto)')

def kind_of(line):
    if line.strip() == '':
        return 'blank'
    m = PARA_START.match(line)
    if m:
        k = m.group(1).strip()
        if k.startswith('synopt'): return 'synopt'
        if k.startswith('p2col'): return 'p2col'
        return {'pstd': 'pstd', 'phang': 'phang', 'phang2': 'phang2', 'pmore': 'pmore', 'pmore2': 'pmore2'}.get(k, 'p')
    m = SINGLE.match(line)
    if m:
        k = m.group(1)
        if k == 'title:': return 'title'
        if k == 'dlgtab:': return 'dlgtab'
        return 'meta'
    return 'raw'

def parse(text):
    lines = text.split('\n')
    # keep track of the final newline: split gives '' as the last element if the file ends with \n
    blocks = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if i == n - 1 and line == '' and text.endswith('\n'):
            break  # the trailing newline belongs to the previous line
        k = kind_of(line)
        if k in ('pstd', 'phang', 'phang2', 'pmore', 'pmore2', 'synopt', 'p2col', 'p'):
            j = i
            buf = [line]
            if '{p_end}' not in line:
                # continue until a line with {p_end}, stopping before a blank line or a new block start
                while j + 1 < n:
                    nxt = lines[j + 1]
                    nk = kind_of(nxt)
                    if nk == 'blank' or PARA_START.match(nxt) or SINGLE.match(nxt):
                        break
                    j += 1
                    buf.append(nxt)
                    if '{p_end}' in nxt:
                        break
            blocks.append((k, i, j, '\n'.join(buf)))
            i = j + 1
        else:
            blocks.append((k, i, i, line))
            i += 1
    out = []
    for idx, (k, a, b, s) in enumerate(blocks):
        smcl = s + '\n'
        out.append({
            'idx': idx, 'kind': k, 'line': a + 1, 'lines': b - a + 1,
            'hash': hashlib.sha1(smcl.encode('utf-8')).hexdigest()[:12],
            'smcl': smcl,
        })
    return out

def main():
    src, dst = sys.argv[1], sys.argv[2]
    ver = ''
    if '--ver' in sys.argv:
        ver = sys.argv[sys.argv.index('--ver') + 1]
    text = open(src, encoding='utf-8').read()
    blocks = parse(text)
    joined = ''.join(b['smcl'] for b in blocks)
    assert joined == text, 'round trip failed: %d vs %d bytes' % (len(joined), len(text))
    editable = sum(1 for b in blocks if b['kind'] not in ('blank', 'meta'))
    meta = {'ver': ver, 'file': src.split('/')[-1], 'nblocks': len(blocks), 'editable': editable,
            'sha1': hashlib.sha1(text.encode('utf-8')).hexdigest()}
    json.dump({'meta': meta, 'blocks': blocks}, open(dst, 'w'), ensure_ascii=False)
    from collections import Counter
    print('round trip OK;', len(blocks), 'blocks;', editable, 'editable;', dict(Counter(b['kind'] for b in blocks)))
    multi = [b for b in blocks if b['lines'] > 1]
    print('multi-line blocks:', len(multi), 'longest:', max(b['lines'] for b in blocks))
    noend = [b for b in blocks if b['kind'] not in ('blank', 'meta', 'title', 'dlgtab', 'raw') and '{p_end}' not in b['smcl']]
    print('paragraph blocks without {p_end}:', len(noend), [(b['line'], b['kind']) for b in noend][:12])
    raws = [b for b in blocks if b['kind'] == 'raw']
    print('raw blocks:', len(raws), [(b['line'], b['smcl'][:50]) for b in raws][:12])

if __name__ == '__main__':
    main()
