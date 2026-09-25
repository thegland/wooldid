#!/usr/bin/env python3
"""Build the help editor page: parse the sthlp into blocks and inject the seed.
Usage: build.py <sthlp> <ver-tag> <out.html>"""
import json, subprocess, sys, os
here = os.path.dirname(os.path.abspath(__file__))
sthlp, ver, out = sys.argv[1], sys.argv[2], sys.argv[3]
blocks_path = os.path.join(here, 'blocks.json')
subprocess.run([sys.executable, os.path.join(here, 'parse_sthlp.py'), sthlp, blocks_path, '--ver', ver], check=True)
seed = open(blocks_path, encoding='utf-8').read()
seed = seed.replace('</', '<\\/')
tpl = open(os.path.join(here, 'help_editor.template.html'), encoding='utf-8').read()
assert tpl.count('__SEED__') == 1
html = tpl.replace('__SEED__', seed)
open(out, 'w', encoding='utf-8').write(html)
print('wrote', out, len(html), 'bytes')
