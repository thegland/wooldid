# helper: load test/compare.py's parser without running its main()
import re
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
SRC=_REPO+'/test/compare.py'
src=open(SRC).read()
src=re.sub(r'\nmain\(\)\s*$', '\n', src)
src=re.sub(r'\nif __name__ == .__main__.:\s*\n\s*main\(\)\s*$', '\n', src)
ns={}
exec(compile(src, SRC, 'exec'), ns)
load=ns['load']; hex2f=ns['hex2f']
RO=_REPO+'/test/suite/out/run_old/'
RN=_REPO+'/test/suite/out/run_new/'
def both(cfg): return load(RO+cfg+'.csv'), load(RN+cfg+'.csv')
