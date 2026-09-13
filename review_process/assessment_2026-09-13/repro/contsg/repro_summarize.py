# Summarize repro_old.log / repro_new.log: re-join Stata's 80-column wrapped lines ("> " continuations),
# print every tagged value side by side, list the decisive cautions, and run test/compare.py for
# bit-identity of the reproduction's dumps against the suite's 2026-09-11 dumps.
import re, subprocess, os
import os, sys
_REPO = os.environ.get("WOOLDID_REPO", "")
_W    = os.environ.get("WOOLDID_SCRATCH", "")
if not _REPO: raise SystemExit("set WOOLDID_REPO=<path>  (see repro/README.md)")
if not _W: raise SystemExit("set WOOLDID_SCRATCH=<path>  (see repro/README.md)")
W=_W
R=_REPO
def joined(path):
    out=[]
    for line in open(path,encoding='utf-8',errors='replace'):
        line=line.rstrip('\n')
        if line.startswith('> ') and out: out[-1]+=line[2:]
        else: out.append(line)
    return out
TAGRE=re.compile(r'^([A-Z0-9]+_[A-Za-z0-9_]+)=(.*)$')
def tags(lines):
    d={}
    for l in lines:
        m=TAGRE.match(l)
        if m and not l.startswith('. '): d[m.group(1)]=m.group(2)
    return d
old=joined(f'{W}/repro_old.log'); new=joined(f'{W}/repro_new.log')
to=tags(old); tn=tags(new)
print('=== BUILD FINGERPRINTS ===')
for name,L in (('old',old),('new',new)):
    for l in L:
        if l.startswith('BUILD_DIR=') or ('wooldid.ado' in l and l.strip().startswith('/')): print(f'  {name}: {l.strip()}')
print('\n=== TAGGED VALUES (old | new) ===')
keys=sorted(set(to)|set(tn), key=lambda k:(k.split('_')[0],k))
for k in keys:
    print(f'{k}\n    OLD: {to.get(k,"-")}\n    NEW: {tn.get(k,"-")}')
print('\n=== SENTINELS ===')
for s in ['REPRO_AMESE_R1_DONE','REPRO_R4W_DONE','REPRO_SIA4_DONE','REPRO_SGLAB_DONE','REPRO_BS8_DONE','REPRO_A17SGV_DONE','REPRO_SIA4GP_DONE','REPRO_SI7_DONE','REPRO_QRSGJT_DONE','REPRO_ALL_OLD_DONE','REPRO_ALL_NEW_DONE']:
    print(f'  {s}: old={sum(1 for l in old if l.strip()==s)} new={sum(1 for l in new if l.strip()==s)}')
print('\n=== DECISIVE CAUTIONS (per config block) ===')
PAT=['completely within a single cluster','NOT IDENTIFIED in this design','not estimable for subgroup','all-zero variance-covariance',
     'standard errors were not identified','rank-deficient in the contrast direction','continuous-treatment event study is skipped',
     'zero variance at a subgroup margin','collapsed','fullprecisionslow','invalid name','conformability','no horizon of the']
def blocks(L):
    cur=None; B={}
    for l in L:
        m=re.match(r'^\. \* ======== ([A-Z0-9]+): config',l)
        if m: cur=m.group(1); B[cur]=[]
        elif cur: B[cur].append(l)
    return B
for name,L in (('OLD',old),('NEW',new)):
    for blk,lines in blocks(L).items():
        hits=[l for l in lines if (l.startswith('Caution') or l.startswith('Error') or l.startswith('Note')) and any(p in l for p in PAT)]
        print(f'  [{name} {blk}] {len(hits)} decisive caution line(s)')
        for h in hits: print(f'      {h[:330]}')
print('\n=== BIT-IDENTITY vs the suite dumps (test/compare.py --rtol 0) ===')
for leg in ('old','new'):
    r=subprocess.run(['python3',f'{R}/test/compare.py',f'{W}/suite_{leg}',f'{W}/{leg}_dumps','--rtol','0'],capture_output=True,text=True)
    print(f'--- {leg}: my {leg}_dumps vs test/suite/out/run_{leg} ---'); print(r.stdout.strip()[:3000])
