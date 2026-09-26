import json,re,sys
T='/Users/hegland/.claude/jobs/ab6fda17/tmp/trim'
b={x['idx']:x for x in json.load(open(T+'/blocks_f556138.json'))}
final={i:b[i]['smcl'] for i in b}
for g in ['G1','G2','G3','G4']:
    p={x['idx']:x for x in json.load(open(f'{T}/proposals_{g}.json'))['proposals']}
    for v in json.load(open(f'{T}/verdicts_{g}.json'))['verdicts']:
        i=v['idx']
        if v['verdict']=='accept': final[i]= '' if p[i]['action']=='delete' else p[i]['new_smcl']
        elif v['verdict']=='fixed_replace': final[i]=v['new_smcl']
def show(i): return re.sub(r'\s+',' ',final[i]).strip()
if __name__=='__main__':
    for a in sys.argv[1:]: print(f'[{a}]', show(int(a))); print()
