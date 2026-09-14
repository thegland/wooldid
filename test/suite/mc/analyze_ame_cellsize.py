#!/usr/bin/env python3
"""Judge mc_ame_cellsize.do: is the continuous adversarial line test correctly sized, and at
what cell size?

wooldid.sthlp:745 says the continuous line test requires cells of "on the order of fifty"
observations. mc_ame_cellsize.do varies exactly that quantity. This scores each run.

Judged the way the house analyzers judge: the rejection COUNT against the central ~95%
region of Binomial(n, alpha), rounded outward. Both p-values here come from an ANALYTIC arm,
so they are continuous - there is no lattice to attain and nominal alpha is the right target.
The att column is the CONTROL: it comes from the same fit, so if both families were low the
DGP or the harness would be at fault rather than the ame variance path.

A Kolmogorov-Smirnov distance from Uniform(0,1) is reported alongside the tail counts because
the two can disagree: at fifty observations per cell the ame p passes KS while still rejecting
about twice as often as nominal in the far-left tail, which is the part a user relies on.

usage: analyze_ame_cellsize.py <csv>      e.g. out3/cellsize/probe_ame_u50.csv
"""
import sys, csv, math

def binom_region(n, p, mass=0.95):
    # central region: smallest [lo,hi] about the mean with coverage >= mass, rounded outward
    pmf = [math.comb(n, k) * p**k * (1-p)**(n-k) for k in range(n+1)]
    mean = int(round(n*p)); lo = hi = mean; tot = pmf[mean]
    while tot < mass:
        dn = pmf[lo-1] if lo > 0 else -1
        up = pmf[hi+1] if hi < n else -1
        if up >= dn: hi += 1; tot += up
        else:        lo -= 1; tot += dn
    return lo, hi, tot

path = sys.argv[1] if len(sys.argv) > 1 else 'probe_ame.csv'
rows = [r for r in csv.DictReader(open(path))]
print(f'rows read: {len(rows)}  from {path}')

bad = [r for r in rows if r['rc'] != '0']
absent = [r for r in rows if r.get('ameabsent') not in ('0', '')]
print(f'P0 integrity: rc!=0 on {len(bad)} rows; ame p absent on {len(absent)} rows')

ps, patts = [], []
for r in rows:
    for col, acc in (('p_ame', ps), ('p_att', patts)):
        v = r[col]
        if v not in ('.', ''):
            try: acc.append(float(v))
            except ValueError: pass
oob = [p for p in ps if not (0.0 <= p <= 1.0)]
print(f'ame p in [0,1]: {len(ps)-len(oob)} of {len(ps)}' + (f'  OUT OF RANGE: {oob[:5]}' if oob else ''))

MINN = 30   # the house minimum-n gate; below it a cell reports INSUFFICIENT_N, never a pass
for label, vals in (('e(adversarial_pval_ame)', ps), ('e(adversarial_pval_att)  [control]', patts)):
    print(f'\n== {label}   n = {len(vals)}')
    if len(vals) < MINN:
        print(f'   INSUFFICIENT_N (n = {len(vals)} < {MINN}) - not scored')
        continue
    for alpha in (0.05, 0.10):
        k = sum(1 for p in vals if p <= alpha)
        lo, hi, cov = binom_region(len(vals), alpha)
        verdict = 'PASS' if lo <= k <= hi else 'FAIL'
        rate = k/len(vals)
        # exact two-sided binomial p for the observed count
        n = len(vals)
        pmf = [math.comb(n, j) * alpha**j * (1-alpha)**(n-j) for j in range(n+1)]
        pk = pmf[k]
        pval = sum(x for x in pmf if x <= pk*(1+1e-12))
        print(f'   alpha {alpha:.2f}: {k:4d} / {n} = {rate:.4f}   band [{lo}, {hi}] '
              f'(cov {cov:.3f})   {verdict}   exact two-sided p = {pval:.3g}')
    q = sorted(vals)
    def qt(f): return q[min(len(q)-1, int(f*len(q)))]
    print(f'   quantiles  min {q[0]:.4f}  q25 {qt(.25):.4f}  med {qt(.5):.4f}  '
          f'q75 {qt(.75):.4f}  max {q[-1]:.4f}')
    # a correctly sized continuous p is Uniform(0,1): report the KS distance too
    n = len(q)
    ks = max(max(abs((i+1)/n - p), abs(i/n - p)) for i, p in enumerate(q))
    crit = 1.358/math.sqrt(n)   # 5% two-sided asymptotic KS critical value
    print(f'   KS distance from Uniform(0,1): {ks:.4f}  (5% critical {crit:.4f})  '
          f'{"UNIFORM not rejected" if ks < crit else "UNIFORMITY REJECTED"}')
