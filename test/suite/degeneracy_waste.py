#!/usr/bin/env python3
"""Find assessment configs that spend their budget on NUMERICALLY DETERMINED surfaces.

The degeneracy scanner's deadness test asks "is this cell missing or a sentinel?". A row whose
se is 2.33e-10 against an estimate of exactly 1 is none of those - it is a perfectly finite
number - so it reads as HEALTHY and the config is credited with testing something. It is not:
the contrast is pinned by the design, and the published se/t/p are floating-point residue.

This finds that class, using the package's OWN convention for it: wooldid_fcdegx
(wooldid.ado:13287) flags a variance at or below 1e-9 of the largest. Same tolerance here.

Then it cross-references INTENT, because a fixture built to exercise a degenerate surface is
doing its job and must not be "fixed":
  - degeneracy_expected.csv  DATASET_DESIGN / DATA_EDGE / PINNED_DEGENERATE / HIDDEN(UNRULED)
  - data_classes.csv         pathological configs
Anything numerically determined but NOT accounted for there is a candidate for a config edit.

usage: degeneracy_waste.py <suitedir> [<dumpdir>]     (default dumpdir: <suitedir>/out/run_new)

RESULT OF THE 2026-09-15 RUN, recorded so the next reader starts from the answer rather than
re-deriving it. Over 497 clean-rc configs on the closing dumps (prov c45ef13):

  * NUMERICALLY DETERMINED rows: 3 configs, all the ext_cov10_rcs_* family (subgroups(group)
    with fe(group_time), which absorbs the subgroup dimension so the pairwise contrast is pinned
    by the design). Ruled closed 2026-09-15: no package guard and no config change, because of
    the 187 configs using subgroups() exactly 4 are absorbed, so 183 carry genuine pairwise
    coverage and these cost nothing.
  * ZERO LIVE RESULTS ROWS (the chunk's "yields only degenerate results when it was intended to
    test other things"): 15 of 497, and NONE of them is wasted coverage --
      2  fastnose configs, where having no SEs is the entire point (scanner: INFO)
      3  PINNED_DEGENERATE by the 2026-09-13 ruling
      4  pathological/stress fixtures named for the behaviour they trigger
      6  lattice(all 3) or contreatcontrols(untreatedZero), which collapse the whole VCE
    The six are NOT caused by cluster(i) - cont_basic, cont_sg and cont_ss use contreat with
    cluster(i) and have live inference. They match the already-ticketed D-6 "lattice(all 3)
    all-zero VCV" candidate, and their surfaces are live-covered elsewhere: lattice 42 of 56
    configs live, makecfxplots 5 with cfx_live>0, untreatedZero 1.
  * THE ONE THIN SPOT, worth watching rather than fixing: contreatcontrols(untreatedZero) has a
    SINGLE config with live inference (ext_cont_ext_untreatedzero). If it ever goes dead the
    option loses live-inference coverage entirely and nothing here would flag it, because this
    tool reports degeneracy, not single-point-of-coverage.
"""
import sys, os, csv, glob, collections
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else '.')
import degeneracy_scan as D

SUITE = sys.argv[1] if len(sys.argv) > 1 else '.'
DUMPS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(SUITE, 'out', 'run_new')
TOL = 1e-9

def se_col(meta):
    cols = meta['cols']
    idx = {c: i + 1 for i, c in enumerate(cols)}
    est = idx.get('estimate', idx.get('b'))
    se = idx.get('se')
    if se is None and 'se_analytic' in idx: se = idx['se_analytic']
    return est, se

exp, _ = D.load_expected(os.path.join(SUITE, 'degeneracy_expected.csv'))
classes = D.load_classes(os.path.join(SUITE, 'data_classes.csv'))
cfgs = [os.path.join(SUITE, 'configs.csv')]
if os.path.isfile(os.path.join(SUITE, 'configs_ext.csv')): cfgs.append(os.path.join(SUITE, 'configs_ext.csv'))
rows = D.load_configs(cfgs)
byname = {r['config']: r for r in rows}

found = collections.OrderedDict()
for r in rows:
    cfg = r['config']
    p = os.path.join(DUMPS, cfg + '.csv')
    if not os.path.exists(p): continue
    try: d = D.load_dump(p)
    except Exception: continue
    if d.get('rc') != '0': continue
    hits = []
    for m, meta in d['matnames'].items():
        if not D.RESULT_MATS.match(m): continue
        est, se = se_col(meta)
        if est is None or se is None: continue
        cells = d['mat'][m]
        ses = []
        for row in range(1, meta['R'] + 1):
            v = cells.get((row, se), '.')
            if D.is_missing(v) or D.is_dead(v): continue
            try: ses.append(abs(float(v)))
            except Exception: pass
        if not ses: continue
        smax = max(ses)
        if smax <= 0: continue
        for row in range(1, meta['R'] + 1):
            lab = meta['rows'][row - 1] if row - 1 < len(meta['rows']) else f'r{row}'
            sv = cells.get((row, se), '.'); ev = cells.get((row, est), '.')
            if D.is_missing(sv) or D.is_dead(sv): continue
            try: s = abs(float(sv)); e = float(ev)
            except Exception: continue
            if s <= 0: continue
            # the package's own rule, applied two ways
            if s <= TOL * smax or (e != 0 and s / abs(e) <= TOL):
                hits.append((m, lab, e, s))
    if hits: found[cfg] = hits

print(f"scanned {len(rows)} configs in {DUMPS}")
print(f"configs with at least one NUMERICALLY DETERMINED row (se <= {TOL:g} of the matrix max, "
      f"or se/|estimate| <= {TOL:g}): {len(found)}\n")
for cfg, hits in found.items():
    e = exp.get(cfg)
    klass = e['class'] if e else '(no expectation row)'
    cite = (e['citation'][:52] if e else '')
    dcls = D.config_class(byname[cfg], classes) if cfg in byname else '?'
    mats = sorted({h[0] for h in hits})
    intended = klass in ('DATASET_DESIGN', 'DATA_EDGE') or cite.startswith('PINNED_DEGENERATE') \
               or dcls == 'pathological'
    tag = 'INTENDED' if intended else '>>> UNINTENDED'
    print(f"{tag:15s} {cfg:38s} rows={len(hits):3d} mats={','.join(mats)[:34]:34s} "
          f"class={klass:16s} data={dcls}")
    for m, lab, ev, sv in hits[:3]:
        print(f"{'':15s}   {m}[{lab}] estimate={ev:.6g} se={sv:.3g}")
    if len(hits) > 3: print(f"{'':15s}   ... {len(hits)-3} more")
