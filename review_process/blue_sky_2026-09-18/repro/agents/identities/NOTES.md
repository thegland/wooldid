# recon/identities — working notes (agent: IDENT)

Pinned build: .../bluesky/build/wooldid.ado (17021 lines). All line refs are THAT file.

## Code facts established (read, not guessed)

### Estimand weights (wooldid.ado:1570-1650)
Built AFTER estimation, on `esamp` = e(sample) (ado:1510-1511).
- `w` = user weight expression, or a tempvar==1 when none (ado:632-640); with NO weight,
  margins is handed `[iw=...]`; with aw/pw it is `[aweight=...]`/`[pweight=...]`; a user
  `fw` is REFUSED (ado:628) but the bootstrap rep path may pass one and it is demoted to
  aweight for margins while the REGRESSION keeps fweight (ado:619-632).
- `customweighttext` = `* cwm` if customweightmultiplier() given; folded into ALL FOUR
  weights but NOT into `w` itself.
- `posttreatment` (ado:1601-1605): 1 = treated unit, post; 2 = treated unit, pre and
  reported; 3 = never-treated in sample; 0 = treated-unit pre rows that are the pooled
  reference; 4 = out of sample.
- `w_per_it` = total(esamp*w) by (i,t)          [USER i and USER t, not _i/_t]
- `w_per_i`  = total(w*esamp) by (posttreatment,i)
- `t_per_i`  = total(tag(i,t)) by (posttreatment,i)   [# distinct i-t cells in the block]
- w_att      = w*cwm
- w_att_it   = w*cwm / w_per_it            -> each (i,t) CELL carries total weight 1
- w_att_i    = w*cwm / w_per_i             -> each UNIT carries total weight 1 in its block
- w_att_itime= w*cwm / (w_per_it*t_per_i)  -> each unit total 1, each of its cells equal
Missing -> 0 in all three.

### Margins subpops (all use `noestimcheck`)
- att family main  (ado:2030): subpop(ifstatement & targetflag==1 & t >= ttre), no over()
- att family pre   (ado:2089): subpop(... targetflag==1 & t < ttre & _t != 0)
- ES horizons      (ado:2153): subpop(... targetflag==1 & _t>=horizonstarter & _t<=horizonend),
                     over(_t), weight = w_att for item==att, ELSE w_att_it (ado:2144)
                     and computed ONCE per horizitem (cap confirm matrix, ado:2146-2149)
- subgroup levels  (ado:2819/3089): same subpop + over(subgroups); weights are the
  FULL-SAMPLE w_`item' (NOT renormalized inside the subgroup)

### targetflag / _t
`_t = t - ttre` then shifted to positive; reference rows -> 0; `targetflag = _t>0 & _t<.`
(ado:1009). Under the POOLED reference targetflag==1 includes PRE-window horizons, which is
why the att subpop adds `& t>=ttre` (ado:816, 12082-12180 in wooldid_reltime).
relyr = horiz - min_t - 1. Under pooled: horizonstarter=1, espre=min_t.

### gteffects (ado:1745-1780 + mata wooldid_gtx ado:13464-13540)
`e(gteffects)` cols: i onset caltime relyr estimate se p lb ub inwindow nobs w_att
w_att_it w_att_i w_att_itime. estimate = the RAW stage-1 coefficient on
`<g>._i#<h>._t#1.targetflag`; se = sqrt(vecdiag(e(V))) of the stage-1 fit; the four weight
columns are the CELL SUMS of the estimand weights over rows with e(sample) & targetflag==1.
Binary only (contreat refused, ado:4696).

### e(b)/e(V) (ado:4044-4050)
e(b) = column 1 of the pooled results matrix; e(V) = diag(se^2) — DIAGONAL ONLY.
Rows named main_att, main_att_it, ..., pretreatment_*, ES_att_post0..., ES_att_it_...,
and `_SG<k>_` variants. e(V_es_*) captured separately from margins r(V) (ado:2160-2175).

### fastnose (ado:1973-1990)
`marginsdof` -> `nose`. EXCEPTION: with subgroups() AND an ES window, `sgmarginsdof` keeps
`df(df_r)` (ado:1981) — so the subgroup-by-horizon call still computes a VCE in fast mode.
`f2bonly` = fastnose & no permute -> e(b_es_*) posted without e(V_es_*).

### replay (ado:27-53 vs live ado:4491-4500)
Live run: display part(1); wooldid_ptp (prints Roth itself); display part(2).
Replay:   display part(1); display part(roth); display part(2); part(perm); part(boot).
=> the Roth block is printed by TWO DIFFERENT code paths. Prime divergence candidate.

### permexact (mata ado:13640-13790)
Enumerates all distinct multiset permutations of the cluster-level onset vector within
permstrata; M = prod over strata of multinomial coefficients; ceiling 10000. Deterministic:
no RNG. Note the code's own comment: M (log-factorial formula) and cols(A) (recursive
enumerator) are computed by INDEPENDENT algorithms and nothing reconciles them; both bucket
`.` as one never-treated level but uniqrows keeps EXTENDED missings (.a,.b) distinct ->
mismatch lead (hand to the permute prober).

### rseed
`set seed rseed` inside wooldid_permrun (ado:7471) and wooldid_bootrun (ado:9015),
i.e. AFTER the observed fit. So same rseed+reps => reproducible.
wooldid_esband saves/restores c(rngstate) and uses a FIXED internal seed (ado:10892,10937)
=> simultaneousbands must be deterministic even with no rseed.

### poisson scale remap (ado:5031-5045)
Under `poisson` WITHOUT `poisexpresults`: semielasticity(eydx)->dydx, ctelas eydx->dydx,
eyex->dyex. So `poisson semielasticity` must print the SAME NUMBERS as plain `poisson`.

### anticipation (ado:743-755)
`gen double ttre = ttre_asgiven - k if ttre_asgiven - k <= maxyearobserved+oos & ifstatement`
Comment claims bit-identity with the hand workaround `gen ttre2 = ttre - k`.
But wooldid_validate's window screen shifts DIFFERENTLY (ado:5347-5357: wdvtlo = wdvtmin + k).

### contreat model (ado:1117, 1476)
`i._i#i._t#i.targetflag##c.dose` -> cell intercept delta_gt AND cell#dose slope beta_gt.
dydx(targetflag) at a row = delta_gt + beta_gt*x_r ; dydx(contreat) = beta_gt ;
dyex = beta_gt*x_r ; eydx = (.)/yhat_r ; eyex = (.)*x_r/yhat_r.
=> elasticity variants are NOT deterministic transforms of each other unless x (or yhat)
is constant over the margins subpop.

## Status
- [x] read weight construction, margins sites, gteffects, e(b)/e(V), fastnose, replay,
      permexact, rseed, poisson remap, anticipation, contreat model
- [ ] write /recon/identities.md

## RESULTS (2026-09-16)
- agg.do/.log  : I1 gteffects aggregation VERIFIED for all 4 estimands x {unweighted, aw, aw+cwm},
                 post and pre blocks. max reldiff 1.2e-15. Checker mata fn `wdaggchk` lives there.
- mech.do/.log : e(b_es_att)/e(V_es_att) post fine with gteffects on (rc 0 — the earlier
                 "missing" was MY bug: `capture` suppresses output).
                 fastnose vs full points: max reldiff 2.4e-16.
                 single cohort balanced -> all 4 estimands BIT-IDENTICAL.
                 single cohort unbalanced -> att==att_it bitwise, att_i==att_itime bitwise.
- sg1.do/.log  : subgroups(one) with one==1 -> reldiff 0 on all 22 rows vs no-subgroups run;
                 every _SG1_ row == its overall twin.
- esmat.do/.log: representative `ereturn list` (names of every posted matrix/scalar).

## DELIVERED
/private/tmp/.../scratchpad/bluesky/recon/identities.md  (20 ranked identities + a
"NOT an identity" list + environment gotchas). DONE.
