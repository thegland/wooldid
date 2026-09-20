# Fix proposals — mechanism `PTP` (`pretrendpower()`, Roth 2022 diagnostics)

Pinned build `.../bluesky/build/wooldid.ado` (sha c1d8979). All line numbers are that build's.
Brief: `fixprop/mech/PTP.md` — two findings, **PTP-1** (high) and **BX-N3** (high).

Everything below was patched into **my own copy** of the ado
(`probes/prop-PTP/patched/wooldid.ado`) and RUN. Nothing was written inside the repo.

| what | file |
|---|---|
| ideal-model verification of all three legs, 200k draws | `probes/prop-PTP/a1_ideal.py` |
| re-analysis of the prober's + refuter's per-replication `W` already on disk | `probes/prop-PTP/a2_repanalysis.py` |
| pinned vs patched, 5 arms, `%21x` | `probes/prop-PTP/p1_beforeafter.do` / `.log` |
| BX-N3 geometry + the new caution + replay | `probes/prop-PTP/p2_bxn3.do` / `.log` |
| runtime cost of everything added | `probes/prop-PTP/p3_cost.do` / `.log` |
| the PS-64 cross-entry-point identity, before and after | `probes/prop-PTP/p4_ps64.do` / `.log` |
| Monte Carlo: does the patched MDT deliver its claimed power | `probes/prop-PTP/p5_mc.do` / `.log` |
| offline rejection analysis of those 750 replications | `probes/prop-PTP/a3_mc.py` |
| the hardened BX-N3 variant (Leg 1b) and the mixed-geometry arms | `probes/prop-PTP/patched2/` , `p6_mixed.do` / `.log` |
| re-run of the ideal-model check, reproduced digit-for-digit | `probes/prop-PTP/a1_ideal_rerun.txt` |
| live-vs-replay on the Hotelling branch, 3 arms | `probes/prop-PTP/p7_replay.do` / `.log` |
| which VCE entry points the gate covers, 5 arms | `probes/prop-PTP/p8_gate.do` / `.log` , `p9_gate2.do` / `.log` |
| working notes | `probes/prop-PTP/NOTES.md` |

**Verdict in one line each.** **PTP-1** is real, is arithmetic rather than a few-cluster artefact,
and has an exact closed-form fix that I patched, ran, and measured to deliver its claim
(size 0.084 -> 0.056, power at the reported MDT 0.696 -> 0.804 against a claim of 0.80); it
**moves reported numbers** and needs an author ruling. **BX-N3** needs no number to move: the
covariance is the correct plug-in for structurally-zero residuals and cannot be repaired, so the
fix is a cause-naming sentence that costs nothing on any clean run.

---

# PROPOSAL 1 — `wooldid_ptp` refers an ESTIMATED covariance to KNOWN-covariance distributions
*covers:* **PTP-1** (and the refuter's SKP-NEW-1, which is the same defect correctly diagnosed).
One fix site: `wooldid_ptp`, `ado:10341-10700`, with a mandatory mirror in the replay
renderer at `ado:11455-11500`.

## What is actually wrong (settled, not re-litigated)

`wooldid_ptp` computes `S = e(V_es_*)`'s pre-block, then

| quantity | line | formula | exact only if |
|---|---|---|---|
| joint pre-trend p | 10519 | `Ftail(keff, df_r, W/keff)`, `W = b'S^-1 b` | `S` is the TRUE covariance |
| `mdt_wald` cutoff | 10559 | `invF(keff, df_r, .95)` + `npnF` | same |
| slope CI | 10512-10514 | `slope ± invttail(df_r,.025)/sqrt(m'S^-1 m)` | same |

With `S` estimated carrying `nu` Wishart degrees of freedom, the exact reference for
`b'S^-1 b` is **Hotelling's T-squared**: `W*(nu-k+1)/(nu*k) ~ F(k, nu-k+1)`.

**Exact size of wooldid's own "5% level" cutoff** (`a1_ideal.py` PART A, closed form):

```
  k \ nu :     15      23      39      79     279
  k = 1  : 0.0500  0.0500  0.0500  0.0500  0.0500
  k = 2  : 0.0611  0.0570  0.0540  0.0519  0.0505
  k = 3  : 0.0784  0.0673  0.0597  0.0546  0.0513
  k = 4  : 0.1034  0.0813  0.0670  0.0579  0.0521
  k = 6  : 0.1840  0.1225  0.0868  0.0663  0.0543
```

A 16-state DiD with `cluster(st)` and `esprelength(3)` sits on the `k=3, nu=15` cell.
`esprelength(6)` on the same panel sits on `k=6, nu=15`: **0.1840**.

## THE BRANCH CONDITION — the one place the refuter is wrong, and it decides the edit

The refuter's summary says the damage is "present at EVERY `df_r` and has nothing to do with
how the covariance was clustered". **That is false on the iid branch, and the code must not
be changed there.** Under iid OLS the pre-block is `s^2 * (a fixed matrix)`, i.e. a SCALAR
multiple of Sigma, not a Wishart. Then

* `W/k = [chi2(k)/k] / [chi2(df_r)/df_r] = F(k, df_r)` **exactly**, and
* the GLS slope does not depend on `s^2` at all, so `(slope-theta)/se_hat ~ t(df_r)` **exactly**.

That is why the prober's Design A (iid, no cluster, `df_r`=279) came back dead on: size
0.0367, power 0.8000 at `mdt_wald(80)` and 0.5000 at `mdt_wald(50)`, slope CI coverage 0.9600.
**Applying a Hotelling scaling there would make a currently-exact test conservative.**

So the whole proposal is gated on the covariance being a cluster sandwich. `e(N_clust)` is
posted at `ado:4281` exactly when a cluster-robust VCE was in force and is visible inside
`wooldid_ptp` (which runs after the post). It is the entire gate:

```stata
  local ptphot = 0
  if `ptpdfok' & !missing(e(N_clust)) & "`e(boot_method)'" == "" local ptphot 1
```

Plain `robust` (no `cluster()`) is a sandwich too, but its Wishart df is of order `N`, where
the correction is numerically nil; it is left out so that the re-bless has one gate and one
story. Under `bootstrap()` the covariance is built from `e(boot_N_V)` replications, for which
`nu = df_r` is the wrong idealization (that is the refuter's SKP-NEW-2, filed uncertain), so
the bootstrap path is excluded and left exactly as it is.

**Which entry points the gate actually covers — measured, not inferred** (`p8_gate.log`,
`p9_gate2.log`, patched build, same panel, `k = 3`). `ado:1721` takes `Nclus` from the
*estimation command's* own `e(N_clust)`, which is why a clustered VCE reached by any route trips
the gate and an unclustered one cannot:

| how the VCE was asked for | `e(N_clust)` | gate | printed reference | `e(roth_slope_dfci_att)` |
|---|---|---|---|---|
| `cluster(st)` | 16 | **ON** | `F(3, 13)`, p = 0.1484 | 13 |
| `customvce(cluster st)` | 16 | **ON** | `F(3, 13)`, p = 0.1504 | 13 |
| `unconditionalse cluster(st)` | 16 | **ON** | `F(3, 13)`, p = 0.3672 | 13 |
| `unconditionalse robust` | `.` | **off** | `F(3, 143)`, p = 0.5423 | 143 (= `df_r`) |

(`unconditionalse` alone is refused at parse — "requires specifying either robust or clustered
standard errors" — so there is no fifth case.)

## The edit

### Leg A — the reference distribution (joint test p, and the noncentrality `mdt_wald` inverts)

**Site 1**, `ado:10369-10370`, immediately after `ptpniscrit`. **Add:**

```stata
  * Hotelling gate. Under iid OLS the pre-block is s2 times a FIXED matrix, so W/k is
  * exactly F(k, df_r) and nothing below applies. Under a cluster-robust sandwich the
  * block is an ESTIMATED covariance carrying ~G-1 Wishart degrees of freedom, and the
  * known-Sigma references are liberal. e(N_clust) is posted (ado:4281) exactly when a
  * cluster-robust VCE was in force; it is the whole gate.
  local ptphot = 0
  if `ptpdfok' & !missing(e(N_clust)) & "`e(boot_method)'" == "" local ptphot 1
```

**Site 2**, `ado:10508-10520`. **Now:**

```stata
    matrix `ptpwork' = `ptpm''*`ptpSi'*`ptpbp'
    local ptpslope = `ptpwork'[1, 1]/`ptpc2'
    local ptpse = 1/sqrt(`ptpc2')
    local ptpll = `ptpslope' - `ptpcrit'*`ptpse'
    local ptpul = `ptpslope' + `ptpcrit'*`ptpse'
    matrix `ptpwork' = `ptpbp''*`ptpSi'*`ptpbp'
    local ptpW = `ptpwork'[1, 1]
    if `ptpdfok' local ptpprep = Ftail(`ptpkeff', `dfr', `ptpW'/`ptpkeff')
    else local ptpprep = chi2tail(`ptpkeff', `ptpW')
```

**Becomes** (Leg A and Leg C together; Leg C is the `ptpseci`/`ptpdfci` block):

```stata
    matrix `ptpwork' = `ptpm''*`ptpSi'*`ptpbp'
    local ptpslope = `ptpwork'[1, 1]/`ptpc2'
    local ptpse = 1/sqrt(`ptpc2')
    matrix `ptpwork' = `ptpbp''*`ptpSi'*`ptpbp'
    local ptpW = `ptpwork'[1, 1]
    * Hotelling degrees of freedom; undefined once keff reaches df_r, so fall back there
    local ptphdf = `dfr' - `ptpkeff' + 1
    local ptphotf = `ptphot'
    if `ptphotf' & `ptphdf' < 1 local ptphotf 0
    * slope CI: EXACT pivot (Anderson step-down decomposition of Hotelling's T2) -
    *   (slope-theta)^2 * c2 * (df_r-keff+1) / (df_r + Q)  ~  F(1, df_r-keff+1),
    * Q = W - slope^2*c2 being the residual quadratic form after the m direction is
    * removed, which is free of theta and so inverts to a closed-form interval.
    * Collapses EXACTLY to the shipped formula at keff = 1 (Q = 0, df_r-keff+1 = df_r).
    local ptpseci = `ptpse'
    local ptpcritci = `ptpcrit'
    local ptpdfci = cond(`ptpdfok', real("`dfr'"), .)
    if `ptphotf' {
      local ptpQ = `ptpW' - `ptpslope'*`ptpslope'*`ptpc2'
      if `ptpQ' < 0 local ptpQ 0
      local ptpseci = `ptpse'*sqrt((`dfr' + `ptpQ')/`ptphdf')
      local ptpcritci = invttail(`ptphdf', (100 - `level')/200)
      local ptpdfci = `ptphdf'
    }
    local ptpll = `ptpslope' - `ptpcritci'*`ptpseci'
    local ptpul = `ptpslope' + `ptpcritci'*`ptpseci'
    if `ptphotf' local ptpprep = Ftail(`ptpkeff', `ptphdf', `ptpW'*`ptphdf'/(`dfr'*`ptpkeff'))
    else if `ptpdfok' local ptpprep = Ftail(`ptpkeff', `dfr', `ptpW'/`ptpkeff')
    else local ptpprep = chi2tail(`ptpkeff', `ptpW')
```

**Site 3**, `ado:10559-10560` (the noncentrality). **Now:**

```stata
      if `ptpdfok' local ptplam = npnF(`ptpkeff', `dfr', invF(`ptpkeff', `dfr', .95), 1 - `ptplev'/100)
      else local ptplam = npnchi2(`ptpkeff', invchi2(`ptpkeff', .95), 1 - `ptplev'/100)
```

**Becomes:**

```stata
      if `ptphotf' local ptplam = npnF(`ptpkeff', `ptphdf', invF(`ptpkeff', `ptphdf', .95), 1 - `ptplev'/100)
      else if `ptpdfok' local ptplam = npnF(`ptpkeff', `dfr', invF(`ptpkeff', `dfr', .95), 1 - `ptplev'/100)
      else local ptplam = npnchi2(`ptpkeff', invchi2(`ptpkeff', .95), 1 - `ptplev'/100)
```

**Site 4**, `ado:10620` (the displayed line must name the reference it now uses):

```stata
    if `ptphotf' dis as result "    joint pre-trend test: F(`ptpkeff', `ptphdf') = " %-8.0g `=`ptpW'*`ptphdf'/(`dfr'*`ptpkeff')' ", p = " %-6.4f `ptpprep'
    else if `ptpdfok' dis as result "    joint pre-trend test: F(`ptpkeff', `dfr') = " %-8.0g `=`ptpW'/`ptpkeff'' ", p = " %-6.4f `ptpprep'
    else dis as result "    joint pre-trend test: chi2(`ptpkeff') = " %-8.0g `ptpW' ", p = " %-6.4f `ptpprep'
```

### Leg B — the plug-in `m'S^-1 m` makes `mdt_wald` systematically optimistic

Separable from Leg A and separately rulable. `c2_hat/c2_true` is distributed `df_r/chi2(df_r-k+1)`
(a standard Wishart result), so `sqrt(lambda/c2_hat)` is a **downward-biased** estimate of the
design's detectable trend, with the exact law `mdt_hat/mdt_true ~ sqrt(chi2(df_r-k+1)/df_r)`.
At `k=3, nu=15` its median is 0.907 and its 5-95 range is [0.63, 1.22].

**Verified, ideal model, 200,000 draws** (`a1_ideal.py` PART B, re-run on 2026-09-19 and
reproduced digit-for-digit; median of the reported `mdt_wald(80)` divided by the design's true
80%-power trend — 1.000 is the target):

| k, nu | shipped | + Leg A only | + Leg A and Leg B |
|---|---|---|---|
| k=3, nu=15 | 0.8888 | 0.9073 | **1.0003** |
| k=4, nu=15 | 0.8366 | 0.8695 | **1.0000** |
| k=6, nu=39 | 0.9142 | 0.9245 | **1.0000** |
| k=3, nu=79 | 0.9823 | 0.9829 | **1.0000** |
| k=3, nu=279 | 0.9951 | 0.9951 | 0.9999 |
| k=1, nu=15 | 0.9768 | 0.9768 | 0.9990 |

Read the `k=1` row: Leg A cannot move it (the Hotelling reference *is* the shipped one at `k=1`),
and the column is still 2.3% optimistic there — that is the whole of Leg B, and it is why the two
legs are separable and separately rulable. Read the middle column against the right one: **Leg A
alone leaves 9-13% of the MDT optimism in place at `nu = 15`**. The pair is what delivers the
number the column claims to be.

**Site 5**, `ado:10561-10564`. **Now:**

```stata
      if !missing(`ptplam') {
        matrix `ptpmat'[`ptplr', 2] = sqrt(`ptplam'/`ptpc2')
        if !missing(`ptphbar') matrix `ptpmat'[`ptplr', 4] = sqrt(`ptplam'/`ptpc2')*(`ptphbar' - `ptpref')
      }
```

**Becomes:**

```stata
      * the plug-in m'S^-1 m is INFLATED (c2hat/c2true ~ df_r/chi2(df_r-keff+1)), so the
      * shipped sqrt(lambda/c2hat) is a downward-biased estimate of the design's
      * detectable trend. Rescale to median-unbiased.
      local ptpc2d = `ptpc2'
      if `ptphotf' local ptpc2d = `ptpc2'*invchi2(`ptphdf', 0.5)/`dfr'
      if !missing(`ptplam') {
        matrix `ptpmat'[`ptplr', 2] = sqrt(`ptplam'/`ptpc2d')
        if !missing(`ptphbar') matrix `ptpmat'[`ptplr', 4] = sqrt(`ptplam'/`ptpc2d')*(`ptphbar' - `ptpref')
      }
```

### Leg C — the slope CI (an EXACT pivot; this is new, the refuter said none existed)

Included in Site 2 above. The refuter's write-up concludes "no closed form for
GLS-with-estimated-S". There is one, and it is cheap. Write `Q = W - slope^2*c2`, the residual
quadratic form after the `m` direction is projected out. `Q` does not involve `theta`, and
Anderson's step-down decomposition of Hotelling's T-squared gives

```
   (slope - theta)^2 * c2 * (nu - k + 1) / (nu + Q)   ~   F(1, nu - k + 1)
```

exactly, which inverts to

```
   slope  ±  invttail(nu-k+1, alpha/2) * (1/sqrt(c2)) * sqrt((nu + Q)/(nu - k + 1))
```

**Verified, 200,000 draws of the ideal model** (`a1_ideal.py`, coverage at nominal 0.95):

| k, nu | shipped | proposed |
|---|---|---|
| k=3, nu=15 | 0.9111 | **0.9507** |
| k=4, nu=15 | 0.8842 | **0.9507** |
| k=6, nu=39 | 0.9139 | **0.9509** |
| k=3, nu=79 | 0.9434 | **0.9494** |
| k=3, nu=279 | 0.9476 | **0.9493** |
| k=1, nu=15 | 0.9501 | 0.9501 (unchanged, as it must be) |

It is exact at every `(k, nu)`, not merely better, and it is algebraically identical to the
shipped formula at `k = 1`.

### Ereturn / display companions (part of Leg C, mandatory)

`e(roth_slope_se_*)` must carry the adjusted standard error, or `CI = slope ± t*se` — which
users check by hand and which the replay renderer at `ado:11487` reconstructs from `e()` —
stops holding. So:

```stata
    cap ereturn scalar roth_slope_se_`ptpty' = `ptpseci'
    cap ereturn scalar roth_slope_dfci_`ptpty' = `ptpdfci'
    cap ereturn hidden scalar roth_slope_seplugin_`ptpty' = `ptpse'
```

`roth_slope_dfci_*` is new and public (a user running `lincom`-style arithmetic needs the df);
`roth_slope_seplugin_*` is hidden and keeps the shipped number reachable for any comparison.

### MANDATORY MIRROR — the replay renderer

`ado:11488` rebuilds the F line from `e(roth_W_*)` and `e(df_r)` with the **shipped** formula.
If only the live path changes, `wooldid` typed alone after the fit prints a different line from
the one the fit printed. Replay must get the same gate:

```stata
    local ptphotr = (`ptpdfok' & !missing(e(N_clust)) & "`e(boot_method)'" == "" & `dfr' - `ptpkeff' + 1 >= 1)
    if `ptphotr' dis as result "    joint pre-trend test: F(`ptpkeff', `=`dfr'-`ptpkeff'+1') = " %-8.0g `=`ptpW'*(`dfr'-`ptpkeff'+1)/(`dfr'*`ptpkeff')' ", p = " %-6.4f `ptpprep'
    else if `ptpdfok' dis as result "    joint pre-trend test: F(`ptpkeff', `dfr') = " %-8.0g `=`ptpW'/`ptpkeff'' ", p = " %-6.4f `ptpprep'
    else dis as result "    joint pre-trend test: chi2(`ptpkeff') = " %-8.0g `ptpW' ", p = " %-6.4f `ptpprep'
```

The slope line and the power matrix need no change once the ereturns carry the adjusted numbers,
and this is worth stating precisely because it is the obvious place for a live-vs-replay
divergence: `ado:11486` prints `ptpll`/`ptpul` read **straight out of** `e(roth_slope_ll_*)` /
`e(roth_slope_ul_*)`. It does **not** rebuild the interval as `slope ± invttail(e(df_r), a/2)*se`,
so there is no second critical value to keep in step. (`e(roth_slope_dfci_*)` is therefore posted
for the *user's* arithmetic, not for the renderer's.)

**Tested — and this is the branch the earlier draft of this document cited the wrong log for.**
`p2_bxn3.log`'s replay arm is `robust`, so `ptphotr = 0` there and the Hotelling replay line never
ran. `p7_replay.do` / `.log` exercises it directly, three arms on the patched build, live line
against replayed line:

| arm | live | replay |
|---|---|---|
| `cluster(st)`, `k=3` (gate ON) | `slope = -.2521666 (se .3056978; 95% CI -.9125866 to .4082535)` / `F(3, 13) = 2.10927, p = 0.1484` | **character-identical** |
| `cluster(st)`, `k=1` (gate on, `hdf = df_r`) | `F(1, 15) = 6.18344, p = 0.0252` | **character-identical** |
| no `cluster()`, `k=3` (gate OFF) | `F(3, 132) = 1.02025, p = 0.3859` | **character-identical** |

Note the middle row: at `k = 1` the Hotelling df is `15 - 1 + 1 = 15`, so the branch prints the
same `F(1, 15)` label the shipped code prints — the `k = 1` identity holds in the renderer too,
not only in the estimator.

### Help edits that go with it (`sthlp:789-791`)

1. "`mdt_wald` powers the joint pre-trend test, **an F test with the regression's residual
   degrees of freedom**" — becomes "an F test referred to `F(keff, df_r)` under iid errors and,
   under `cluster()`, to the Hotelling small-sample reference `F(keff, df_r-keff+1)`".
2. "(the minimum detectable trend itself **depends only on the design's precision, never on
   the estimated coefficients**)" — true of the *estimand*, false of the *number printed*.
   It is an estimate of the design's precision, formed from the estimated covariance, with a
   23% coefficient of variation at `G = 16, k = 3`. This sentence should say so explicitly;
   it is the sentence both the prober and the refuter singled out.

## MEASURED before / after

### (i) On my own 16-state panel, one fit, `%21x` (`p1_beforeafter.log`)

`wooldid y coh t ttre, fe(st) cluster(st) esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(50 80)`
— `df_r = 15`, `e(N_clust) = 16`, `keff = 3`, **0 cautions print on this geometry**
(independently confirming the prober's central claim that nothing warns the user here).

| | pinned | patched | ratio |
|---|---|---|---|
| `roth_slope_att` | -0.2521665900889164 | -0.2521665900889164 | **bit-identical** |
| `roth_W_att` | 7.30131420641624 | 7.30131420641624 | **bit-identical** |
| `roth_slope_se_att` | 0.2394292468629917 | 0.3056978488599751 | x1.2768 |
| CI | [-0.76250, +0.25816] | [-0.91259, +0.40825] | |
| `roth_prep_att` (joint p) | 0.1052063928135773 | 0.1484382479273842 | x1.4109 |
| `mdt_wald(80)` | 0.9018458960326692 | 1.0149413237942049 | x1.1254 |
| `bias_wald(80)` | 2.946029927040053 | 3.3154749910610697 | x1.1254 |
| `roth_slope_dfci_att` | (did not exist) | 13 | new |

At `esprelength(4)` (`k=4`) on the same panel: joint p x2.2499, `mdt_wald(80)` x1.1953,
slope se x1.5339.

*(One number in this table needs a footnote so it is not read as a contradiction of (iv): this
single fit's shipped `mdt_wald(80)` is 0.902, while the median over the 250 fits in (iv) — drawn
from the **same** generator — is 0.553. That spread is not an inconsistency, it is the finding's
own measurement restated: the reported MDT has a 23% coefficient of variation on this design and
ranges 0.296-0.925 over repeated fits. It is the strongest single argument for the `sthlp:791`
correction below.)*

### (ii) The three no-change controls — all three are exact, not approximate

| arm | what it proves | result |
|---|---|---|
| `esprelength(1)`, `cluster(st)` (`k=1`) | Hotelling collapses to `F(1,df_r)` and `Q = 0` | `slope`, `se`, `ll`, `ul`, `prep`, `W`, `nissize` ALL **bit-identical** |
| no `robust`, no `cluster()` (iid, `df_r`=132) | the gate holds | every scalar and every `mdt`/`bias` cell **bit-identical** |
| `robust`, no `cluster()` (`df_r`=132) | the gate holds | every scalar and every `mdt`/`bias` cell **bit-identical** |

The only thing that moves on those three arms is the new `e(roth_slope_dfci_*)` scalar, which
did not exist before. At `k=1` the one number that does move is `mdt_wald`, and that is Leg B,
not Leg A — the `c2` plug-in bias is real at `k=1` too (x1.0228 there).

### (iii) Re-analysis of the 1,100 null replications already on disk (`a2_repanalysis.py`)

The PTP prober and the PTP refuter each recorded per-replication `W`. Applying the two cutoffs
offline (`W` cutoff `9.86215` -> `11.80569` at `k=3, nu=15`) needs no new Stata at all:

| null stream | n | size now | size proposed |
|---|---|---|---|
| prober, pooled (E0 400 + repro 200) | 600 | 0.1050 [.082,.132] | **0.0683 [.049,.092]** |
| refuter, pooled (s3 250 + s4 250) | 500 | 0.1200 [.093,.152] | **0.0800 [.058,.107]** |

and the reported column, recomputed exactly (max relative error vs wooldid's own `mdt80`
column: 1.22e-15):

```
   mdt_wald(80), median over the prober's 400 Design-E fits
       shipped                        0.54630
       + Leg A (Hotelling)            0.55763
       + Leg A and Leg B              0.61481
   the refuter's honest 80% point (their MC covariance)   0.63007
   the prober's honest 80% point (size-corrected interp)  0.657
```

**Be honest about what is left.** The proposal removes about half to two thirds of the measured
excess size and about 80% of the measured MDT optimism. The residual (size ~0.068-0.080 against
0.05; MDT 0.615 against 0.630) is genuine few-cluster CRVE heterogeneity on a cohort-level
contrast identified off four states per cohort. **That residual is the already-ruled item**
(`known_bydesign` S1: analytic SEs understated at few clusters, no plug-in correction) and this
proposal does not claim to remove it. What it does remove is the part that is arithmetic and
that is present even in the best possible case.

### (iv) Monte Carlo on the patched build — does the corrected column deliver its claim?

`p5_mc.do` / `p5_mc.log` / `a3_mc.py`, seed 20260919 — an **independent seed stream** from both
the prober's and the refuter's. 250 replications per arm on the Design-E geometry (16 states,
4 onset cohorts, `i` = cohort, `fe(st)`, `cluster(st)`, AR(1) within state, `esprelength(3)`,
so `k = 3`, `df_r = 15`, `G = 16`). Every replication runs the **patched** build and records
`e(roth_W_att)`, the new hidden `e(roth_slope_seplugin_att)`, `e(roth_slope_se_att)` and the
`mdt_wald(80)` cell. Three arms: `delta = 0`; `delta = 0.55256`, the median **shipped**
`mdt_wald(80)` on this design; `delta = 0.62185`, the median **patched** one. `W` is recorded
per replication, so all three cutoffs are applied offline to the same draws.

`W` cutoffs at `k=3, nu=15`: shipped `k*invF(3,15,.95)` = **9.8622**; Hotelling
`nu*k*invF(3,13,.95)/(nu-k+1)` = **11.8057**; and the arm's own empirical 95th percentile of
`W` under the null = **11.9320**. *The honest cutoff the data demands and the cutoff the
proposal computes from first principles agree to 1.1%*, while the shipped one is 21% too small.

| arm (n = 250) | delta | shipped cutoff | Hotelling cutoff (proposed) | empirical size-corrected |
|---|---|---|---|---|
| null | 0 | **0.0840** [.0556,.1250] | **0.0560** [.0336,.0918] | 0.0520 |
| at the **shipped** mdt80 | 0.5526 | 0.7600 [.7034,.8088] | 0.6960 [.6364,.7497] | 0.6920 |
| at the **patched** mdt80 | 0.6219 | 0.8720 [.8249,.9079] | **0.8040** [.7503,.8484] | 0.7960 |

Read the diagonal: **the patched `mdt_wald(80)` is detected by the patched test with power
0.8040, CP95 [.750,.848], against a claim of 0.80** — and 0.796 after size correction. The
shipped `mdt_wald(80)` is detected by an honestly-sized test with power 0.692-0.696, not 0.80.
The null row is the same statement from the other side: the shipped reference rejects at 0.0840
with a CP95 interval that excludes 0.05, the proposed reference at 0.0560 with one that contains
it.

This stream's shipped size (0.0840) is milder than the prober's (0.1050) and the refuter's
(0.1200); the correction is measured on its own draws and does not depend on which stream the
reader believes.

## Does it move a reported number? YES — and exactly which

Under `cluster()` (or `customvce(cluster ...)`) with `pretrendpower()` and `keff >= 2`:

| moves | by |
|---|---|
| `e(roth_prep_*)` and the printed `joint pre-trend test:` line (value AND its df label) | large; x1.41 at `k=3,G=16`, x2.25 at `k=4,G=16` |
| `e(roth_slope_se_*)`, `e(roth_slope_ll_*)`, `e(roth_slope_ul_*)` and the printed slope line | x1.28 (`k=3`), x1.53 (`k=4`) on the se |
| `e(roth_power_*)` columns `mdt_wald` and `bias_wald` | x1.125 (`k=3`), x1.195 (`k=4`); Leg B alone is x1.023 at `k=1` |
| NEW: `e(roth_slope_dfci_*)`, hidden `e(roth_slope_seplugin_*)` | — |

**Does not move:** `e(roth_slope_*)`, `e(roth_W_*)`, `e(roth_hbar_*)`, `e(roth_npre_*)`,
`e(roth_keff_*)`, `e(roth_nissize_*)`, `mdt_nis`, `bias_nis`, and every number on every run
without `cluster()`. Nothing outside `pretrendpower()` moves.

**Why `mdt_nis` / `roth_nissize` are deliberately left alone**, since PTP-1's claim covers them
too (`mdt_nis(80)` achieving 0.7133, realised NIS size 0.1125 against a reported 0.0781). That
shortfall is a **different defect with a different fix**: the NIS leg evaluates a *normal*
probability at a *t* critical value (`ptpniscrit = invttail(df_r, .025)` handed to
`mvnormalcvqp`), which the original prober filed separately as their PTP-2 and which the pinned
build's own external check confirms — wooldid's `mdt_nis` sits at exactly
`invttail(279,.025)/invnormal(.975) = 1.0043` times Roth's `pretrends` answer on the same `b`/`V`.
The two consistent repairs are (a) use `invnormal(.975)`, making wooldid's NIS leg agree with
Roth's package to machine precision but abandoning the small-sample critical value the rest of
the table uses, or (b) compute a genuine multivariate-**t** rectangle probability, which Stata
has no native function for and which would need a Mata simulation inside a diagnostic that is
supposed to be cheap. Those are different tradeoffs from Leg A/B/C and belong to a separate
ruling, so the Hotelling gate is written not to touch them and `roth_nissize_*` comes back
bit-identical on every arm of `p1_beforeafter.log`. **Scoped out, not overlooked.**

## Blast radius — the PS-64 identity breaks, and the author must rule on it

PS-64 (2026-09-07/08) is explicit: "one statistic, one rule, in BOTH" entry points, and
`ado:1408` records that `wooldid_ptp` computes its joint pre-trend test "from the identical
analytic `e(V_es_*)`" as `jointtests`. **Measured** (`p4_ps64.log`, `jointtests` +
`pretrendpower(80)` on the same fit):

```
  PINNED    e(joint_horiz) = +1.aeece6076e5ccX-004   e(roth_prep_att) = +1.aeece6076e5ccX-004   BIT-IDENTICAL
  PATCHED   e(joint_horiz) = +1.aeece6076e5ccX-004   e(roth_prep_att) = +1.3000646296d86X-003   0.10521 vs 0.14844
```

`jointtests`' `joint_horiz` comes from `margins ..., contrast(overjoint)` (`ado:2268`), i.e.
Stata's own `F(q, df_r)` convention; `wooldid_ptp` forms the same quadratic form by hand.
Three ways to land this, in order of size:

* **(A) `pretrendpower()` only — what I patched and what I recommend.** The two entry points
  then differ, and the reason is statable in one sentence: `pretrendpower()`'s whole product is
  an *inverted design quantity* (`mdt_wald`, `bias_wald`, the slope interval), where a
  systematic one-directional error is the whole of the output; `jointtests` prints one p-value
  under Stata's own post-cluster convention, which no comparator package corrects either. PS-64
  would need amending to "one statistic, two references, and the reason is stated".
* **(B) Both.** Re-refer `joint_horiz`, `joint_horiz_i`, `joint_Choriz`, `joint_Choriz_i` and
  their `_s#` subgroup twins through the same scaling (the inputs `r(chi2)`/`r(F)`, `r(df)` are
  available right where `r(p)` is read). Coherent, and correct for exactly the same reason —
  but it moves **every** joint-test p on **every** clustered run, plus the `adversarial`
  statistic at `ado:14219`, which shares the convention by design. That is a large re-bless and
  a separate ruling.
* **(C) Nothing; document.** `sthlp:789`'s F-test sentence becomes conditional and
  `sthlp:791`'s "depends only on the design's precision" is corrected. Costs nothing, moves
  nothing, and leaves a printed "5% level" that is 7.8% (k=3) or 18.4% (k=6) at `G=16`.

I recommend **A**, with the help edits from C shipped alongside it whichever way A/B goes.

## Cost at runtime — free

Every added operation is scalar arithmetic on quantities already in hand: no data pass, no
matrix of larger dimension, no extra `mata` call. `p3_cost.log`: 10,000 repetitions of
`Ftail` + `npnF` + `invF` + `invttail` + `invchi2` take **0.35 s**, i.e. 35 microseconds per
evaluation set. `wooldid_ptp` evaluates them at most (4 families) x (number of requested power
levels) times, so the whole proposal costs well under a millisecond on any dataset of any size.

## Alternatives rejected

* **A caution keyed on `e(N_clust)`** (the prober's implied remedy, from their `G=80`
  "discriminating control"). Rejected on the refuter's evidence, which I re-derived: the error
  is a function of `(k, df_r)` and is 0.1034 at `k=4` with *the same* 16 clusters, while at
  `k=1` it is exactly 0.05 at every `G`. A cluster-count caution would fire where the test is
  exact and stay silent where it is worst. It also leaves the wrong number on the screen.
* **A cluster wild bootstrap for the joint pre-trend test.** Correct and strictly better, but
  `pretrendpower()` runs *after* the ereturn post with no estimation command available to it,
  and the measured cost (the prober priced a bootstrap-per-replication MC at ~19 hours) puts
  it out of scope for a column that is meant to be a cheap design diagnostic.
* **Reporting the MDT as an interval** rather than median-unbiasing it. Honest — its CV is 23%
  at `G=16` — but it adds two columns to a five-column table and a second ruling. Better served
  by the `sthlp:791` sentence saying the number is an estimate.
* **Leg B on the iid branch too.** Statistically right (`c2hat/c2true ~ nu/chi2(nu)` there, so
  the factor is `sqrt(nu/invchi2(nu,0.5))`), and it would be 1.0228 at `nu=15` and 1.0009 at
  `nu=279`. Rejected for the *shipping* version only because it puts a second gate and a second
  re-bless story into the same edit; recorded here so the author can take it if he wants the
  statistically complete version.
* **Leaving `e(roth_slope_se_*)` as the plug-in and adding an `_adj` twin.** Rejected: the
  replay renderer rebuilds the printed slope line from `e(roth_slope_se_*)`, so the shipped
  scalar has to be the one the interval is built from or live and replay disagree.

## The test that would prove it correct — four, three of which need no oracle

1. **The `k = 1` identity (no oracle, no simulation).** At one pre-period coefficient the
   Hotelling reference collapses to `F(1, df_r)` and `Q = 0` identically, so
   `roth_prep`, `roth_slope_se`, `roth_slope_ll/ul` and `roth_W` must be **bit-identical**
   before and after. Run. `p1_beforeafter.log` arm C: all seven scalars bit-identical.
2. **The gate identity (no oracle).** Any run without `cluster()` must be bit-identical
   throughout. Run. `p1_beforeafter.log` arms D (iid) and E (robust): every scalar and every
   matrix cell bit-identical.
3. **Coverage/size under the ideal model (closed form + simulation, no wooldid needed).**
   `b ~ N(0,Sigma)`, `S ~ Wishart(Sigma,nu)/nu` independent — the best case for any covariance
   estimator. The corrected joint test must have size 0.05 and the corrected slope interval
   coverage 0.95 at every `(k, nu)`; the shipped ones must not. Run, 200,000 draws,
   `a1_ideal.py`: size 0.0785 -> 0.0502 and coverage 0.9111 -> 0.9507 at `k=3,nu=15`;
   0.1031 -> 0.0495 and 0.8842 -> 0.9507 at `k=4,nu=15`; unchanged at `k=1`.
4. **The self-validating Monte Carlo (the decisive one; oracle = the design itself).** Take the
   column's own claim as the test: plant a linear pre-trend of exactly the size the patched
   `mdt_wald(80)` reports, on the geometry that fit it, and count rejections. It must come back
   at 0.80. Run, 250 replications on an independent seed stream, section (iv) above:
   **0.8040 CP95 [.7503,.8484]**, with the shipped column's own MDT coming back at 0.696. This
   test needs no distributional theory at all and would catch an error in any of the three legs.


---

# PROPOSAL 2 — the rank-deficiency caution does not name its cause; name it, for free
*covers:* **BX-N3**. One fix site: `wooldid_ptp`'s rank-deficiency branch, `ado:10505`, with the
same mandatory mirror in the replay renderer (`ado:11479`). Two optional legs follow it.

## What is actually wrong

With `i` the unit — `wooldid y unit t onset`, the default way the command is called — every
treated (`i`, relative-period) cell is a **single observation** carrying its own dummy, so that
observation's residual is **exactly zero by construction**. Both `robust` and any `cluster()`
sandwich form `sum_g X_g' e_g e_g' X_g`; the treated rows contribute **nothing at all**. The
delta-method covariance `margins` hands back is then built from control-cell variation alone.
Where two event-study horizons load on the control cells through weight vectors that are
proportional, the estimated covariance reports `corr = 1` and a variance ratio equal to the
square of that scale — although the point estimates, which *do* load on the treated rows, are
not in that ratio.

That is BX-N3's measurement exactly: `corr(pre7, pre4) = 1`, `Var(pre7)/Var(pre4) = 4.000`, with
`b_pre7 = -.2489` and `b_pre4 = -.0806`; and the bootstrap gives the direction the analytic
matrix calls a structural zero an ordinary variance (`v'S_boot v = 0.110`, 71% of the block's
mean bootstrap variance).

**My own reproduction** (`p2_bxn3.do` / `.log`, built from scratch: 30 states x 2 counties,
onset waves at t = 6 and t = 10, one never-treated wave, 14 periods, `robust`,
`esfixedbaseperiod esprelength(8) espostlength(3) pretrendpower(80)`):

| | `i` = unit (the default call) | discriminator: `i` = cohort, `fe(unit)` |
|---|---|---|
| max treated (`i`, relyr) cell | **1** | 20 |
| pre-block rank | **7 of 8** | 8 of 8 |
| smallest eigenvalue | **6.765e-17** | 3.11e-2 |
| largest off-diagonal \|corr\| | **0.9999999999999992** (rows 1, 5) | 0.631 |
| what wooldid says about the cause | **nothing** | (n/a) |

The pinned build prints its rank-deficiency caution and stops there. The user is told the
*consequence* ("a generalized inverse retains only the independent subset") and nothing about
the *cause*, and the one screen in the program that does name a cause (`ado:1413-1444`) keys on
`i` being nested in `cluster()` — which is false here, because there is no `cluster()` at all.

## Why the covariance cannot be repaired, only disclosed

Given residuals that are structurally zero, the sandwich is the **correct plug-in**; it is the
residuals that are uninformative, not the arithmetic. No HC variant rescues it: HC2 and HC3
divide by `1 - h_ii`, and `h_ii = 1` **exactly** on a singleton treated cell, so both are
undefined rather than corrected. The honest remedies are a coarser `i` with `fe()`, or
`bootstrap()`/`permute()`, and **wooldid already has all three**. So the whole deliverable here
is a sentence that tells the user which one they need — which is exactly what the finding asks
for ("either the caution should name singleton treated cells, or...").

## The edit

**Site**, `ado:10505`. **Now** (one line):

```stata
    if `ptpkeff' < `ptpk' dis as text "Caution: the `ptpty' pre-period variance-covariance block is rank-deficient (rank `ptpkeff' of `ptpk'). ... interpret with care."
```

**Becomes:**

```stata
    * BX-N3: the deficiency is almost always a SINGLETON-CELL design (sthlp's "Design
    * conditions for every consumer of e(V_es_*)", condition Two). Priced to be free on a
    * clean run: the gtools pass below is reached ONLY on a block that is ALREADY
    * rank-deficient, and only on the live analytic path, where the worker's tempvars are
    * in hand (the replay call at ado:449 passes hba()/hbi() and no tvar/ttrevar/ivar; the
    * bootstrap worker call exits at ado:10421 under hbaronly). Nothing here moves a number
    if `ptpkeff' < `ptpk' {
      local ptpsingmax .
      if "`hba'" == "" & "`hbi'" == "" & "`tvar'" != "" & "`ttrevar'" != "" & "`ivar'" != "" {
        tempvar ptpsr ptpsn
        cap qui gen double `ptpsr' = `tvar' - `ttrevar' if e(sample) & `ttrevar' != .
        cap qui gegen long `ptpsn' = count(`ptpsr'), by(`ivar' `ptpsr')
        cap qui sum `ptpsn' if `ptpsr' < ., meanonly
        if _rc == 0 local ptpsingmax = r(max)
        cap drop `ptpsr' `ptpsn'
      }
      dis as text "Caution: the `ptpty' pre-period variance-covariance block is rank-deficient (rank `ptpkeff' of `ptpk'). ... interpret with care."   /* UNCHANGED text */
      if `ptpsingmax' == 1 dis as text "  The usual cause is a singleton-cell design, and this run is one: every treated (i, relative-period) cell holds exactly one observation, so the treated residuals are identically zero and this covariance is built from control-cell variation alone - two pre-period horizons can then be reported as exactly proportional when their estimates are not. It is not caused by the clustering, and it is present under robust. A coarser i (cohorts, or treatment-onset groups) with unit fixed effects supplied through fe() removes it, and bootstrap() or permute() give inference that does not read this matrix. See the Design conditions for every consumer of e(V_es_*) in the help, condition Two."
      if `ptpsingmax' > 1 & `ptpsingmax' < . dis as text "  Treated (i, relative-period) cells here hold up to `ptpsingmax' observations, so the singleton-cell degeneracy the help describes under the Design conditions for every consumer of e(V_es_*) is not the cause; the deficiency is in the estimated covariance itself. bootstrap() or permute() give inference that does not read this matrix."
      cap ereturn hidden scalar roth_cellmax_`ptpty' = `ptpsingmax'
    }
```

**Replay mirror**, after `ado:11479` — the renderer has no data access, so it reads the scalar
the live run posted:

```stata
    if `ptpkeff' < `ptpk' & e(roth_cellmax_`ptpty') == 1 dis as text "  The usual cause is a singleton-cell design, and this run is one: ... condition Two."
```

**One refinement I recommend folding in** (`ptpsingdone`): the family loop runs over up to four
families, so a run with two rank-deficient blocks pays the pass twice. Hoisting it behind a
`if "`ptpsingdone'" == "" { ... local ptpsingdone 1 }` sentinel makes it once per run. The
version I actually **ran** computes it per rank-deficient family; the measured costs below are
for one pass.

## Cost at runtime — free on every clean run, and small on a dirty one

Three gates stand in front of the one data pass: the block must already be rank-deficient
(`ptpkeff < ptpk`), the run must be the live analytic one, and the user must have asked for
`pretrendpower()` at all. **Measured** (`p3_cost.log`, `p6_mixed.log`):

| | 1,000,000 rows | 5,000,000 rows |
|---|---|---|
| `gen` + `gegen count(), by(i relyr)` + `sum` | **0.171 s** | **0.495 s** |
| the same plus the two extra `count` passes of the Leg-1b variant below | — | 0.583 s |

and **zero** whenever the block is full rank — confirmed, not asserted: arms 2, 3 and 4 of
`p6_mixed.log` (`keff = 8 of 8`) post `e(roth_cellmax_att) = .`, so the branch was never entered.

## Does it move a reported number? **NO**

`p2_bxn3.log`, pinned vs patched on identical data: `keff = 7`, `eigmin = 6.765e-17`,
`eigmax = 1.0660415`, `max |corr| = .9999999999999992` — **character-identical** on both builds,
as are the estimate tables and the whole `Pre-Trend Power Diagnostics` block. The only `e()`
addition is a **hidden** scalar. Nothing needs a re-bless. (The Hotelling gate of Proposal 1 is
off on this arm — it is `robust`, not `cluster()` — so the two proposals were exercised
independently on the same patched build.)

## Blast radius

Confined to the rank-deficient branch of `wooldid_ptp` plus its replay mirror. The risks and
what was done about them: **(a)** a data command inside `wooldid_ptp`, which runs after the
`ereturn post` — handled by using only `gen`/`gegen`/`sum` on tempvars, every one under `cap`,
and dropping them, so a failure leaves `ptpsingmax` missing and no sentence prints;
**(b)** replay divergence — handled by the mirror and the posted scalar, **tested**
(`p2_bxn3.log`, `REPLAY_PATCHED_A` prints the identical two-sentence caution);
**(c)** the bootstrap worker call — it exits at `ado:10421` under `hbaronly` and never reaches
this code, so no bootstrap run pays anything.

## Alternatives rejected

* **Refuse `pretrendpower()` on singleton-cell designs.** Directly contradicts PS-64 (author
  ruling 2026-09-08: compute with a named caution rather than refuse, because "a hard error sent
  users into a documented dead end"). It would also strand the documented `honestdid` handoff,
  which consumes the same matrix with no option flag wooldid could screen on.
* **Repair the covariance.** There is nothing to repair it with: HC2/HC3 need `1 - h_ii` and
  `h_ii = 1` exactly. The repair that exists is `bootstrap()`, already shipped, already the
  sentence's recommendation.
* **Widen the `ado:1413` screen's text instead.** That screen's predicate is about clustering
  (`i` nested in `cluster()`); BX-N3 fires under plain `robust` with no `cluster()` at all, so
  widening the text cannot reach it without widening the predicate — which is Leg 3 below.
* **Say it unconditionally** ("the usual cause is a singleton-cell design"), with no data pass.
  Costs literally nothing and is the fallback if the author refuses the pass on principle. I
  rejected it because it would be *wrong text* on the runs where the cells are not singletons,
  and the whole complaint in BX-N3 is a caution that names the wrong cause.
* **Leg 1b — report the SHARE of treated observations in singleton cells, not just the max.**
  Written and run (`patched2/`, `p6_mixed.do`): it replaces `max == 1` with
  `share == 1 / 0 < share < 1 / share == 0` and prints the percentage. It costs two extra
  `count` passes (0.583 s vs 0.495 s at 5M rows). I do **not** recommend it as the default,
  because I could not construct the case it protects against: across four geometries
  (`i` = unit; `i` = cohort + `fe(unit)`; and a *mixed* `i`, unit-level for half the states and
  cohort-level for the other half, with and without `fe(unit)`), **the only rank-deficient block
  was the all-singleton one** — the other three came back `keff = 8 of 8` and never entered the
  branch. Recorded here because the tested version's second sentence ("...is not the cause")
  would be wrong if such a design does exist, and the author may prefer the hardened wording for
  that reason alone.

## Leg 2 — the bootstrap caution at `ado:6914` (text only, zero cost, recommended)

BX-N3's third consequence is that under `boottstat` every replication inherits the degeneracy, so
"the bootstrapped `joint_horiz` joint test is unavailable" fires at every `reps()`. That message
already names rank deficiency but not the design that causes it, and the diagnostic pass is not
available there (the estimation frame is not the live one). The zero-risk edit is a **pointer,
not a diagnosis**:

```stata
  local wdbc`wdbncaut' the bootstrapped `ename' joint test is unavailable (no usable replication statistic - rank-deficient replication covariances, or a horizon missing from the replications; a unit-level i is the common cause - see the Design conditions for every consumer of e(V_es_*) in the help)
```

No number moves, no data is touched. **Untested at runtime** (I did not run a bootstrap arm).

## Leg 3 — extend the `ado:1413` screen to singleton-cell designs (offered for ruling, UNTESTED)

The help's condition Two says that on a singleton-cell design the event-study standard errors
"look entirely ordinary and the object is not measuring what it appears to" **even when the block
is full rank**. Nothing says that at runtime. The screen at `ado:1413-1444` already fires under
exactly the right trigger — `jointtests` or `pretrendpower()`, i.e. the users who asked for a
number computed from this matrix — and already pays two `gegen tag` passes over the full sample
per cluster variable, so one more pass of the same order is in budget. The arm would be: compute
the max treated (`i`, relyr) cell size on `ifstatement`; if it is 1, print a caution naming
condition Two and the three remedies.

I am **not** recommending it as written, for two reasons the author should weigh rather than me:
`i` = unit is the *default* way the command is called, so this would newly caution a large share
of ordinary `pretrendpower()`/`jointtests` runs; and it pays the pass on every such run, not only
on the already-broken ones. It belongs in the same ruling as Proposal 1's option (A)/(B)/(C),
because it is the same question — how loudly should a design-level caveat be said.

## The test that would prove it correct — three run, one honest gap

1. **The discriminator (run, no oracle).** The same data, the same outcome, the same VCE, with
   `i` coarsened to the onset cohort and unit effects moved into `fe()`: the block must become
   full rank and the sentence must not print. `p2_bxn3.log` arm B: `keff = 8 of 8`,
   `eigmin = 3.11e-2`, max `|corr| = .631`, `e(roth_cellmax_att) = .`, nothing printed.
2. **The no-number-moves identity (run, no oracle).** Pinned vs patched on identical data:
   every reported scalar, every eigenvalue and every correlation character-identical.
   `p2_bxn3.log`, `SENTINEL_..._PINNED` vs `SENTINEL_..._PATCHED`.
3. **Replay (run, no oracle).** `wooldid` typed alone after the fit must print the identical
   two-sentence caution from `e()` alone. `p2_bxn3.log`, `REPLAY_PATCHED_A`.
4. **The cell count against an independent counter (max: run; share: NOT established).** The
   `bysort i relyr: gen _N` reference in `p2_bxn3.do` agrees with the ado's `gegen count()` on
   the maximum (1 and 20 on the two geometries). My reference program for the *share* had an
   `r(N)`-ordering bug (`sum, meanonly` overwrote `r(N)` before it was read), so the Leg-1b
   share is **not** independently confirmed — the ado's own ordering is correct, but I am not
   claiming a verified number for it.

---

# What I did not test, and what would settle it

* **Proposal 1 option (B)** — re-referring `jointtests`' own `joint_horiz` family through the
  same Hotelling scaling. Correct for the same reason; I did not patch it, because it moves every
  joint-test p on every clustered run and is a separate re-bless.
* **The `poisson` / `ppmlhdfe` branch.** `ptpdfok` is 0 on some poisson paths, so the reference
  falls back to chi-squared and the gate is off. Whether a clustered poisson run posts a usable
  `df_r` (and therefore gets the Hotelling reference) is not measured here; the prober filed a
  separate finding (their PTP-3) on that branch.
* **`bootstrap()` + `pretrendpower()`.** Deliberately excluded from the gate by
  `"`e(boot_method)'" == ""`, because the covariance there carries `e(boot_N_V)` replications and
  `nu = df_r` is the wrong idealization. What the right `nu` is there is open.
* **Leg 2 and Leg 3 of Proposal 2** — written, not run.
