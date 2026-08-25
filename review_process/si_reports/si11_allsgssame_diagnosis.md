# SI-11 diagnosis — the continuous SG-equality ("all SGs have the same ame") joint tests

**Date:** 2026-08-25. **Author directive:** extensive diagnosis incl. econometric causes; "do you
get this result even if clustering is turned off?" **Builds compared:** pre-Group-A
(`git show 4077c5c:wooldidinstall/wooldid.ado`, byte-identical chain to the f4aeb60 vintage used
for the T-M2 evidence bank) vs current HEAD (`wooldidinstall/wooldid.ado`, post-AMESE R1/R3/R4).
All probes in this directory (`w1_*`/`w2_*` do-files + logs, `instrument.py`, `instr_*.ado`,
datasets); nothing in the repo was modified.

## Verdict in one paragraph

SI-11 is **not a code defect**. The chain code (wooldid_sgchain -> `test` -> rank guard ->
contrast(overjoint) fallback, wooldid.ado L2874-2941) is correct and symmetric with the working
att twin, including the continuous-vs-binary `r(p)[1,1]` vs `[1,2]` reads (no T-M2-class whiff).
The missing values are the **family-A (G1) cluster-in-cell covariance degeneracy**: under
cohort-level `cluster(i)` with `contreat`, the posted margins covariance of the per-subgroup AMEs
is **exactly rank 1** (every entry of e(V) identical — measured: all entries 1.916e-07 on
`cont_sg`, all .03971533 on the rich panel), so the cross-SG equality contrast has *exactly zero*
estimated variance. What the user then sees is a floating-point accident of rank detection:
full drop -> "." (the SI-11 symptom), partial drop -> garbage p≈5e-80, no drop -> garbage
p≈1e-109. With clustering off (robust / no vce) or with non-contained clustering
(cluster(unit)) the tests come back with real, sane values on BOTH builds. At HEAD the
cluster(i) cases are intentionally skipped by the AMESE R3 gate with an explanatory caution —
which is the correct refusal of exactly this pathology. **Recommendation: close SI-11 as
diagnosed-and-already-fixed-by-Group-A**, with two small optional hardenings listed at the end.

## 1. Reproduction matrix

Data: `cont` = frozen battery `test/data/dta_cont.dta` + `sg=mod(i,2)+1` (2 SGs, att-only
default, exact battery `cont_sg` config with only the vce token swapped);
`rich2` = 15 cohorts (12 staggered + 3 never-treated) x 8 units x 16 periods, dose varying
within every cohort-period cell, sg=mod(i,3)+1, **post ame differs by sg in truth**
((0.3+0.3*sg)*dose) + sg-differentiated dose pre-trends;
`richo` = the T-M2 rich DGP verbatim (equal post ame 0.6*dose across sgs — post-side null TRUE).
Cells show `e(joint_ame_allsgssame)` / `e(joint_pre_ame_allsgssame)`; "absent" = scalar never
posted, "." = posted missing.

| config | pre-Group-A 4077c5c | current HEAD |
|---|---|---|
| cont x cluster(i) | **. / .** (posted missing — the SI-11 baseline symptom) | **absent / absent** + R1 containment caution |
| cont x robust | .3658 / .5741 | identical |
| cont x no vce | .4262 / .6455 | identical |
| rich2 x cluster(i) | 1.03e-109 / 1.86e-98 (**garbage**, F≈2.6e16) | absent / absent + caution |
| rich2 x robust | 6.3e-41 / .0046 (real; null false by design) | identical |
| rich2 x no vce | 8.0e-10 / .3752 | identical |
| rich2 x cluster(unit) | 2.2e-30 / .0032 (predicate does NOT fire) | identical |
| richo x cluster(i) | **. / 5.2e-80** (matches the banked tm2 rich_stock pins exactly) | (not run; predicate fires — absent) |
| richo x robust (negative control) | .6669 / .0046 — test correctly finds nothing post | identical |
| cont x unconditionalse cluster(i) | .1915 / .1268 (real, sane) | .1915 / .1268 — chain runs UNGATED (R1 never fires under unconditionalse, by design) and is healthy |

att-side twins posted real values in every cell of the matrix, on both builds.
The unconditionalse row is the one named carve-out where the R1 predicate does not fire yet the
chain still runs at HEAD: probed because A.6 lists it as an R1 non-fire. Empirically clean —
vce(unconditional) is a different variance construction (reg-based, covariate-sampling term)
and the G1 annihilation does not carry over on the battery config. (The rich2 x unconditionalse
leg refuses r(198) by wooldid's own validation — "unconditional standard errors may be
unreliable with non-att estimand types" — unrelated to SI-11. Whether vce(unconditional) is
VALID under containment is a family-A question, but it shows no SI-11 symptom.)
Identity check: the instrumented builds' e() values equal the stock runs' values exactly in
every compared cell (display-only instrumentation, T-M2 protocol).

## 2. The author's question, answered directly

**No — the missing values do NOT persist when clustering is turned off.** Under `robust`, under
no vce option at all, and under unit-level `cluster(uid)` (clusters not containing treated
cohort-period cells), the SG-equality ame joint tests return real, plausible p-values on both
the old and current builds — including a clean non-rejection (p=.667) on the negative-control
DGP whose post-period ame's are truly equal, and strong rejections (6e-41) where they truly
differ. The failure is specific to cohort-level clustering that triggers cell containment,
i.e. exactly the AMESE/family-A G1 predicate. It is an econometric degeneracy, not a code path
that loses the value.

## 3. Mechanism trace (instrumented pre-Group-A build, cluster(i) cells)

Chain sites: current L2874-2893 (post) / L2920-2938 (pre); pre-GA L2731-2749 / L2768-2787.
The trace markers (SI11A/B/F/C in `w2_prega.log`) show, per step:

1. **margins pwcompare(effects) post succeeds.** Posted e(b) is correct:
   `intensity:1bn.sg intensity:2.sg` (cont) / `dose:1bn.sg dose:2.sg dose:3.sg` (rich).
2. **wooldid_sgchain builds the right chain**: `[intensity]1bn.sg = [intensity]2.sg`
   (sgchainn=2), `[dose]1bn.sg = [dose]2.sg = [dose]3.sg` (sgchainn=3). No name mismatch —
   the T-M2 class is ruled out here.
3. **The posted e(V) is exactly rank 1** — the G1 signature. cont post-side:
   V = [[1.916e-07, 1.916e-07],[1.916e-07, 1.916e-07]]; richo pre-side: all nine entries
   .03971533. Every SG's ame carries the identical pooled-slope cluster variance and the
   identical covariance: the equality contrast's estimated variance is **exactly 0**.
4. **`test` reacts to the zero-variance constraint by dropping it**: r(drop)=1, and then
   - cont post & pre, richo post: ALL constraints dropped -> r(df)=0, r(F)=., r(p)=.
   - richo pre: 1 of 2 dropped -> r(df)=1, and the surviving noise-rank constraint yields
     F=2.37e12, p=5.2e-80.
   - rich2 (all four sites): the last-bit noise in V keeps numerical full rank -> r(df)=2,
     F=2.6e16, p=1e-109 — the guard PASSES and the garbage p posts directly.
5. **The rank guard** (`r(df) == sgchainn-1`, L2880/L2926) correctly fails on the dropped
   cases -> **fallback margins contrast(overjoint)** runs (L2889/L2935):
   - r(p) comes back **1x1** in every continuous case (confirming the per-shape index
     convention: the ame-side `[1,1]` read is correct; binary r(p) is 1x2, read `[1,2]`).
   - full-drop cases: margins' own joint contrast is rank-0 -> r(p)[1,1] = **"."** while the
     same r(table) shows the degenerate pairwise row (contrast .0023253, se 1.76e-10,
     t=1.3e7): the local is set to ".", ereturned as a missing scalar (L3671-3672), displayed
     as "p-value ...: ." (L5658-5659).
   - partial-drop case (richo pre): overjoint also reports the residual-rank p=5.2e-80,
     which posts and prints as if real.

So the "missing" e() scalar is the path: rank-1 V -> test drops all constraints -> fallback
overjoint p = ".". The equally-diseased siblings post 5.2e-80 / 1e-109 instead — same
invalidity, unluckier rank detection. On the current build under cluster(i), `amescrub`
(AMESE R3, L2854) empties `sgpost` before any of this: the chain never runs, the locals stay
empty, the scalars are absent, and the R1 caution ("...treated i-t cell(s) lie completely
within a single cluster, preventing valid analytic estimation of SEs for AME estimates...")
explains why.

## 4. Code-vs-econometric verdict per failure (four-bucket taxonomy)

| failure cell | bucket | detail |
|---|---|---|
| cluster(i) contained, pre-GA "." | **(iii) not estimable by the requested plug-in method — a correct refusal delivered silently and by accident** | The equality contrast of ame's genuinely has zero estimated variance under the G1 rank-1 V; no valid plug-in CRVE test exists on this V. The "." is the honest answer, but it arrived with no caution, and only when the numerics happened to expose the rank deficiency. |
| cluster(i) contained, pre-GA garbage p (5.2e-80 / 1e-109) | **(iii), presented as if (i) — silent wrong output, the worse sibling** | Same disease; rank detection failed (rich2) or half-failed (richo pre); an astronomically small p posts with no flag. This is the family-A "silent-wrong" class that motivated R1/R3. |
| cluster(i) contained, HEAD absent + caution | **correct refusal** (Group A working as designed) | Not a defect. In principle a valid test of the SG-equality null exists via non-plug-in routes (cluster bootstrap / IM-style dispersion = the skipped R2; impvalues is barred with subgroups at L3943), so "absent" is currently the ceiling of what wooldid offers. |
| robust / no vce / cluster(unit) | **estimable; values post on both builds** | The test exists and behaves (rejects a false null, accepts a true one). Caveat inherited from G5/NF-3 (not SI-11): small-cell leverage bias makes robust/unit-cluster joint tests oversized (SE ratios 0.68-0.97); that lives in SI-A3/R3's "elsewhere" leg. |

## 5. att-vs-ame asymmetry findings

Line-by-line, the two chains are symmetric: same sgpost gating, same sgchain builder (att adds
`half` to skip the 0b.targetflag base-outcome half of the posted e(b) — verified correct in the
trace), same rank guard, same fallback, and the per-shape r(p) reads ([1,2] binary / [1,1]
continuous) are BOTH correct. The divergence is purely in the posted covariance:

- att site, cont x cluster(i): e(V) = [[2.010e-06, 2.027e-06],[2.027e-06, 2.046e-06]] —
  near-rank-1 but not exactly (contrast variance ≈ 2e-9 > 0 in the stored bits) -> r(drop)=0,
  r(df)=1, p=5.09e-19 posts. att cells lack the absorbed-slope channel (G1a: understated, not
  annihilated), so the att equality tests always survive numerically — though under cluster(i)
  they carry the G1a ~0.66-0.75 SE understatement and are themselves suspect (open SI-A1;
  their F≈77000 on this config has the same degenerate flavor).
- ame site: the within-cell channel is annihilated bit-exactly -> exact rank 1 -> drops/garbage.

So the att twin "working" is not evidence the ame chain is coded differently — it is the same
code fed a non-degenerate (merely biased) V.

## 6. Suggested fixes with blast radius (no changes made)

**F0 (recommended): close SI-11 as diagnosed; no estimator change.** The pre-GA missingness is
family-A G1; the Group-A R1/R3 gate is the fix and is already in place, with the caution and
the standing-issues A-section as documentation. Optionally add one sentence to the sthlp
SG-equality/jointtests entries: under cohort-level clustering with contreat these tests are
withheld (containment caution names the alternatives). Blast radius: none (doc only).

**F1 (optional hardening): stop the fallback from posting residual-rank p's.** At the four
chain sites (wooldid.ado L2879-2891, L2925-2938 and binary twins L2634-2647, L2678-2692):
when `r(drop)==1` after `test `sgchain'`, the current code falls back to contrast(overjoint)
and posts whatever it returns — which on a partially-detected degenerate V is a residual-rank
garbage p (measured 5.2e-80). Change: treat `r(drop)==1` as "not estimable" (leave the local
empty + one-time caution) and reserve the overjoint fallback for genuine test failures (r(111)
etc.). NOTE the E1 design intent was bit-compatibility with the pre-E1 contrast engine (which
posted the same garbage), so this is an author-ruling behavior change. Blast radius at HEAD:
effectively nil on the battery (the reachable degenerate configs are already amescrub-gated;
no current baseline pins these scalars — `grep ame_allsgssame test/baseline/*.csv` is empty),
but it protects any future non-containment degeneracy. The one HEAD config that runs the
chain under cluster(i) — the unconditionalse carve-out — tested clean (see matrix). Do NOT
extend to a rank/eigenvalue guard on e(V) itself — family-A ruled that out (G2 spoofed full
rank).

**F2 (optional tidy): posted-missing "." vs absent inconsistency + one uncapped read.**
(a) At the ereturn site (L3671-3672 and att twin L3665-3666) a local equal to "." posts as a
missing scalar and prints as "p-value ...: ." with no explanation (the pre-GA symptom; still
reachable at HEAD in principle via a fallback "." on a non-gated config). Suppressing
"."-valued locals (`!= "."` guard) plus a caution unifies the two missing modes. No current
baseline pins a "." there (the tm2 evidence bank pins are pre-R3 and deliberately frozen).
(b) L2890: `local joint_C`item'_sdiff = r(p)[1,1]` after the post-side fallback is uncapped —
if margins ever returns no r(p) this aborts wooldid with a raw error; the pre twin (L2936)
caps the same read. Symmetrize with `cap` + empty-on-failure. Latent robustness only; no
observed trigger. Blast radius of F2: none on current baselines.

**Explicitly not recommended:** any attempt to "make the number appear" under cluster(i)
containment with the plug-in V — the only honest routes are the skipped R2 replacement
estimator (cluster bootstrap / IM-style dispersion of per-cohort ame's) or finer clustering,
exactly as the R1 caution already tells the user.

## Probe inventory (this directory)

- `mkdata.do` -> `rich2.dta`, `rich_orig.dta`; `dta_cont.dta` copied from test/data.
- `w1_core.do` + `w1_prega.do`/`w1_current.do` -> `w1_prega.log`, `w1_current.log`
  (stock 8-config VCE matrix, both builds; SI11ERET bracket dumps).
- `instrument.py` -> `instr_prega.ado`, `instr_current.ado` (display-only SI11A/B/F/C dumps
  at the four chain sites; placement verified by diff; e() identity vs stock verified).
- `w2_core.do` + `w2_prega.do`/`w2_current.do` -> `w2_prega.log` (the mechanism trace:
  e(b)/e(V)/chain/r(df)/r(drop)/fallback r(p) at every step), `w2_current.log` (HEAD:
  ame chain skipped under cluster(i) with sgpost emptied; healthy full traces under robust).
- `w3_uncond.do` -> `w3_uncond.log` (the unconditionalse carve-out probe, both builds).
- `old_prega.ado` = `git show 4077c5c:wooldidinstall/wooldid.ado`.
