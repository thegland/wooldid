# LENS 2 (realism & impact) — tiebreaker on fp-contreat-B1

## What is NOT in dispute (I accept it, both prior agents measured it independently)
The engine split is real: on contreat()+lattice() panels the wooldid DEFAULT (reghdfe)
path returns max|e(V)| = 0 while plain regress on the identical fit (identical N, rss,
point estimates) returns a finite robust covariance. Prober 4+ designs; skeptic 6/6 on
an independent fixture. I am not re-litigating existence. My job is severity and remedy.

## The question my lens actually has to answer
1. How likely is a real applied economist to hit it?
2. What EXACTLY do they lose — is anything valid lost, or only a blank?
3. Would they notice? (loud vs silent)
4. Is the prober's proposed remedy ("adopt the regress covariance, wooldid already
   knows how to drive regress") SAFE?

## The thing neither prior agent connected  (my central observation)
The skeptic's OWN B2 Monte Carlo, run on the EXACT design where B1 fires
(R7 part b, 12 draws, lattice(all 3), customvce(robust) = the regress engine):
    ame mean -0.342, MC SD 2.519, range [-6.091, +2.906];  TRUTH 0.400
    reported se: mean 0.341, range [0.305, 0.434]      -> SE is 7.4x TOO SMALL
So on the ame row, the covariance that B1 wants wooldid to adopt is NOT a valid standard
error. B1's remedy delivers B2's defect. Under this lens (silent wrong number >> loud
failure) that inverts the sign of the recommendation on that row.
The open question is whether the SAME holds for the att row, which is cell-intercept
driven and may be perfectly well identified — that is the casualty worth rescuing.

## Probes
- L1_impact.do  : 12-draw MC, 3 arms on the same draw, att AND ame tracked separately.
    A = lattice(all 3)  (dose-quantile lattice, prober/skeptic design)
    B = lattice(reg3)   (substantive lattice INDEPENDENT of dose, wide within-bin dose
                         variation — the most natural applied use of lattice())
    N = no lattice      (calibration baseline; note the skeptic's no-lattice SE/SD was
                         0.75, not 1.0, because i=state leaves the county FE unabsorbed,
                         so compare A and B against N, not against 1.0)
  Records SE/SD and 95% coverage per row per arm, plus zero-V incidence on the default
  path for A and B.
- L2_side.do    : what the user SEES (caution text), what e(V) CONTAINS (zeros would mean
  a downstream test/lincom silently reports t=inf/p=0; missings mean loud), a downstream
  test/lincom, boottstat on the zero-V spec (studentized t = (b_r-b)/0 — number or
  refusal?), and incidence across realistic specs incl. cluster(cty) (FINER than i, so
  the documented R1 containment predicate cannot fire).

## Results — L1_impact.do (12 draws, my own fixture, my own DGP truth)

SE vs realized Monte Carlo sampling SD, and 95% coverage. All three arms use
customvce(robust) (the regress engine), i.e. exactly the covariance fp-contreat-B1
proposes wooldid should adopt on the default path.

| row | arm                  | mean    | MC SD   | mean SE | SE/SD | cover95 |
|-----|----------------------|---------|---------|---------|-------|---------|
| att | A lattice(all 3)     | 0.85649 | 0.12465 | 0.16695 | 1.339 | 12/12   |
| att | B lattice(reg3)      | 0.85063 | 0.12468 | 0.16565 | 1.329 | 12/12   |
| att | N no lattice         | 0.85362 | 0.12396 | 0.15974 | 1.289 | 12/12   |
| ame | A lattice(all 3)     | 0.75096 | 1.18272 | 0.32450 | 0.274 |  5/12   |
| ame | B lattice(reg3)      | 0.63188 | 0.55366 | 0.15573 | 0.281 |  5/12   |
| ame | N no lattice         | 0.40053 | 0.13291 | 0.11370 | 0.855 | 11/12   |

true att (per draw) mean 0.89852 sd 0.03454 ;  true ame = 0.400 exactly.

zero-V incidence on the wooldid DEFAULT (reghdfe) path:
    A lattice(all 3)  12/12   (deterministic, matches prober + skeptic)
    B lattice(reg3)    8/12   (INTERMITTENT — same specification, different sample,
                               SEs appear on 4 draws and vanish on 8)

### What this settles — the two rows are NOT the same casualty
- **att row: a genuine loss.** Valid in every arm — unbiased-ish, SE/SD ~1.3
  (conservative), 12/12 coverage, and the lattice does not degrade it at all
  (0.12465 / 0.12468 / 0.12396 MC SD across the three arms — identical). The default
  path throws this away. This is the real user cost.
- **ame row: NOT a loss — the blank is protective.** Under ANY lattice (dose-quantile
  A or substantive independent B) the regress SE is 3.6x too small and the nominal 95%
  interval covers 42%. Without a lattice the same engine is fine (SE/SD 0.855, 11/12).
  So the covariance B1 wants adopted is, on the ame row, a wrong number. Adopting it
  wholesale would convert a loud blank into a silent 3.6x-understated SE.
- Calibration check: the no-lattice arm's att SE/SD is 1.289, so the lattice arms' 1.33
  is not a lattice artifact; and the ame collapse from 0.855 to 0.27-0.28 IS.

### Consequence for the remedy
The safe remedy is the NARROW one: reinstate the att-family rescue (the AMESE R4 refit
removed 2026-09-08, ado:1546-1565) and keep the ame side blank with its caution. The
premise recorded for that removal — "cell-level containment has already made the ame
inference invalid" — is factually wrong here (no cluster() at all, so no containment),
which is the prober's point and my numbers support it; but the CONCLUSION it reached
about the ame happens to be right for a different reason (lattice fragility, 42%
coverage). So: prober right about the att, wrong to generalize it to the ame.


## Final state (all probes complete)
L1 (12-draw MC, 3 arms), L2 (cautions, e(V), test/lincom, boottstat, incidence S5/S6),
L4 (isolated bootstrap legs), L5 (14 draws: what the default path PRINTS when it does not
zero out, + the bootstrap remedy). L3 was killed as redundant after L2 completed S5/S6.
VERDICT: PARTIAL, severity medium (down from high), realistic true.
The single most important number I added: on the 5/14 draws where the DEFAULT engine does
NOT zero out under lattice(reg3), it prints an ame SE of 0.113 against a realized sampling
SD of 1.741 and covers 2/5. The zero-V blank is the only thing preventing that number from
being printed on the other 9. B1's proposed remedy would print it on all 14.
