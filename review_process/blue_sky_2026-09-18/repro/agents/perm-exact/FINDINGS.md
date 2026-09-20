# perm-exact — findings
Agent: `permexact` against an independent brute-force enumeration (EXACT ORACLE method).
Build: `/private/tmp/.../scratchpad/bluesky/build/wooldid.ado`, asserted by `which wooldid`
in every probe log ("version 1.4 2026-08-19").

---

## FINDING PX-1 — `permstepdown`'s p_raw / p_singlestep / p_stepdown are documented and
## printed as ADD-ONE shares, but under `permexact` the kernel computes plain exact shares

**Seeded?** No — found by me.

**Claim.** Under a complete `permexact` enumeration the step-down kernel
(`wooldid_permstep`, ado:16290-16320) uses the plain share `nex/n`, but the runtime footer
(ado:11700) and the help (sthlp:514) both describe that row unconditionally as the add-one
share. The two differ by a factor of about two at the floor: `1/M` versus `1/(M+1)` — and
with one exceedance, `1/M` versus `2/(M+1)`.

**Code.** ado:16302-16320, three times:
```
out[2, h] = (exk | memb ? nex/n : (1 + nex)/(n + 1))   // p_raw
out[3, h] = (exk | memb ? nex/n : (1 + nex)/(n + 1))   // p_singlestep
pj        = (exk | memb ? nex/n : (1 + nex)/(n + 1))   // p_stepdown
```
`exk = (M < . & n == M)` is exactly "this is a complete enumeration".

**Documentation.**
- ado:11700 (the footer actually printed by a `permexact` run):
  `"... p_raw is the two-sided |t| share (add-one) in the step-down's own metric ..."` —
  no exactness branch, unlike its sibling at ado:11681 which *is* exactness-aware
  (`if "e(perm_exact)" == "1" ... "the p-values are exact shares with no Monte-Carlo error"`).
- sthlp:514: `"The stored kernel's p_raw row (and the footer's first column) is the
  two-sided |t| share, add-one, in the step-down's own studentized metric"`.

**Observed** (`p10_stepdown.do`, M = 30 enumerated assignments, 30/30 completed):
```
e(perm_stepdown)[4,3]
                   ES_att_post0       ES_att_post1       ES_att_post2
       t_obs  2.193871675920451  2.446479886586404  2.835961821942529
       p_raw  0.033333333333333  0.033333333333333  0.033333333333333
p_singlestep  0.033333333333333  0.033333333333333  0.033333333333333
  p_stepdown  0.033333333333333  0.033333333333333  0.033333333333333
```
0.0333333... = 1/30 = the plain share. The add-one value the footer and help promise is
(1 + 1)/(30 + 1) = 0.0645161. The same run prints
`"... p_raw is the two-sided |t| share (add-one) in the step-down's own metric ..."`.

**Expected.** Either the label is exactness-aware (as the coefficient-p footer at ado:11681
already is, and as sthlp:503 and sthlp:452 already are for the equality and joint tests:
*"with the same add-one form for sampled draws and plain shares under permexact"*), or the
kernel uses add-one. The computed number is the right one for an enumeration — the package's
own general statement is at sthlp:654, *"permexact needs no add-one while sampled draws do"* —
so the defect is in the two descriptions, not in the arithmetic.


**The 0.05 straddle.** This is not only a resolution-reading issue. In the repro, with one
exceedance, the plain share is 1/30 = 0.0333 and the add-one share is 2/31 = 0.0645 — on
opposite sides of 0.05. A reader who takes the stated convention at face value reads a
family-wise-adjusted p of 0.0333 as a number that convention cannot produce. Every
step-down number the run prints is a multiple of 1/30 and none of 1/31: the pre-treatment
family gives p_raw 0.900 = 27/30 and 0.400 = 12/30, p_singlestep 0.966666… = 29/30 and
0.700 = 21/30 (`e(perm_stepdown_pre)` in the same log).

**Why it matters.** The stated convention is what a reader uses to judge resolution. A user
told the row is add-one concludes its floor is 1/(n+1) and that it cannot equal the table's
convention; under `permexact` both conclusions are wrong. It also makes the sentence
"compare it with the step-down p rather than with the table" misleading: under `permexact`
the step-down's p_raw and the table's `permabsp`/`permrawp` p are the same share.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact/p10_stepdown.do`
(self-contained: builds its own 6-unit panel, prints 1/M and 2/(M+1) beside the matrix).
Run: `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p10_stepdown.do`

`realistic` true · `severity` low · `confidence` certain · `surface` permstepdown + permexact


---

## FINDING PX-2 — the `permexact` add-one caution describes a branch the code does not take:
## at complete-enumeration-but-not-a-member the share is plain, and the p can be 0

**Seeded?** No — found by me. **Reachability from the command line: UNPROVEN (see below).**

**Claim.** `wooldid_permtail` selects the plain share on `exact | memb` (ado:15992), where
`exact` is COMPLETENESS of the enumeration; but `wooldid_permstat` flags an element into
`wdp_addone` on `!memb` ALONE (ado:16147), and the caution the caller then prints (ado:8297)
asserts the add-one form is in use:

> "Caution: under permexact the observed assignment's own replication did not carry the
>  following elements, so their reference set does not contain the observed statistic and
>  their p-values use the add-one form (1 + count)/(n + 1) rather than an exact share: ..."

For an element with `exact = 1` and `memb = 0` that sentence is false in both halves that
matter: the code takes the plain share, and it takes it over a reference set that (by the
caution's own admission) does not contain the observed statistic — so the add-one floor of
1/(n+1) is not applied and the p can be exactly 0. Every sibling kernel has the same shape
(`wooldid_permqf` ado:16221-16249, `wooldid_permstep` ado:16293-16320, `wooldid_permpw`).

**Observed** (`p12_permtail_unit.do`, the three kernels extracted verbatim from the pinned
build and called directly; 30 complete draws, observed statistic one ulp below all of them):
```
PERMTAIL exact=1 memb=1 : p = 0                       (unreachable: memb=1 means sobs is a draw)
PERMTAIL exact=1 memb=0 : p = 0          <-- the element wdp_addone names
PERMTAIL exact=0 memb=0 : p = .064516129 <-- what the caution PROMISES, 2/(n+1)
p_lo / p_hi at exact=1 memb=0 : 0  1         (at exact=0 memb=0 : .032258065  1)
permside(left) exact=1 memb=0 : 0            (add-one branch would give .032258065)
Monte-Carlo se of the share at exact=1 memb=0 : .   (left missing, i.e. treated as exact)
```

**Expected.** Either the share branch keys on membership alone (`memb` rather than
`exact | memb`) for a flagged element, or the caution says what actually happens. As it
stands the run prints a message promising a floor of 1/(n+1) while reporting 0 — and
`known_bydesign` §2 ranks a p of 0 on a non-`permrawp` permutation p as reportable, which is
precisely the value this branch produces.

**Reachability.** I could NOT drive wooldid into `exact = 1, memb = 0` from the command
line. `memb` was TRUE for every one of 54 + 30 + 16 + 1 permuted elements across six
designs (p01/p04/p05/p06), because the rep-0 comparator and the identity replication are
the same fit, and the two estimation paths the review brief names as not bit-reproducible
run-to-run (`unconditionalse`, `customvce(cluster ...)`) are REFUSED by `permute()`
(p07b_edges.do, E10-E12, r(198)). So this is a demonstrated contradiction between a branch
and the message that documents it, of unknown reachability — filed `uncertain` for that
reason, not because the branch's behaviour is in doubt.

**Repro.** `/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact/p12_permtail_unit.do`
(self-contained: extracts the kernels itself, prints all three branches side by side).
Run: `"/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do p12_permtail_unit.do`

`realistic` false (no command-line path found) · `severity` low · `confidence` uncertain
· `surface` permexact + wooldid_permtail / wdp_addone caution

---

## NOT findings — what the oracle certified CLEAN (detail in NOTES.md)

The seeded lead (M and A produced by two unreconciled algorithms) does **not** reproduce on
any reachable design, and every documented p-value convention matched an independent
brute-force enumeration bitwise. See NOTES.md for the full list and for the one latent
inconsistency (`wdp_addone`'s caution text versus permtail's `exact | memb` branch) that I
could not reach from the command line.
