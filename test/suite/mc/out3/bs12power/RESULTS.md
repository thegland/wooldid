# AB12/bt power — the bootstrapped adversarial statistic under a real violation

Author ruling 2026-09-14: "measure power first, then decide." The size ladder had found
`AB12/bt` rejecting **0 of 200** under the null at both 5% and 10%, where the analytic arm on
the same geometry rejects 72.5%. That establishes only that the bootstrap does not reject under
the null — a conservative test and an inert test are indistinguishable there. This measures the
other side.

## Design

`mc_adv_bs12_power.do`. AB12's geometry, VCE and options exactly — `cluster(i)` with
`safety(off)`, the containment rung — so the ONLY difference from the size campaign is the plant.

The plant is a **step**: effect 0 in every pre period and TAU in every post period. This matters.
The size ladder's own plant is `SLOPE*(relyr+1)`, which is zero at the reference period, so it is
a line THROUGH the reference period and BS-12's one-parameter null contains every such line —
raising SLOPE there is a second null, not power. A step is not a line through the reference
period unless TAU is 0, so it is a genuine violation. TAU 0 reproduces the null arm and was run
first as the harness control.

40 reps, B = 199, TAU = 0.5, 4 shards. Both arms recorded per rep on the same dataset, so the
comparison is paired.

## Result

    rc = 0 on all 80 rows; advabsent 0; ncomp = 199 and nfail = 0 on every bt row.

| arm | n | reject @ .05 | reject @ .10 | min p | median p | max p |
|---|---|---|---|---|---|---|
| `an` analytic | 40 | **35/40 = 0.875** | 37/40 = 0.925 | 1.71e-09 | 0.0033 | 0.231 |
| `bt` bootstrapped | 40 | **0/40 = 0.000** | 1/40 = 0.025 | **0.080** | 0.560 | 0.865 |

**Paired over the same 40 datasets: the analytic arm rejects on 35, the bootstrapped arm rejects
on 0, and the two never agree on a rejection.** On those 35 the analytic p reaches 1.7e-09.

## Reading

The bootstrapped adversarial statistic on a containment geometry is **inert, not conservative**.
It does not reject under the null (0/200, which looked like a repair of the analytic arm's 72.5%
over-rejection) and it does not reject under a violation the analytic arm detects at p = 1e-09.
A test that rejects nothing is not a corrected test; it carries no information about the
hypothesis at all, and a user who sees a non-significant `AB12/bt` p learns nothing.

**Scope.** One geometry (containment, `cluster(i)`), one plant shape, one TAU, 40 reps. The
conclusion needs no more than this: 0 rejections against a paired 35 is not a power level that
becomes adequate with more draws. What is NOT established is whether the same inertness affects
the bootstrapped adversarial p on non-containment geometries, which the size ladder does not
cover because `bt` runs only where `cluster()` nests `i`.
