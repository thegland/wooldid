# LENS 0 (independent reproduction) -- FCM-1: VERDICT CONFIRMED, severity high

Build asserted: .../scratchpad/bluesky/build/wooldid.ado, *! version 1.4 2026-08-19.
Repro: L0b_fcm1.do / L0b_fcm1.log  (first attempt L0_fcm1.do died on MY OWN bad line
`local mm : all macros`, r(198), after arm A; arm A's numbers agree to every digit).

## The four arms (reps 250, rseed 31415, cluster(distr), att att_i att_it)
100 schools in 20 districts, T=10, i=school, cluster(district) COARSER than i.
The district label is attached AFTER y is generated, so it reaches the command ONLY
through cluster(): no relabelling can move a point estimate. Confirmed: b(att) is
+1.e2240256c1c85X-001 and b(att_i) is +1.de30ff7f19451X-001 in ALL FOUR arms, and
N=1000, N_clust=20, completed=250, failed=0 in all four.

arm                         PREDNEST str  se(att)    se(att_it)  se(att_i)
A  baseline                     1     0   0.182595   0.182595    0.182402
B  1 school moves mid-panel     0     1   0.184071   0.142393    0.141058
D  10 schools move mid-panel    0    10   0.179999   0.140252    0.137478
C  the SAME 10 schools and the
   SAME targets, but for ALL t  1     0   0.173922   0.173922    0.172964

## What each row buys
1. ONE straddler out of 100 cuts se(att_i) by 22.7% and se(att_it) by 22.0%.
   se(att) -- the observation-weighted estimand, which does not take this leg -- goes
   UP 0.8%, which is the sensible response to moving rows between clusters.
2. The next NINE straddlers buy 2.5% more (D/B = 0.975). The rule is a cliff, not a dose.
3. THE CONTROL NOBODY RAN: arm C moves the SAME ten schools to the SAME ten target
   districts, differing from D only in whether the move is recorded for the whole panel
   or from 2015 on. C keeps PREDNEST=1 and se(att_i) stays at 0.9483 of baseline; D flips
   PREDNEST to 0 and falls to 0.7537 of baseline. D/C = 0.795 on identical data movements.
   And arm C is the LARGEST perturbation of the four by the analytic clustered SE
   (se_ana(att_i): A 0.044869, B 0.048985, D 0.047367, C 0.049376 -- C is +10.0%, the
   biggest). The arm that moves the cluster structure most moves the bootstrap SE least.
   The driver is the global flag, not the data change. This is not open to another reading.
4. In arm B the printed table carries att and att_it with BIT-IDENTICAL point estimates
   (+1.e2240256c1c85X-001) and standard errors 0.184071 and 0.142393 -- a 22.6% gap
   between two rows of one table that are the same number, with nothing on screen to
   explain it. Arm A prints them identical to the last digit.

## Which arm is calibrated -- from my own data, no oracle needed
SD over the 20 districts of the district mean treatment effect = 0.894869, so the
between-district component of a pairs-cluster SE is 0.894869/sqrt(20) = 0.2001.
Nesting arm 0.1824 sits just under it; keep-i arm 0.1411 is 29% below it. The nesting
arm is the calibrated one and keep-i is anti-conservative, exactly as claimed.
Structural reason, which needs no calibration: with cluster COARSER than i every cell
lies wholly inside one cluster, so duplicating a cluster cannot move any cell mean;
under keep-i the duplicate cells are identical and the equal-cell weighting counts them
once, so a drawn cluster's only effect is presence/absence. Variance can only be removed.

## The claim's ARITHMETIC IS WRONG (skeptic's correction 1 independently confirmed)
ratio_oracle.py simulates the two replication statistics directly, no Stata, 400k draws:
  NEST  replicate = sum(M_c th_c)/sum(M_c), M ~ Multinomial(G,1/G)  -> Var = s^2/G
  KEEP  replicate = sum(I_c th_c)/sum(I_c), I_c = 1{M_c>=1}, p=1-1/e
        = mean of a random SUBSET of ~pG of the G clusters -> Var = s^2(1-p)/(pG)
  SE ratio -> sqrt((1-p)/p) = 0.7629;  simulated G=20 0.7803, G=24 0.7767, G=200 0.7639.
My measured B/A = 0.7733 and D/C = 0.7948 sit on that curve at G=20.
The claim's sqrt(Var(Bernoulli(1-1/e))/Var(Poisson(1))) = sqrt(0.2325) = 0.4822 drops the
renormalisation by the realized number of included clusters and matches nothing measured
by anybody (0.482 would imply a 34.5% size, not the ~16% the claim itself quotes; the
claim computed its 16% from its own MEASURED 0.713, not from its stated 0.482).
Correct implied size of a nominal 5% two-sided test at my measured 0.773: 13.0%.

## "Silently" -- confirmed
Zero "Caution" lines in any of the four bootstrap logs. No `ereturn local/scalar/matrix`
anywhere in the pinned ado has a name containing nest, reid or regime. The only user-
visible trace of `wdb_nest` in the whole program is the A10 caution at ado:406-412, and
it is gated on `wdb_nest == 1` -- it fires only in the regime that is NOT the problem.

## Classification -- I side with the skeptic, with one addition
sthlp:532 says the keep-i rule applies "where the cluster is finer than i, OR CUTS ACROSS
IT", and one straddler does make the two partitions cross-cutting, so the code does what
the help's words say. It is a design defect, not a doc-vs-code divergence. The addition:
the help states the rule as a property of the design ("where the cluster nests i"), never
says the test is all-or-nothing over units, and never warns that one unit out of a
hundred decides it for all hundred -- so a user cannot predict which rule their
99%-nesting panel will get. That is a real documentation shortfall on top of the defect.
Worth noting too that the code's own comment at ado:8980-8986 names this exact harm --
a coarser cluster "would silently take the keep-i rule and collapse att_i's composition
variance" -- as the thing the `if smp` placement exists to prevent. The author saw the
failure mode and guarded the sample-concept route to it, not the one-straddler route.

## NOT measured
att_itime (the claim's third estimand). It is a cell-weighted estimand and takes the same
physical leg per sthlp:532, and att_it behaved identically to att_i here, so I expect the
same; I did not run it. Stated as inference, not measurement.
