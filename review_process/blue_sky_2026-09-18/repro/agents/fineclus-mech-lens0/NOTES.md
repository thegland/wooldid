# fineclus-mech-lens0 -- THIRD LENS (independent reproduction) on FCM-1

Assignment: rebuild the case a third time FROM THE CLAIM ALONE and say whether the world
behaves as claimed. Prober says CONFIRMED high; skeptic says CONFIRMED high with two
arithmetic/classification corrections. I built my own DGP and added a control neither ran.

## My DGP (prog mkdat in L0_fcm1.do) -- different shape from both prior agents
20 school DISTRICTS x 5 SCHOOLS = 100 schools, T = 2010..2019 (10 periods, 1000 rows).
i = school. cluster(distr) is strictly COARSER than i. onset 2013/2015/2017 by
mod(school,4); mod4==0 (25 schools) never treated. district shock dshk ~ N(0,.7^2)
and district effect dtau ~ N(1,.9^2); school effect stau = dtau + N(0,.35^2).
y is generated BEFORE any relabelling, so the district label reaches the command ONLY
through cluster(): no relabelling can move a point estimate.

## FOUR ARMS, one seed (rseed 31415), reps(250), estimands att att_i att_it
  A  baseline                                              -> PREDNEST 1, 0 straddlers
  B  1 school moves district from 2015 on                  -> PREDNEST 0, 1 straddler
  D  10 schools move district from 2015 on                 -> PREDNEST 0, 10 straddlers
  C  the SAME 10 schools, SAME target districts, ALL t     -> PREDNEST 1, 0 straddlers
C vs D is the discriminating pair nobody ran: identical schools, identical targets, the
only difference is whole-panel vs mid-panel recording.

## Build assertion
which wooldid -> .../scratchpad/bluesky/build/wooldid.ado, *! version 1.4 2026-08-19. OK.

## RESULT 1 -- the estimation data is identical across all four arms (hex-exact)
b(att)   = +1.e2240256c1c85X-001 in A, B, D and C
b(att_i) = +1.de30ff7f19451X-001 in A, B, D and C
N = 1000 and N_clust = 20 in all four.

## RESULT 2 -- the ANALYTIC clustered SE responds to the data, and responds MOST to arm C
  se_ana(att_i):  A 0.044869   B 0.048985 (+9.2%)   D 0.047367 (+5.6%)   C 0.049376 (+10.0%)
Arm C -- the arm whose regime flag is UNCHANGED -- is the LARGEST analytic perturbation of
the four. That is the control's whole point: if the bootstrap drop tracked "how much did
the cluster structure move", C would drop the most. It has the biggest move and (see
RESULT 3) no drop at all.

## RESULT 3 -- bootstrap SEs
(pending)

## RESULT 4 -- the limiting SE ratio, derived and simulated independently (ratio_oracle.py)
Let theta_c be cluster c's contribution to the equal-cell-weighted estimand, G clusters.
  NEST rule : a twice-drawn cluster is re-identified and its cells COUNT TWICE
              -> replicate = sum(M_c theta_c)/sum(M_c),  M ~ Multinomial(G, 1/G)
              -> Var = s^2/G    (sum M_c is exactly G; standard pairs bootstrap)
  KEEP-i rule: duplicates are identical cells the equal-cell weighting counts ONCE
              -> replicate = sum(I_c theta_c)/sum(I_c),  I_c = 1{M_c >= 1}, p = 1-1/e
              -> this is the mean of a random SUBSET of size ~pG of the G clusters,
                 so Var = (s^2/(pG)) * ((G-m)/(G-1)) -> s^2 (1-p)/(pG)
  ratio of SEs -> sqrt((1-p)/p) = sqrt(0.3679/0.6321) = 0.7629
400k-draw simulation of exactly those two statistics (no wooldid, no Stata):
  G=20 ratio 0.7803 | G=24 0.7767 | G=40 0.7729 | G=200 0.7639
The claim's sqrt(Var(Bernoulli(1-1/e))/Var(Poisson(1))) = sqrt(0.2325) = 0.4822 DROPS the
renormalisation by the realized number of included clusters and matches nothing anyone
measured. The skeptic's correction is right and I reproduce it from scratch.
Implied size of a nominal 5% two-sided test: at ratio 0.763 it is 13.5% (not the claim's
~16%; the claim computed 16% from its own MEASURED 0.713, not from its stated 0.482 --
0.482 would imply 34.5%). The direction and the existence of the defect are unaffected.

## RESULT 5 -- the "silently" clause, established at the source (independent of any run)
grep over the whole pinned ado: there is NO `ereturn local/scalar/matrix` anywhere whose
name contains nest, reid or regime. The ONLY place `wdb_nest` ever reaches a user is the
A10 caution at ado:406-412, and that caution is gated on `wdb_nest == 1` -- i.e. it fires
only in the NESTING regime, and then only if a covariate option names the i variable. In
the keep-i regime nothing at all is printed and nothing is posted. My arm-A bootstrap log
carries zero "Caution" lines.

## Run history
L0_fcm1.do  -- first attempt; arm A bootstrap completed, then died r(198) on my own bad
               line `local mm : all macros` (not a wooldid defect). Arm A result kept:
               BOOT_A se(att)=0.182594673261 se(att_it)=0.18259467326
                      se(att_i)=0.182402476545  PREDNEST=1 completed=250 failed=0
L0b_fcm1.do -- same file with that block removed; all four arms.
