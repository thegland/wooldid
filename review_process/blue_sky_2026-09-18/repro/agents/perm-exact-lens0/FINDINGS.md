# PERMEX-PX-2 — third lens (INDEPENDENT REPRODUCTION). VERDICT: PARTIAL, severity low, not realistic.

Pinned build md5 444892481be8d50ec962a102d6abef51. All line numbers below are EXACT on it
(the prober's 15992 / 16147 are off by a few lines; the sites are the ones they meant).

## A. THE CODE-LEVEL CONTRADICTION: CONFIRMED, and it is FOUR-WAY, not two-way.
Four sites are supposed to describe one class of element. Three of them key on the same
predicate; one does not.

  ado:15982  `if (exact | memb) { plo = m0/n; phi = m1/n } else { (1+m)/(n+1) }`
             -> the plain share is taken on `exact | memb`.
  ado:16029  `if (!exact & !memb) { out[16..18] = binomial se + Clopper-Pearson }`
             -> the Monte-Carlo se / CI companions are filled on the exact complement.
             Its own comment (ado:16026-16028): "The only element under permexact that still
             gets a binomial SE is one whose observed row is genuinely absent - THE SAME
             ELEMENT THAT TAKES THE ADD-ONE FORM."
  ado:16130  `if (M < . & !memb & sobs < . & sum(s :< .) >= 2) nonm = nonm + <name>`
             -> the wdp_addone flag keys on `!memb` ALONE. It is missing the `!exk &` conjunct
             that 15982 and 16029 both carry.  ado:16127 sets exk = COMPLETENESS.
  ado:8297   the caution printed for that list: "...their p-values use the add-one form
             (1 + count)/(n + 1) rather than an exact share".
  ado:11681  the exact-run FOOTER, which the prober and skeptic both missed:
             "e(perm_p_se) and e(perm_p_ci) are left missing (AN ELEMENT NAMED BELOW AS TAKING
             THE ADD-ONE FORM CARRIES BOTH)."  Gated on 16029, i.e. on `!exact & !memb`.

So for an element with exk=1 and memb=0 the run would (i) take the plain share, (ii) leave
e(perm_p_se)/e(perm_p_ci) missing, and (iii) print two separate user-visible statements saying
the opposite of both. The one-line fix is `!exk &` at ado:16130.

The same `exact | memb` OR sits in wooldid_permqf (16211-16247), wooldid_permstep (16283-16320),
wooldid_permpw (16161-16172) and wooldid_permcellwire (16495-16497); none of those feeds
wdp_addone, so they carry the latent state with no message at all.

## B. REACHABILITY: NOT reached. Closed further than the skeptic closed it.
exk=1 & memb=0 needs the identity replication to complete for the element and to carry a
statistic that is not bitwise equal to the rep-0 comparator's.

Closed by READING (each verified on the pinned build):
 - permexact only exists under permute(shuffle): ado:5924 refuses permexact under
   permute(probabilities), which is the one generator that returns before setting wdp_obsrow.
   So orow is never 0 under permexact.
 - ado:13716-13721 scans every enumerated column for the identity; wooldid_permx_msp is a
   correct distinct-multiset-permutation enumerator so cols(A)==M, hence 1 <= orow <= M and the
   `orow > rows(Am)` escape in wooldid_permmemb (ado:15940) cannot fire either.
 - ado:782 refuses a ttre not constant within i; ado:5979 refuses a cluster that mixes onset
   values or mixes treated and never-treated rows. Together these make the identity assignment
   write back EXACTLY the observed onset on every design-sample row - there is no data
   difference for the identity replication to inherit.
 - Both estimation paths the review brief flags as "not bit-reproducible run to run" are
   unreachable under permute(): ado:5827 refuses permute()+unconditionalse, and permute()
   REQUIRES cluster() (measured, L0_px2.log) while ado:4668 refuses customvce() together with
   cluster(). This is the route I most expected to work and it is shut.

Closed by MEASUREMENT (mine):
 - L0_px2.do: two back-to-back identical wooldid calls, e(b) compared bit for bit, on ten
   paths: plain, fastnose, cluster(i), robust, all four estimands, aweights, poisson/ppmlhdfe,
   customvce(cluster i). NDIFF=0 everywhere.
 - L0_px2b.do PART 1: the same comparison with 24 intervening fits on REASSIGNED onsets between
   the two identity fits - the actual shape of the replication loop - on cluster(i), the
   20-column all-estimands+event-study spec, poisson, contreat(dose), and fastnose. NDIFF=0
   everywhere, including the 20-column vector.
   => a wooldid fit is a bitwise-deterministic function of its input on every path reachable
   under permute(). That is precisely the property memb=1 needs.
 - L0_px2b.do PART 2: nine permexact designs (M=105, and M=3 under permstrata) with a runtime
   detector - an element in e(perm_addone) but NOT in e(perm_incomplete) is exactly exk=1 &
   memb=0. Designs: base; 20-element all-estimands + esprelength(3) espostlength(3) +
   jointtests; poisson; aweights; contreat+contreatcontrols(basic); permtstat; unbalanced
   panel; permstrata; a panel with a block of missing y. Every one: ADDONE=[], INCOMP=[],
   FALSECAUTION=[], no p equal to 0, FAIL=0, BLANK=0.

## C. ON THE PROBER'S EVIDENCE
Their OBSERVED is a unit test of kernels lifted out with sed and fed inputs (30 complete draws,
an observed statistic one ulp below all of them, exact=1, memb=0) that no wooldid run produces.
It is sound evidence about wooldid_permtail and no evidence about the command. I agree with the
skeptic on that. I do not agree that it leaves nothing: it correctly shows what the branch WOULD
do, and the four-site inconsistency at ado:16130 is real whether or not it is reachable.

## D. WHAT WOULD CHANGE MY MIND
A single path under permute() where two identical wooldid fits differ in any bit, or any way to
get orow outside 1..M, or a route that lets the identity replication see different data. I
attacked all three and found none. Note also that even if it did fire, the damage is bounded:
under a COMPLETE enumeration the identity replication IS in the reference set, one ulp away from
sobs, so the one-sided SHARE is understated by exactly 1/M. Under the default equal-tailed rule
p = 2*min(plo, phi) that is 2/M on the reported p (in my designs, 2/105 reported as exactly 0);
under permside()/permabsp/permrawp it is 1/M.

## E. ONE OPEN EDGE I DID NOT CLOSE
Nobody - prober, skeptic, or me - has ever made e(perm_addone) fire at all. There is a route I
identified and did not spend a run on: an element whose value is MISSING in the identity
replication but present in others (the ame `blank` mechanism at ado:8283, or a degenerate
observed se under permtstat). That gives exk=0 AND memb=0, i.e. the add-one form really is taken
and the caution at ado:8297 fires TRUTHFULLY. So I cannot say whether the caution is dead code
or has a reachable correct domain; I can only say that the over-broad term it adds - exk=1 with
memb=0 - is unreachable. The fix (`!exk &` at ado:16130) is the same either way, and under the
ame/blank route it changes nothing.
