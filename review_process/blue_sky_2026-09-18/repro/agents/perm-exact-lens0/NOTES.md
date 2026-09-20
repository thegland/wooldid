# PERMEX-PX-2 — third lens (independent reproduction). WORKING NOTES

## What I verified by reading the pinned build myself (md5 444892481be8d50ec962a102d6abef51)
- ado:15992 `if (exact | memb) { plo = m0/n; phi = m1/n } else { add-one }` — plain share on the OR.
- ado:16127 `exk = (M < . & sum(s :< .) == M)` — COMPLETENESS of the enumeration, not membership.
- ado:16130 `if (M < . & !memb & sobs < . & sum(s :< .) >= 2) nonm = nonm + ...` — the wdp_addone
  flag keys on `!memb` ALONE; it does NOT require `!exk`.
- ado:8297 caution asserts "their p-values use the add-one form (1 + count)/(n + 1)".
- ado:16028 (permtail's own comment) says "The only element under permexact that still gets a
  binomial SE is one whose observed row is genuinely absent - THE SAME ELEMENT THAT TAKES THE
  ADD-ONE FORM."  Rows 16-18 are gated on `!exact & !memb`.  So the add-one branch, the binomial-SE
  branch and the caution are all supposed to describe ONE class; the caution's gate is the odd one out.
  => The textual/code contradiction is CONFIRMED, independently, and it is three-way
  (16130 vs 15992 vs 16018-16028), not two-way.

## Reachability — what I closed by reading
- permexact is only reachable under permute(shuffle) (ado:5924 refuses permexact with
  permute(probabilities)), so orow is set by the enumerator, never 0.
- ado:13716-13721: the identity column is searched over ALL cols(A); wooldid_permx_msp is a correct
  distinct-multiset-permutation enumerator, so cols(A)==M and 1 <= orow <= M always.
- ado:782 refuses a ttre that is not constant within i.
- ado:5979 refuses a cluster that mixes onset values or mixes treated/never-treated rows.
  => the identity assignment writes back EXACTLY the observed onset on every design-sample row.
- ado:5827 refuses permute()+unconditionalse; ado:4668 refuses customvce() together with cluster(),
  and permute() REQUIRES cluster() (measured, run 1) — so both of the review brief's
  "not bit-reproducible run to run" estimation paths are unreachable under permute().

## Measured, run 1 (L0_px2.do / .log) — determinism of the fit
Two back-to-back identical wooldid calls, e(b) compared bit for bit:
  plain 0 diffs, fastnose 0, cluster(i) 0, robust 0, all-four-estimands 0,
  aweights 0, poisson/ppmlhdfe 0, customvce(cluster i) 0.
=> a wooldid fit is a bitwise-deterministic function of its input on every path I could reach.
Also measured: permute() requires cluster(); customvce+cluster is refused, so
customvce()+permute() is unreachable.

## Run 2 (L0_px2b.do) — in flight
Adds (a) INTERLEAVED determinism (identity fit, 24 shuffled-onset fits, identity fit again),
(b) 9 permexact designs with a runtime detector: an element in e(perm_addone) but NOT in
e(perm_incomplete) is exactly exk=1 & memb=0.
