# SKEPTIC notes -- refuting the "boot-reproduce" findings (BR-1..BR-4)

Build asserted pinned: probe_br1.log:59 shows .../scratchpad/bluesky/build/wooldid.ado
I read the original prober's do-files but ran NONE of them. Every number below is from
my own do-files in this directory on my own generated data.

## VERDICTS (all four reproduced independently -> none refuted)
- BR-1 permsave comma   CONFIRMED  low      probe_br1.do/.log
- BR-2 reps caution     CONFIRMED  medium   probe_br2.do/.log
- BR-3 unclustered legs CONFIRMED  low      probe_br3.do/.log
- BR-4 smp vs e(sample) CONFIRMED  high     probe_br4.do/.log + probe_br4c.do/.log
  (one sentence of their EXPECTED corrected: the percentile family is NOT immune)

## NEW (mine, seeded=false)
- SK-1 the RF-6 comma refusal's printed remedy is itself wrong  (probe_br2.log)
- SK-2 permute() shares the universe defect; 18 onset_ columns for 12 estimation
       clusters, 24/40 replications assign an onset to a data-less unit
       (probe_permuniv.do/.log)
- SK-3 sort-order reproducibility of the two legs  (probe_sort.do/.log) -- RUNNING
- (probe_save.do queued: e(bootdist) vs bootsave file vs rep-0 row)

## NEGATIVE RESULTS (probed hard, clean)
- e(boot_rngstate) replay works EXACTLY: mreldif 0 on the draw matrix, bit-identical
  e(V), after `set rngstate` (probe_rng.log:244)
- same rseed() twice: mreldif 0 (probe_rng.log:379)
- bootsave() characteristics present and right: _dta[wooldid_engine]=bootstrap,
  [wooldid_seed]=31415, [wooldid_rngstate] 5011 chars (probe_rng.log:460-468)
- e(boot_seed) correctly missing on a no-rseed run

## NOT INDEPENDENTLY VERIFIED
- BR-4's sub-claim that the count leg stops matching Stata's own `bootstrap` prefix
  replication for replication. Mechanism is clear (the prefix drops observations
  outside e(sample) before resampling unless `nodrop`), and it is the least
  consequential part of BR-4, so I spent the budget on the df/caution/realism legs
  instead. Flagged rather than assumed.

## KEY MECHANISM (source, pinned)
ado:682   qui replace `ifvar' = 0 if `i' >= . | `t' >= . | `y' >= . | `w' >= . | `w' <= 0
ado:704   same for controls()/ccc_absorb()/contreat()/interactivecontrols()/non-i cluster()
ado:615   c_local wdp_if `"`if'"'            <- the RAW if is what is forwarded
ado:8917  bootrun   smp   = ((`ifr') == 1)   <- markout never applied
ado:5957  permctl   insmp = ((`ifr') == 1) & !missing(`clus')
ado:7479  permrun   smp   = ((`ifr') == 1)
ado:663-670  the comment that claims "the estimation sample, the permutation design
             and the bootstrap frame agree" -- the invariant that is violated
ado:10233 ereturn scalar df_r = `wdbdf'  (= nclus - 1, nclus counted on smp)
ado:6456  `if "`clus'" != "" {`  -- opens the block the reps() caution is trapped in
ado:9301  count path: bsample, cluster() idcluster()  ("plain bsample diverges - A17")
ado:9345  reid path unclustered: plain `bsample'      <- the diverging form
ado:5866  "permsave(): the name rules are bootsave()'s" -- they are not (no comma check)
ado:6287  subinstr(",", " ") builds the wrong remedy in the bootsave comma message

## FINAL STATE (all probes complete)
probe_sort:  COUNTCL 1.641e-15 / REIDCL 1.745e-15 (order-invariant)
             COUNTUN 0.48457418 / REIDUN 0.77309209 (order-DEPENDENT)
probe_save:  e(bootdist) == bootsave file to <=1 ulp; _rep 0 row == e(b); chars present
probe_br4c:  e(boot_failed)=0 in all four routes (no failure inflation at this geometry)
probe_permuniv: ARM A 18 onset_ columns for e(N_clust)=12; onset_14 non-missing in 24/40
DEDUPE FLAG: another agent is running r02_phantom.do -- SK-2/BR-4 may overlap with theirs.
