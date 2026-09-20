# SKEPTIC pass on the e() contract findings (EC-1..EC-8)

Build: scratchpad/bluesky/build via `adopath ++`. All data generated here, own seeds.
Rule followed: every verdict rests on MY do-file, not the prober's.

## Probes
- p1_level.do / .log   -> EC-2   DONE
- p2_bootrank.do/.log  -> EC-1   RUNNING (relaunched after a brace bug)
- p3_wt_sg.do/.log     -> EC-6, EC-7  DONE
- p4_replay.do/.log    -> EC-3, EC-5, EC-8  RUNNING
- p5_esV.do/.log       -> EC-4   QUEUED

## Verdicts so far

### EC-2  CONFIRMED, and the failing set is ~40x larger than reported
p1_level.log. Mechanism reproduced end to end:
  `local x = 99.9` -> |99.90000000000001| ; `set level` on that text -> r(198).
wooldid itself: level(90.1)=r198, level(91.1)=r198, level(94.9)=r198, level(99.9)=r198,
level(99.99)=r198, while 90/92.5/97.5/99/99.5/68.3/93.7 all rc 0.
regress and `set level` accept ALL of them (LITERAL block, rc 0 for 12/12).
EXHAUSTIVE: of the 90 one-decimal levels 90.0-99.9, **40 fail** (90.1 90.4 90.6 90.9
91.1 91.4 91.6 91.9 ... 99.9); of the 1000 two-decimal levels in [90,100), **320 fail**.
The prober's 3-of-12 grid badly understates the breadth.
Session route reproduced: `set level 99.9` (rc 0) then wooldid with NO level() -> r(198);
e() keeps the previous run (e(level)=92.5, cols lb_92p5ci/ub_92p5ci).

### EC-6  CONFIRMED, bit-exactly
p3_wt_sg.log. Same panel, att + ES:
  unweighted se hex +1.583ad4d457ba4X-003
  [aw=pw]        se +1.57f3b3c758339X-003   e(vce)=ols
  [pw=pw]        se +1.00814b4a093acX-002   e(vce)=ols   e(vcetype) absent
  [aw=pw] robust se +1.00814b4a093acX-002   e(vce)=robust  <-- BIT-IDENTICAL to [pw=]
  [pw=pw] robust se +1.00814b4a093acX-002
  [aw=pw] cluster == [pw=pw] cluster (hex +1.2eeae120fd10dX-002)
e(wtype)/e(wexp) never posted on any weighted run (grep of the ado: zero `ereturn local
wtype|wexp`). Control: `reg y i.t [pw=pw]` posts wtype=pweight wexp="= pw" vce=robust.

### EC-7  CONFIRMED, one line
p3_wt_sg.log: lincom main_att_s1 -> r(111); lincom main_SG1_att -> rc 0;
pretreatment_att_s1 r(111) / pretreatment_SG1_att rc 0; ES_att_post0_s1 r(111) /
ES_att_SG1_post0 rc 0. e(numpre_att_s1) EXISTS (so _s# is right everywhere else).
grep -c '_SG' wooldid.sthlp = 0. sthlp:1307 is the defective sentence.

### EC-1  (partial data) rank(e(V)) == rank(e(V_boot)) == rank(centered e(bootdist))
at reps 5 (4) and reps 8 (7), i.e. exactly reps-1. So wooldid is NOT dropping
information when it posts; the rank is the honest rank of the replication cloud.
Waiting on reps 40/60 on IRREGULAR data to see whether "saturates at 10" survives.
