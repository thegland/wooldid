# SKEPTIC pass on FW-01 / FW-02 / FW-03 (fweights / weight surface)

Build: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build (adopath ++)

## Zero-cost code reads done first (no Stata)
- ado:555 `syntax ... [aw fw pw / ]` -- the `/` means `exp` holds the RAW expression text.
- ado:607 `local w `exp''  -> confirmed: the raw text is carried.
- ado:1619 `gegen `w_per_it' = total(`esamp' * `w')` -- so the DENOMINATOR is also
  text-pasted, not just the numerator the seeded finding names. `total(esamp*wf+wm)`.
- ado:1630/1637/1644/1650 numerator sites as claimed.
- ado:628-632 fweight refusal; ado:633 demotes fweight->aweight in wdrepmode only AFTER
  `weightinserter' was built at 619, so in replication mode the REGRESSION really does get
  an fweight.
- ado:1448 `if "`robust'"=="robust" local vcespecifier vce(robust)` -- under [pw] with no
  vce option `vcespecifier' is EMPTY, so reghdfe is handed [pw=w] bare and makes its own
  robust V. Confirms FW-02's mechanism.
- ado:4291-4301 the vce-label ladder: no pweight branch. Confirms FW-02's mechanism.
- ado:4676 unconditionalse gate: `& "`weight'" != "pweight"` -- sthlp:268's claim that
  pweights satisfy the robust requirement IS implemented. (Killed my candidate new finding
  N1 for free.)
- ado:4623 the aweight caution under unconditionalse is aweight-only; consistent.
- sthlp syntax line says `[aw/pw=weights]` -- fw is NOT documented as accepted.

## Probes
(see per-probe sections below, appended as they run)

## P1 (p1_exprgrid.do/.log) -- FW-01(a) attribution
Build asserted: which wooldid -> the pinned build path, version 1.4 2026-08-19.
Grid of 9 weight expressions x {wooldid default, wooldid nofastcollapse, wooldid
customvce(robust), reghdfe direct, regress direct}:
  expr          wd-default  wd-nofc  wd-customvce  reghdfe  regress
  wtot           0           0        0             0        0
  2*k            111         111      0             111      0
  k/2            198         198      0             198      0
  k*2            111         111      0             111      0
  int(k)         101         101      0             101      0
  wf+wm          198         198      0             198      0
  (wf+wm)        198         198      0             198      0
  1/pscore       198         198      0             198      0
  cond(k>2,k,1)  101         101      0             101      0
=> The default-path refusal is INHERITED FROM reghdfe verbatim (same rc, same message:
   "+ invalid name", "variable 2*k not found", "1 invalid name"). nofastcollapse does not
   change it, so the collapse is NOT the fix site. wooldid's own paths (customvce/regress)
   accept every expression. `regress` accepts all nine.
   REALISM NOTE: [pw=1/pscore] (inverse-probability weighting) is the single most common
   weight expression in applied work and it dies r(198) "1 invalid name" on the DEFAULT path.

## P2 (p2_mechanism.do/.log) -- FW-01(b) mechanism, my own asymmetric test
[aw=wf+zero] (zero==0 everywhere): att, att_it, att_i, att_itime ALL BITWISE IDENTICAL
  (+1.76fb4bd29be62X+000 = 1.464771975439270), SEs identical too.
[aw=wf]     : 1.464771975, 1.458344351, 1.415271428, 1.417598312 -- four distinct values.
[aw=zero+wm] == [aw=wm] BITWISE on all four estimands and all four SEs.
=> exactly the asymmetric signature of verbatim text pasting. CONFIRMED.
[aw=wf+wm] vs [aw=wtot] (wtot==wf+wm asserted), customvce(robust):
  att      1.470710677532787 vs 1.470710677532787  (exact)
  att_it   1.468148949125649 vs 1.461823093540426  reldiff 4.33e-03
  att_i    1.468059193713792 vs 1.425732957679875  reldiff 2.97e-02   <-- 3%
  att_itime1.468041600513423 vs 1.420357846860070  reldiff 3.36e-02   <-- 3.4%
  (larger than the seeded finding's max 5.1e-03)
customweightmultiplier(cwm) + [aw=wf+wm] vs [aw=wtot]: att 1.464679018 vs 1.455847603,
  reldiff 6.07e-03 -- so even att breaks under cwm, as claimed.
unconditionalse robust safety(off): att exact; att_it 1.468148949 vs 1.461823094 and
  SE(att_it) 0.143077097 vs 0.126915925 -- a 12.7% SE gap.

## P3 (p3_pw_scale_collapse.do/.log)
FW-02 CONFIRMED and mechanism verified: [pw=pop] e(V) vs [aw=pop] robust e(V):
  mreldif = EXACTLY 0.  vs [aw=pop] iid: mreldif 4.5e-3.  e(vce)=[ols], e(vcetype)=[.]
  under pw.  reghdfe direct with [pw]: robust/Robust.  regress direct with [pw]: robust/Robust.
  Header prints "Standard Errors: " (blank).
NEW (negative): weight-SCALE invariance holds. [aw=pop] vs [aw=pop1000] (=1000*pop as a
  pre-made variable): all four estimands and SEs agree to 1-2 ulp (e.g. att
  +1.473577ab01cc8 vs +1.473577ab01cc6). Same under [pw]. ado:1655's claim stands.
NEW (negative): default cell collapse vs nofastcollapse with a weight that VARIES WITHIN
  the cohort-by-period cell (i=state, 8 counties per state-year, pop lognormal): all four
  estimands and all four SEs agree to <=1 ulp under iid, robust AND cluster(i).

## P4 (p4_fw_collapsed.do/.log) -- FW-03
Plain-regress oracle, 8717 microdata rows -> 360 (state,year) cells:
  micro cluster(state)  se ratio 1.000  (b reldif 2e-14)
  [fw=n] cells cluster  se ratio 1.000000000  <-- fw DOES reproduce microdata cluster SE
  [aw=n] cells cluster  se ratio 1.060009592
  micro iid             1.000
  [fw=n] cells iid      0.192690664  <-- fw is CATASTROPHICALLY wrong under iid
  [aw=n] cells iid      1.006426015  <-- aw is nearly exact under iid
  => the aw/pw ratio 1.060009592 equals EXACTLY
     sqrt( [(360-1)/(360-42)] / [(8717-1)/(8717-42)] ) = 1.060009...
     i.e. the ENTIRE gap is Stata's (N-1)/(N-k) small-sample multiplier. Meat and bread
     are identical. This is generic Stata aweight behaviour, not anything wooldid does.
wooldid legs: [fw=n] -> "Error: wooldid does not accept fweights. Specify aweights or
  pweights." r(198). [aw=n]==[pw=n] cells bitwise. att b matches microdata to ~1e-14.
  att SE cells/micro under cluster(state) ~ 1.259 (bigger than the plain-regress 1.060
  because reghdfe's K is much larger relative to N=360).
esttab check: e(vcetype)=="Robust" on a clustered regress, yet `esttab, se` prints
  "Standard errors in parentheses", NOT "Robust standard errors". So FW-02's stated
  downstream consequence (estout/esttab printing a false "Robust" footer) is NOT what
  esttab does by default -- the label is wrong, but the named harm is not demonstrated.

## P5 (p5_newprobes.do/.log)
(B) sthlp:710 poisson weight-type gates: [pw=]+poisson rc 0; [aw=]+poisson r(198)
    "Error: poisson only accepts pweights, not aweights."; safety(off) does not lift it.
    MATCHES the help. Negative result.
(C) zero survey weights (1760 of 8717 rows with w==0): wooldid reports N=6957 and
    e(N_dropped)=1760, and the run is BITWISE identical to the hand-restricted
    `if w0>0' run. regress [aw=w0] also reports N=6957. Transparent. Negative result.
(A) my df-multiplier prediction used e(df_r), which wooldid overrides to G-1=29 on both
    legs, so that leg of the test was uninformative. The plain-regress proof in P4 stands.

## P6 (p6_pw_engines.do/.log)
[pw] + unconditionalse WITHOUT robust/cluster RUNS (ado:4676's pweight exemption is real);
  [aw] in the same call is refused r(198). Matches sthlp:268. Negative result.
[pw] + bootstrap(standard) reps(25) cluster(state): runs, e(vce)=bootstrap.
[pw] + permute(shuffle): runs.
NEW (positive, mine): unconditionalse + [pw] + `robust' and unconditionalse + [pw] WITHOUT
  `robust' produce the SAME e(V) bitwise (se +1.246441c7f884cX-002 both) but are LABELLED
  DIFFERENTLY: e(vce)=[robust]/e(vcetype)=[Robust] in the first, e(vce)=[unconditional] in
  the second. ado:4291's ladder tests `robust' before `unconditionalse', so on the
  unconditionalse path e(vce) names the STAGE-1 regression's vce, not what e(V) is -- the
  same ado:4283 contract FW-02 is about. Cosmetic, but it is the same defect class.

## P7 (p7_esttab.do/.log) -- redo of the esttab check, single-model tables
The p4 esttab test put a Robust model and an iid model in ONE table, which could have
suppressed a vcetype-keyed footer. Redone with three SEPARATE single-model tables:
  esttab rclu  (e(vcetype)=Robust, cluster) -> "Standard errors in parentheses"
  esttab rrob  (e(vcetype)=Robust, robust)  -> "Standard errors in parentheses"
  esttab riid  (e(vcetype)=.,      iid)     -> "Standard errors in parentheses"
So esttab's DEFAULT footer does not read e(vcetype). But `estout rrob, stats(vcetype N)`
DOES print a "vcetype  Robust" row on request. Conclusion unchanged: FW-02's label claim is
true, its automatic-publication harm is not; a user who asks for the vce row or reads
ereturn list does get the false "ols". FW-02 stays PARTIAL / low.
