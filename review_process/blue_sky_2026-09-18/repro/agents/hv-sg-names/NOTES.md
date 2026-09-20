# hv-sg-names — verification notes

## Claim under test
1. sthlp:1307: "Under subgroups() each subgroup's rows take the suffix _s#."
   Code (ado:4015) injects an INFIX `_SG<value>_` instead.
2. sthlp:1305-1306: "for the other estimand types the same names with att replaced by
   att_it, att_i, att_itime or ame" implies ES_att_i_* / ES_att_itime_* exist.
   Code (ado:6828-6835 comment) says att_it / att_i / att_itime SHARE one ES equation
   named ES_att_it.

## Static evidence gathered (before running)
- ado:4015  `local rownamespool_sginj = subinstr("`rownamessg`sgcase''",":","_SG`sgcase'_",.)`
- ado:10910 `local wdspx = cond("`wdsg'" == "", "ES_`wdsb'_", "ES_`wdsb'_SG`wdsg'_")`
- ado:11061 same for e(V) prefix
- ado:14638 Mata: `return("ES_" + regexs(1) + "_SG" + regexs(2) + "_" + nm)`
- ado:6828-6835 AUTHOR COMMENT: "the event-study row pool gives a SINGLE equation name,
  ES_att_it, to att_it, att_i and att_itime alike ... Adding ES_att_i / ES_att_itime to
  the family list finds no columns"
- ado:11024 AUTHOR COMMENT: "the row equations of e(wdidsgresults_3) ... are NOT suffixed
  (the _SG3_ injection touches only the stacked e(b)/e(V) names)"
- sthlp:1315 says the ES matrices are `e(b_es_att_s{sgcase})` -> a THIRD scheme
- **`grep -c '_SG' wooldid.sthlp` == 0**: the true e(b) convention appears NOWHERE in
  the help file.

## Run log
- p1_names.do run 1: crashed r(198) on my own bad `local mn : r(dummy)`. Fixed.
  Partial result already confirmed:
  RUN A_att (att only, esprelength(2) espostlength(2)):
    main_att, pretreatment_att, ES_att_pre2, ES_att_pre1, ES_att_post0, ES_att_post1
  (matches the help exactly for the no-subgroup att case)

## Result summary (all probes complete)
BOTH halves of the recon claim CONFIRMED, with the mechanism as stated.

p1_names.do  -> e(b) names across estimands x subgroups; consequence leg (test/lincom/esttab/coefplot)
p2_consumers.do -> bootstrap tables + bootsave .dta varnames; negative/non-integer sg values
p3_perm_and_mats.do -> permute tables; full e() matrix inventory; sg values 3/17
p4_ame_and_silent.do -> ame family names on an identified dose; event_plot; silent-consumer sweep
p5_semantics_and_overflow.do -> att_i-alone semantics; 7/8-digit subgroup-value name overflow
p6_sharpen.do / p6b_coefplot.do -> ES_att vs ES_att_it magnitude; rigorous coefplot drop() count

Key numbers:
- C_sgA e(b) = main_SG1_att ... ES_att_it_SG2_post1 (48 cols). NEVER `_s#`.
- grep -c '_SG' wooldid.sthlp == 0
- test main_att_s1 = main_att_s2 -> r(111);  test main_SG1_att = main_SG2_att -> rc 0
- lincom ES_att_i_post0 -> r(111); lincom ES_att_it_post0 -> rc 0
- coefplot drop(*_s1*): 12 of 12 coefs still plotted, rc 0, no message
- ES_att_post0 = 0.65035 vs ES_att_it_post0 = 1.23917 on the unequal-cohort design

Self-caught non-findings (my own DGP mistakes, corrected):
- p2/p3 run O: a unit-constant dose under i=unit is absorbed -> documented amevoid path,
  NOT a missing-ame-ES defect. Redone with an obs-level dose in p4 RUN P: all ame names present.
- p5 first attempt: never-treated y was missing because coheff was left `.`; fixed.
