# datashapes recon — NOTES (final)

DELIVERABLE (written): .../bluesky/recon/data_shapes.md
GENERATORS (canonical): .../bluesky/recon/gen/{gen_core,gen_shapes,gen_mutate,gen_truth}.do
                        identical copy in ./gen/  (shasums matched after every edit)

## status: COMPLETE
- gen_core.do    engine, ~38 knobs (gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0
                 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub
                 gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype
                 gwspread gwzero gwneg gwmiss gymiss gxmiss gpois)
- gen_shapes.do  61 named shapes + GOVR "name value ..." override mechanism
- gen_mutate.do  15 mutators
- gen_truth.do   design truth of the four estimands
- runs: run_build.do (61/61 build, 15/15 mutate), run_smoke.do (57 shapes x 1 wooldid call),
        run_extra.do (onepre/onsetgap/dose+cc/i=cohort/espre/fcol/awpw), run_het.do,
        run_het2.do (exact recovery unbal + hetero_unbal), run_final.do (i=cohort harness,
        four-distinct-estimand config)

## the measured results that matter
- EXACT RECOVERY works: canon_exact, unbal (7.03e-16), hetero_unbal (1.6e-15),
  rcs+gtauhet (3.9e-16 on all four), hetero_unbal at i=cohort with gaisd 0 (1.9e-16).
  => the "t == mint forced reference" rule does NOT bias points on an unbalanced panel.
- BEST PROBE CONFIG: SHAPE rcs, GOVR "gysd 0 gtauhet 0.4" -> att 1.135377607,
  att_it 1.230907782, att_i 1.423516370, att_itime 1.422383070, all exact.
- fcol vs nofastcollapse on rcs: BIT-IDENTICAL, with and without within-cell controls.
- [aw=w] vs [pw=w] point estimates: BIT-IDENTICAL.
- espre(3) espost(3) moves att by 3.086e-02 relative (documented obliquely in the escomp entry).
- admin (971 units, i=unit): r(907) maxvar, raw Stata message. i=cohort: rc 0.
- poisson (120 units, i=unit): r(322) perfect separation. i=cohort: rc 0.
- ed_tscale (t %td vs ttre a year): r(459) whose message does NOT name the scale mismatch.
- ed_badw: "Error: weights must not be negative." r(411). wzero: silent 111-row drop.
- sgvary (subgroup varying within unit): r(459), informative.

## ado mechanisms cited in the deliverable
weights 1594-1650; att subpop :2030; pre subpop :2089; _i = user's i :927/:12181; fit :1462;
reltime forced reference :12079, nopre drop :12046, BS-33 min_t :12170 vs abs() twin :889;
fcol key :1902, eligibility :5068-5125; markout :670; if-rule :667; integer grid :5341;
float 2^24 :4967; F13 :775-790; 32-char name skips :6700 :8097 :8814 :8821 :9990 :10026 :10911;
replay gate :26-52.
